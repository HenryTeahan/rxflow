"""Bond movement filtering and reuse of existing endpoint conformers."""
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from rdkit import Chem

from hpc_workflow import connect_database, filter_endpoint_conformers
from rxflow_endpoints import connectivity_changes, mapped_molecule


def geometry(smiles, distance, shift=0.):
    molecule = mapped_molecule(smiles)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.SetAtomPosition(0, (shift, 0., 0.))
    conformer.SetAtomPosition(1, (shift + distance, 0., 0.))
    molecule.AddConformer(conformer)
    return molecule


class BondMovementTests(unittest.TestCase):
    def test_movement_is_checked_without_a_perceived_graph_change(self):
        with patch('rxembed.pipeline.metrics.connectivity', return_value=([], [])):
            for before, after in ((1., 2.6), (2.5, 4.1), (5., 3.4)):
                with self.subTest(before=before, after=after):
                    changes = connectivity_changes(geometry('[C:1][C:2]', after),
                                                   geometry('[C:1][C:2]', before))
                    self.assertEqual(len(changes['distance_changes']), 1)
                    self.assertAlmostEqual(changes['distance_changes'][0]['change_angstrom'], after - before)

    def test_exactly_one_point_five_is_allowed_and_greater_is_rejected(self):
        with patch('rxembed.pipeline.metrics.connectivity', return_value=([], [])):
            before = geometry('[C:1][C:2]', 1.5)
            self.assertFalse(any(connectivity_changes(geometry('[C:1][C:2]', 3.), before).values()))
            self.assertTrue(connectivity_changes(geometry('[C:1][C:2]', 3.001), before)['distance_changes'])

    def test_rigid_translation_does_not_change_bond_lengths(self):
        with patch('rxembed.pipeline.metrics.connectivity', return_value=([], [])):
            self.assertFalse(any(connectivity_changes(geometry('[C:1][C:2]', 1.5, shift=100.),
                                                      geometry('[C:1][C:2]', 1.5)).values()))

    def test_lost_coordination_is_detected_even_below_the_movement_limit(self):
        # A 1.25 A movement can still release an initially bound donor.
        changes = connectivity_changes(geometry('[Pd+2:1]<-[NH3:2]', 3.4),
                                       geometry('[Pd+2:1]<-[NH3:2]', 2.15))
        self.assertFalse(changes['distance_changes'])
        self.assertEqual(changes['coordination'][0]['left'], [1])

    def test_intact_metal_bond_above_one_point_five_and_two_is_allowed(self):
        changes = connectivity_changes(geometry('[Pd+2:1]<-[NH3:2]', 2.24),
                                       geometry('[Pd+2:1]<-[NH3:2]', 2.15))
        self.assertFalse(any(changes.values()), changes)


class SavedConformerFilterTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.database = self.root / 'workflow.sqlite'
        self.smiles = '[Pd+2:1]<-[NH3:2]'
        with closing(connect_database(self.database)) as db:
            for job_id, ligand_id, distance, energy in ((1, '11', 3.99, -100.),
                                                       (2, '11', 2.24, -99.99),
                                                       (3, '12', 3.99, -100.)):
                initial, optimized = f'input_{job_id}.xyz', f'optimized_{job_id}.xyz'
                (self.root / initial).write_text(Chem.MolToXYZBlock(geometry(self.smiles, 2.15)))
                (self.root / optimized).write_text(Chem.MolToXYZBlock(geometry(self.smiles, distance)))
                db.execute("""INSERT INTO jobs (id,ligand_id,pair_index,seating,endpoint,conformer,
                    status,stage,input_xyz,xyz_file,smiles,charge,settings_json,electronic_energy_eh,
                    gibbs_free_energy_eh,retained) VALUES (?,?,1,4,'reactant',?,'succeeded','complete',
                    ?,?,?,2,'{}',?,?,?)""", (job_id, ligand_id, job_id, initial, optimized, self.smiles,
                                           energy, energy + .01, int(job_id == 1)))

    def rows(self):
        with closing(connect_database(self.database)) as db:
            return {row['id']: dict(row) for row in db.execute('SELECT * FROM jobs')}

    def test_dry_run_is_read_only(self):
        before = hashlib.sha256(self.database.read_bytes()).hexdigest()
        report = filter_endpoint_conformers(self.root, ligand_id='11', dry_run=True)
        self.assertEqual(report['rejected_ids'], [1])
        self.assertEqual(report['selections'][0]['selected_job_id'], 2)
        self.assertEqual(hashlib.sha256(self.database.read_bytes()).hexdigest(), before)
        self.assertFalse((self.root / 'endpoint_filter_audits').exists())

    def test_reuses_passing_conformer_preserves_history_and_unrelated_ligands(self):
        before = self.rows()
        report = filter_endpoint_conformers(self.root, ligand_id='11')
        after = self.rows()
        self.assertEqual(report['rejected_ids'], [1])
        self.assertEqual(report['queued_for_rerun'], [])
        self.assertEqual(after[1]['status'], 'failed')
        self.assertEqual(after[1]['stage'], 'connectivity_filter')
        self.assertEqual(after[1]['gibbs_free_energy_eh'], before[1]['gibbs_free_energy_eh'])
        self.assertEqual(after[2]['status'], 'succeeded')
        self.assertEqual(after[2]['retained'], 1)
        self.assertEqual(after[3], before[3])
        with closing(sqlite3.connect(report['backup_file'])) as backup:
            self.assertEqual(backup.execute('SELECT status FROM jobs WHERE id=1').fetchone()[0], 'succeeded')
        self.assertEqual(json.loads(Path(report['report_file']).read_text())['rejected_ids'], [1])

    def test_no_valid_alternative_prepares_original_embedding_for_rerun(self):
        with closing(connect_database(self.database)) as db:
            db.execute("UPDATE jobs SET status='failed' WHERE id=2")
        report = filter_endpoint_conformers(self.root, ligand_id='11')
        after = self.rows()
        self.assertEqual(report['queued_for_rerun'], [1])
        self.assertEqual(after[1]['status'], 'pending')
        self.assertEqual(after[1]['stage'], 'opt')
        self.assertEqual(after[1]['input_xyz'], 'input_1.xyz')
        self.assertIsNone(after[1]['gibbs_free_energy_eh'])
        self.assertIsNone(after[1]['xyz_file'])


if __name__ == '__main__':
    unittest.main()
