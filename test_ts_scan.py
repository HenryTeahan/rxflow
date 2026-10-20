"""Check physical peak selection and correspondence with each candidate's own pair."""

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.io import write, read

from hpc_workflow import connect_database, reset_failed_transition_states, ts_recovery_status
from ts_scan import best_scored_candidate, interior_peak, scan_inputs, score_ts_scan


class ScanTests(unittest.TestCase):
    def frames(self, count=5):
        return [Atoms('CHN', positions=[[0, 0, 0], [1 + i*.1, 0, 0], [3, 0, 0]])
                for i in range(count)]

    def test_global_endpoint_maximum_does_not_hide_the_interior_local_peak(self):
        peak, motion = interior_peak(self.frames(), [10, 0, 2, 1, -1], 0, 1, 2)
        self.assertEqual(peak, 2)
        self.assertGreater(motion['C_H_change_angstrom'], 0)
        self.assertLess(motion['H_acceptor_change_angstrom'], 0)

    def test_highest_interior_frame_is_not_automatically_a_local_peak(self):
        with self.assertRaisesRegex(ValueError, 'No interior'):
            interior_peak(self.frames(), [10, 8, 6, 4, 2], 0, 1, 2)

    def test_highest_of_several_valid_local_peaks_is_selected(self):
        peak, _ = interior_peak(self.frames(7), [0, 1, 0, 3, 0, 2, 0], 0, 1, 2)
        self.assertEqual(peak, 3)

    def test_bad_or_reversed_motion_and_partial_scores_are_rejected(self):
        frames = self.frames()
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            interior_peak(frames, [0, 1, np.nan, 1, 0], 0, 1, 2)
        with self.assertRaisesRegex(ValueError, 'No interior'):
            interior_peak(list(reversed(frames)), [0, 0, 2, 0, 0], 0, 1, 2)

    def test_original_candidate_uses_its_ranked_pair_instead_of_group_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = base/'seeds/ligand_01/ts_guess'
            source.mkdir(parents=True)
            files = {}
            for key in ('profile_xyz', 'reactant_xyz', 'product_xyz'):
                path = source/f'candidate_2_{key}.xyz'
                path.write_text('candidate two')
                files[key] = str(path)
            (source/'pair_001_seating_04_selection.json').write_text(json.dumps(
                dict(candidates=[dict(rank=1), dict(rank=2, **files)])))
            job = dict(ligand_id='1', pair_index=1, seating=4, conformer=2,
                       input_xyz='ligand_01/ts_guess/pair_001_seating_04_conformer_002_input.xyz')
            found = scan_inputs(job, base/'results', base/'seeds')
            self.assertEqual({k: str(p) for k, p in found.items()}, files)

    def test_scoring_preserves_atom_order_files_and_energy_only_context(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frames = self.frames()
            for name, atoms in [('input.xyz', frames[1]), ('path.xyz', frames),
                                ('reactant.xyz', frames[0]), ('product.xyz', frames[-1])]:
                write(root/name, atoms, format='xyz')
            job = dict(id=1, attempt=1, input_xyz='input.xyz', charge=0, multiplicity=1,
                       smiles='[C:1]([H:2])#[N:3]', donor=0, hydrogen=1, acceptor=2,
                       ts_guess_json=json.dumps(dict(profile_xyz='path.xyz',
                           reactant_xyz='reactant.xyz', product_xyz='product.xyz')))
            original = (root/'input.xyz').read_bytes()
            result = SimpleNamespace(energies_ev=(10., 0., 2., 1., -1.), failures={},
                                     n_evaluations=5, n_gradients=0)
            with patch('rxts.scoring.score', return_value=result) as scorer:
                report = score_ts_scan(job, root, dict(cpus=1))
            self.assertEqual(scorer.call_args.kwargs,
                             dict(method='gxtb', solvent=None, charge=0, multiplicity=1, n_cores=1))
            self.assertEqual(report['selected_frame'], 2)
            self.assertTrue(np.allclose(read(root/report['guess_xyz']).positions, frames[2].positions))
            self.assertEqual((root/'input.xyz').read_bytes(), original)
            self.assertTrue((root/report['report_file']).is_file())

    def test_failed_scoring_retains_partial_evidence_without_publishing_guess(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            job = dict(id=1, attempt=1, input_xyz='input.xyz')
            with self.assertRaisesRegex(ValueError, 'Missing rxts'):
                score_ts_scan(job, root, {})
            report = json.loads((root/'job_1_run_01_gxtb_scan/scoring.json').read_text())
            self.assertFalse(report['valid'])
            self.assertNotIn('guess_xyz', report)

    def test_group_neb_chooses_lowest_energy_scored_candidate_with_matching_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            conn = connect_database(':memory:')
            self.addCleanup(conn.close)
            for index, energy in [(1, -10), (2, -20)]:
                paths = {key: f'{index}_{key}.xyz' for key in ('guess_xyz', 'reactant_xyz', 'product_xyz')}
                for file in paths.values(): (root/file).write_text('candidate')
                report = dict(valid=True, selected_energy_ev=energy, **paths)
                conn.execute("""INSERT INTO jobs (id,ligand_id,pair_index,seating,endpoint,conformer,
                    status,ts_strategy,input_xyz,smiles,charge,settings_json,ts_scan_json)
                    VALUES (?, '1',1,1,'TS',?,'failed','direct','unused','unused',0,'{}',?)""",
                    (index, index, json.dumps(report)))
            chosen = best_scored_candidate(conn, dict(ligand_id='1', pair_index=1, seating=1), root)
            self.assertEqual(chosen['source_job_id'], 2)
            self.assertEqual(chosen['reactant_xyz'], '2_reactant_xyz.xyz')


class ResetTests(unittest.TestCase):
    def test_reset_archives_history_and_preserves_success_active_and_sampling_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            db = connect_database(root/'workflow.sqlite')
            cases = [('failed', 'optts', 'direct'), ('failed', 'neb', 'neb'),
                     ('failed', 'sampling', None), ('succeeded', 'complete', 'direct'),
                     ('running', 'optts', 'direct'), ('failed', 'ts_guess', 'direct')]
            for index, (status, stage, strategy) in enumerate(cases, 1):
                db.execute("""INSERT INTO jobs (id,ligand_id,pair_index,seating,endpoint,conformer,
                    status,stage,attempt,input_xyz,smiles,charge,settings_json,ts_strategy,
                    gibbs_free_energy_eh) VALUES (?,'1',1,1,'TS',?,?,?,4,'seed.xyz','unused',0,'{}',?,-1)""",
                    (index, index, status, stage, strategy))
            before = [dict(r) for r in db.execute('SELECT * FROM jobs ORDER BY id')]
            db.close()
            dry = reset_failed_transition_states(root, dry_run=True)
            self.assertEqual(dry['ids'], [1, 2, 6])
            self.assertFalse((root/'ts_workflow_resets').exists())
            report = reset_failed_transition_states(root)
            self.assertEqual(report['ids'], [1, 2, 6])
            self.assertEqual(report['sampling_failures'], [3])
            self.assertEqual(report['active_ids'], [5])
            archive = Path(report['backup_file']).parent
            saved = json.loads((archive/'previous_jobs.json').read_text())
            self.assertEqual(saved, [before[i] for i in (0, 1, 5)])
            db = connect_database(root/'workflow.sqlite')
            after = [dict(r) for r in db.execute('SELECT * FROM jobs ORDER BY id')]
            db.close()
            self.assertEqual(after[2:5], before[2:5])
            self.assertEqual([(r['status'], r['stage'], r['attempt']) for r in after[:2]],
                             [('pending', 'opt', 4), ('pending', 'neb', 4)])
            self.assertEqual(after[5]['stage'], 'ts_guess')
            status = ts_recovery_status(root)
            self.assertEqual(status['pending'], 3)
            self.assertEqual(status['groups_with_TS_before'], 1)
            self.assertEqual(status['groups_with_TS_now'], 1)
            self.assertEqual(status['newly_covered_groups'], [])
            self.assertEqual(reset_failed_transition_states(root)['reset'], 0)


if __name__ == '__main__':
    unittest.main()
