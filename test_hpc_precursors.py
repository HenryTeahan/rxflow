"""Reference graph, fresh sampling, queue isolation, and reaction-diagram linking tests."""

from copy import deepcopy
from contextlib import closing
import csv
from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from rdkit import Chem
import rxembed

from hpc_precursors import (acid_base_plan, balanced_diagram_references, connect_precursors,
                            bound_reactant_plan, prepare_bound_reactants,
                            handle_precursors,
                            precursor_graph, precursor_links, precursor_plan, reference_root,
                            precursor_status, prepare_precursors, sample_precursor, submit_precursors)
from hpc_workflow import claim_job, main, run_job
from rxflow_endpoints import mapped_molecule


class PrecursorTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.root = self.base / 'results'
        self.root.mkdir()
        self.config = dict(paths=dict(ligands_csv='ligands.csv', source='absent_seeds', results='results'),
                           sampling=dict(seed=12345),
                           calculation=dict(cpus=1, memory_mb=2000, opt_attempts=2, ts_attempts=2,
                                            neb_fallback=True, solvent='h2o', timeout_seconds=86400),
                           executables=dict(orca_exe='orca', external_exe='oet_gxtb_alpb', xtb_exe='xtb'))
        with (self.base / 'ligands.csv').open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ID', 'ligand_smiles', 'central_metal', 'metal_geom', 'denticity', 'substrate'])
            writer.writerow(['11', '[O-]C1=NC(CC2=NC(=CC=C2)C2=NC=CC=C2)=CC=C1',
                             'Pd+2', 'spl', 3, 'COC(=O)C1=CC=NC=C1'])
        with closing(sqlite3.connect(self.root / 'workflow.sqlite')) as c, c:
            c.execute('''CREATE TABLE jobs (ligand_id TEXT,binding_permutation INTEGER,binding_mode INTEGER,
                coordination_isomer INTEGER,seating INTEGER,donor_maps TEXT,multiplicity INTEGER,
                pair_index INTEGER,substrate_site INTEGER,resonance INTEGER,substrate TEXT)''')
            for pair, site, resonance in [(1, 6, 0), (2, 6, 1), (3, 7, 0), (4, 7, 1)]:
                c.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                          ('11', 1, 4, 1, 4, '[3,7,13]', 1, pair, site, resonance,
                           'COC(=O)C1=CC=NC=C1'))
        self.plan = precursor_plan(self.config, self.base, self.root)
        native, isomer = precursor_graph(self.plan[0])
        substrate = Chem.AddHs(Chem.MolFromSmiles('COC(=O)C1=CC=NC=C1'))
        with closing(sqlite3.connect(self.root / 'workflow.sqlite')) as db, db:
            db.execute("ALTER TABLE jobs ADD COLUMN endpoint TEXT DEFAULT 'reactant'")
            db.execute('ALTER TABLE jobs ADD COLUMN smiles TEXT')
            for pair, site in ((1, 6), (2, 6), (3, 7), (4, 7)):
                donor = native.GetNumAtoms() + site - 1
                editable = Chem.RWMol(Chem.CombineMols(native, substrate))
                editable.AddBond(donor, native.GetNumAtoms() - 1, Chem.BondType.DATIVE)
                mol = editable.GetMol()
                Chem.SanitizeMol(mol)
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIdx() + 1)
                slots = {slot: atom for slot, atom in enumerate(isomer.centres[0].vertices) if atom is not None}
                slots[3] = donor
                smiles = rxembed.cxsmiles(rxembed.Isomer(mol, 'spl', slots))
                db.execute('UPDATE jobs SET smiles=? WHERE pair_index=?', (smiles, pair))

    def fresh_sample(self, spec, *, n, seed):
        molecule, source = precursor_graph(spec)
        self.assertEqual(molecule.GetNumConformers(), 0)
        geometries = []
        for i in range(n):
            geometry = Chem.Mol(molecule)
            conformer = Chem.Conformer(geometry.GetNumAtoms())
            for atom in geometry.GetAtoms():
                conformer.SetAtomPosition(atom.GetIdx(), (atom.GetIdx() * 0.1, i * 0.1, 0))
            geometry.AddConformer(conformer)
            geometries.append(geometry)
        return molecule, geometries

    def prepare(self):
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample) as embed:
            summary = prepare_precursors(self.config, self.base, self.root)
        return summary, embed

    def test_native_protonation_charges_and_binding_arrangement(self):
        self.assertEqual(len(self.plan), 3)
        deprot, prot, substrate = [precursor_graph(spec) for spec in self.plan]
        self.assertEqual([Chem.GetFormalCharge(item[0]) for item in (deprot, prot, substrate)], [1, 2, 0])
        self.assertEqual(deprot[0].GetNumConformers(), 0)
        self.assertEqual(prot[0].GetNumAtoms(), deprot[0].GetNumAtoms() + 1)
        self.assertEqual(prot[1].centres[0].vertices, deprot[1].centres[0].vertices)
        self.assertEqual(prot[0].GetAtomWithIdx(0).GetFormalCharge(), 0)
        self.assertEqual(prot[0].GetAtomWithIdx(0).GetSymbol(), 'O')
        self.assertTrue(any(a.GetAtomicNum() == 1 for a in prot[0].GetAtomWithIdx(0).GetNeighbors()))
        self.assertEqual(substrate[0].GetNumConformers(), 0)
        self.assertFalse(any(a.GetAtomicNum() in (28, 46) for a in substrate[0].GetAtoms()))

    def test_fresh_embed_receives_n10_without_coordinates(self):
        spec = self.plan[1]
        expected, geometries = self.fresh_sample(spec, n=10, seed=12345)
        class Ensemble:
            ids = list(range(10))
            def restored_mol(self, cid):
                return geometries[cid]
        with patch('hpc_precursors.rxembed.embed', return_value=Ensemble()) as embed:
            graph, results = sample_precursor(spec)
        self.assertEqual(len(results), 10)
        self.assertEqual(embed.call_args.kwargs, dict(n=10, seed=12345, threads=1))
        self.assertEqual(graph.GetNumConformers(), 0)
        self.assertEqual(embed.call_args.args[0].mol.GetNumConformers(), 0)
        self.assertEqual(Chem.GetFormalCharge(results[0]), 2)

    def test_embedded_charge_mismatch_is_rejected(self):
        spec = self.plan[1]
        expected, wrong = self.fresh_sample(self.plan[0], n=1, seed=12345)
        class Ensemble:
            ids = [0]
            def restored_mol(self, cid):
                return wrong[cid]
        with patch('hpc_precursors.rxembed.embed', return_value=Ensemble()):
            with self.assertRaisesRegex(ValueError, 'wrong total charge'):
                sample_precursor(spec)

    def test_preparation_is_isolated_idempotent_and_linked_to_every_pair(self):
        original = hashlib.sha256((self.root / 'workflow.sqlite').read_bytes()).digest()
        summary, embed = self.prepare()
        self.assertEqual(summary['species'], 3)
        self.assertEqual(summary['added'], 30)
        self.assertEqual(embed.call_count, 3)
        self.assertEqual(len(list((self.root / 'ligand_11' / 'complex').glob('*_input.xyz'))), 20)
        self.assertEqual(len(list((self.root / 'substrate').glob('*_input.xyz'))), 10)
        self.assertEqual(hashlib.sha256((self.root / 'workflow.sqlite').read_bytes()).digest(), original)
        with patch('hpc_precursors.sample_precursor', side_effect=AssertionError('Must not reembed')):
            again = prepare_precursors(self.config, self.base, self.root)
        self.assertEqual(again['added'], 0)
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            self.assertEqual(c.execute('SELECT count(*) FROM jobs').fetchone()[0], 30)
            links = c.execute('SELECT * FROM precursor_links ORDER BY pair_index').fetchall()
            self.assertEqual(len(links), 4)
            self.assertEqual({r['resonance'] for r in links}, {0, 1})
            self.assertEqual({r['protonation_site'] for r in links}, {1})
            energies = c.execute('SELECT * FROM precursor_diagram_references').fetchall()
            self.assertTrue(all(r['deprotonated_complex_plus_substrate_eh'] is None for r in energies))

    def test_partial_embedding_records_missing_seeds_and_can_repair_them(self):
        def partial(spec, *, n, seed):
            mol, geometries = self.fresh_sample(spec, n=n, seed=seed)
            return mol, geometries[:9]
        with patch('hpc_precursors.sample_precursor', side_effect=partial):
            result = prepare_precursors(self.config, self.base, self.root)
        self.assertEqual(len(result['errors']), 3)
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            self.assertEqual(c.execute("SELECT count(*) FROM jobs WHERE stage='sampling'").fetchone()[0], 3)
        repaired, embed = self.prepare()
        self.assertEqual(repaired['repaired'], 3)
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            self.assertEqual(c.execute("SELECT count(*) FROM jobs WHERE status='pending'").fetchone()[0], 30)

    def test_lowest_gibbs_ignores_retained_filter_and_does_not_merge_resonance_routes(self):
        self.prepare()
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            for role, energy in [('deprotonated', -100.), ('protonated', -101.), ('substrate', -20.)]:
                c.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=?,retained=1 WHERE endpoint=? AND conformer=1", (energy, role))
            c.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=-102,retained=0 WHERE endpoint='protonated' AND conformer=2")
            c.execute("UPDATE jobs SET status='failed',gibbs_free_energy_eh=-1000 WHERE endpoint='protonated' AND conformer=3")
            result = precursor_status(c)
            self.assertEqual(len(result['diagram_references']), 4)
            for r in result['diagram_references']:
                self.assertEqual(r['protonated_complex_gibbs_eh'], -102.)
                self.assertEqual(r['deprotonated_complex_plus_substrate_eh'], -120.)
                self.assertEqual(r['protonated_complex_plus_substrate_eh'], -122.)

    def test_changed_protocol_is_refused(self):
        self.prepare()
        changed = deepcopy(self.config)
        changed['calculation']['solvent'] = 'thf'
        with self.assertRaisesRegex(ValueError, 'chemistry/settings changed'):
            prepare_precursors(changed, self.base, self.root)
        with self.assertRaisesRegex(ValueError, 'chemistry/settings changed'):
            prepare_precursors(self.config, self.base, self.root, conformers=5)

    def test_worker_reuses_minimum_optimization_and_frequency_validation(self):
        self.prepare()
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            job = claim_job(c)
            xyz = (self.root / job['input_xyz']).read_text()
            results = [dict(out_file='opt.out', xyz=xyz, converged=True),
                       dict(out_file='freq.out', electronic_energy_eh=-100., gibbs_free_energy_eh=-99.,
                            frequencies=[0.]*6+[10.], imaginary_count=0)]
            with patch('hpc_workflow.run_orca', side_effect=results) as orca, \
                    patch('hpc_workflow.connectivity_changes', return_value={}):
                run_job(c, job, self.root, {})
            self.assertEqual([call.args[3] for call in orca.call_args_list], ['opt', 'freq'])
            self.assertEqual(job['charge'], 1)
            self.assertEqual(c.execute('SELECT status FROM jobs WHERE id=?', (job['id'],)).fetchone()[0], 'succeeded')
            self.assertEqual(len(precursor_status(c)['lowest_gibbs']), 1)

    def test_submission_has_prepare_dependency_and_workers_use_same_config(self):
        config_path = self.base / 'hpc.toml'
        with patch('hpc_precursors.subprocess.check_output', side_effect=['100\n', '101\n']) as submit:
            result = submit_precursors(config_path, self.config, workers=7)
        self.assertEqual((result['prepare'], result['workers']), ('100', '101'))
        self.assertEqual(result['species'], 3)
        first, second = [call.args[0] for call in submit.call_args_list]
        self.assertIn('--cpus-per-task=1', first)
        self.assertIn('--dependency=afterok:100', second)
        self.assertIn('--array=0-6%7', second)
        self.assertIn('precursors prepare --conformers 10', first[-1])
        self.assertIn(str(config_path), second[-1])
        self.assertIn('precursors worker', second[-1])

    def test_cli_dispatches_precursors_to_separate_workflow(self):
        config_path = self.base / 'hpc.toml'
        config_path.write_text('')
        with patch('hpc_workflow.tomllib.loads', return_value=self.config), \
                patch('hpc_precursors.handle_precursors') as handle:
            main(['--config', str(config_path), 'precursors', 'status'])
        args, config, base, root, parser = handle.call_args.args
        self.assertEqual((args.command, args.precursor_command), ('precursors', 'status'))
        self.assertEqual(args.conformers, 10)
        self.assertEqual(root, self.root.resolve())

    def test_substrate_change_cannot_link_to_unrelated_reaction(self):
        changed = deepcopy(self.plan)
        changed[-1]['substrate'] = 'C'
        with self.assertRaisesRegex(ValueError, 'absent from the configured CSV'):
            precursor_links(changed, self.root / 'workflow.sqlite')

    def test_missing_database_or_invalid_complex_charge_is_refused(self):
        bad = deepcopy(self.plan[0])
        bad['charge'] = 2
        with self.assertRaisesRegex(ValueError, 'charge mismatch'):
            precursor_graph(bad)
        charged = dict(self.plan[-1], substrate='[O-]c1ccccc1')
        with self.assertRaisesRegex(ValueError, 'neutral'):
            from hpc_precursors import canonical_substrate
            canonical_substrate(charged['substrate'])
        with self.assertRaisesRegex(ValueError, 'reaction workflow first'):
            precursor_plan(self.config, self.base, self.base / 'missing')

    def test_acid_base_states_balance_atoms_and_charge(self):
        acid, acetate = [precursor_graph(spec)[0] for spec in acid_base_plan()]
        deprot, prot, substrate = [precursor_graph(spec)[0] for spec in self.plan]
        self.assertEqual([Chem.GetFormalCharge(m) for m in (acid, acetate)], [0, -1])
        def atoms(*molecules):
            return Counter(a.GetAtomicNum() for m in molecules for a in m.GetAtoms())
        self.assertEqual(atoms(prot, substrate, acetate), atoms(deprot, substrate, acid))
        self.assertEqual(Chem.GetFormalCharge(prot) + Chem.GetFormalCharge(acetate), 1)
        self.assertEqual(Chem.GetFormalCharge(deprot) + Chem.GetFormalCharge(acid), 1)

    def test_acid_base_real_sampling_keeps_requested_conformers(self):
        for spec in acid_base_plan():
            with self.subTest(spec=spec['reference_name']):
                molecule, conformers = sample_precursor(spec, n=10, seed=12345)
                self.assertEqual(len(conformers), 10)
                self.assertEqual(molecule.GetNumConformers(), 0)
                self.assertTrue(all(Chem.GetFormalCharge(m) == spec['charge'] for m in conformers))

    def test_shared_reference_preparation_and_balanced_diagram_totals(self):
        self.prepare()
        root = reference_root(self.config, self.base)
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample):
            summary = prepare_precursors(self.config, self.base, root, reference=True)
        self.assertEqual(summary['added'], 20)
        self.assertEqual(len(list((root / 'acetic_acid').glob('*_input.xyz'))), 10)
        self.assertEqual(len(list((root / 'acetate').glob('*_input.xyz'))), 10)
        with closing(connect_precursors(root / 'references.sqlite')) as ref:
            for name, energy in [('acetic_acid', -30.), ('acetate', -29.)]:
                ref.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=? WHERE ligand_id=? AND conformer=1", (energy, name))
        with closing(sqlite3.connect(self.root / 'workflow.sqlite')) as reaction, reaction:
            reaction.execute('ALTER TABLE jobs ADD COLUMN status TEXT')
            reaction.execute('ALTER TABLE jobs ADD COLUMN gibbs_free_energy_eh REAL')
            reaction.execute("UPDATE jobs SET endpoint='reactant',status='succeeded',gibbs_free_energy_eh=-121")
        with closing(connect_precursors(self.root / 'precursors.sqlite')) as c:
            for name, energy in [('deprotonated', -100.), ('protonated', -101.), ('substrate', -20.)]:
                c.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=? WHERE endpoint=? AND conformer=1", (energy, name))
            rows = balanced_diagram_references(c, root / 'references.sqlite', self.root / 'workflow.sqlite')
            self.assertEqual(len(rows), 4)
            for row in rows:
                self.assertEqual(row['protonated_start_gibbs_eh'], -150.)
                self.assertEqual(row['deprotonated_start_gibbs_eh'], -150.)
                self.assertEqual(row['reactant_with_acid_gibbs_eh'], -151.)
                self.assertIsNone(row['product_with_acid_gibbs_eh'])
                self.assertIsNone(row['TS_with_acid_gibbs_eh'])
            with closing(connect_precursors(root / 'references.sqlite')) as ref:
                protocol = json.loads(ref.execute("SELECT value FROM precursor_metadata WHERE key='protocol'").fetchone()[0])
                protocol['settings']['solvent'] = 'thf'
                ref.execute("UPDATE precursor_metadata SET value=? WHERE key='protocol'", (json.dumps(protocol),))
            with self.assertRaisesRegex(ValueError, 'settings differ'):
                balanced_diagram_references(c, root / 'references.sqlite', self.root / 'workflow.sqlite')

    def test_reference_slurm_commands_use_global_directory_and_acid_base_action(self):
        with patch('hpc_precursors.subprocess.check_output', side_effect=['100\n', '101\n']) as submit:
            result = submit_precursors(self.base / 'hpc.toml', self.config, workers=2, reference=True)
        self.assertEqual(result['species'], 4)
        self.assertIn('acid-base prepare --conformers 10', submit.call_args_list[0].args[0][-1])
        self.assertIn('acid-base worker', submit.call_args_list[1].args[0][-1])

    def test_bound_states_share_resonance_routes_and_keep_native_protonation_and_substrate(self):
        plan = bound_reactant_plan(self.config, self.base, self.root)
        self.assertEqual(len(plan), 2)
        self.assertEqual([[p['pair_index'] for p in spec['reaction_pairs']] for spec in plan], [[1, 2], [3, 4]])
        acid, acetate = [precursor_graph(spec)[0] for spec in acid_base_plan()]
        inventory = lambda *mols: Counter(atom.GetAtomicNum() for mol in mols for atom in mol.GetAtoms())
        for spec in plan:
            mol, source = precursor_graph(spec)
            reactant = mapped_molecule(spec['reactant_smiles'])
            self.assertEqual(Chem.GetFormalCharge(mol), 2)
            self.assertEqual(mol.GetNumConformers(), 0)
            self.assertEqual(source.mol.GetNumConformers(), 0)
            self.assertEqual(mol.GetNumAtoms(), reactant.GetNumAtoms() + 1)
            self.assertEqual(mol.GetAtomWithIdx(mol.GetNumAtoms() - 1).GetNeighbors()[0].GetIdx(), spec['protonation_site'] - 1)
            dative = lambda m: [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()
                                if b.GetBondType() == Chem.BondType.DATIVE]
            self.assertCountEqual(dative(mol), dative(reactant))
            self.assertEqual(inventory(mol, acetate), inventory(reactant, acid))

    def test_bound_protonation_uses_native_csv_site_despite_shifted_reactant_resonance(self):
        spec = bound_reactant_plan(self.config, self.base, self.root)[0]
        original = mapped_molecule(spec['reactant_smiles'])
        shifted = next((form for form in Chem.ResonanceMolSupplier(original)
                        if form.GetAtomWithIdx(spec['protonation_site'] - 1).GetFormalCharge() >= 0), None)
        self.assertIsNotNone(shifted)
        spec = dict(spec, reactant_smiles=Chem.MolToCXSmiles(shifted))
        mol, source = precursor_graph(spec)
        proton = mol.GetAtomWithIdx(mol.GetNumAtoms() - 1)
        self.assertEqual(proton.GetNeighbors()[0].GetIdx(), spec['protonation_site'] - 1)
        self.assertEqual(Chem.GetFormalCharge(mol), 2)

    def test_bound_preparation_is_additive_and_maps_two_resonances_to_one_state(self):
        self.prepare()
        before = hashlib.sha256((self.root / 'precursors.sqlite').read_bytes()).digest()
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample) as sample:
            summary = prepare_bound_reactants(self.config, self.base, self.root)
        self.assertEqual(summary['species'], 2)
        self.assertEqual(summary['added'], 20)
        self.assertEqual(sample.call_count, 2)
        self.assertTrue(all(call.kwargs['n'] == 10 for call in sample.call_args_list))
        self.assertEqual(hashlib.sha256((self.root / 'precursors.sqlite').read_bytes()).digest(), before)
        with closing(connect_precursors(self.root / 'acid_base.sqlite', bound=True)) as db:
            self.assertEqual([tuple(row) for row in db.execute('SELECT pair_index,bound_pair_index FROM bound_reactant_links ORDER BY pair_index')],
                             [(1, 1), (2, 1), (3, 3), (4, 3)])
            self.assertEqual(db.execute('SELECT DISTINCT charge FROM jobs').fetchone()[0], 2)
            db.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=-100 WHERE id=1")
        with patch('hpc_precursors.sample_precursor', side_effect=AssertionError('Must not reembed')):
            self.assertEqual(prepare_bound_reactants(self.config, self.base, self.root)['added'], 0)
        with closing(connect_precursors(self.root / 'acid_base.sqlite', bound=True)) as db:
            self.assertEqual(tuple(db.execute('SELECT status,gibbs_free_energy_eh FROM jobs WHERE id=1').fetchone()), ('succeeded', -100.))

    def test_acid_base_prepare_reuses_completed_global_references_and_submit_can_repeat(self):
        refs = reference_root(self.config, self.base)
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample):
            prepare_precursors(self.config, self.base, refs, reference=True)
        with closing(connect_precursors(refs / 'references.sqlite')) as db:
            db.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=-30")
        args = SimpleNamespace(command='acid-base', precursor_command='prepare', conformers=10)
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample) as sample, patch('builtins.print'):
            handle_precursors(args, self.config, self.base, self.root, argparse.ArgumentParser())
        self.assertEqual(sample.call_count, 2)
        with closing(connect_precursors(refs / 'references.sqlite')) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM jobs WHERE status='succeeded'").fetchone()[0], 20)
        with patch('hpc_precursors.subprocess.check_output', side_effect=['100\n', '101\n']):
            result = submit_precursors(self.base / 'hpc.toml', self.config, reference=True)
        self.assertEqual(result['prepare'], '100')

    def test_protonated_bound_uses_existing_minimum_optimization_and_frequency(self):
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample):
            prepare_bound_reactants(self.config, self.base, self.root)
        with closing(connect_precursors(self.root / 'acid_base.sqlite', bound=True)) as db:
            job = claim_job(db)
            xyz = (self.root / job['input_xyz']).read_text()
            results = [dict(out_file='opt.out', xyz=xyz, converged=True),
                       dict(out_file='freq.out', electronic_energy_eh=-100., gibbs_free_energy_eh=-99.,
                            frequencies=[0.]*6+[10.], imaginary_count=0)]
            with patch('hpc_workflow.run_orca', side_effect=results) as orca, \
                    patch('hpc_workflow.connectivity_changes', return_value={}):
                run_job(db, job, self.root, {})
            self.assertEqual(job['endpoint'], 'protonated_bound')
            self.assertEqual(job['charge'], 2)
            self.assertEqual([call.args[3] for call in orca.call_args_list], ['opt', 'freq'])
            self.assertEqual(db.execute('SELECT status FROM jobs WHERE id=?', (job['id'],)).fetchone()[0], 'succeeded')

    def test_acid_base_worker_drains_shared_references_and_its_own_bound_queue(self):
        refs = reference_root(self.config, self.base)
        with patch('hpc_precursors.sample_precursor', side_effect=self.fresh_sample):
            prepare_precursors(self.config, self.base, refs, reference=True)
            prepare_bound_reactants(self.config, self.base, self.root)
        calls = []
        def finish(db, job, root, executables):
            calls.append((job['endpoint'], root))
            db.execute("UPDATE jobs SET status='succeeded',stage='complete' WHERE id=?", (job['id'],))
        args = SimpleNamespace(command='acid-base', precursor_command='worker')
        with patch('hpc_precursors.shutil.which', return_value='/usr/bin/true'), \
                patch('hpc_precursors.run_job', side_effect=finish), patch('builtins.print'):
            handle_precursors(args, self.config, self.base, self.root, argparse.ArgumentParser())
        self.assertEqual(Counter(calls), Counter({('substrate', refs): 20, ('protonated_bound', self.root): 20}))
        self.assertEqual(reference_root(self.config, self.base), (self.base / 'solution_references').resolve())


if __name__ == '__main__':
    unittest.main()
