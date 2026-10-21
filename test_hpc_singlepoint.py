"""Single-point provenance, queueing, energy replacement and ORCA validation."""

from contextlib import closing
from concurrent.futures import ThreadPoolExecutor
import argparse
import io
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hpc_precursors import connect_precursors, load_diagram_inputs
from hpc_singlepoint import (METHOD, calculation_settings, claim_singlepoint, connect_singlepoints,
                            handle_singlepoints, mpi_executable, plan_singlepoints, prepare_singlepoints, read_jobs, run_singlepoint,
                            selected_jobs, submit_singlepoints)
from hpc_workflow import connect_database
from orca_hpc import collect_singlepoint, orca_input
from rxflow import HARTREE_TO_KCAL

XYZ = "3\nsource geometry\nO 0 0 0\nH 0.9572 0 0\nH -0.239 0.927 0\n"
SMILES = "[O:1]([H:2])[H:3]"


class SinglepointTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name).resolve()
        self.root = self.base / 'results'
        self.refs = self.base / 'solution_references'
        self.root.mkdir()
        self.refs.mkdir()
        self.config = dict(paths=dict(results='results'), calculation={}, executables=dict(orca_exe='orca'))
        self.mpi_bin = self.base / 'mpi' / 'bin'
        self.mpi_bin.mkdir(parents=True)
        (self.mpi_bin.parent / 'lib').mkdir()
        self.mpirun = self.mpi_bin / 'mpirun'
        self.mpirun.write_text('#!/bin/sh\nprintf "SCF CONVERGED AFTER 7 CYCLES\\nFINAL SINGLE POINT ENERGY -100.5\\nORCA TERMINATED NORMALLY\\n"\n')
        self.mpirun.chmod(0o755)
        self.config['executables']['mpirun_exe'] = str(self.mpirun)
        self.reaction = self.root / 'workflow.sqlite'
        self.precursor = self.root / 'precursors.sqlite'
        self.reference = self.refs / 'references.sqlite'
        for path, connector in [(self.reaction, connect_database), (self.precursor, connect_precursors),
                                (self.reference, connect_precursors)]:
            with closing(connector(path)):
                pass
        self.add(self.reaction, 1, 'reactant', -10.1, -10.2, retained=0)
        self.add(self.reaction, 2, 'reactant', -10.05, -10.2, retained=1, conformer=2)
        self.add(self.reaction, 3, 'TS', -10.0, -10.15)
        self.add(self.reaction, 4, 'product', -10.12, -10.25)
        self.add(self.reaction, 5, 'reactant', -99., -99., status='failed', conformer=3)
        self.add(self.precursor, 1, 'deprotonated', -9., -9.1, seating=0)
        self.add(self.precursor, 2, 'protonated', -9.4, -9.55, seating=1)
        self.add(self.precursor, 3, 'substrate', -1., -1.1, ligand='substrate_001', pair=0, seating=0)
        self.add(self.reference, 1, 'substrate', -5., -5.1, ligand='acetic_acid', pair=0, seating=0)
        self.add(self.reference, 2, 'substrate', -4.5, -4.7, ligand='acetate', pair=0, seating=0, charge=-1)
        with closing(connect_precursors(self.precursor)) as db:
            db.execute("INSERT INTO precursor_links VALUES ('1',1,4,1,6,0,1,'substrate_001')")

    def add(self, path, identifier, role, gibbs, electronic, *, retained=1, status='succeeded',
            conformer=1, ligand='1', pair=1, seating=4, charge=0):
        xyz_path = Path('fixtures') / f'{path.stem}_{identifier}.xyz'
        target = path.parent / xyz_path
        target.parent.mkdir(exist_ok=True)
        target.write_text(XYZ)
        connector = connect_database if path == self.reaction else connect_precursors
        with closing(connector(path)) as db:
            db.execute('''INSERT INTO jobs (id,ligand_id,pair_index,seating,endpoint,conformer,status,
                attempt,input_xyz,xyz_file,smiles,charge,multiplicity,settings_json,
                electronic_energy_eh,gibbs_free_energy_eh,retained,binding_permutation)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                (identifier, ligand, pair, seating, role, conformer, status, 1, str(xyz_path), str(xyz_path),
                 SMILES, charge, 1, '{}', electronic, gibbs, retained, pair))

    def prepare(self):
        return prepare_singlepoints(self.config, self.base, self.root)

    def finish(self, *, omit_reference=None):
        for path in (self.root / 'singlepoints.sqlite', self.refs / 'singlepoints.sqlite'):
            with closing(connect_singlepoints(path)) as db:
                for row in db.execute('SELECT * FROM jobs').fetchall():
                    if row['ligand_id'] == omit_reference:
                        continue
                    energy = -1000. - row['source_job_id']
                    gibbs = row['gxtb_gibbs_free_energy_eh'] - row['gxtb_electronic_energy_eh'] + energy
                    db.execute("UPDATE jobs SET status='succeeded',electronic_energy_eh=?,corrected_gibbs_free_energy_eh=? WHERE id=?",
                               (energy, gibbs, row['id']))

    def load(self, mode='auto'):
        return load_diagram_inputs(self.reaction, self.precursor, self.reference, energy_mode=mode)

    def test_selected_geometry_ignores_retention_and_failed_rows_and_preparation_is_idempotent(self):
        before = {path: hashlib.sha256(path.read_bytes()).digest() for path in (self.reaction, self.precursor, self.reference)}
        plan = plan_singlepoints(self.config, self.base, self.root)
        self.assertEqual(plan['selected_states'], 8)
        result = self.prepare()
        self.assertEqual(result['added'], 8)
        self.assertFalse(result['errors'])
        self.assertEqual(self.prepare()['existing'], 8)
        rows = read_jobs(self.root / 'singlepoints.sqlite')
        selected = [row['source_job_id'] for row in rows if row['source_kind'] == 'reaction']
        self.assertEqual(selected, [1, 3, 4])
        for path in before:
            self.assertEqual(hashlib.sha256(path.read_bytes()).digest(), before[path])
        for root in (self.root, self.refs):
            for inp in root.glob('**/singlepoint/*_input.inp'):
                text = inp.read_text()
                self.assertIn('! PBE0 def2-TZVP D4 SP', text)
                self.assertNotIn('TightSCF', text)
                self.assertNotIn('MaxIter', text)
                self.assertIn('SMDsolvent "water"', text)
                self.assertNotIn('ExtOpt', text)
                self.assertNotIn('NumFreq', text)
                self.assertNotIn('%geom', text)

    def test_ts_window_uses_numfreq_gibbs_and_includes_boundary_per_group(self):
        for identifier, relative, status in ((6, 1., 'succeeded'), (7, 3., 'succeeded'),
                                             (8, 3.001, 'succeeded'), (9, -100., 'failed')):
            self.add(self.reaction, identifier, 'TS', -10. + relative / HARTREE_TO_KCAL,
                     -100. - identifier, conformer=identifier, retained=0, status=status)
        # A much lower TS in another pathway must not eliminate this group's seeds.
        self.add(self.reaction, 10, 'TS', -20., -21., pair=2)
        selected = selected_jobs(read_jobs(self.reaction))
        ts = [row for row in selected if row['endpoint'] == 'TS']
        self.assertEqual([row['id'] for row in ts], [3, 6, 7, 10])
        plan = self.prepare()
        self.assertEqual(plan['ts_energy_window_kcal'], 3.)
        self.assertEqual(plan['added'], 11)
        queued = [row['source_job_id'] for row in read_jobs(self.root / 'singlepoints.sqlite')
                  if row['endpoint'] == 'TS']
        self.assertEqual(queued, [3, 6, 7, 10])
        self.assertEqual(self.prepare()['added'], 0)

    def test_dft_waits_for_entire_ts_window_and_can_change_selected_ts(self):
        self.add(self.reaction, 6, 'TS', -10. + 2. / HARTREE_TO_KCAL, -10.15,
                 conformer=-2, retained=0)
        self.prepare()
        self.finish()
        with closing(connect_singlepoints(self.root / 'singlepoints.sqlite')) as db:
            db.execute("UPDATE jobs SET status='pending' WHERE source_kind='reaction' AND source_job_id=6")
        partial = self.load()
        self.assertEqual(partial['singlepoint_summary']['ready'], 8)
        self.assertEqual(partial['singlepoint_summary']['selected_states'], 9)
        self.assertEqual(partial['singlepoint_summary']['active_mode'], 'gxtb')
        self.finish()
        result = self.load()
        self.assertEqual(result['singlepoint_summary']['active_mode'], 'dft')
        ts = [row for row in result['jobs'] if row['endpoint'] == 'TS']
        self.assertEqual({row['id'] for row in ts if row['singlepoint_selected']}, {3, 6})
        best = min(ts, key=lambda row: row['gibbs_free_energy_eh'])
        self.assertEqual(best['id'], 6)
        self.assertAlmostEqual(result['diagram_references'][0]['TS_with_acid_gibbs_eh'],
                               best['gibbs_free_energy_eh'] + result['reference_minima']['acetic_acid']['gibbs_free_energy_eh'])

    def test_new_lower_ts_skips_pending_outside_window_and_can_reactivate_it(self):
        self.prepare()
        self.add(self.reaction, 6, 'TS', -10. - 4. / HARTREE_TO_KCAL, -10.15, conformer=-2)
        prepared = self.prepare()
        self.assertEqual(prepared['skipped'], 1)
        rows = [row for row in read_jobs(self.root / 'singlepoints.sqlite') if row['endpoint'] == 'TS']
        self.assertEqual({row['source_job_id']: row['status'] for row in rows}, {3: 'skipped', 6: 'pending'})
        with closing(connect_database(self.reaction)) as db:
            db.execute("UPDATE jobs SET status='failed' WHERE id=6")
        prepared = self.prepare()
        self.assertEqual(prepared['reactivated'], 1)
        rows = [row for row in read_jobs(self.root / 'singlepoints.sqlite') if row['endpoint'] == 'TS']
        self.assertEqual({row['source_job_id']: row['status'] for row in rows}, {3: 'pending', 6: 'skipped'})

    def test_completed_dft_is_preserved_when_ts_moves_outside_window(self):
        self.prepare()
        self.finish()
        before = next(row for row in read_jobs(self.root / 'singlepoints.sqlite')
                      if row['endpoint'] == 'TS')
        self.add(self.reaction, 6, 'TS', -10. - 4. / HARTREE_TO_KCAL, -10.15, conformer=-2)
        self.prepare()
        after = next(row for row in read_jobs(self.root / 'singlepoints.sqlite') if row['id'] == before['id'])
        self.assertEqual(after, before)

    def test_singlepoint_input_has_exact_requested_method_and_preserves_charge_and_spin(self):
        text = orca_input(XYZ, 'singlepoint', -1, 1, calculation_settings(self.config))
        self.assertIn('* xyz -1 1', text)
        self.assertIn('%pal nprocs 4 end', text)
        self.assertIn('%maxcore 1500', text)
        self.assertIn('smd true', text)
        with self.assertRaisesRegex(ValueError, 'positive integers'):
            calculation_settings(dict(singlepoint=dict(cpus=0)))

    def test_encounter_reuses_identical_electronic_singlepoint_with_its_own_thermal_correction(self):
        self.prepare()
        self.finish()
        source = self.root / 'encounters.sqlite'
        with closing(connect_precursors(source, encounter=True)):
            pass
        self.add(source, 1, 'encounter', -11., -11.3, seating=0)
        result = self.prepare()
        self.assertEqual(result['added'], 1)
        self.assertEqual(result['reused'], 1)
        row = next(r for r in read_jobs(self.root / 'singlepoints.sqlite') if r['source_kind'] == 'encounter')
        old = next(r for r in read_jobs(self.root / 'singlepoints.sqlite') if r['id'] == row['reused_from_id'])
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual(row['electronic_energy_eh'], old['electronic_energy_eh'])
        self.assertAlmostEqual(row['corrected_gibbs_free_energy_eh'], .3 + old['electronic_energy_eh'])
        data = self.load()
        self.assertEqual(data['singlepoint_summary']['active_mode'], 'dft')
        link = data['diagram_references'][0]
        self.assertEqual(link['encounter_job_id'], 1)
        self.assertAlmostEqual(link['encounter_gibbs_eh'], row['corrected_gibbs_free_energy_eh'])
        self.assertAlmostEqual(link['encounter_with_acid_gibbs_eh'],
                               row['corrected_gibbs_free_energy_eh'] + link['acetic_acid_gibbs_eh'])

    def test_singlepoint_reuse_requires_matching_geometry_charge_and_spin(self):
        self.prepare()
        self.finish()
        source = self.root / 'encounters.sqlite'
        with closing(connect_precursors(source, encounter=True)):
            pass
        self.add(source, 1, 'encounter', -11., -11.3, seating=0, charge=1)
        self.add(source, 2, 'encounter', -11., -11.3, seating=0, pair=2)
        self.add(source, 3, 'encounter', -11., -11.3, seating=0, pair=3)
        with closing(connect_precursors(source, encounter=True)) as db:
            db.execute('UPDATE jobs SET multiplicity=3 WHERE id=2')
        geometry = next(r for r in read_jobs(source) if r['id'] == 3)
        (self.root / geometry['xyz_file']).write_text(XYZ.replace('0.9572', '0.9580'))
        result = self.prepare()
        self.assertEqual(result['reused'], 0)
        rows = [r for r in read_jobs(self.root / 'singlepoints.sqlite') if r['source_kind'] == 'encounter']
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(r['status'] == 'pending' and r['reused_from_id'] is None for r in rows))

    def test_shared_references_are_not_duplicated_by_another_run(self):
        self.prepare()
        other = self.base / 'other'
        other.mkdir()
        with closing(connect_database(other / 'workflow.sqlite')):
            pass
        result = prepare_singlepoints(self.config, self.base, other)
        self.assertEqual(result['added'], 0)
        self.assertEqual(result['existing'], 2)
        self.assertEqual(len(read_jobs(self.refs / 'singlepoints.sqlite')), 2)

    def test_new_bound_states_queue_only_missing_singlepoints_and_share_resonance_links(self):
        self.prepare()
        self.finish()
        source = self.root / 'acid_base.sqlite'
        with closing(connect_precursors(source, bound=True)):
            pass
        for identifier, site, conformer, gibbs, electronic in ((1, 1, 1, -11., -11.1),
                                                              (2, 1, 2, -10., -10.1),
                                                              (3, 3, 1, -10.9, -11.)):
            self.add(source, identifier, 'protonated_bound', gibbs, electronic, conformer=conformer, charge=2)
            with closing(connect_precursors(source, bound=True)) as db:
                db.execute('UPDATE jobs SET protonation_site=?,substrate_site=6 WHERE id=?', (site, identifier))
        with closing(connect_precursors(self.precursor)) as db:
            db.execute("INSERT INTO precursor_links VALUES ('1',2,4,1,6,1,1,'substrate_001')")
        prepared = self.prepare()
        self.assertEqual(prepared['added'], 2)
        self.assertEqual(prepared['existing'], 8)
        queued = [row for row in read_jobs(self.root / 'singlepoints.sqlite') if row['source_kind'] == 'bound_reactant']
        self.assertEqual({row['source_job_id'] for row in queued}, {1, 3})
        self.assertEqual({row['charge'] for row in queued}, {2})
        self.assertEqual(self.load()['singlepoint_summary']['active_mode'], 'gxtb')
        partial = self.load('dft')
        self.assertTrue(all(row['protonated_bound_with_acetate_gibbs_eh'] is None for row in partial['diagram_references']))
        self.finish()
        corrected = self.load()
        self.assertEqual(corrected['singlepoint_summary']['ready'], 10)
        self.assertEqual(corrected['singlepoint_summary']['active_mode'], 'dft')
        links = corrected['diagram_references']
        self.assertEqual(len(links), 2)
        self.assertEqual(links[0]['protonated_bound_job_id'], links[1]['protonated_bound_job_id'])
        for link in links:
            self.assertAlmostEqual(link['protonated_bound_gibbs_eh'], -11. - (-11.1) - 1001.)
            self.assertAlmostEqual(link['protonated_bound_with_acetate_gibbs_eh'],
                                   link['protonated_bound_gibbs_eh'] + link['acetate_gibbs_eh'])
        self.assertEqual(self.prepare()['added'], 0)

    def test_concurrent_preparation_does_not_duplicate_or_truncate_shared_inputs(self):
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda _: self.prepare(), range(2)))
        self.assertEqual(sum(result['added'] for result in results), 8)
        self.assertEqual(len(read_jobs(self.refs / 'singlepoints.sqlite')), 2)
        for root in (self.root, self.refs):
            for job in read_jobs(root / 'singlepoints.sqlite'):
                self.assertEqual((root / job['input_xyz']).read_text(), XYZ)

    def test_auto_switches_only_after_every_selected_state_and_reference_is_complete(self):
        self.prepare()
        self.finish(omit_reference='acetate')
        partial = self.load()
        self.assertEqual(partial['singlepoint_summary']['ready'], 7)
        self.assertEqual(partial['singlepoint_summary']['active_mode'], 'gxtb')
        self.assertEqual(partial['jobs'][0]['gibbs_free_energy_eh'], -10.1)
        self.finish()
        corrected = self.load()
        self.assertEqual(corrected['singlepoint_summary']['active_mode'], 'dft')
        selected = next(row for row in corrected['jobs'] if row['id'] == 1)
        self.assertAlmostEqual(selected['gibbs_free_energy_eh'], -10.1 - (-10.2) - 1001.)
        self.assertEqual(selected['gxtb_gibbs_free_energy_eh'], -10.1)
        self.assertEqual(selected['gxtb_electronic_energy_eh'], -10.2)
        self.assertIsNone(next(row for row in corrected['jobs'] if row['id'] == 2)['gibbs_free_energy_eh'])
        link = corrected['diagram_references'][0]
        acid = corrected['reference_minima']['acetic_acid']['gibbs_free_energy_eh']
        self.assertAlmostEqual(link['reactant_with_acid_gibbs_eh'], selected['gibbs_free_energy_eh'] + acid)
        self.assertAlmostEqual(link['protonated_start_gibbs_eh'],
                               link['protonated_complex_gibbs_eh'] + link['substrate_gibbs_eh'] + link['acetate_gibbs_eh'])
        self.assertAlmostEqual(link['deprotonated_start_gibbs_eh'],
                               link['deprotonated_complex_gibbs_eh'] + link['substrate_gibbs_eh'] + acid)
        raw = self.load('gxtb')
        self.assertEqual(raw['jobs'][0]['gibbs_free_energy_eh'], -10.1)

    def test_partial_dft_mode_leaves_missing_reference_blank_instead_of_mixing_models(self):
        self.prepare()
        self.finish(omit_reference='acetate')
        result = self.load('dft')
        self.assertEqual(result['singlepoint_summary']['active_mode'], 'dft')
        self.assertNotIn('acetate', result['reference_minima'])
        self.assertIsNone(result['diagram_references'][0]['protonated_start_gibbs_eh'])
        self.assertIsNotNone(result['diagram_references'][0]['deprotonated_start_gibbs_eh'])
        self.assertIsNone(next(row for row in result['reference_jobs'] if row['ligand_id'] == 'acetate')['gibbs_free_energy_eh'])

    def test_changed_source_geometry_invalidates_completed_correction(self):
        self.prepare()
        self.finish()
        row = next(row for row in read_jobs(self.reaction) if row['id'] == 1)
        path = self.root / row['xyz_file']
        path.write_text(XYZ.replace('0.9572', '0.9580'))
        result = self.load()
        self.assertEqual(result['singlepoint_summary']['ready'], 7)
        self.assertEqual(result['singlepoint_summary']['active_mode'], 'gxtb')
        self.assertEqual(self.prepare()['added'], 1)
        self.assertEqual(len(read_jobs(self.root / 'singlepoints.sqlite')), 7)

    def test_queue_claim_and_real_subprocess_store_matched_thermal_correction(self):
        self.prepare()
        executable = self.base / 'fake-orca'
        executable.write_text(
            '#!/bin/sh\n'
            '[ "$2" = "--use-hwthread-cpus" ] || exit 20\n'
            'case "$PWD" in "$SCRATCH"/rxflow-orca-*) ;; *) exit 21 ;; esac\n'
            'case "$LD_LIBRARY_PATH" in *"/mpi/lib:"*) ;; *) exit 22 ;; esac\n'
            'mpirun\n')
        executable.chmod(0o755)
        scratch = self.base / 'node-scratch'
        with closing(connect_singlepoints(self.root / 'singlepoints.sqlite')) as db:
            job = claim_singlepoint(db)
            self.assertEqual(job['attempt'], 1)
            with patch.dict('os.environ', PATH='/usr/bin:/bin', LD_LIBRARY_PATH='', SCRATCH=str(scratch)):
                run_singlepoint(db, job, self.root, str(executable), str(self.mpirun))
            result = dict(db.execute('SELECT * FROM jobs WHERE id=?', (job['id'],)).fetchone())
            self.assertEqual(result['status'], 'succeeded')
            self.assertEqual(result['electronic_energy_eh'], -100.5)
            self.assertAlmostEqual(result['corrected_gibbs_free_energy_eh'], -100.4)
            self.assertIn('ORCA TERMINATED NORMALLY', (self.root / result['out_file']).read_text())
            self.assertEqual(db.execute('SELECT status FROM jobs WHERE id=?', (job['id'],)).fetchone()[0], 'succeeded')
        self.assertEqual(list(scratch.iterdir()), [])

    def test_missing_configured_mpi_does_not_block_submission(self):
        self.config['executables']['mpirun_exe'] = str(self.base / 'absent' / 'mpirun')
        with patch('hpc_singlepoint.subprocess.check_output', return_value='123\n') as submit:
            result = submit_singlepoints(self.base / 'hpc.toml', self.config)
        submit.assert_called_once()
        self.assertEqual(result['pending'], 8)
        self.assertEqual(result['workers'], '123')

    def test_worker_falls_back_to_inherited_path_when_configured_mpi_is_absent(self):
        self.prepare()
        executable = self.base / 'fake-orca'
        executable.write_text(
            '#!/bin/sh\n'
            '[ "$2" = "--use-hwthread-cpus" ] || exit 20\n'
            'mpirun\n')
        executable.chmod(0o755)
        self.config['executables']['singlepoint_orca_exe'] = str(executable)
        self.config['executables']['mpirun_exe'] = str(self.base / 'absent-mpi' / 'bin' / 'mpirun')
        with patch.dict('os.environ', PATH=f'{self.mpi_bin}:/usr/bin:/bin', LD_LIBRARY_PATH='',
                        SCRATCH=str(self.base / 'node-scratch')):
            handle_singlepoints(SimpleNamespace(precursor_command='worker'), self.config,
                                self.base, self.root, argparse.ArgumentParser())
        rows = read_jobs(self.root / 'singlepoints.sqlite') + read_jobs(self.refs / 'singlepoints.sqlite')
        self.assertEqual(len(rows), 8)
        for row in rows:
            self.assertEqual(row['status'], 'succeeded', row['error'])
            self.assertEqual(row['electronic_energy_eh'], -100.5)
            self.assertAlmostEqual(row['corrected_gibbs_free_energy_eh'],
                                   row['gxtb_gibbs_free_energy_eh'] - row['gxtb_electronic_energy_eh'] - 100.5)

    def test_worker_uses_fallback_installation_when_preferred_mpi_is_absent(self):
        self.prepare()
        executable = self.base / 'fake-orca'
        executable.write_text('#!/bin/sh\nmpirun\n')
        executable.chmod(0o755)
        self.config['executables']['singlepoint_orca_exe'] = str(executable)
        self.config['executables']['mpirun_exe'] = str(self.base / 'absent-mpi' / 'bin' / 'mpirun')
        self.config['executables']['mpirun_fallback_exe'] = str(self.mpirun)
        with patch.dict('os.environ', PATH='/usr/bin:/bin', LD_LIBRARY_PATH='',
                        SCRATCH=str(self.base / 'node-scratch')):
            handle_singlepoints(SimpleNamespace(precursor_command='worker'), self.config,
                                self.base, self.root, argparse.ArgumentParser())
        rows = read_jobs(self.root / 'singlepoints.sqlite') + read_jobs(self.refs / 'singlepoints.sqlite')
        self.assertEqual(len(rows), 8)
        for row in rows:
            self.assertEqual(row['status'], 'succeeded', row['error'])
            self.assertEqual(row['electronic_energy_eh'], -100.5)

    def test_worker_with_no_mpi_leaves_all_calculations_unclaimed(self):
        self.prepare()
        executable = self.base / 'fake-orca'
        executable.write_text('#!/bin/sh\nexit 0\n')
        executable.chmod(0o755)
        self.config['executables']['singlepoint_orca_exe'] = str(executable)
        self.config['executables']['mpirun_exe'] = str(self.base / 'absent-mpi' / 'bin' / 'mpirun')
        with patch.dict('os.environ', PATH=str(self.base / 'empty-path')), \
                patch('sys.stderr', new_callable=io.StringIO) as stderr, \
                patch('hpc_singlepoint.run_singlepoint') as run:
            with self.assertRaises(SystemExit):
                handle_singlepoints(SimpleNamespace(precursor_command='worker'), self.config,
                                    self.base, self.root, argparse.ArgumentParser())
            self.assertIn('jobs were left pending', stderr.getvalue())
        run.assert_not_called()
        for row in read_jobs(self.root / 'singlepoints.sqlite') + read_jobs(self.refs / 'singlepoints.sqlite'):
            self.assertEqual(row['status'], 'pending')
            self.assertEqual(row['attempt'], 0)

    def test_worker_uses_singlepoint_orca_without_requiring_the_gxtb_executable(self):
        self.prepare()
        executable = self.base / 'dft-orca'
        executable.write_text('#!/bin/sh\nexit 0\n')
        executable.chmod(0o755)
        self.config['executables']['singlepoint_orca_exe'] = str(executable)
        self.config['executables']['orca_exe'] = str(self.base / 'missing-gxtb-orca')
        with patch('hpc_singlepoint.run_singlepoint') as run:
            handle_singlepoints(SimpleNamespace(precursor_command='worker'), self.config,
                                self.base, self.root, argparse.ArgumentParser())
        self.assertEqual(run.call_count, 8)
        for call in run.call_args_list:
            self.assertEqual(call.args[3], str(executable))
            self.assertEqual(call.args[4], str(self.mpirun))

    def test_tampered_frozen_geometry_fails_without_running_orca(self):
        self.prepare()
        with closing(connect_singlepoints(self.root / 'singlepoints.sqlite')) as db:
            job = claim_singlepoint(db)
            (self.root / job['input_xyz']).write_text(XYZ + '\n')
            with patch('hpc_singlepoint.run_orca') as run:
                run_singlepoint(db, job, self.root, 'orca')
            run.assert_not_called()
            self.assertEqual(db.execute('SELECT status FROM jobs WHERE id=?', (job['id'],)).fetchone()[0], 'failed')

    def test_parser_requires_normal_termination_and_scf_convergence(self):
        good = 'SCF CONVERGED AFTER 4 CYCLES\nFINAL SINGLE POINT ENERGY -1.23D+02\nORCA TERMINATED NORMALLY'
        self.assertEqual(collect_singlepoint(good)['electronic_energy_eh'], -123.)
        for bad in [good.replace('SCF CONVERGED AFTER', 'SCF NOT CONVERGED AFTER'),
                    good.replace('ORCA TERMINATED NORMALLY', ''),
                    good.replace('-1.23D+02', '1E999')]:
            with self.assertRaises(ValueError):
                collect_singlepoint(bad)

    def test_slurm_submission_uses_dft_resources_and_only_the_singlepoint_worker(self):
        self.prepare()
        with patch('hpc_singlepoint.subprocess.check_output', return_value='123;cluster\n') as submit:
            result = submit_singlepoints(self.base / 'hpc.toml', self.config, workers=2, resume=True)
        self.assertEqual(result['workers'], '123')
        args = submit.call_args.args[0]
        self.assertIn('--nodes=1', args)
        self.assertIn('--ntasks=4', args)
        self.assertIn('--cpus-per-task=1', args)
        self.assertIn('--no-requeue', args)
        self.assertIn('--mem=8000', args)
        self.assertIn('--array=0-1%2', args)
        self.assertIn('singlepoint worker', args[-1])


if __name__ == '__main__':
    unittest.main()
