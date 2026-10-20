"""TS preoptimization, OptTS retries and preoptimized NEB endpoints."""

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hpc_workflow import claim_job, connect_database, regenerate_ts_guess, run_job
from orca_hpc import orca_input, ts_preoptimization_atoms


class TSWorkflowTests(unittest.TestCase):
    def run_candidate(self, results, attempts=2, strategy="direct", mode_checks=None,
                      recovery=None, neb_fallback=True, changes=None, follow_neb=False, refresh_guess=False,
                      scan_recovery=False, scan_error=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / "ligand_01" / "ts_guess"
            folder.mkdir(parents=True)
            (folder / "input.xyz").write_text("initial geometry")
            connection = connect_database(":memory:")
            self.addCleanup(connection.close)
            connection.execute(
                """INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                   status,attempt,input_xyz,smiles,charge,settings_json,ts_strategy)
                   VALUES ('1',1,1,'TS',1,'running',1,?, 'unused',1,?,?)""",
                ("ligand_01/ts_guess/input.xyz", json.dumps({"ts_attempts": attempts,
                 "opt_attempts": 2, "neb_fallback": neb_fallback, "ts_scan_recovery": scan_recovery}), strategy))
            job = dict(connection.execute("SELECT * FROM jobs").fetchone())
            job["ts_mode_recovery"] = recovery
            if refresh_guess:
                job["stage"] = "ts_guess"
            for role in ("reactant", "product"):
                reactive = root / "ligand_01" / role / "input.xyz_reactive.xyz"
                reactive.parent.mkdir(parents=True)
                reactive.write_text(f"original {role}")
                connection.execute(
                    """INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                       input_xyz,smiles,charge,settings_json) VALUES ('1',1,1,?,1,?,?,1,'{}')""",
                    (role, str(reactive.relative_to(root)), f"{role} graph"))
            if follow_neb:
                connection.execute("""INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                    input_xyz,smiles,charge,settings_json,ts_strategy)
                    VALUES ('1',1,1,'TS',0,?,'unused',1,?,'neb')""",
                    (job["input_xyz"], job["settings_json"]))
            valid_mode = dict(valid=True, fingerprint="mock", reason="All required TS atoms participate")

            def score(*args):
                if scan_error:
                    raise ValueError(scan_error)
                (folder / 'scored_guess.xyz').write_text('scored interior peak')
                return dict(valid=True, guess_xyz='ligand_01/ts_guess/scored_guess.xyz',
                            selected_energy_ev=-10., selected_frame=7,
                            reactant_xyz='ligand_01/reactant/input.xyz_reactive.xyz',
                            product_xyz='ligand_01/product/input.xyz_reactive.xyz')

            def refresh(connection, original_job, results_root, source, sampling):
                fresh = folder / "fresh_rxts_guess.xyz"
                fresh.write_text("fresh rxts geometry")
                return {**original_job, "input_xyz": str(fresh.relative_to(root)),
                        "ts_guess_json": json.dumps({f"{role}_xyz": f"ligand_01/{role}/input.xyz_reactive.xyz"
                                                     for role in ("reactant", "product")})}

            with patch("hpc_workflow.run_orca", side_effect=results) as calculate, \
                    patch("hpc_workflow.check_ts_mode", side_effect=mode_checks, return_value=valid_mode), \
                    patch("hpc_workflow.molecule_with_xyz"), \
                    patch("hpc_workflow.connectivity_changes", return_value=changes or {}), \
                    patch("hpc_workflow.score_ts_scan", side_effect=score), \
                    patch("hpc_workflow.regenerate_ts_guess", side_effect=refresh):
                run_job(connection, job, root, {})
                if follow_neb:
                    self.assertEqual(connection.execute("SELECT status FROM jobs WHERE id=?", (job["id"],)).fetchone()[0],
                                     "failed")
                    connection.execute("UPDATE jobs SET status='succeeded' WHERE endpoint!='TS'")
                    job = claim_job(connection)
                    self.assertIsNotNone(job)
                    self.assertEqual(job["ts_strategy"], "neb")
                    run_job(connection, job, root, {})
            row = dict(connection.execute("SELECT * FROM jobs WHERE id=?", (job["id"],)).fetchone())
            for role in ("reactant", "product"):
                reactive = root / "ligand_01" / role / "input.xyz_reactive.xyz"
                self.assertEqual(reactive.read_text(), f"original {role}")
            optimized = (root / row["xyz_file"]).read_text() if row["xyz_file"] else None
            return row, calculate.call_args_list, optimized

    @staticmethod
    def optimization(xyz, *, converged=False, max_cycles=True):
        return dict(out_file="optimization.out", xyz=xyz, converged=converged,
                    max_cycles_reached=max_cycles)

    @staticmethod
    def frequency(imaginary_count=1):
        return dict(out_file="frequency.out", electronic_energy_eh=-1.0,
                    gibbs_free_energy_eh=-0.9, frequencies=[-500.0] * imaginary_count,
                    imaginary_count=imaginary_count)

    def neb_results(self):
        return [self.optimization("optimized reactant", converged=True),
                self.optimization("optimized product", converged=True),
                self.optimization("NEB TS", converged=True), self.frequency()]

    def test_first_step_preoptimizes_before_optts_and_stops_after_success(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("preoptimized guess", converged=True),
            self.optimization("converged TS", converged=True), self.frequency()])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["ts_preopt", "optts", "freq"])
        self.assertEqual([c.args[4] for c in calls],
                         ["initial geometry", "preoptimized guess", "converged TS"])
        self.assertEqual(optimized, "converged TS")

    def test_reset_candidate_regenerates_guess_before_any_orca_stage(self):
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized guess", converged=True),
            self.optimization("TS", converged=True), self.frequency()], refresh_guess=True)
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual(calls[0].args[4], 'fresh rxts geometry')
        self.assertEqual([c.args[3] for c in calls], ['ts_preopt', 'optts', 'freq'])

    def test_reset_neb_uses_matching_pair_from_regenerated_guess(self):
        row, calls, _ = self.run_candidate(self.neb_results(), strategy='neb', refresh_guess=True)
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual([c.args[3] for c in calls], ['neb_reactant_opt', 'neb_product_opt', 'neb', 'freq'])
        self.assertEqual(calls[2].args[4], 'fresh rxts geometry')
        self.assertEqual(calls[2].kwargs['neb_endpoints'], ['optimized reactant', 'optimized product'])

    def test_cycle_limited_preoptimization_uses_last_geometry_then_two_optts(self):
        row, calls, _ = self.run_candidate([
            self.optimization("50-cycle geometry"), self.optimization("partial TS"),
            self.optimization("converged TS", converged=True), self.frequency()])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["ts_preopt", "optts", "optts", "freq"])
        self.assertEqual([c.args[4] for c in calls],
                         ["initial geometry", "50-cycle geometry", "partial TS", "converged TS"])
        self.assertEqual([c.kwargs["attempt"] for c in calls], [1, 1, 2, 1])

    def test_optts_budget_is_preserved_without_second_ts_recovery_cycle(self):
        for attempts in (1, 2, 3):
            with self.subTest(attempts=attempts):
                row, calls, _ = self.run_candidate([
                    self.optimization("preoptimized", converged=True),
                    *[self.optimization(f"TS attempt {i}") for i in range(attempts)]], attempts=attempts)
                self.assertEqual(row["status"], "failed")
                self.assertEqual(row["error"], "OptTS attempts did not converge")
                self.assertEqual([c.args[3] for c in calls], ["ts_preopt"] + ["optts"] * attempts)
                self.assertIsNone(row["gibbs_free_energy_eh"])

    def test_preoptimization_aborts_without_geometry_or_on_non_cycle_failure(self):
        for result in (self.optimization(None), self.optimization("partial", max_cycles=False)):
            with self.subTest(result=result):
                row, calls, _ = self.run_candidate([result])
                self.assertEqual(row["status"], "failed")
                self.assertIn("preoptimization failed", row["error"])
                self.assertEqual([c.args[3] for c in calls], ["ts_preopt"])

    def test_missing_optts_geometry_reuses_preoptimized_input_for_retry(self):
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True), self.optimization(None),
            self.optimization("TS", converged=True), self.frequency()])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[4] for c in calls[1:3]], ["preoptimized", "preoptimized"])

    def test_audited_candidate_uses_same_new_first_step(self):
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("TS", converged=True), self.frequency()], recovery="relax")
        self.assertEqual(row["status"], "succeeded")
        self.assertIsNone(row["ts_mode_recovery"])
        self.assertEqual([c.args[3] for c in calls], ["ts_preopt", "optts", "freq"])

    def test_wrong_frequency_count_rejects_candidate_without_extra_optts_phase(self):
        for count in (0, 2):
            with self.subTest(count=count):
                row, calls, _ = self.run_candidate([
                    self.optimization("preoptimized", converged=True),
                    self.optimization("TS", converged=True), self.frequency(count)])
                self.assertEqual(row["status"], "failed")
                self.assertEqual(row["imaginary_count"], count)
                self.assertEqual([c.args[3] for c in calls], ["ts_preopt", "optts", "freq"])

    def test_wrong_mode_rejects_direct_candidate_for_shared_neb_fallback(self):
        bad = dict(valid=False, fingerprint="bad", reason="TS imaginary mode is missing required atoms [1]")
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("wrong-mode TS", converged=True), self.frequency()], mode_checks=[bad])
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["ts_mode_status"], "failed")
        self.assertIn("missing required atoms", row["error"])
        self.assertEqual([c.args[3] for c in calls], ["ts_preopt", "optts", "freq"])

    def test_unreadable_mode_or_frequency_is_rejected_for_neb(self):
        for mode_checks, freq in (([IndexError("incomplete mode matrix")], self.frequency()),
                                  (None, ValueError("incomplete frequency table"))):
            with self.subTest(mode_checks=mode_checks):
                row, calls, _ = self.run_candidate([
                    self.optimization("preoptimized", converged=True),
                    self.optimization("TS", converged=True), freq], mode_checks=mode_checks)
                self.assertEqual(row["status"], "failed")
                self.assertEqual(len(calls), 3)
                if mode_checks:
                    self.assertEqual(row["ts_mode_status"], "error")

    def test_historical_failed_optts_preoptimizes_coembedded_endpoints_then_neb(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("TS 1"), self.optimization("TS 2"),
            *self.neb_results()], strategy=None)
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls],
                         ["ts_preopt", "optts", "optts", "neb_reactant_opt", "neb_product_opt", "neb", "freq"])
        self.assertEqual(calls[3].args[4], "original reactant")
        self.assertEqual(calls[4].args[4], "original product")
        self.assertEqual(calls[5].args[4], "initial geometry")
        self.assertEqual(calls[5].kwargs["neb_endpoints"], ["optimized reactant", "optimized product"])
        self.assertEqual(optimized, "NEB TS")

    def test_wrong_mode_uses_same_preoptimized_endpoint_neb_fallback(self):
        bad = dict(valid=False, fingerprint="bad", reason="TS imaginary mode is missing required atoms [1]")
        good = dict(valid=True, fingerprint="good", reason="All required TS atoms participate")
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("wrong-mode TS", converged=True), self.frequency(),
            *self.neb_results()], strategy=None, mode_checks=[bad, good])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual(row["ts_mode_fingerprint"], "good")
        self.assertEqual([c.args[3] for c in calls],
                         ["ts_preopt", "optts", "freq", "neb_reactant_opt", "neb_product_opt", "neb", "freq"])
        self.assertEqual([c.kwargs["attempt"] for c in calls if c.args[3] == "freq"], [1, 2])

    def test_wrong_direct_mode_releases_shared_neb_with_preoptimized_endpoints(self):
        bad = dict(valid=False, fingerprint="bad", reason="TS imaginary mode is missing required atoms [1]")
        good = dict(valid=True, fingerprint="good", reason="All required TS atoms participate")
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("wrong-mode TS", converged=True), self.frequency(),
            *self.neb_results()], follow_neb=True, mode_checks=[bad, good])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual(row["ts_strategy"], "neb")
        self.assertEqual(row["ts_mode_status"], "passed")
        self.assertEqual([c.args[3] for c in calls],
                         ["ts_preopt", "optts", "freq", "neb_reactant_opt", "neb_product_opt", "neb", "freq"])
        self.assertEqual(calls[-2].kwargs["neb_endpoints"], ["optimized reactant", "optimized product"])

    def test_neb_preoptimizes_each_original_endpoint_using_its_own_graph(self):
        row, calls, optimized = self.run_candidate(self.neb_results(), strategy="neb")
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["neb_reactant_opt", "neb_product_opt", "neb", "freq"])
        for call, role in zip(calls, ("reactant", "product")):
            self.assertEqual(call.args[0]["endpoint"], role)
            self.assertEqual(call.args[0]["smiles"], f"{role} graph")
            self.assertEqual(call.args[4], f"original {role}")
        self.assertEqual(calls[2].kwargs["neb_endpoints"], ["optimized reactant", "optimized product"])
        self.assertEqual(optimized, "NEB TS")

    def test_neb_endpoint_optimization_restarts_once_from_its_last_geometry(self):
        row, calls, _ = self.run_candidate([
            self.optimization("partial reactant"), *self.neb_results()], strategy="neb")
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls[:3]],
                         ["neb_reactant_opt", "neb_reactant_opt", "neb_product_opt"])
        self.assertEqual(calls[1].args[4], "partial reactant")
        self.assertEqual([c.kwargs["attempt"] for c in calls[:2]], [1, 2])

    def test_unconverged_endpoint_stops_before_neb(self):
        row, calls, _ = self.run_candidate([
            self.optimization("reactant 1"), self.optimization("reactant 2")], strategy="neb")
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["error"], "ORCA endpoint optimization did not converge")
        self.assertEqual([c.args[3] for c in calls], ["neb_reactant_opt", "neb_reactant_opt"])

    def test_endpoint_connectivity_check_is_kept_before_neb(self):
        row, calls, _ = self.run_candidate([
            self.optimization("changed reactant", converged=True)], strategy="neb", changes={"bonds": [1]})
        self.assertEqual(row["status"], "failed")
        self.assertIn("Graph change", row["error"])
        self.assertEqual([c.args[3] for c in calls], ["neb_reactant_opt"])

    def test_disabled_neb_does_not_prepare_endpoints(self):
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("TS 1"), self.optimization("TS 2")], strategy=None, neb_fallback=False)
        self.assertEqual(row["status"], "failed")
        self.assertEqual(len(calls), 3)

    def test_audited_neb_uses_preoptimization_and_two_optts_before_fallback(self):
        row, calls, _ = self.run_candidate([
            self.optimization("preoptimized", converged=True),
            self.optimization("TS 1"), self.optimization("TS 2"),
            *self.neb_results()], strategy="neb", recovery="relax")
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls],
                         ["ts_preopt", "optts", "optts", "neb_reactant_opt", "neb_product_opt", "neb", "freq"])

    def test_neb_wrong_mode_is_not_accepted(self):
        bad = dict(valid=False, fingerprint="bad", reason="TS imaginary mode is missing required atoms [1]")
        row, calls, _ = self.run_candidate(self.neb_results(), strategy="neb", mode_checks=[bad])
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["ts_mode_status"], "failed")
        self.assertEqual(calls[-1].args[3], "freq")

    def test_gxtb_recovery_runs_after_two_failures_and_stops_on_valid_ts(self):
        row, calls, _ = self.run_candidate([
            self.optimization('preopt', converged=True), self.optimization('failed 1'),
            self.optimization('failed 2'), self.optimization('scored preopt', converged=True),
            self.optimization('failed 3'), self.optimization('scored TS', converged=True),
            self.frequency()], scan_recovery=True)
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual([c.args[3] for c in calls],
                         ['ts_preopt', 'optts', 'optts', 'ts_preopt', 'optts', 'optts', 'freq'])
        self.assertEqual([c.kwargs['attempt'] for c in calls], [1, 1, 2, 2, 3, 4, 2])
        self.assertEqual(calls[3].args[4], 'scored interior peak')
        self.assertEqual(json.loads(row['ts_scan_json'])['selected_frame'], 7)

    def test_gxtb_is_not_used_when_initial_optts_succeeds(self):
        row, calls, _ = self.run_candidate([
            self.optimization('preopt', converged=True), self.optimization('TS', converged=True),
            self.frequency()], scan_recovery=True, scan_error='must not be called')
        self.assertEqual(row['status'], 'succeeded')
        self.assertIsNone(row['ts_scan_json'])

    def test_scoring_failure_does_not_invent_a_guess_or_launch_neb(self):
        row, calls, _ = self.run_candidate([
            self.optimization('preopt', converged=True), self.optimization('failed 1'),
            self.optimization('failed 2')], strategy=None, scan_recovery=True,
            scan_error='No interior g-xTB local maximum')
        self.assertEqual(row['status'], 'failed')
        self.assertIn('No interior', row['error'])
        self.assertEqual(len(calls), 3)
        self.assertFalse(json.loads(row['ts_scan_json'])['valid'])

    def test_shared_neb_reuses_failed_candidates_scored_guess_and_pair(self):
        row, calls, _ = self.run_candidate([
            self.optimization('preopt', converged=True), self.optimization('failed 1'),
            self.optimization('failed 2'), self.optimization('scored preopt', converged=True),
            self.optimization('failed 3'), self.optimization('failed 4'),
            *self.neb_results()], scan_recovery=True, follow_neb=True)
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual(calls[-2].args[4], 'scored interior peak')
        self.assertEqual(calls[-2].kwargs['neb_endpoints'], ['optimized reactant', 'optimized product'])
        self.assertEqual(json.loads(row['ts_scan_json'])['source_job_id'], 1)

    def test_failed_frequency_validation_uses_scored_optts_before_neb(self):
        row, calls, _ = self.run_candidate([
            self.optimization('preopt', converged=True), self.optimization('bad TS', converged=True),
            self.frequency(2), self.optimization('scored preopt', converged=True),
            self.optimization('good TS', converged=True), self.frequency()], scan_recovery=True)
        self.assertEqual(row['status'], 'succeeded')
        self.assertEqual(row['imaginary_count'], 1)
        self.assertEqual(calls[3].args[4], 'scored interior peak')


class GuessRegenerationTests(unittest.TestCase):
    SMILES = '[C:1]([H:2])#[N:3]'
    XYZ = '3\nfixture\nC 0 0 0\nH 1.1 0 0\nN 2.7 0 0\n'

    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)/'results'
        self.source = Path(folder.name)/'seeds'
        (self.root/'ligand_01/ts_guess').mkdir(parents=True)
        (self.source/'ligand_01').mkdir(parents=True)
        self.original = self.root/'ligand_01/ts_guess/original_input.xyz'
        self.original.write_text(self.XYZ)
        with (self.source/'ligand_01/reaction_pairs.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['reactant_smiles'])
            writer.writeheader()
            writer.writerow(dict(reactant_smiles=self.SMILES))
        self.db = connect_database(':memory:')
        self.addCleanup(self.db.close)
        self.db.execute("""INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,status,stage,
            attempt,input_xyz,smiles,charge,multiplicity,settings_json,ts_strategy)
            VALUES ('1',1,1,'TS',1,'running','ts_guess',2,?,?,0,1,'{}','direct')""",
            (str(self.original.relative_to(self.root)), self.SMILES))
        self.job = dict(self.db.execute('SELECT * FROM jobs').fetchone())

    def generate(self, left, right, output, **kwargs):
        for role, name in (('ts_guess', 'reaction.xyz'), ('reactant', 'reaction_reactive.xyz'),
                           ('product', 'reaction_reactive.xyz')):
            (output/role).mkdir(parents=True, exist_ok=True)
            (output/role/name).write_text(self.XYZ)
        (output/'ts_guess/reaction_selection.json').write_text('{}')
        return dict(selected_ts_rmsd=.2, seed=kwargs['seed'], ts_mode=dict(donor=0, hydrogen=1, acceptor=2))

    def test_existing_job_regenerates_original_seating_and_keeps_old_input(self):
        with patch('hpc_workflow.create_endpoint_isomers', return_value=[(0, 'reactant isomer', 'product isomer')]), \
                patch('hpc_workflow.generate_ts_guess', side_effect=self.generate) as generate:
            job = regenerate_ts_guess(self.db, self.job, self.root, self.source,
                                       dict(seed=42, coembed_conformers=10))
        self.assertEqual(generate.call_args.args[:2], ('reactant isomer', 'product isomer'))
        self.assertEqual(generate.call_args.kwargs['n'], 10)
        self.assertNotEqual(generate.call_args.kwargs['seed'], 42)
        self.assertNotEqual(job['input_xyz'], self.job['input_xyz'])
        self.assertEqual(job['attempt'], self.job['attempt'])
        self.assertEqual(self.original.read_text(), self.XYZ)
        metadata = json.loads(job['ts_guess_json'])
        self.assertTrue(metadata['valid'])
        for key in ('reactant_xyz', 'product_xyz', 'selection_file'):
            self.assertTrue((self.root/metadata[key]).is_file())
        self.assertEqual(self.db.execute('SELECT input_xyz FROM jobs').fetchone()[0], job['input_xyz'])

    def test_exhausted_guess_generation_fails_before_any_orca_calculation(self):
        with patch('hpc_workflow.create_endpoint_isomers', return_value=[(0, 'left', 'right')]), \
                patch('hpc_workflow.generate_ts_guess', side_effect=ValueError('No TS guesses survived after 3 attempts')), \
                patch('hpc_workflow.run_orca') as orca:
            run_job(self.db, self.job, self.root, {}, source=self.source)
        orca.assert_not_called()
        row = self.db.execute('SELECT * FROM jobs').fetchone()
        self.assertEqual((row['status'], row['stage']), ('failed', 'ts_guess'))
        self.assertEqual(row['input_xyz'], self.job['input_xyz'])
        self.assertIn('3 attempts', row['error'])

    def test_source_resonance_notation_does_not_discard_the_original_binding(self):
        original = '[O-:1][C:2](=[N+:3]([H:5])[H:6])[C:4]([H:7])([H:8])[H:9]'
        stored = '[O:1]=[C:2]([N:3]([H:5])[H:6])[C:4]([H:7])([H:8])[H:9]'
        with (self.source/'ligand_01/reaction_pairs.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['reactant_smiles'])
            writer.writeheader()
            writer.writerow(dict(reactant_smiles=original))
        job = {**self.job, 'smiles': stored}
        with patch('hpc_workflow.create_endpoint_isomers', return_value=[(0, 'left', 'right')]), \
                patch('hpc_workflow.generate_ts_guess', side_effect=RuntimeError('reached regeneration')):
            with self.assertRaisesRegex(RuntimeError, 'reached regeneration'):
                regenerate_ts_guess(self.db, job, self.root, self.source, {})


class TSPreoptimizationInputTests(unittest.TestCase):
    SMILES = '[O:1]=[C:2]([H:3])[N:4]([H:5])([H:6])->[Pd+2:7]<-[C-:8]([H:9])([H:10])[H:11]'
    SETTINGS = dict(xtb_exe="xtb", external_exe="oet_gxtb_alpb", solvent="h2o", cpus=1, memory_mb=2000)
    XYZ = '3\n\nC 0 0 0\nH 1 0 0\nN 2 0 0\n'

    def test_freezes_o_c_n_metal_path_plus_donor_and_transferring_h_only(self):
        atoms = ts_preoptimization_atoms(dict(smiles=self.SMILES, donor=7, hydrogen=8, acceptor=0))
        self.assertEqual(atoms, (0, 1, 3, 6, 7, 8))
        text = orca_input(self.XYZ, 'ts_preopt', 1, 1, self.SETTINGS, atoms)
        self.assertIn('%geom\n MaxIter 50\n Constraints\n', text)
        for index in atoms:
            self.assertIn(f'{{C {index} C}}', text)
        for index in (2, 4, 5, 9, 10):
            self.assertNotIn(f'{{C {index} C}}', text)

    def test_directly_bound_acceptor_does_not_freeze_unrelated_ligand_atoms(self):
        self.assertEqual(ts_preoptimization_atoms(dict(smiles=self.SMILES, donor=7, hydrogen=8, acceptor=3)),
                         (3, 6, 7, 8))

    def test_missing_metal_or_disconnected_acceptor_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'central coordination metal'):
            ts_preoptimization_atoms(dict(smiles='[C:1]([H:2])#[N:3]', donor=0, hydrogen=1, acceptor=2))
        smiles = '[O-:1].[N:2]([H:3])([H:4])([H:5])->[Pd+2:6]<-[C-:7]([H:8])([H:9])[H:10]'
        with self.assertRaisesRegex(ValueError, 'no graph path'):
            ts_preoptimization_atoms(dict(smiles=smiles, donor=6, hydrogen=7, acceptor=0))

    def test_only_ts_preoptimization_has_cycle_limit_and_all_stages_keep_alpb(self):
        for stage in ('ts_preopt', 'optts', 'neb_reactant_opt', 'neb_product_opt', 'neb', 'opt', 'freq'):
            with self.subTest(stage=stage):
                text = orca_input(self.XYZ, stage, 1, 1, self.SETTINGS, atoms=(0, 1, 2))
                self.assertIn('ProgExt "oet_gxtb_alpb"', text)
                self.assertIn('Ext_Params "-e xtb -s h2o"', text)
                self.assertEqual('MaxIter' in text, stage == 'ts_preopt')
                self.assertEqual('Constraints' in text, stage == 'ts_preopt')


if __name__ == '__main__':
    unittest.main()
