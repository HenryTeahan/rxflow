"""Regression checks for direct OptTS recovery without running ORCA."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hpc_workflow import connect_database, run_job
from orca_hpc import orca_input


class DirectOptTSRetryTests(unittest.TestCase):
    def run_candidate(self, results, attempts=2, strategy="direct"):
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
                ("ligand_01/ts_guess/input.xyz", json.dumps({"ts_attempts": attempts, "neb_fallback": True}), strategy),
            )
            job = dict(connection.execute("SELECT * FROM jobs").fetchone())
            if strategy == "neb":
                for role in ("reactant", "product"):
                    reactive = root / "ligand_01" / role / "input.xyz_reactive.xyz"
                    reactive.parent.mkdir(parents=True)
                    reactive.write_text(f"original {role}")
                    connection.execute(
                        """INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                           input_xyz,smiles,charge,settings_json) VALUES ('1',1,1,?,1,?,'unused',1,'{}')""",
                        (role, str(reactive.relative_to(root))),
                    )
            with patch("hpc_workflow.run_orca", side_effect=results) as calculate, \
                    patch("hpc_workflow.molecule_with_xyz"):
                run_job(connection, job, root, {})
            row = dict(connection.execute("SELECT * FROM jobs WHERE id=?", (job["id"],)).fetchone())
            if strategy == "neb":
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

    def test_initial_optts_continues_then_runs_frequency_without_recovery(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("intermediate geometry"),
            self.optimization("converged geometry", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["optts", "optts", "freq"])
        self.assertEqual([c.args[4] for c in calls],
                         ["initial geometry", "intermediate geometry", "converged geometry"])
        self.assertEqual([c.kwargs["attempt"] for c in calls], [1, 2, 1])
        self.assertEqual(optimized, "converged geometry")

    def test_complete_recovery_sequence_restarts_from_seed_and_preserves_attempt_numbers(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("failed direct 1"), self.optimization("failed direct 2"),
            self.optimization("partial relaxation"),
            self.optimization("relaxed geometry", converged=True),
            self.optimization("partial recovered TS"),
            self.optimization("recovered TS", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls],
                         ["optts", "optts", "relax", "relax", "optts", "optts", "freq"])
        self.assertEqual([c.args[4] for c in calls],
                         ["initial geometry", "failed direct 1", "initial geometry", "partial relaxation",
                          "relaxed geometry", "partial recovered TS", "recovered TS"])
        self.assertEqual([c.kwargs["attempt"] for c in calls], [1, 2, 1, 2, 3, 4, 2])
        self.assertEqual(optimized, "recovered TS")

    def test_every_phase_respects_configured_budget_and_stops_when_relaxation_fails(self):
        for attempts in (1, 2, 3):
            with self.subTest(attempts=attempts):
                row, calls, optimized = self.run_candidate(
                    [self.optimization(f"direct {i}") for i in range(attempts)]
                    + [self.optimization(f"relaxation {i}") for i in range(attempts)], attempts)
                self.assertEqual(row["status"], "failed")
                self.assertEqual(row["stage"], "relax")
                self.assertEqual(row["error"], "Constrained TS relaxation did not converge")
                self.assertEqual([c.args[3] for c in calls],
                                 ["optts"] * attempts + ["relax"] * attempts)
                self.assertEqual([c.kwargs["attempt"] for c in calls],
                                 list(range(1, attempts + 1)) * 2)
                self.assertEqual(calls[attempts].args[4], "initial geometry")
                self.assertIsNone(optimized)

    def test_failure_after_recovered_optts_does_not_run_frequency_or_neb(self):
        row, calls, _ = self.run_candidate([
            self.optimization("direct 1"), self.optimization("direct 2"),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered 1"), self.optimization("recovered 2"),
        ])
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["stage"], "optts")
        self.assertEqual(row["error"], "OptTS after constrained TS relaxation did not converge")
        self.assertEqual([c.args[3] for c in calls], ["optts", "optts", "relax", "optts", "optts"])
        self.assertEqual([c.kwargs["attempt"] for c in calls[-2:]], [3, 4])
        self.assertIsNone(row["gibbs_free_energy_eh"])

    def test_stops_immediately_after_success_at_each_phase(self):
        row, calls, _ = self.run_candidate([
            self.optimization("direct", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["optts", "freq"])
        row, calls, _ = self.run_candidate([
            self.optimization("direct 1"), self.optimization("direct 2"),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["optts", "optts", "relax", "optts", "freq"])

    def test_wrong_initial_frequency_count_triggers_recovery_from_original_guess(self):
        for imaginary_count in (0, 2):
            with self.subTest(imaginary_count=imaginary_count):
                row, calls, optimized = self.run_candidate([
                    self.optimization("direct TS", converged=True), self.frequency(imaginary_count),
                    self.optimization("relaxed", converged=True),
                    self.optimization("recovered", converged=True), self.frequency(),
                ])
                self.assertEqual(row["status"], "succeeded")
                self.assertEqual(row["imaginary_count"], 1)
                self.assertEqual([c.args[3] for c in calls], ["optts", "freq", "relax", "optts", "freq"])
                self.assertEqual(calls[2].args[4], "initial geometry")
                self.assertEqual([c.kwargs["attempt"] for c in calls if c.args[3] == "freq"], [1, 2])
                self.assertEqual(optimized, "recovered")

    def test_malformed_initial_numfreq_triggers_recovery(self):
        row, calls, _ = self.run_candidate([
            self.optimization("direct TS", converged=True),
            ValueError("ORCA NumFreq did not terminate normally"),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["optts", "freq", "relax", "optts", "freq"])
        self.assertEqual(calls[2].args[4], "initial geometry")

    def test_invalid_final_frequency_fails_without_another_recovery_cycle(self):
        row, calls, _ = self.run_candidate([
            self.optimization("direct 1"), self.optimization("direct 2"),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered", converged=True), self.frequency(2),
        ])
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["stage"], "freq")
        self.assertEqual(row["error"], "Expected 1 imaginary frequencies, found 2")
        self.assertEqual(row["imaginary_count"], 2)
        self.assertEqual(len(calls), 5)

    def test_malformed_final_numfreq_fails_without_another_recovery_cycle(self):
        row, calls, _ = self.run_candidate([
            self.optimization("direct 1"), self.optimization("direct 2"),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered", converged=True), ValueError("incomplete frequency table"),
        ])
        self.assertEqual(row["status"], "failed")
        self.assertEqual(row["stage"], "freq")
        self.assertEqual(row["error"], "incomplete frequency table")
        self.assertEqual(len(calls), 5)

    def test_missing_geometry_repeats_previous_input_then_recovers_from_seed(self):
        row, calls, _ = self.run_candidate([
            self.optimization(None), self.optimization(None),
            self.optimization("relaxed", converged=True),
            self.optimization("recovered", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[4] for c in calls[:3]], ["initial geometry"] * 3)

    def test_non_cycle_limit_nonconvergence_also_gets_a_repeat(self):
        row, calls, _ = self.run_candidate([
            self.optimization("partial", max_cycles=False),
            self.optimization("direct TS", converged=True), self.frequency(),
        ])
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["optts", "optts", "freq"])
        self.assertEqual(calls[1].args[4], "partial")

    def test_failed_recovery_clears_geometry_and_energies_from_initial_frequency(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("invalid direct TS", converged=True), self.frequency(2),
            self.optimization("relax 1"), self.optimization("relax 2"),
        ])
        self.assertEqual(row["status"], "failed")
        self.assertIsNone(row["xyz_file"])
        self.assertIsNone(row["gibbs_free_energy_eh"])
        self.assertIsNone(row["imaginary_count"])
        self.assertIsNone(optimized)
        self.assertEqual([c.args[3] for c in calls], ["optts", "freq", "relax", "relax"])

    def test_constrained_input_freezes_only_reaction_centre_and_optts_releases_it(self):
        xyz = "3\n\nC 0 0 0\nH 1 0 0\nN 2 0 0\n"
        settings = dict(xtb_exe="xtb", external_exe="oet_gxtb_alpb", solvent="h2o", cpus=1, memory_mb=2000)
        relaxation = orca_input(xyz, "relax", 1, 1, settings, atoms=(0, 1, 2))
        self.assertIn("! ExtOpt Opt\n", relaxation)
        self.assertIn("Constraints\n  {C 0 C}\n  {C 1 C}\n  {C 2 C}\n end", relaxation)
        optts = orca_input(xyz, "optts", 1, 1, settings, atoms=(0, 1, 2))
        self.assertNotIn("Constraints", optts)
        self.assertIn("Calc_Hess true", optts)
        self.assertIn("NumHess true", optts)

    def test_neb_still_runs_frequency_without_overwriting_its_endpoints(self):
        row, calls, optimized = self.run_candidate([
            self.optimization("NEB TS", converged=True), self.frequency(),
        ], strategy="neb")
        self.assertEqual(row["status"], "succeeded")
        self.assertEqual([c.args[3] for c in calls], ["neb", "freq"])
        self.assertEqual(calls[0].kwargs["neb_endpoints"], ["original reactant", "original product"])
        self.assertEqual(optimized, "NEB TS")
        self.assertEqual(row["xyz_file"], "ligand_01/ts_guess/job_1_run_01_optimized.xyz")


if __name__ == "__main__":
    unittest.main()
