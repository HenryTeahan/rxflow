"""RMSD pruning and the once-per-group TS conformational sampling stage."""

import csv
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from conformer_sampling import prune_conformers
from hpc_sampling import sample_ligand, sample_transition_state
from hpc_workflow import (claim_job, connect_database, prepare_jobs, reset_failed_transition_states,
                          run_job, select_transition_state)


def molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


class PruningTests(unittest.TestCase):
    def test_rigid_duplicates_merge_without_modifying_coordinates(self):
        original = molecule()
        translated = Chem.Mol(original)
        conf = translated.GetConformer()
        for index, pos in enumerate(conf.GetPositions()):
            conf.SetAtomPosition(index, pos + [10., -3., 2.])
        different = Chem.Mol(original)
        different.GetConformer().SetAtomPosition(0, [10., 20., 30.])
        before = [mol.GetConformer().GetPositions().copy() for mol in (original, translated, different)]
        kept = prune_conformers([original, translated, different])
        self.assertEqual(len(kept), 2)
        self.assertIn(different, kept)
        for mol, expected in zip((original, translated, different), before):
            np.testing.assert_array_equal(mol.GetConformer().GetPositions(), expected)

    def test_only_rmsd_is_passed_to_prism(self):
        geometries = [molecule(), molecule()]
        with patch("prism_pruner.pruner.prune_by_rmsd", return_value=(None, [True, False])) as prism:
            self.assertEqual(prune_conformers(geometries, 0.3), geometries[:1])
        self.assertEqual(prism.call_args.kwargs, dict(max_rmsd=0.3, max_dev=float("inf")))

    def test_distinct_stereoisomers_are_never_merged(self):
        left = Chem.AddHs(Chem.MolFromSmiles("C[C@H](F)Cl"))
        right = Chem.AddHs(Chem.MolFromSmiles("C[C@@H](F)Cl"))
        AllChem.EmbedMolecule(left, randomSeed=12)
        right.AddConformer(Chem.Conformer(left.GetConformer()))
        self.assertEqual(prune_conformers([left, right]), [left, right])

    def test_real_reembedding_preserves_fixed_ts_atoms(self):
        mol = molecule()
        job = dict(smiles=Chem.MolToSmiles(mol), charge=0, donor=0, hydrogen=4, acceptor=1)
        geometries, report = sample_transition_state(job, Chem.MolToXYZBlock(mol), n=4)
        self.assertGreater(len(geometries), 0)
        self.assertLessEqual(len(geometries), 4)
        self.assertEqual(report["fixed_atoms"], [0, 4, 1])
        reference = Chem.MolFromXYZBlock(Chem.MolToXYZBlock(mol)).GetConformer().GetPositions()[[0, 4, 1]]
        for geometry in geometries:
            np.testing.assert_allclose(geometry.GetConformer().GetPositions()[[0, 4, 1]], reference, atol=1e-6)


class TSSamplingTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.db = connect_database(self.root / "workflow.sqlite")
        self.addCleanup(self.db.close)
        self.mol = molecule()
        (self.root / "optimized.xyz").write_text(Chem.MolToXYZBlock(self.mol))

    def add_ts(self, conformer, energy=-10., status="succeeded", mode="passed", strategy="direct"):
        cursor = self.db.execute("""INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
            input_xyz,xyz_file,smiles,charge,settings_json,donor,hydrogen,acceptor,
            status,stage,gibbs_free_energy_eh,imaginary_count,ts_mode_status,ts_strategy)
            VALUES ('1',1,1,'TS',?,'optimized.xyz','optimized.xyz',?,0,?,0,4,1,?,'complete',?,1,?,?)""",
            (conformer, Chem.MolToSmiles(self.mol), json.dumps(dict(ts_attempts=2, neb_fallback=True)),
             status, energy, mode, strategy))
        return cursor.lastrowid

    def test_sampling_waits_and_selects_lowest_valid_ts_once(self):
        self.add_ts(1, -10.)
        best = self.add_ts(2, -11.)
        active = self.add_ts(3, -12., status="running")
        self.assertIsNone(claim_job(self.db, {}))
        self.db.execute("UPDATE jobs SET status='failed' WHERE id=?", (active,))
        task = claim_job(self.db, {})
        self.assertEqual(task["ts_strategy"], "sampling")
        self.assertEqual(task["ts_parent_id"], best)
        self.assertEqual(json.loads(task["ts_sampling_json"])["requested"], 20)
        self.assertIsNone(claim_job(self.db, {}))
        report = dict(requested=20, embedded=20, retained=2, fixed_atoms=[0, 4, 1])
        with patch("hpc_sampling.sample_transition_state", return_value=([self.mol, self.mol], report)) as embed:
            run_job(self.db, task, self.root, {})
        self.assertEqual(embed.call_args.kwargs["n"], 20)
        children = self.db.execute("SELECT * FROM jobs WHERE ts_strategy='sampled'").fetchall()
        self.assertEqual(len(children), 2)
        self.assertEqual({row["ts_parent_id"] for row in children}, {best})
        self.assertEqual({row["conformer"] for row in children}, {-2, -3})
        self.assertTrue(all((self.root / row["input_xyz"]).is_file() for row in children))
        self.db.execute("UPDATE jobs SET status='succeeded',gibbs_free_energy_eh=-20,imaginary_count=1,"
                        "ts_mode_status='passed' WHERE ts_strategy='sampled'")
        select_transition_state(self.db, task)
        retained, = self.db.execute("SELECT * FROM jobs WHERE retained=1").fetchall()
        self.assertEqual(retained["ts_strategy"], "sampled")
        self.assertIsNone(claim_job(self.db, {}))
        self.assertEqual(self.db.execute("SELECT count(*) FROM jobs WHERE ts_strategy='sampling'").fetchone()[0], 1)

    def test_all_failed_or_invalid_groups_do_not_sample(self):
        self.add_ts(1, status="failed")
        self.add_ts(2, mode="failed")
        self.assertIsNone(claim_job(self.db, {}))

    def test_neb_success_can_seed_sampling(self):
        self.add_ts(1, status="failed")
        parent = self.add_ts(0, strategy="neb")
        task = claim_job(self.db, {})
        self.assertEqual(task["ts_parent_id"], parent)

    def test_competing_workers_claim_only_one_sampling_task(self):
        self.add_ts(1)

        def claim(_):
            db = connect_database(self.root / "workflow.sqlite")
            try:
                return claim_job(db, {})
            finally:
                db.close()

        with ThreadPoolExecutor(max_workers=2) as pool:
            claims = list(pool.map(claim, range(2)))
        self.assertEqual(sum(task is not None for task in claims), 1)
        self.assertEqual(self.db.execute("SELECT count(*) FROM jobs WHERE ts_strategy='sampling'").fetchone()[0], 1)

    def test_failure_preserves_parent_and_retry_reuses_sampling_protocol(self):
        parent = self.add_ts(1)
        task = claim_job(self.db, {"ts_conformers": 20, "prune_rmsd": 0.3})
        with patch("hpc_sampling.sample_transition_state", side_effect=ValueError("embedding failed")):
            run_job(self.db, task, self.root, {})
        self.assertIsNone(claim_job(self.db, {}))
        self.assertEqual(self.db.execute("SELECT status FROM jobs WHERE id=?", (parent,)).fetchone()[0], "succeeded")
        reset_failed_transition_states(self.root)
        retried = claim_job(self.db, {})
        self.assertEqual(retried["id"], task["id"])
        self.assertEqual(retried["stage"], "ts_sampling")
        self.assertEqual(retried["ts_sampling_json"], task["ts_sampling_json"])

    def test_sampled_ts_uses_preopt_optts_frequency_and_no_scan_recovery(self):
        self.add_ts(1)
        task = claim_job(self.db, {})
        with patch("hpc_sampling.sample_transition_state", return_value=([self.mol], dict(retained=1))):
            run_job(self.db, task, self.root, {})
        child = claim_job(self.db, {})
        optimization = dict(converged=True, xyz=Chem.MolToXYZBlock(self.mol), out_file="opt.out")
        frequency = dict(out_file="freq.out", electronic_energy_eh=-11., gibbs_free_energy_eh=-10.5,
                         imaginary_count=1, frequencies=[-500.])
        with patch("hpc_workflow.run_orca", side_effect=[optimization, optimization, frequency]) as orca, \
                patch("hpc_workflow.check_ts_mode", return_value=dict(valid=True, fingerprint="check")), \
                patch("hpc_workflow.score_ts_scan", side_effect=AssertionError("No rxts scan for sampled TS")):
            run_job(self.db, child, self.root, {})
        self.assertEqual([call.args[3] for call in orca.call_args_list], ["ts_preopt", "optts", "freq"])
        self.assertEqual(self.db.execute("SELECT status FROM jobs WHERE id=?", (child["id"],)).fetchone()[0], "succeeded")

    def test_sampled_failure_keeps_parent_and_does_not_start_scan_or_neb(self):
        parent = self.add_ts(1)
        select_transition_state(self.db, dict(ligand_id="1", pair_index=1, seating=1))
        task = claim_job(self.db, {})
        with patch("hpc_sampling.sample_transition_state", return_value=([self.mol], dict(retained=1))):
            run_job(self.db, task, self.root, {})
        child = claim_job(self.db, {})
        preopt = dict(converged=True, xyz=Chem.MolToXYZBlock(self.mol), out_file="opt.out")
        failed = {**preopt, "converged": False}
        with patch("hpc_workflow.run_orca", side_effect=[preopt, failed, failed]) as orca, \
                patch("hpc_workflow.score_ts_scan", side_effect=AssertionError("No saved scan")):
            run_job(self.db, child, self.root, {})
        self.assertEqual([call.args[3] for call in orca.call_args_list], ["ts_preopt", "optts", "optts"])
        parent_row = self.db.execute("SELECT * FROM jobs WHERE id=?", (parent,)).fetchone()
        self.assertEqual((parent_row["status"], parent_row["retained"]), ("succeeded", 1))


class EndpointSamplingTests(unittest.TestCase):
    def test_pruned_counts_import_without_missing_conformer_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = base / "source" / "ligand_01"
            root = base / "results"
            source.mkdir(parents=True)
            root.mkdir()
            mol = molecule()
            smiles = Chem.MolToSmiles(mol)

            def screen(*args, **kwargs):
                (source / "input.json").write_text(json.dumps(dict(ID="1", substrate="C")))
                with (source / "reaction_pairs.csv").open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["mode", "isomer", "donors", "substrate_site", "resonance"])
                    writer.writeheader()
                    writer.writerow(dict(mode=0, isomer=0, donors="[]", substrate_site=1, resonance=0))
                for role in ("reactant", "product", "ts_guess"):
                    (source / role).mkdir()
                return [{}]

            ensemble = SimpleNamespace(ids=[0, 1], restored_mol=lambda cid: Chem.Mol(mol))
            with patch("hpc_sampling.generate_reaction_pairs", side_effect=screen), \
                    patch("hpc_sampling.create_endpoint_isomers", return_value=[(0, mol, mol)]), \
                    patch("hpc_sampling.rxembed.cxsmiles", return_value=smiles), \
                    patch("hpc_sampling.rxembed.embed", return_value=ensemble) as embed, \
                    patch("hpc_sampling.generate_ts_guess"):
                summary = sample_ligand(dict(ID="1", substrate="C"), source, "unused")
            self.assertEqual(summary["endpoint_conformers"], 20)
            self.assertEqual(summary["errors"], [])
            self.assertTrue(all(call.kwargs["n"] == 20 for call in embed.call_args_list))
            self.assertEqual([record["retained"] for record in summary["endpoint_sampling"]], [1, 1])
            db = connect_database(":memory:")
            self.addCleanup(db.close)
            prepare_jobs(db, source.parent, root, dict(neb_fallback=False))
            endpoints = db.execute("SELECT * FROM jobs WHERE endpoint!='TS'").fetchall()
            self.assertEqual(len(endpoints), 2)
            self.assertEqual({row["status"] for row in endpoints}, {"pending"})


if __name__ == "__main__":
    unittest.main()
