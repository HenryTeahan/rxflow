"""Append scored tridentate N,N,N ligand-3 isomers omitted by the 20 kcal/mol screen."""

import csv
import json
from pathlib import Path
import tomllib

from rdkit import Chem

from hpc_sampling import sample_reaction_pairs
from rxflow import (create_permutations, read_ligands,
                    reaction_pairs_for_coordination_mode)


ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = [
    ("ED_dF_ligands.csv", "ED_dF_hpc_seeds", "hpc_ED_dF.toml"),
    ("ED_pyr_ligands.csv", "ED_pyr_hpc_seeds", "hpc_ED_pyr.toml"),
    ("Ni_ligands.csv", "Ni_hpc_seeds", "hpc_ni.toml"),
    ("Pd_TM_ligands.csv", "hpc_tm_seeds", "hpc.toml"),
]


def main():
    for ligand_csv, seed_name, config_name in SYSTEMS:
        config = tomllib.loads((ROOT / config_name).read_text())
        sampling = config.get("sampling", {})
        ligand = next(row for row in read_ligands(ROOT / ligand_csv) if row["ID"] == "3")
        seed = ROOT / seed_name / "ligand_03"
        modes_path = seed / "coordination_modes.json"
        modes = json.loads(modes_path.read_text())
        ligand_graph = Chem.AddHs(Chem.MolFromSmiles(ligand["ligand_smiles"]))
        nnn_modes = [record for record in modes if len(record["donors"]) == 3 and all(
            ligand_graph.GetAtomWithIdx(index).GetAtomicNum() == 7 for index in record["donors"])]
        if not nnn_modes:
            print(f"{seed_name}: no scored tridentate N,N,N modes")
            continue

        pair_path = seed / "reaction_pairs.csv"
        with pair_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            existing_pairs = list(reader)
        existing = {(int(row["mode"]), tuple(json.loads(row["donors"])), int(row["isomer"]))
                    for row in existing_pairs}
        mode_molecules = {tuple(donors): molecule for donors, molecule in
                          create_permutations(ligand["ligand_smiles"], ligand["metal"], ligand["denticity"])}
        additions = []
        for record in nnn_modes:
            record["retained"] = True
            identity = (record["mode"], tuple(record["donors"]), record["isomer"])
            if identity in existing:
                continue
            additions.extend(reaction_pairs_for_coordination_mode(
                ligand, record["mode"], record["donors"], record["isomer"],
                mode_molecules[tuple(record["donors"])], record["coordination_smiles"]))

        if additions:
            with pair_path.open("a", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writerows(additions)
            summary = sample_reaction_pairs(
                ligand, additions, seed, first_pair_index=len(existing_pairs) + 1,
                endpoint_conformers=sampling.get("endpoint_conformers", 20),
                coembed_conformers=sampling.get("coembed_conformers", 10),
                ts_candidates=sampling.get("ts_candidates", 5),
                ts_conformers=sampling.get("ts_conformers", 20),
                prune_rmsd=sampling.get("prune_rmsd", 0.5),
                seed=sampling.get("seed", 12345))
            sampling_path = seed / "sampling.json"
            previous = json.loads(sampling_path.read_text()) if sampling_path.exists() else {}
            previous["ligand_id"] = ligand["ID"]
            previous["groups"] = previous.get("groups", 0) + summary["groups"]
            previous["errors"] = previous.get("errors", []) + summary["errors"]
            previous["endpoint_sampling"] = previous.get("endpoint_sampling", []) + summary["endpoint_sampling"]
            previous.update({key: value for key, value in summary.items()
                             if key not in ("groups", "errors", "endpoint_sampling")})
            previous["status"] = "completed_with_errors" if previous["errors"] else "complete"
            sampling_path.write_text(json.dumps(previous, indent=2) + "\n")
        modes_path.write_text(json.dumps(modes, indent=2) + "\n")
        print(f"{seed_name}: NNN isomers={len(nnn_modes)}, appended pairs={len(additions)}")


if __name__ == "__main__":
    main()
