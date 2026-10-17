"""Create HPC calculation seeds directly from one CSV ligand row."""

import json
from pathlib import Path

from rdkit import Chem
import rxembed

from rxflow import generate_reaction_pairs
from rxflow_endpoints import create_endpoint_isomers, generate_ts_guess


def sample_ligand(row, output, executable, *, endpoint_conformers=3, coembed_conformers=10,
                  ts_candidates=5, seed=12345):
    """Screen binding modes, embed endpoints, and coembed TS guesses without endpoint optimization."""
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Sampling already exists in {output}; choose a fresh source directory")
    pairs = generate_reaction_pairs(row, output, executable=executable, seed=seed)
    errors = []
    groups = 0
    for index, pair in enumerate(pairs, 1):
        for seating, reactant, product in create_endpoint_isomers(pair):
            prefix = f"pair_{index:03d}_seating_{seating + 1:02d}"
            (output / f"{prefix}_endpoints.json").write_text(json.dumps({
                "reactant": rxembed.cxsmiles(reactant), "product": rxembed.cxsmiles(product),
                "substrate": row["substrate"]}, indent=2) + "\n")
            for role, isomer in (("reactant", reactant), ("product", product)):
                try:
                    ensemble = rxembed.embed(isomer, n=endpoint_conformers, seed=seed, threads=1)
                    if not ensemble.ids:
                        raise ValueError(f"No {role} conformers embedded")
                    for conformer, cid in enumerate(ensemble.ids, 1):
                        path = output / role / f"{prefix}_conformer_{conformer:03d}_input.xyz"
                        path.write_text(Chem.MolToXYZBlock(ensemble.restored_mol(cid)))
                    if len(ensemble.ids) < endpoint_conformers:
                        raise ValueError(f"Only {len(ensemble.ids)} of {endpoint_conformers} {role} conformers embedded")
                except Exception as exc:
                    errors.append({"pair": index, "seating": seating + 1, "endpoint": role, "error": str(exc)})
            try:
                generate_ts_guess(reactant, product, output, prefix=prefix, n=coembed_conformers,
                                  keep=ts_candidates, seed=seed)
                groups += 1
            except Exception as exc:
                errors.append({"pair": index, "seating": seating + 1, "endpoint": "TS", "error": str(exc)})
    summary = dict(ligand_id=row["ID"], groups=groups, errors=errors,
                   endpoint_conformers=endpoint_conformers, coembed_conformers=coembed_conformers,
                   ts_candidates=ts_candidates, seed=seed)
    if not groups and not errors:
        errors.append({"error": "No reaction/seating groups produced"})
    summary["status"] = "completed_with_errors" if errors else "complete"
    (output / "sampling.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
