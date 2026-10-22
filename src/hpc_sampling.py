"""Create HPC calculation seeds directly from one CSV ligand row."""

import json
from pathlib import Path

from rdkit import Chem
import rxembed

from conformer_sampling import DEFAULT_PRUNE_RMSD, embedded_conformers
from rxflow import generate_reaction_pairs
from rxflow_endpoints import create_endpoint_isomers, generate_ts_guess


def sample_ligand(row, output, executable, *, endpoint_conformers=20, coembed_conformers=10,
                  ts_candidates=5, ts_conformers=20, prune_rmsd=DEFAULT_PRUNE_RMSD, seed=12345):
    """Screen binding modes, embed endpoints, and coembed TS guesses without endpoint optimization."""
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Sampling already exists in {output}; choose a fresh source directory")
    pairs = generate_reaction_pairs(row, output, executable=executable, seed=seed, prune_rmsd=prune_rmsd)
    summary = sample_reaction_pairs(
        row, pairs, output, first_pair_index=1, endpoint_conformers=endpoint_conformers,
        coembed_conformers=coembed_conformers, ts_candidates=ts_candidates,
        ts_conformers=ts_conformers, prune_rmsd=prune_rmsd, seed=seed)
    summary["ligand_id"] = row["ID"]
    (output / "sampling.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def sample_reaction_pairs(row, pairs, output, *, first_pair_index, endpoint_conformers=20,
                          coembed_conformers=10, ts_candidates=5, ts_conformers=20,
                          prune_rmsd=DEFAULT_PRUNE_RMSD, seed=12345):
    """Sample endpoint and TS seeds for appended reaction-pair rows in an existing seed tree."""
    output = Path(output)
    errors = []
    endpoint_sampling = []
    groups = 0
    for index, pair in enumerate(pairs, first_pair_index):
        for seating, reactant, product in create_endpoint_isomers(pair):
            prefix = f"pair_{index:03d}_seating_{seating + 1:02d}"
            descriptor = output / f"{prefix}_endpoints.json"
            if descriptor.exists():
                raise FileExistsError(f"Endpoint seeds already exist for pair {index}, seating {seating + 1}")
            descriptor.write_text(json.dumps({
                "reactant": rxembed.cxsmiles(reactant), "product": rxembed.cxsmiles(product),
                "substrate": row["substrate"]}, indent=2) + "\n")
            for role, isomer in (("reactant", reactant), ("product", product)):
                record = dict(pair=index, seating=seating + 1, endpoint=role,
                              requested=endpoint_conformers, embedded=0, retained=0)
                endpoint_sampling.append(record)
                try:
                    ensemble = rxembed.embed(isomer, n=endpoint_conformers, seed=seed, threads=1)
                    ensembles = [ensemble] if hasattr(ensemble, "restored_mol") else list(ensemble)
                    record["embedded"] = sum(len(item.ids) for item in ensembles)
                    conformers = [mol for batch in embedded_conformers(ensemble, prune_rmsd) for mol in batch]
                    record["retained"] = len(conformers)
                    if not conformers:
                        raise ValueError(f"No {role} conformers embedded")
                    for conformer, mol in enumerate(conformers, 1):
                        path = output / role / f"{prefix}_conformer_{conformer:03d}_input.xyz"
                        path.write_text(Chem.MolToXYZBlock(mol))
                except Exception as exc:
                    errors.append({"pair": index, "seating": seating + 1, "endpoint": role, "error": str(exc)})
            try:
                generate_ts_guess(reactant, product, output, prefix=prefix, n=coembed_conformers,
                                  keep=ts_candidates, seed=seed)
                groups += 1
            except Exception as exc:
                errors.append({"pair": index, "seating": seating + 1, "endpoint": "TS", "error": str(exc)})
    summary = dict(groups=groups, errors=errors,
                   endpoint_conformers=endpoint_conformers, coembed_conformers=coembed_conformers,
                   ts_candidates=ts_candidates, ts_conformers=ts_conformers, prune_rmsd=prune_rmsd,
                   endpoint_sampling=endpoint_sampling, seed=seed)
    if not groups and not errors:
        errors.append({"error": "No reaction/seating groups produced"})
    summary["status"] = "completed_with_errors" if errors else "complete"
    return summary


def sample_transition_state(job, xyz, *, n=20, seed=12345, max_rmsd=DEFAULT_PRUNE_RMSD):
    """Reembed the optimized TS with the three mapped transfer atoms fixed."""
    import numpy as np
    from orca_hpc import molecule_with_xyz

    molecule = molecule_with_xyz(job["smiles"], xyz)
    fixed = [job.get(name) for name in ("donor", "hydrogen", "acceptor")]
    if (any(type(i) is not int or not 0 <= i < molecule.GetNumAtoms() for i in fixed)
            or len(set(fixed)) != 3):
        raise ValueError("TS sampling requires distinct donor, hydrogen and acceptor indices")
    embedded = rxembed.embed(molecule, fix=fixed, n=n, seed=seed, threads=1, charge=job["charge"])
    batches = embedded_conformers(embedded, max_rmsd)
    conformers = [mol for batch in batches for mol in batch]
    expected = [(a.GetAtomicNum(), a.GetAtomMapNum()) for a in molecule.GetAtoms()]
    reference = molecule.GetConformer().GetPositions()[fixed]
    for mol in conformers:
        if [(a.GetAtomicNum(), a.GetAtomMapNum()) for a in mol.GetAtoms()] != expected:
            raise ValueError("TS embedding changed atom inventory/map order")
        if Chem.GetFormalCharge(mol) != job["charge"]:
            raise ValueError("TS embedding changed the molecular charge")
        # rxembed preserves the fixed core's geometry in each seed's frame.
        # Rigidly fit the entire conformer back to the parent's Cartesian frame.
        positions = mol.GetConformer().GetPositions()
        moving = positions[fixed]
        left, _, right = np.linalg.svd((moving - moving.mean(0)).T @ (reference - reference.mean(0)))
        if np.linalg.det(left @ right) < 0:
            left[:, -1] *= -1  # Rotate without reflecting the stereochemistry.
        fitted = (positions - moving.mean(0)) @ (left @ right) + reference.mean(0)
        for index, position in enumerate(fitted):
            mol.GetConformer().SetAtomPosition(index, position)
        if not np.allclose(mol.GetConformer().GetPositions()[fixed], reference, atol=1e-6, rtol=0):
            raise ValueError("TS embedding moved a fixed transfer atom")
    ensembles = [embedded] if hasattr(embedded, "restored_mol") else list(embedded)
    report = dict(requested=n, embedded=sum(len(ensemble.ids) for ensemble in ensembles),
                  retained=len(conformers), seed=seed, fixed_atoms=fixed,
                  index_base=0, prune_rmsd=max_rmsd)
    return conformers, report
