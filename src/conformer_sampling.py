"""RMSD-only pruning of embedded geometries, without changing coordinates."""

import math

import numpy as np
from rdkit import Chem


DEFAULT_PRUNE_RMSD = 0.5


def prune_conformers(geometries, max_rmsd=DEFAULT_PRUNE_RMSD):
    """Use prism's heavy-atom RMSD mask, separately for each molecular identity."""
    from prism_pruner.pruner import prune_by_rmsd

    if not math.isfinite(max_rmsd) or max_rmsd <= 0:
        raise ValueError("Pruning RMSD must be finite and positive")
    geometries = list(geometries)
    groups = {}
    for index, mol in enumerate(geometries):
        # Keep defined stereo/coordination graphs separate, including atom-map order.
        identity = (Chem.MolToSmiles(mol), tuple(a.GetAtomMapNum() for a in mol.GetAtoms()))
        groups.setdefault(identity, []).append(index)
    kept = set()
    for indices in groups.values():
        if len(indices) == 1:
            kept.update(indices)
            continue
        structures = np.array([geometries[i].GetConformer().GetPositions() for i in indices])
        atoms = np.array([a.GetSymbol() for a in geometries[indices[0]].GetAtoms()])
        # Disable prism's extra maximum-atom-displacement criterion. No energies
        # are supplied: only aligned RMSD decides which conformers are duplicates.
        _, mask = prune_by_rmsd(structures, atoms, max_rmsd=max_rmsd, max_dev=float("inf"))
        kept.update(i for i, retain in zip(indices, mask, strict=True) if retain)
    return [mol for i, mol in enumerate(geometries) if i in kept]


def embedded_conformers(embedded, max_rmsd=DEFAULT_PRUNE_RMSD):
    """Prune each rxembed stereo/coordination ensemble independently."""
    ensembles = [embedded] if hasattr(embedded, "restored_mol") else list(embedded)
    return [prune_conformers([ensemble.restored_mol(cid) for cid in ensemble.ids], max_rmsd)
            for ensemble in ensembles]
