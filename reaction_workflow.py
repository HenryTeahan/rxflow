#!/usr/bin/env python3
"""Enumerate metal-mediated H-transfer candidates from deprotonated ligand SMILES.

The stages follow reaction_smarts.ipynb: three-donor metal binding modes, square-
planar embedding, g-xTB single-point filtering, substrate coordination, reaction
SMARTS over resonance forms, reactant/product conformer ranking, and one
VeloxChem MM TS guess per retained ligand binding isomer.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import subprocess
import tempfile
from contextlib import nullcontext
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata as package_metadata
from itertools import combinations
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdDistGeom
from rdkit.Geometry import Point3D

import rxembed
from conformer_sampling import prune_conformers
from workflow_store import WorkflowStore, digest, import_ligand


LOG = logging.getLogger("reaction_workflow")
HARTREE_TO_KCAL = 627.5094740631
DEFAULT_SUBSTRATE = "[H]C1=CC(OC)=C(OC)C(OC)=C1"
DEFAULT_METAL_SMILES = "[Pd+2]"
ENERGY_RE = re.compile(
    r"TOTAL ENERGY\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)\s+Eh"
)
ORIGINAL_INDEX = "_workflow_atom_index"

CSV_COLUMNS = [
    "ligand_id", "ligand_smiles", "metal_smiles", "donor_indices", "ligand_isomer",
    "ligand_conformer", "ligand_energy_eh", "substrate_atom_index",
    "reactant_smiles", "reactant_isomer", "reactant_conformer",
    "reactant_energy_eh", "product_smiles", "resonance_form",
    "acceptor_atom_index", "acceptor_element", "acceptor_metal_bound",
    "product_isomer", "product_conformer", "product_energy_eh",
    "reactant_xyz", "product_xyz", "ts_guess_xyz",
    "ts_reactant_input_xyz", "ts_product_input_xyz",
]
LIGAND_COLUMNS = [
    "donor_indices", "ligand_isomer", "ligand_conformer", "energy_eh",
    "relative_kcal_mol", "formal_charge", "retained",
]
FAILURE_COLUMNS = ["ligand_id", "stage", "candidate", "error"]


@dataclass
class Geometry:
    mol: Chem.Mol
    energy_eh: float
    isomer: int
    conformer: int
    xyz_path: Path | None = None
    database_id: int | None = None


class SamplingRecorder:
    """Attach sampling operations to a chemical state and its parent geometry."""

    def __init__(self, args, key, role, graph, donors, parent_geometry=None, **metadata):
        self.args = args
        self.store = args._store
        self.ligand_id = args._run_ligand_id
        self.key = f"{args._sampling_scope}:{key}"
        self.parent_geometry = parent_geometry
        self.donors = donors
        parent_state = None
        if parent_geometry is not None:
            row = self.store.connection.execute(
                "SELECT i.state_id FROM geometry g JOIN isomer i ON i.id = g.isomer_id WHERE g.id = ?",
                (parent_geometry,),
            ).fetchone()
            parent_state = row[0] if row else None
        graph_copy = Chem.Mol(graph)
        graph_copy.RemoveAllConformers()
        self.state_id = self.store.state(
            self.ligand_id, self.key, role, graph_smiles=smiles(graph),
            graph_block=Chem.MolToMolBlock(graph_copy, forceV3000=True),
            formal_charge=Chem.GetFormalCharge(graph), metal_smiles=args.metal_smiles,
            donors=donors, parent_state_id=parent_state, atom_index_base=0,
            atom_count=graph.GetNumAtoms(), **metadata,
        )

    def operation(self, kind, suffix, parameters=None, **targets):
        return self.store.operation(self.ligand_id, f"{self.key}:{kind}:{suffix}", kind,
                                    parameters=parameters, state_id=self.state_id, **targets)

    def record_geometry(self, mol, isomer_index, conformer_index, seed,
                        source="rxembed.embed", generated_by=None):
        local_index = None if source == "rxembed.minimize" else isomer_index
        isomer_id = self.store.isomer(self.state_id, local_index,
                                      "local_relaxation" if local_index is None else "rxembed.metal",
                                      None if local_index is None else "SPL")
        key = f"{self.key}:{source}:{isomer_index}:{conformer_index}"
        return self.store.geometry(
            self.ligand_id, key, source, isomer_id=isomer_id, conformer_index=conformer_index,
            seed=seed, parent_geometry_id=self.parent_geometry,
            xyz=Chem.MolToXYZBlock(mol), attempt_id=generated_by,
        )

    def score(self, mol, isomer_index, conformer_index, seed, source="rxembed.embed",
              generated_by=None, geometry_id=None, scratch=None):
        if geometry_id is None:
            geometry_id = self.record_geometry(mol, isomer_index, conformer_index, seed, source, generated_by)
        key = f"{self.key}:{source}:{isomer_index}:{conformer_index}"
        parameters = {"method": "g-xTB", "charge": Chem.GetFormalCharge(mol),
                      "threads": self.args.xtb_threads, "executable": str(self.args.xtb.resolve()),
                      "spin_setting": "engine_default", "solvent": None}
        with self.operation("gxtb_sp", key, parameters, geometry_id=geometry_id) as attempt:
            calculation = self.store.calculation(
                attempt, geometry_id, "g-xTB", charge=parameters["charge"], settings=parameters,
                engine_version=getattr(self.args, "_xtb_version", None),
            )
            energy = gxtb_energy(mol, self.args.xtb, scratch or self.args.output_dir / "scratch", self.args.xtb_threads,
                                 result_callback=lambda result: self.store.logs(
                                     attempt, result.stdout, result.stderr, result.returncode))
            self.store.energy(calculation, energy)
        return Geometry(mol, energy, isomer_index, conformer_index, database_id=geometry_id)


def sampling_recorder(args, key, role, graph, donors, parent_geometry=None, **metadata):
    if getattr(args, "_store", None) is None:
        return None
    return SamplingRecorder(args, key, role, graph, donors, parent_geometry, **metadata)


def record_cutoff(args, candidates, retained, cutoff):
    if getattr(args, "_store", None) is None or not candidates:
        return
    minimum = min(candidates, key=lambda item: item.geometry.energy_eh)
    retained_ids = {item.geometry.database_id for item in retained}
    for item in candidates:
        args._store.select_geometry(
            item.geometry.database_id,
            "retained" if item.geometry.database_id in retained_ids else "rejected_by_cutoff",
            policy="relative_electronic_energy", cutoff_kcal_mol=cutoff,
            reference_geometry_id=minimum.geometry.database_id,
            relative_kcal_mol=(item.geometry.energy_eh - minimum.geometry.energy_eh) * HARTREE_TO_KCAL,
        )


@dataclass
class LigandMode:
    key: str
    donors: tuple[int, ...]
    ligand_smiles: str
    geometry: Geometry


@dataclass
class ReactantMode:
    key: str
    ligand: LigandMode
    substrate_atom_index: int
    graph: Chem.Mol
    geometry: Geometry


@dataclass
class ProductMode:
    reactant: ReactantMode
    smiles: str
    resonance_form: int
    acceptor_atom_index: int
    acceptor_element: str
    acceptor_metal_bound: bool
    geometry: Geometry


@dataclass(frozen=True)
class MetalSpec:
    smiles: str
    atomic_number: int
    symbol: str


def parse_metal_smiles(value: str) -> MetalSpec:
    mol = Chem.MolFromSmiles(value)
    if mol is None or mol.GetNumAtoms() != 1 or mol.GetNumBonds() != 0:
        raise ValueError("--metal-smiles must be a single bracketed atom, e.g. [Ni+2]")
    atom = mol.GetAtomWithIdx(0)
    canonical = Chem.MolToSmiles(mol)
    if (not canonical.startswith("[") or not canonical.endswith("]")
        or atom.GetAtomicNum() == 0 or atom.GetAtomMapNum()
        or atom.GetIsotope() or atom.GetTotalNumHs()):
        raise ValueError("--metal-smiles must be an unmapped bracketed atom without H or isotope")
    return MetalSpec(canonical, atom.GetAtomicNum(), atom.GetSymbol())


def reaction_smarts(metal: MetalSpec) -> str:
    mapped_metal = f"{metal.smiles[:-1]}:3]"
    return (
        f"([c:1](-[#1:4])->{mapped_metal}.[#7,#8,#16;-1:2])"
        ">>"
        f"([c-:1]->{mapped_metal}.[*+0:2]-[#1:4])"
    )


def install_rdkit_rxembed_compatibility() -> None:
    """Bridge rxembed's older bounds arguments to RDKit 2026's API.

    The installed rxembed calls three RDKit options removed in RDKit 2026.
    This process-local adapter lets its current embedding path run. The old
    set14bounds=False distinction is unavailable in this RDKit build.
    """
    old = rdDistGeom.GetMoleculeBoundsMatrix
    try:
        old(Chem.MolFromSmiles("CC"), set14bounds=False)
        return
    except (TypeError, ValueError):
        pass

    def bounds_compat(mol, *args, **kwargs):
        kwargs.pop("set14bounds", None)
        kwargs.pop("embedParams", None)
        return old(mol, *args, **kwargs)

    rdDistGeom.GetMoleculeBoundsMatrix = bounds_compat
    if not hasattr(rdDistGeom.KDG(), "useLegacyImplementation"):
        rdDistGeom.EmbedParameters.useLegacyImplementation = property(
            lambda self: False, lambda self, value: None
        )
    LOG.warning("Using RDKit 2026 compatibility adapter for installed rxembed")


def tag_atom_order(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetIntProp(ORIGINAL_INDEX, atom.GetIdx())


def restore_atom_order(mol: Chem.Mol, reference: Chem.Mol) -> Chem.Mol:
    """Use persistent atom properties to undo any embedding reordering."""
    count = reference.GetNumAtoms()
    if mol.GetNumAtoms() != count:
        raise ValueError(f"Atom count changed during embedding: {count} -> {mol.GetNumAtoms()}")
    by_original = {}
    for atom in mol.GetAtoms():
        if not atom.HasProp(ORIGINAL_INDEX):
            raise ValueError("Embedding removed atom identity properties")
        original = atom.GetIntProp(ORIGINAL_INDEX)
        if original in by_original:
            raise ValueError(f"Duplicate atom identity {original} after embedding")
        by_original[original] = atom.GetIdx()
    if set(by_original) != set(range(count)):
        raise ValueError("Embedding changed atom identities")
    ordered = Chem.RenumberAtoms(mol, [by_original[i] for i in range(count)])
    if any(
        ordered.GetAtomWithIdx(i).GetAtomicNum() != reference.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(count)
    ):
        raise ValueError("Atom elements changed during embedding")
    return ordered


def one_conformer(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    single = Chem.Mol(mol)
    conformer = Chem.Conformer(single.GetConformer(conf_id))
    single.RemoveAllConformers()
    single.AddConformer(conformer, assignId=True)
    return single


def smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), canonical=True)


def gxtb_energy(mol: Chem.Mol, xtb: Path, scratch: Path, threads: int, *, result_callback=None) -> float:
    """Score one geometry with the same g-xTB single point used in the notebook."""
    charge = Chem.GetFormalCharge(mol)
    with tempfile.TemporaryDirectory(prefix="gxtb-", dir=scratch) as work:
        work_path = Path(work)
        xyz = work_path / "candidate.xyz"
        xyz.write_text(Chem.MolToXYZBlock(mol))
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)
        run = subprocess.run(
            [str(xtb), xyz.name, "--gxtb", "--charge", str(charge), "--sp"],
            cwd=work_path, env=env, capture_output=True, text=True, check=False,
        )
        if result_callback is not None:
            result_callback(run)
        if run.returncode:
            raise RuntimeError(f"g-xTB exit {run.returncode}: {(run.stderr or run.stdout)[-600:]}")
        matches = ENERGY_RE.findall(run.stdout)
        if not matches:
            raise RuntimeError("g-xTB output has no TOTAL ENERGY in Eh")
        energy = float(matches[-1].replace("D", "E").replace("d", "e"))
        if not math.isfinite(energy):
            raise RuntimeError("g-xTB returned a non-finite energy")
        return energy


def within_cutoff(items: list, cutoff_kcal: float) -> list:
    if not items:
        return []
    minimum = min(item.geometry.energy_eh for item in items)
    return [
        item for item in items
        if (item.geometry.energy_eh - minimum) * HARTREE_TO_KCAL <= cutoff_kcal + 1e-8
    ]


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    """Replace a CSV only after its complete contents have been written."""
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", prefix=f".{path.name}.", suffix=".tmp",
        dir=path.parent, delete=False,
    ) as handle:
        temporary = Path(handle.name)
        try:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_results_csv(path: Path) -> list[dict]:
    rows = read_csv(path)
    for row in rows:
        if "acceptor_pd_bound" in row:
            row["acceptor_metal_bound"] = row.pop("acceptor_pd_bound")
            row.setdefault("metal_smiles", DEFAULT_METAL_SMILES)
    return rows


def run_options(args) -> dict:
    return {
        "smiles_column": args.smiles_column,
        "id_column": args.id_column,
        "substrate": args.substrate,
        "metal_smiles": args.metal_smiles,
        "xtb": str(args.xtb.resolve()),
        "conformers": args.conformers,
        "cutoff_kcal": args.cutoff_kcal,
        "n_donors": args.n_donors,
        "seed": args.seed,
        "xtb_threads": args.xtb_threads,
        "ts_mm_steps": args.ts_mm_steps,
        "ts_lambda_points": args.ts_lambda_points,
        "ts_snapshots": args.ts_snapshots,
        "max_pairs_per_ligand": args.max_pairs_per_ligand,
    }


def ligand_directory(args, row: dict, row_number: int) -> Path:
    ligand_id = str(row.get(args.id_column) or row_number)
    return args.output_dir / f"ligand_{row_number:04d}_{safe_name(ligand_id)}"


def load_checkpoint(path: Path, row: dict, row_number: int, args):
    state_path = path / "status.json"
    if not state_path.is_file():
        return None
    state = json.loads(state_path.read_text())
    ligand_id = str(row.get(args.id_column) or row_number)
    saved_options = dict(state.get("options") or {})
    saved_options.setdefault("metal_smiles", DEFAULT_METAL_SMILES)
    if (state.get("status") != "complete"
        or state.get("row_number") != row_number
        or state.get("ligand_id") != ligand_id
        or state.get("ligand_smiles") != str(row.get(args.smiles_column) or "").strip()
        or saved_options != run_options(args)):
        raise ValueError(f"Checkpoint does not match this input or run configuration: {state_path}")
    results_path = path / "ligand_results.csv"
    failures_path = path / "ligand_failures.csv"
    if not results_path.is_file() or not failures_path.is_file():
        raise ValueError(f"Completed checkpoint is missing its CSV files: {path}")
    return read_results_csv(results_path), read_csv(failures_path)


def save_checkpoint(path: Path, row: dict, row_number: int, args,
                    results: list[dict], failures: list[dict]) -> None:
    write_csv(path / "ligand_results.csv", CSV_COLUMNS, results)
    write_csv(path / "ligand_failures.csv", FAILURE_COLUMNS, failures)
    state = {
        "status": "complete",
        "row_number": row_number,
        "ligand_id": str(row.get(args.id_column) or row_number),
        "ligand_smiles": str(row.get(args.smiles_column) or "").strip(),
        "options": run_options(args),
    }
    if getattr(args, "_store", None):
        state["database"] = str(args._store.path)
        state["database_run_id"] = args._run_id
        state["database_run_ligand_id"] = args._run_ligand_id
    with tempfile.NamedTemporaryFile(
        mode="w", prefix=".status.", suffix=".tmp", dir=path, delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(state, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, path / "status.json")


def embedded_geometries(
    graph: Chem.Mol,
    *,
    n_conformers: int,
    seed: int,
    xtb: Path,
    scratch: Path,
    xtb_threads: int,
    template: tuple[Chem.Mol, dict[int, int]] | None = None,
    recorder: SamplingRecorder | None = None,
):
    """Yield every successfully scored conformer of every square-planar isomer."""
    tag_atom_order(graph)
    with recorder.operation("enumerate_isomers", "all", {"geometry": "SPL"}) if recorder else nullcontext():
        isomers = list(rxembed.metal(graph, geometry="SPL"))
        if recorder:
            for index in range(len(isomers)):
                recorder.store.isomer(recorder.state_id, index)
    for isomer_index, isomer in enumerate(isomers):
        try:
            isomer_id = recorder.store.isomer(recorder.state_id, isomer_index) if recorder else None
            parameters = {"n_conformers": n_conformers, "seed": seed + isomer_index,
                          "template_parent_geometry_id": recorder.parent_geometry if recorder else None,
                          "template_atom_map": template[1] if template else None}
            with recorder.operation("embed", str(isomer_index), parameters, isomer_id=isomer_id) if recorder else nullcontext() as attempt:
                ensemble = rxembed.embed(
                    isomer, n=n_conformers, seed=seed + isomer_index,
                    template=template,
                )
                result = restore_atom_order(ensemble.mol, graph)
                if result.GetNumConformers() == 0:
                    raise ValueError("Embedding returned no conformers")
                conformers = prune_conformers([one_conformer(result, conf.GetId())
                                              for conf in result.GetConformers()])
                generated = []
                for conformer_index, candidate in enumerate(conformers):
                    geometry_id = recorder.record_geometry(candidate, isomer_index, conformer_index,
                                                            seed + isomer_index, generated_by=attempt) if recorder else None
                    generated.append((conformer_index, candidate, geometry_id))
        except (RuntimeError, ValueError, TimeoutError) as exc:
            LOG.warning("Embedding isomer %s failed: %s", isomer_index, exc)
            continue
        for conformer_index, candidate, geometry_id in generated:
            try:
                if recorder:
                    scored = recorder.score(candidate, isomer_index, conformer_index,
                                             seed + isomer_index, geometry_id=geometry_id, scratch=scratch)
                else:
                    scored = Geometry(candidate, gxtb_energy(candidate, xtb, scratch, xtb_threads),
                                      isomer_index, conformer_index)
            except RuntimeError as exc:
                LOG.warning("g-xTB isomer %s conformer %s: %s", isomer_index, conformer_index, exc)
                continue
            yield scored


def ligand_complexes(ligand_smiles: str, n_donors: int, metal: MetalSpec):
    ligand = Chem.MolFromSmiles(ligand_smiles)
    if ligand is None:
        raise ValueError(f"Invalid ligand SMILES: {ligand_smiles!r}")
    ligand = Chem.AddHs(ligand)
    donor_indices = [
        atom.GetIdx() for atom in ligand.GetAtoms()
        if atom.GetAtomicNum() not in (1, 6)
    ]
    metal_mol = Chem.MolFromSmiles(metal.smiles)
    for donors in combinations(donor_indices, n_donors):
        editable = Chem.RWMol(Chem.CombineMols(ligand, metal_mol))
        metal_index = ligand.GetNumAtoms()
        for donor_index in donors:
            editable.AddBond(donor_index, metal_index, Chem.BondType.DATIVE)
        complex_mol = editable.GetMol()
        complex_mol.UpdatePropertyCache(strict=False)
        yield donors, complex_mol


def substrate_sites(substrate: Chem.Mol) -> list[int]:
    """One representative per symmetry-distinct aromatic C-H site."""
    ranks = Chem.CanonicalRankAtoms(substrate, breakTies=False)
    unique = {}
    for atom in substrate.GetAtoms():
        if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
            continue
        if any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()):
            unique.setdefault(ranks[atom.GetIdx()], atom.GetIdx())
    return list(unique.values())


def build_reactant(ligand_geometry: Chem.Mol, substrate: Chem.Mol, site: int,
                   metal: MetalSpec) -> Chem.Mol:
    ligand_count = ligand_geometry.GetNumAtoms()
    combined = Chem.CombineMols(ligand_geometry, substrate)
    editable = Chem.RWMol(combined)
    metal_indices = [a.GetIdx() for a in ligand_geometry.GetAtoms()
                     if a.GetAtomicNum() == metal.atomic_number]
    if len(metal_indices) != 1:
        raise ValueError(f"Expected one {metal.symbol}, found {metal_indices}")
    editable.AddBond(ligand_count + site, metal_indices[0], Chem.BondType.DATIVE)
    reactant = editable.GetMol()
    reactant.RemoveAllConformers()
    Chem.SanitizeMol(reactant)
    return reactant


def transfer_site(reactant: Chem.Mol, metal: MetalSpec) -> tuple[int, int]:
    sites = []
    for carbon in reactant.GetAtoms():
        if carbon.GetAtomicNum() != 6:
            continue
        donates_to_metal = any(
            bond.GetBondType() == Chem.BondType.DATIVE
            and bond.GetBeginAtomIdx() == carbon.GetIdx()
            and bond.GetEndAtom().GetAtomicNum() == metal.atomic_number
            for bond in carbon.GetBonds()
        )
        if donates_to_metal:
            sites.extend(
                (neighbor.GetIdx(), carbon.GetIdx())
                for neighbor in carbon.GetNeighbors()
                if neighbor.GetAtomicNum() == 1
            )
    if len(sites) != 1:
        raise ValueError(f"Expected one {metal.symbol}-bound substrate C-H, found {sites}")
    return sites[0]


def products_from_resonance(reactant: Chem.Mol, metal: MetalSpec):
    """Keep distinct protonation products and preserve the reactant atom order."""
    rxn = AllChem.ReactionFromSmarts(reaction_smarts(metal))
    hydrogen_index, carbon_index = transfer_site(reactant, metal)
    seen = set()
    for form_index, form in enumerate(Chem.ResonanceMolSupplier(reactant)):
        for (raw_product,) in rxn.RunReactants((form,)):
            try:
                source_indices = [
                    atom.GetIntProp("react_atom_idx") for atom in raw_product.GetAtoms()
                ]
                if sorted(source_indices) != list(range(reactant.GetNumAtoms())):
                    raise ValueError("Reaction did not preserve every reactant atom")
                acceptor = next(
                    atom for atom in raw_product.GetAtoms()
                    if atom.HasProp("old_mapno") and atom.GetIntProp("old_mapno") == 2
                )
                acceptor_index = acceptor.GetIntProp("react_atom_idx")
                product = Chem.RenumberAtoms(
                    raw_product,
                    [source_indices.index(i) for i in range(reactant.GetNumAtoms())],
                )
                Chem.SanitizeMol(product)
                if product.GetNumAtoms() != reactant.GetNumAtoms():
                    raise ValueError("Product atom count changed")
                if any(
                    product.GetAtomWithIdx(i).GetAtomicNum()
                    != reactant.GetAtomWithIdx(i).GetAtomicNum()
                    for i in range(reactant.GetNumAtoms())
                ):
                    raise ValueError("Product atom elements changed")
                h_neighbors = [a.GetIdx() for a in product.GetAtomWithIdx(hydrogen_index).GetNeighbors()]
                if h_neighbors != [acceptor_index]:
                    raise ValueError(f"Transferring H is bonded to {h_neighbors}, expected {acceptor_index}")
                if carbon_index in h_neighbors:
                    raise ValueError("Substrate C-H bond was retained")
                product_smiles = smiles(product)
            except (ValueError, StopIteration, Chem.MolSanitizeException) as exc:
                LOG.warning("Rejected reaction product from resonance %s: %s", form_index, exc)
                continue
            if product_smiles in seen:
                continue
            seen.add(product_smiles)
            acceptor_atom = form.GetAtomWithIdx(acceptor_index)
            metal_bound = any(
                bond.GetBondType() == Chem.BondType.DATIVE
                and bond.GetBeginAtomIdx() == acceptor_index
                and bond.GetEndAtom().GetAtomicNum() == metal.atomic_number
                for bond in acceptor_atom.GetBonds()
            )
            product.RemoveAllConformers()
            yield (product, product_smiles, form_index, acceptor_index,
                   acceptor_atom.GetSymbol(), metal_bound)


def product_from_reactant(
    reactant: ReactantMode, product_graph: Chem.Mol, acceptor_index: int,
    *, n_conformers: int, seed: int,
    xtb: Path, scratch: Path, xtb_threads: int, metal: MetalSpec,
    recorder: SamplingRecorder | None = None,
) -> Geometry:
    """Return the lowest-energy product sampled from this reactant conformer."""
    product = Chem.Mol(product_graph)
    source = reactant.geometry.mol
    if product.GetNumAtoms() != source.GetNumAtoms():
        raise ValueError("Reactant/product atom counts differ")
    product.RemoveAllConformers()
    product.AddConformer(Chem.Conformer(source.GetConformer()), assignId=True)
    tag_atom_order(product)
    hydrogen_index, _ = transfer_site(source, metal)
    conf = product.GetConformer()
    acceptor_position = conf.GetAtomPosition(acceptor_index)
    hydrogen_position = conf.GetAtomPosition(hydrogen_index)
    direction = hydrogen_position - acceptor_position
    length = direction.Length()
    if length < 1e-8:
        raise ValueError("Transferring H overlaps the proton acceptor")
    direction *= 1.02 / length
    conf.SetAtomPosition(
        hydrogen_index,
        Point3D(
            acceptor_position.x + direction.x,
            acceptor_position.y + direction.y,
            acceptor_position.z + direction.z,
        ),
    )
    candidates = []
    try:
        with recorder.operation("local_relaxation", "product") if recorder else nullcontext() as attempt:
            relaxed = restore_atom_order(rxembed.minimize(product).mol, product)
            if relaxed.GetNumConformers() != 1:
                raise ValueError("Product relaxation did not return exactly one conformer")
            geometry_id = recorder.record_geometry(relaxed, reactant.geometry.isomer, reactant.geometry.conformer,
                                                    None, "rxembed.minimize", attempt) if recorder else None
        if recorder:
            candidates.append(recorder.score(relaxed, reactant.geometry.isomer, reactant.geometry.conformer,
                                              None, "rxembed.minimize", geometry_id=geometry_id, scratch=scratch))
        else:
            energy = gxtb_energy(relaxed, xtb, scratch, xtb_threads)
            candidates.append(Geometry(relaxed, energy, reactant.geometry.isomer, reactant.geometry.conformer))
    except (RuntimeError, ValueError, TimeoutError) as exc:
        LOG.info("Local product relaxation failed (%s); trying templated embedding", exc)

    # Sample product conformers while holding the parent's non-coordinating
    # heavy-atom scaffold as a template. The product remains associated with
    # this specific reactant, even if local relaxation above failed.
    donor_atoms = {
        bond.GetBeginAtomIdx() for bond in product.GetBonds()
        if bond.GetBondType() == Chem.BondType.DATIVE
    }
    scaffold_map = {
        atom.GetIdx(): atom.GetIdx() for atom in product.GetAtoms()
        if atom.GetAtomicNum() > 1
        and atom.GetAtomicNum() != metal.atomic_number
        and atom.GetIdx() not in donor_atoms
    }
    try:
        candidates.extend(embedded_geometries(
            product_graph, n_conformers=n_conformers, seed=seed,
            xtb=xtb, scratch=scratch, xtb_threads=xtb_threads,
            template=(source, scaffold_map),
            recorder=recorder,
        ))
    except (RuntimeError, ValueError, TimeoutError) as exc:
        if not candidates:
            raise
        LOG.warning("Product conformer embedding failed; retaining local relaxation: %s", exc)
    if not candidates:
        raise ValueError("No product geometry survived relaxation or templated embedding")
    selected = min(candidates, key=lambda candidate: candidate.energy_eh)
    if recorder:
        for candidate in candidates:
            recorder.store.select_geometry(candidate.database_id, "retained" if candidate is selected else "not_selected",
                                           policy="lowest_energy_product_conformer", selected_geometry_id=selected.database_id)
    return selected


def save_geometry(item: Geometry, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(Chem.MolToXYZBlock(item.mol))
    item.xyz_path = path.resolve()


def ts_guess(
    reactant: Chem.Mol, product: Chem.Mol, path: Path,
    *, reactant_path: Path, product_path: Path,
    metal: MetalSpec, mm_steps: int = 1000,
    lambda_points: int = 21, snapshots: int = 10,
) -> None:
    """Align endpoint coordinates and save VeloxChem's MM scan TS guess."""
    import veloxchem as vlx

    reactant_copy = Chem.Mol(reactant)
    product_copy = Chem.Mol(product)
    if reactant_copy.GetNumAtoms() != product_copy.GetNumAtoms():
        raise ValueError("Reactant/product atom counts differ")
    if any(
        reactant_copy.GetAtomWithIdx(i).GetAtomicNum()
        != product_copy.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(reactant_copy.GetNumAtoms())
    ):
        raise ValueError("Reactant/product atom order differs")
    heavy = [(i, i) for i, atom in enumerate(reactant_copy.GetAtoms()) if atom.GetAtomicNum() > 1]
    rdMolAlign.AlignMol(product_copy, reactant_copy, atomMap=heavy)
    charge = Chem.GetFormalCharge(reactant_copy)
    if Chem.GetFormalCharge(product_copy) != charge:
        raise ValueError("Reactant/product charges differ")
    hydrogen_index, carbon_index = transfer_site(reactant_copy, metal)
    reactant_conf = reactant_copy.GetConformer()
    carbon_position = reactant_conf.GetAtomPosition(carbon_index)
    hydrogen_position = reactant_conf.GetAtomPosition(hydrogen_index)
    direction = hydrogen_position - carbon_position
    old_distance = direction.Length()
    if old_distance < 1e-8:
        raise ValueError("Transferring H and substrate carbon overlap")
    scaled = direction * (1.10 / old_distance)
    reactant_conf.SetAtomPosition(
        hydrogen_index,
        Point3D(
            carbon_position.x + scaled.x,
            carbon_position.y + scaled.y,
            carbon_position.z + scaled.z,
        ),
    )
    product_h_neighbors = [
        atom.GetIdx() for atom in product_copy.GetAtomWithIdx(hydrogen_index).GetNeighbors()
    ]
    if len(product_h_neighbors) != 1:
        raise ValueError("Transferred hydrogen has no unique product acceptor")
    acceptor_index = product_h_neighbors[0]
    reactant_endpoint_xyz = Chem.MolToXYZBlock(reactant_copy)
    product_endpoint_xyz = Chem.MolToXYZBlock(product_copy)

    with tempfile.TemporaryDirectory(prefix="velox-ts-", dir=path.parent) as work:
        original_cwd = Path.cwd()
        try:
            os.chdir(work)
            start = vlx.Molecule.read_xyz_string(reactant_endpoint_xyz)
            end = vlx.Molecule.read_xyz_string(product_endpoint_xyz)
            for endpoint in (start, end):
                endpoint.set_charge(charge)
                endpoint.set_multiplicity(1)
                if not endpoint.check_multiplicity():
                    raise ValueError("Invalid endpoint charge/multiplicity")
            guesser = vlx.TransitionStateGuesser()
            guesser.do_qm_scan = False
            guesser.save_results_file = False
            guesser.mm_steps = mm_steps
            guesser.lambda_vector = [i / (lambda_points - 1) for i in range(lambda_points)]
            guesser.conformer_snapshots = snapshots
            # The RDKit product was already renumbered to the reactant atom order.
            # Supplying this identity map prevents the XYZ-based matcher from
            # exchanging equivalent hydrogens or inventing a different reaction.
            mapping = {i + 1: i + 1 for i in range(reactant_copy.GetNumAtoms())}
            result = guesser.find_transition_state(
                start, end, product_mapping=mapping,
            )
            breaking = {tuple(sorted(bond)) for bond in result["breaking_bonds"]}
            forming = {tuple(sorted(bond)) for bond in result["forming_bonds"]}
            expected_break = tuple(sorted((hydrogen_index, carbon_index)))
            expected_form = tuple(sorted((hydrogen_index, acceptor_index)))
            if expected_break not in breaking or expected_form not in forming:
                raise ValueError(
                    f"VeloxChem perceived {breaking=} and {forming=}; expected "
                    f"C-H {expected_break} to break and acceptor-H {expected_form} to form"
                )
            xyz = result.get("max_mm_xyz")
            if not xyz:
                raise RuntimeError("VeloxChem returned no max_mm_xyz")
        finally:
            os.chdir(original_cwd)
    path.write_text(xyz)
    reactant_path.write_text(reactant_endpoint_xyz)
    product_path.write_text(product_endpoint_xyz)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._") or "ligand"


def process_ligand(row: dict, row_number: int, args, scratch: Path, failures: list[dict]):
    ligand_id = str(row.get(args.id_column) or row_number)
    ligand_smiles = str(row.get(args.smiles_column) or "").strip()
    if not ligand_smiles:
        raise ValueError(f"Ligand {ligand_id}: missing {args.smiles_column!r}")
    ligand_dir = ligand_directory(args, row, row_number)
    ligand_dir.mkdir(parents=True, exist_ok=True)
    metal = parse_metal_smiles(args.metal_smiles)
    substrate = Chem.MolFromSmiles(args.substrate)
    if substrate is None:
        raise ValueError(f"Invalid substrate SMILES: {args.substrate!r}")
    substrate = Chem.AddHs(substrate)
    sites = substrate_sites(substrate)
    if not sites:
        raise ValueError("Substrate has no aromatic C-H site")

    def failed(stage, label, exc):
        LOG.warning("Ligand %s %s %s: %s", ligand_id, stage, label, exc)
        failures.append({"ligand_id": ligand_id, "stage": stage,
                         "candidate": str(label), "error": str(exc)})
        if getattr(args, "_store", None):
            task = args._store.create_task(args._run_ligand_id,
                                          f"{args._sampling_scope}:failure:{stage}:{label}", stage,
                                          parameters={"candidate": str(label)})
            attempt = args._store.start_attempt(task)
            args._store.finish_attempt(attempt, "failed", exc)

    # Stage 1: one ligand+metal conformer per donor combination and coordination
    # isomer, followed by the ligand-only g-xTB energy cutoff.
    ligand_modes = []
    for donor_number, (donors, complex_mol) in enumerate(
        ligand_complexes(ligand_smiles, args.n_donors, metal)
    ):
        key = f"d{donor_number}"
        try:
            recorder = sampling_recorder(args, key, "ligand_metal", complex_mol, donors)
            for geom in embedded_geometries(
                complex_mol, n_conformers=1,
                seed=args.seed + donor_number * 1000, xtb=args.xtb,
                scratch=scratch, xtb_threads=args.xtb_threads,
                recorder=recorder,
            ):
                ligand_modes.append(LigandMode(key, donors, ligand_smiles, geom))
        except (RuntimeError, ValueError, TimeoutError) as exc:
            failed("ligand", key, exc)
    scored_ligand_modes = ligand_modes
    ligand_modes = within_cutoff(scored_ligand_modes, args.cutoff_kcal)
    record_cutoff(args, scored_ligand_modes, ligand_modes, args.cutoff_kcal)
    minimum_ligand_energy = min(
        (mode.geometry.energy_eh for mode in scored_ligand_modes), default=None
    )
    retained_keys = {(mode.key, mode.geometry.isomer) for mode in ligand_modes}
    energy_rows = []
    for mode in sorted(scored_ligand_modes, key=lambda item: item.geometry.energy_eh):
        relative = (mode.geometry.energy_eh - minimum_ligand_energy) * HARTREE_TO_KCAL
        retained = (mode.key, mode.geometry.isomer) in retained_keys
        energy_rows.append({
            "donor_indices": ",".join(map(str, mode.donors)),
            "ligand_isomer": mode.geometry.isomer,
            "ligand_conformer": mode.geometry.conformer,
            "energy_eh": mode.geometry.energy_eh,
            "relative_kcal_mol": relative,
            "formal_charge": Chem.GetFormalCharge(mode.geometry.mol),
            "retained": retained,
        })
        if retained:
            LOG.info(
                "Ligand %s retained %s isomer donors=%s isomer=%s charge=%+d: ΔE=%.2f kcal/mol (%.8f Eh)",
                ligand_id, metal.symbol, mode.donors, mode.geometry.isomer,
                Chem.GetFormalCharge(mode.geometry.mol), relative,
                mode.geometry.energy_eh,
            )
    write_csv(ligand_dir / "ligand_metal_isomers.csv", LIGAND_COLUMNS, energy_rows)
    LOG.info(
        "Ligand %s: kept %s/%s scored ligand+%s isomers within %.1f kcal/mol",
        ligand_id, len(ligand_modes), len(scored_ligand_modes), metal.symbol,
        args.cutoff_kcal,
    )

    # Stage 2: both symmetry-distinct substrate sites, five seeds per binding mode.
    reactants = []
    for ligand_mode_index, ligand_mode in enumerate(ligand_modes):
        ligand_geom = ligand_mode.geometry.mol
        for site in sites:
            key = f"l{ligand_mode_index}_s{site}"
            try:
                graph = build_reactant(ligand_geom, substrate, site, metal)
                recorder = sampling_recorder(args, key, "reactant", graph, ligand_mode.donors,
                                             ligand_mode.geometry.database_id, substrate_atom_index=site)
                # rxembed permits at most one pinned coordination-sphere atom.
                # Fix the ligand scaffold while its metal and donor positions relax.
                donor_atoms = {
                    bond.GetBeginAtomIdx() for bond in ligand_geom.GetBonds()
                    if bond.GetBondType() == Chem.BondType.DATIVE
                }
                scaffold_map = {
                    atom.GetIdx(): atom.GetIdx() for atom in ligand_geom.GetAtoms()
                    if atom.GetAtomicNum() > 1
                    and atom.GetAtomicNum() != metal.atomic_number
                    and atom.GetIdx() not in donor_atoms
                }
                template = (ligand_geom, scaffold_map) if scaffold_map else None
                for geom in embedded_geometries(
                    graph, n_conformers=args.conformers,
                    seed=args.seed + 100000 + ligand_mode_index * 1000 + site * 10,
                    xtb=args.xtb, scratch=scratch, xtb_threads=args.xtb_threads,
                    template=template,
                    recorder=recorder,
                ):
                    reactants.append(ReactantMode(key, ligand_mode, site, graph, geom))
            except (RuntimeError, ValueError, TimeoutError, Chem.MolSanitizeException) as exc:
                failed("reactant", key, exc)
    LOG.info("Ligand %s: %s scored reactant geometries", ligand_id, len(reactants))

    # Stage 3: one endpoint pair per donor set / square-planar ligand isomer.
    # Try reactants in energy order. For the lowest one with viable products,
    # sample product conformers and keep its lowest-energy product channel.
    reactants_by_isomer = {}
    for reactant in reactants:
        isomer_key = (reactant.ligand.donors, reactant.ligand.geometry.isomer)
        reactants_by_isomer.setdefault(isomer_key, []).append(reactant)
    for key, modes in reactants_by_isomer.items():
        retained = within_cutoff(modes, args.cutoff_kcal)
        record_cutoff(args, modes, retained, args.cutoff_kcal)
        reactants_by_isomer[key] = retained
    products = []
    graphs = {mode.key: mode.graph for mode in reactants}
    reaction_channels = {
        key: list(products_from_resonance(graph, metal)) for key, graph in graphs.items()
    }
    for isomer_number, (isomer_key, candidates) in enumerate(reactants_by_isomer.items()):
        candidates.sort(key=lambda mode: mode.geometry.energy_eh)
        for reactant_number, reactant in enumerate(candidates):
            corresponding_products = []
            for product_number, info in enumerate(reaction_channels[reactant.key]):
                product_graph, product_smiles, form, acceptor, element, metal_bound = info
                try:
                    recorder = sampling_recorder(
                        args, f"product:{reactant.geometry.database_id}:{product_number}", "product",
                        product_graph, reactant.ligand.donors, reactant.geometry.database_id,
                        substrate_atom_index=reactant.substrate_atom_index, resonance_form=form,
                        acceptor_atom_index=acceptor, acceptor_element=element, acceptor_metal_bound=metal_bound,
                    )
                    geom = product_from_reactant(
                        reactant, product_graph, acceptor,
                        n_conformers=args.conformers,
                        seed=args.seed + 200000 + isomer_number * 100000
                        + reactant_number * 1000 + product_number * 10,
                        xtb=args.xtb, scratch=scratch, xtb_threads=args.xtb_threads,
                        metal=metal,
                        recorder=recorder,
                    )
                    corresponding_products.append(ProductMode(
                        reactant, product_smiles, form, acceptor, element, metal_bound, geom
                    ))
                except (RuntimeError, ValueError, TimeoutError) as exc:
                    failed("product", f"{reactant.key}_c{reactant.geometry.conformer}_p{product_number}", exc)
            if corresponding_products:
                selected = min(
                    corresponding_products, key=lambda mode: mode.geometry.energy_eh
                )
                products.append(selected)
                if getattr(args, "_store", None):
                    for channel in corresponding_products:
                        args._store.select_geometry(channel.geometry.database_id,
                                                    "retained" if channel is selected else "not_selected",
                                                    policy="lowest_energy_product_channel",
                                                    selected_geometry_id=selected.geometry.database_id)
                LOG.info(
                    "Ligand %s isomer %s: selected reactant %.8f Eh and product %.8f Eh",
                    ligand_id, isomer_key, reactant.geometry.energy_eh,
                    selected.geometry.energy_eh,
                )
                break
        else:
            LOG.warning("Ligand %s isomer %s: no valid reactant/product pair", ligand_id, isomer_key)
    LOG.info("Ligand %s: %s ligand-isomer endpoint pairs selected", ligand_id, len(products))

    rows = []
    for pair_number, product in enumerate(products, 1):
        reactant = product.reactant
        store = getattr(args, "_store", None)
        pair_id = None
        task_parameters = {"method": "VeloxChem MM scan", "mm_steps": args.ts_mm_steps,
                           "lambda_points": args.ts_lambda_points, "snapshots": args.ts_snapshots,
                           "charge": Chem.GetFormalCharge(reactant.geometry.mol), "multiplicity": 1}
        task_key = f"{getattr(args, '_sampling_scope', '')}:ts_guess:{pair_number}"
        if store:
            pair_id = store.pair(
                args._run_ligand_id, f"{args._sampling_scope}:pair:{pair_number}",
                reactant.geometry.database_id, product.geometry.database_id, reactant.ligand.geometry.database_id,
                atom_mapping={i: i for i in range(reactant.geometry.mol.GetNumAtoms())},
                reason={"policy": "lowest_reactant_with_viable_product_then_lowest_product",
                        "cutoff_kcal_mol": args.cutoff_kcal}, atom_index_base=0,
            )
            for geometry in (reactant.geometry, product.geometry):
                store.select_geometry(geometry.database_id, "selected_for_pair", pair_id=pair_id)
            store.create_task(args._run_ligand_id, task_key, "ts_guess", parameters=task_parameters, pair_id=pair_id)
        if args.max_pairs_per_ligand is not None and pair_number > args.max_pairs_per_ligand:
            continue
        pair_dir = ligand_dir / f"pair_{pair_number:05d}"
        pair_dir.mkdir()
        reactant_xyz = pair_dir / "reactant.xyz"
        product_xyz = pair_dir / "product.xyz"
        guess_xyz = pair_dir / "ts_guess.xyz"
        ts_reactant_xyz = pair_dir / "ts_reactant_input.xyz"
        ts_product_xyz = pair_dir / "ts_product_input.xyz"
        save_geometry(reactant.geometry, reactant_xyz)
        save_geometry(product.geometry, product_xyz)
        try:
            with store.operation(args._run_ligand_id, task_key, "ts_guess", parameters=task_parameters, pair_id=pair_id) if store else nullcontext() as attempt:
                ts_guess(
                    reactant.geometry.mol, product.geometry.mol, guess_xyz,
                    reactant_path=ts_reactant_xyz, product_path=ts_product_xyz,
                    mm_steps=args.ts_mm_steps,
                    lambda_points=args.ts_lambda_points,
                    snapshots=args.ts_snapshots,
                    metal=metal,
                )
                if store:
                    guess = store.geometry(args._run_ligand_id, f"{args._sampling_scope}:ts_guess:{pair_number}",
                                           "veloxchem_mm_scan", parent_geometry_id=reactant.geometry.database_id,
                                           xyz=guess_xyz.read_bytes(), attempt_id=attempt, pair_id=pair_id)
                    store.pair_outputs(pair_id, guess,
                                       store.file_artifact(ts_reactant_xyz, "ts_endpoint_xyz"),
                                       store.file_artifact(ts_product_xyz, "ts_endpoint_xyz"))
        except (RuntimeError, ValueError, OSError, ImportError) as exc:
            failed("ts", pair_number, exc)
            continue
        rows.append({
            "ligand_id": ligand_id,
            "ligand_smiles": ligand_smiles,
            "metal_smiles": metal.smiles,
            "donor_indices": ",".join(map(str, reactant.ligand.donors)),
            "ligand_isomer": reactant.ligand.geometry.isomer,
            "ligand_conformer": reactant.ligand.geometry.conformer,
            "ligand_energy_eh": reactant.ligand.geometry.energy_eh,
            "substrate_atom_index": reactant.substrate_atom_index,
            "reactant_smiles": smiles(reactant.graph),
            "reactant_isomer": reactant.geometry.isomer,
            "reactant_conformer": reactant.geometry.conformer,
            "reactant_energy_eh": reactant.geometry.energy_eh,
            "product_smiles": product.smiles,
            "resonance_form": product.resonance_form,
            "acceptor_atom_index": product.acceptor_atom_index,
            "acceptor_element": product.acceptor_element,
            "acceptor_metal_bound": product.acceptor_metal_bound,
            "product_isomer": product.geometry.isomer,
            "product_conformer": product.geometry.conformer,
            "product_energy_eh": product.geometry.energy_eh,
            "reactant_xyz": str(reactant_xyz.resolve()),
            "product_xyz": str(product_xyz.resolve()),
            "ts_guess_xyz": str(guess_xyz.resolve()),
            "ts_reactant_input_xyz": str(ts_reactant_xyz.resolve()),
            "ts_product_input_xyz": str(ts_product_xyz.resolve()),
        })
    LOG.info("Ligand %s: %s reactant/product/TS pairs", ligand_id, len(rows))
    return rows


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligands-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--database", type=Path,
                        help="SQLite database path (default: OUTPUT_DIR/workflow.sqlite)")
    parser.add_argument("--smiles-column", default="deprotonated")
    parser.add_argument("--id-column", default="ID")
    parser.add_argument("--substrate", default=DEFAULT_SUBSTRATE)
    parser.add_argument("--metal-smiles", "--metal", default=DEFAULT_METAL_SMILES,
                        help="Single-atom metal SMILES, e.g. '[Ni+2]' (default: [Pd+2])")
    parser.add_argument("--xtb", type=Path, default=Path("/home/henry/opt/xtb-6.7.1/bin/xtb"))
    parser.add_argument("--conformers", type=int, default=5,
                        help="Seeds per reactant/product coordination isomer; ligand+metal uses one")
    parser.add_argument("--cutoff-kcal", type=float, default=20.0)
    parser.add_argument("--n-donors", type=int, default=3)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--xtb-threads", type=int, default=1)
    parser.add_argument("--ts-mm-steps", type=int, default=1000)
    parser.add_argument("--ts-lambda-points", type=int, default=21)
    parser.add_argument("--ts-snapshots", type=int, default=10)
    parser.add_argument("--max-ligands", type=int, default=None,
                        help="Limit ligand rows for a small demonstration run")
    parser.add_argument("--max-pairs-per-ligand", type=int, default=None,
                        help="Cap selected ligand-isomer TS guesses after energy filtering")
    parser.add_argument("--resume-from", type=int, metavar="ROW",
                        help="Resume in an existing output folder from this 1-based input CSV row")
    parser.add_argument("--archive-incomplete", action="store_true",
                        help="On resume, move nonempty unfinished ligand folders aside before retrying")
    args = parser.parse_args(argv)
    try:
        args.metal_smiles = parse_metal_smiles(args.metal_smiles).smiles
    except ValueError as exc:
        parser.error(str(exc))
    if min(args.conformers, args.n_donors, args.xtb_threads,
           args.ts_mm_steps, args.ts_snapshots) < 1 or args.ts_lambda_points < 3:
        parser.error("counts must be positive and ts-lambda-points must be at least 3")
    if not math.isfinite(args.cutoff_kcal) or args.cutoff_kcal < 0:
        parser.error("cutoff-kcal must be non-negative")
    if args.max_ligands is not None and args.max_ligands < 1:
        parser.error("--max-ligands must be positive")
    if args.max_pairs_per_ligand is not None and args.max_pairs_per_ligand < 0:
        parser.error("--max-pairs-per-ligand must be non-negative")
    if args.resume_from is not None and args.resume_from < 1:
        parser.error("--resume-from must be a positive CSV row number")
    if args.archive_incomplete and args.resume_from is None:
        parser.error("--archive-incomplete requires --resume-from")
    if not args.xtb.is_file():
        parser.error(f"xTB executable not found: {args.xtb}")
    if not args.ligands_csv.is_file():
        parser.error(f"Ligand CSV not found: {args.ligands_csv}")
    args.output_dir = args.output_dir.resolve()
    args.database = args.database.resolve() if args.database else args.output_dir / "workflow.sqlite"
    if args.resume_from is None:
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            parser.error(f"Output directory is not empty: {args.output_dir}")
    elif not args.output_dir.is_dir():
        parser.error(f"Resume output directory does not exist: {args.output_dir}")
    return args


def execute_workflow(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scratch = args.output_dir / "scratch"
    scratch.mkdir(exist_ok=args.resume_from is not None)
    with args.ligands_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if args.smiles_column not in (reader.fieldnames or []):
            raise ValueError(f"CSV has no {args.smiles_column!r} column")
        input_rows = list(reader)
    if args.max_ligands is not None:
        input_rows = input_rows[:args.max_ligands]
    if args.resume_from is not None and args.resume_from > len(input_rows):
        raise ValueError(f"--resume-from {args.resume_from} exceeds {len(input_rows)} input rows")
    if args.resume_from is not None:
        ligand_ids = [str(row.get(args.id_column) or number)
                      for number, row in enumerate(input_rows, 1)]
        if len(set(ligand_ids)) != len(ligand_ids):
            raise ValueError("Resuming requires unique ligand IDs in the input CSV")

    old_results_path = args.output_dir / "reaction_pairs.csv"
    old_failures_path = args.output_dir / "failures.csv"
    old_results_by_id = {}
    old_failures_by_id = {}
    if args.resume_from is not None and old_results_path.is_file():
        for result in read_results_csv(old_results_path):
            old_results_by_id.setdefault(result["ligand_id"], []).append(result)
    if args.resume_from is not None and old_failures_path.is_file():
        for failure in read_csv(old_failures_path):
            old_failures_by_id.setdefault(failure["ligand_id"], []).append(failure)

    all_rows = []
    failures = []
    missing_legacy_metadata = False

    def database_ligand(row, number):
        return args._store.ensure_ligand(args._run_id, args._input_hash, number,
                                        str(row.get(args.id_column) or number),
                                        str(row.get(args.smiles_column) or "").strip())

    def restore_database(row, number, directory, results, local_failures):
        ligand = database_ligand(row, number)
        if not args._store.connection.execute("SELECT 1 FROM chemical_state WHERE run_ligand_id = ?", (ligand,)).fetchone():
            import_ligand(args._store, ligand, directory, results, local_failures, args.metal_smiles)
        has_errors = local_failures or args._store.has_errors(args._run_id, ligand)
        args._store.set_ligand_status(ligand, "completed_with_errors" if has_errors else "succeeded" if results else "no_candidates")
    for row_number, row in enumerate(input_rows, 1):
        if args.resume_from is None or row_number >= args.resume_from:
            break
        ligand_id = str(row.get(args.id_column) or row_number)
        ligand_dir = ligand_directory(args, row, row_number)
        checkpoint = load_checkpoint(ligand_dir, row, row_number, args)
        if checkpoint is not None:
            previous_rows, previous_failures = checkpoint
            LOG.info("Resume: ligand %s complete; loaded its checkpoint", ligand_id)
        elif ligand_dir.is_dir() and any(ligand_dir.iterdir()):
            previous_rows = old_results_by_id.get(ligand_id, [])
            previous_failures = old_failures_by_id.get(ligand_id, [])
            if not old_results_path.is_file():
                missing_legacy_metadata = True
                LOG.warning(
                    "Resume: preserving legacy ligand %s folder, but its CSV metadata is unavailable",
                    ligand_id,
                )
            else:
                LOG.info("Resume: loaded legacy ligand %s rows from %s", ligand_id, old_results_path)
        else:
            raise ValueError(
                f"Cannot skip ligand {ligand_id}: no completed checkpoint or legacy output folder"
            )
        all_rows.extend(previous_rows)
        failures.extend(previous_failures)
        restore_database(row, row_number, ligand_dir, previous_rows, previous_failures)

    if missing_legacy_metadata:
        output_csv = args.output_dir / f"reaction_pairs_resumed_from_{args.resume_from:04d}.csv"
        failure_csv = args.output_dir / f"failures_resumed_from_{args.resume_from:04d}.csv"
        LOG.warning("Earlier legacy rows cannot be recovered; resumed rows will be in %s", output_csv)
    else:
        output_csv = old_results_path
        failure_csv = old_failures_path

    def save_aggregate():
        write_csv(output_csv, CSV_COLUMNS, all_rows)
        write_csv(failure_csv, FAILURE_COLUMNS, failures)

    start_row = args.resume_from or 1
    for row_number in range(start_row, len(input_rows) + 1):
        row = input_rows[row_number - 1]
        ligand_dir = ligand_directory(args, row, row_number)
        ligand_id = str(row.get(args.id_column) or row_number)
        if args.resume_from is not None:
            checkpoint = load_checkpoint(ligand_dir, row, row_number, args)
            if checkpoint is not None:
                previous_rows, previous_failures = checkpoint
                all_rows.extend(previous_rows)
                failures.extend(previous_failures)
                restore_database(row, row_number, ligand_dir, previous_rows, previous_failures)
                LOG.info("Resume: ligand %s complete; skipping it", ligand_id)
                save_aggregate()
                continue
            if ligand_dir.exists() and any(ligand_dir.iterdir()):
                if not args.archive_incomplete:
                    raise ValueError(
                        f"Ligand {ligand_id} has unfinished files in {ligand_dir}; "
                        "stop the old run and use --archive-incomplete to preserve them before retrying"
                    )
                stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                archive = ligand_dir.with_name(f"{ligand_dir.name}.interrupted-{stamp}")
                suffix = 1
                while archive.exists():
                    archive = ligand_dir.with_name(
                        f"{ligand_dir.name}.interrupted-{stamp}-{suffix}"
                    )
                    suffix += 1
                ligand_dir.rename(archive)
                LOG.warning("Archived unfinished ligand %s output at %s", ligand_id, archive)
            LOG.info("Resume: processing ligand %s (input row %s)", ligand_id, row_number)
        local_failures = []
        args._run_ligand_id = database_ligand(row, row_number)
        args._store.interrupt_active(args._run_ligand_id)
        args._store.set_ligand_status(args._run_ligand_id, "running")
        try:
            with args._store.operation(args._run_ligand_id, "process_ligand", "process_ligand",
                                       attempt_metadata={"software": args._software}) as attempt:
                args._sampling_scope = f"attempt:{attempt}"
                results = process_ligand(row, row_number, args, scratch, local_failures)
                save_checkpoint(ligand_dir, row, row_number, args, results, local_failures)
                all_rows.extend(results)
            has_errors = local_failures or args._store.has_errors(args._run_id, args._run_ligand_id)
            has_pairs = args._store.connection.execute("SELECT 1 FROM reaction_pair WHERE run_ligand_id = ?",
                                                       (args._run_ligand_id,)).fetchone()
            args._store.set_ligand_status(args._run_ligand_id,
                                         "completed_with_errors" if has_errors else "succeeded" if results or has_pairs else "no_candidates")
        except (RuntimeError, ValueError, OSError) as exc:
            LOG.exception("Ligand %s failed", ligand_id)
            local_failures.append({"ligand_id": ligand_id, "stage": "ligand",
                                   "candidate": "all", "error": str(exc)})
            args._store.set_ligand_status(args._run_ligand_id, "failed")
        failures.extend(local_failures)
        save_aggregate()
    LOG.info("Wrote %s pairs to %s", len(all_rows), output_csv)
    if failures:
        LOG.warning("%s candidates failed; see failures.csv", len(failures))
    return 0 if all_rows else 1


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args(argv)
    install_rdkit_rxembed_compatibility()
    args._input_hash = digest(args.ligands_csv.read_bytes())
    software = {"workflow_sha256": digest(Path(__file__).read_bytes()),
                "store_sha256": digest(Path(__file__).with_name("workflow_store.py").read_bytes())}
    for package in ("rdkit", "rxembed", "veloxchem"):
        try:
            software[package] = package_metadata.version(package)
        except package_metadata.PackageNotFoundError:
            software[package] = "unknown"
    try:
        version = subprocess.run([str(args.xtb), "--version"], capture_output=True, text=True, timeout=5)
        args._xtb_version = (version.stdout or version.stderr).strip() if version.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        args._xtb_version = None
    software["xtb"] = args._xtb_version or "unknown"
    args._software = software
    with WorkflowStore(args.database) as store:
        args._store = store
        try:
            args._run_id = store.create_run(f"workflow:{args.output_dir}", args._input_hash, run_options(args), software)
        except ValueError as exc:
            raise ValueError(f"Checkpoint does not match this input or run configuration: {exc}") from exc
        store.set_run_status(args._run_id, "running")
        try:
            result = execute_workflow(args)
        except BaseException as exc:
            status = "interrupted" if isinstance(exc, (KeyboardInterrupt, SystemExit)) else "failed"
            store.set_run_status(args._run_id, status)
            if getattr(args, "_run_ligand_id", None) is not None:
                row = store.connection.execute("SELECT status FROM run_ligand WHERE id = ?", (args._run_ligand_id,)).fetchone()
                if row and row[0] == "running":
                    store.set_ligand_status(args._run_ligand_id, status)
            raise
        store.set_run_status(args._run_id, "completed_with_errors" if store.has_errors(args._run_id) else "succeeded" if result == 0 else "no_candidates")
        return result


if __name__ == "__main__":
    raise SystemExit(main())
