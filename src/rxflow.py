"""Generate mapped reaction pairs following reaction_smarts.ipynb.

Keep explicit hydrogens and stable atom maps for later coembedding and TS work.
The original workflow is retained in legacy/reaction_workflow.py as a reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from itertools import combinations
from importlib.metadata import distribution

from rdkit import Chem
from rdkit.Chem import AllChem


LOG = logging.getLogger(__name__)
HARTREE_TO_KCAL = 627.5094740631
SUBSTRATE = "[H]C1=CC(OC)=C(OC)C(OC)=C1"
ENERGY = re.compile(r"TOTAL ENERGY\s+([-+\d.eEdD]+)\s+Eh")


def rxembed_revision():
    """Report the installed Git revision, including non-editable Git installations."""
    metadata = json.loads(distribution("rxembed").read_text("direct_url.json") or "{}")
    if metadata.get("vcs_info", {}).get("commit_id"):
        return metadata["vcs_info"]["commit_id"]
    import rxembed

    checkout = Path(rxembed.__file__).resolve().parents[2]
    result = subprocess.run(["git", "-C", str(checkout), "rev-parse", "HEAD"], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def read_ligands(path):
    """Read the supplied CSV; accept metal or central_metal as the metal column."""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        required = {"ID", "ligand_smiles", "substrate", "metal_geom", "denticity"}
        if not required <= columns or not {"metal", "central_metal"} & columns:
            raise ValueError(f"CSV requires {sorted(required)} and metal or central_metal")
        rows = list(reader)
    identifiers = set()
    for row in rows:
        row["ID"] = row["ID"].strip()
        if not row["ID"] or row["ID"] in identifiers:
            raise ValueError("Ligand IDs must be nonempty and unique")
        identifiers.add(row["ID"])
        row["substrate"] = (row["substrate"] or "").strip()
        substrate = Chem.MolFromSmiles(row["substrate"])
        if substrate is None or not substrate.GetNumAtoms():
            raise ValueError(f"Ligand {row['ID']}: invalid substrate SMILES")
        value = (row.get("metal") or row.get("central_metal") or "").strip()
        row["metal"] = value if value.startswith("[") else f"[{value}]"
        metal = Chem.MolFromSmiles(row["metal"])
        if metal is None or metal.GetNumAtoms() != 1:
            raise ValueError(f"Ligand {row['ID']}: expected one metal atom")
        row["metal_geom"] = row["metal_geom"].strip().upper()
        row["denticity"] = int(row["denticity"])
        if row["denticity"] < 1:
            raise ValueError("Denticity must be positive")
    return rows


def create_permutations(ligand_smiles, metal, denticity):
    """Enumerate donor subsets with ligand-to-metal dative bonds, as in the notebook."""
    ligand = Chem.MolFromSmiles(ligand_smiles)
    if ligand is None:
        raise ValueError(f"Invalid ligand SMILES: {ligand_smiles}")
    ligand = Chem.AddHs(ligand)
    donors = [a.GetIdx() for a in ligand.GetAtoms() if a.GetAtomicNum() not in (1, 6)]
    if denticity > len(donors):
        raise ValueError(f"Denticity {denticity} exceeds {len(donors)} donor candidates")
    for selected in combinations(donors, denticity):
        editable = Chem.RWMol(Chem.CombineMols(ligand, Chem.MolFromSmiles(metal)))
        for index in selected:
            editable.AddBond(index, ligand.GetNumAtoms(), Chem.BondType.DATIVE)
        mol = editable.GetMol()
        Chem.SanitizeMol(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        yield selected, mol


def retain_coordination_mode(donors, ligand, relative_kcal_mol, cutoff=20.0):
    """Keep low-energy modes plus tridentate N,N,N arrangements.

    Tridentate nitrogen binding is chemically informative for the exchange
    screen even when the initial gxTB coordination score puts it just outside
    the broad conformer-pruning window.
    """
    if relative_kcal_mol is None or relative_kcal_mol <= cutoff + 1e-8:
        return True
    return len(donors) == 3 and all(ligand.GetAtomWithIdx(index).GetAtomicNum() == 7
                                    for index in donors)


def reaction_pairs_for_coordination_mode(row, mode_index, donors, isomer_index, mol,
                                         coordination_smiles):
    """Build the substrate reaction pairs for one already-screened metal isomer."""
    pairs = []
    for site, reactant in add_substrate(mol, row["substrate"]):
        for resonance, acceptor, product in create_products(reactant, row["metal"]):
            pairs.append({"ligand_id": row["ID"], "mode": mode_index,
                          "donors": json.dumps(donors), "isomer": isomer_index,
                          "metal": row["metal"], "metal_geom": row["metal_geom"],
                          "coordination_smiles": coordination_smiles,
                          "substrate": row["substrate"], "site": site,
                          "substrate_site": int(site.removeprefix("site_")),
                          "resonance": resonance, "acceptor_map": acceptor + 1,
                          "reactant_smiles": Chem.MolToSmiles(reactant),
                          "product_smiles": Chem.MolToSmiles(product)})
    return pairs


def ligand_directory(results, ligand_id):
    """Name ligand folders by CSV ID, with at least two digits for numeric IDs."""
    identifier = str(ligand_id)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", identifier):
        raise ValueError(f"Invalid ligand ID for a directory: {identifier!r}")
    return Path(results) / f"ligand_{identifier.zfill(2)}"


def run_gxtb(mol, directory, executable, *, name="calculation", optimize=False, solvent=None):
    """Run g-xTB in temporary scratch and save named files directly in directory."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise ValueError(f"Invalid calculation name: {name!r}")
    if any(directory.glob(f"{name}_*")):
        raise FileExistsError(f"Calculation {name} already exists in {directory}; choose another name")
    input_name = f"{name}_input.xyz"
    (directory / input_name).write_text(Chem.MolToXYZBlock(mol))
    command = [str(executable), input_name, "--gxtb", "--charge",
               str(Chem.GetFormalCharge(mol)),
               "--opt" if optimize else "--sp"]
    if solvent:
        command.extend(["--alpb", solvent])
    (directory / f"{name}_command.json").write_text(json.dumps(command, indent=2) + "\n")
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    with tempfile.TemporaryDirectory(prefix="rxflow-gxtb-") as scratch:
        work = Path(scratch)
        (work / input_name).write_text(Chem.MolToXYZBlock(mol))
        result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True)
        for path in work.iterdir():
            if path.suffix in (".xyz", ".log") and path.name != input_name:
                shutil.copy2(path, directory / f"{name}_{path.name}")
    (directory / f"{name}_stdout.log").write_text(result.stdout)
    (directory / f"{name}_stderr.log").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"g-xTB exited {result.returncode}; see {directory}: {result.stderr.strip()[-500:]}")
    matches = ENERGY.findall(result.stdout)
    if not matches:
        raise RuntimeError(f"No TOTAL ENERGY in {directory / f'{name}_stdout.log'}")
    energy = float(matches[-1].replace("D", "E").replace("d", "e"))
    if not math.isfinite(energy):
        raise ValueError("g-xTB returned a nonfinite energy")
    geometry = Chem.Mol(mol)
    if optimize:
        if "GEOMETRY OPTIMIZATION CONVERGED" not in result.stdout:
            raise RuntimeError(f"g-xTB optimization did not converge; see {directory}")
        optimized = Chem.MolFromXYZFile(str(directory / f"{name}_xtbopt.xyz"))
        if optimized is None or [a.GetSymbol() for a in optimized.GetAtoms()] != [
            a.GetSymbol() for a in mol.GetAtoms()
        ]:
            raise ValueError("g-xTB changed the atom inventory")
        geometry.RemoveAllConformers()
        geometry.AddConformer(Chem.Conformer(optimized.GetConformer()), assignId=True)
    return energy, geometry


def add_substrate(complex_mol, substrate_smiles=SUBSTRATE):
    """Attach one representative of each symmetry-distinct aromatic C-H site."""
    substrate = Chem.MolFromSmiles(substrate_smiles)
    if substrate is None or not substrate.GetNumAtoms():
        raise ValueError(f"Invalid substrate SMILES: {substrate_smiles}")
    for atom in substrate.GetAtoms():
        atom.SetAtomMapNum(0)  # Input map labels must not break molecular symmetry.
    substrate = Chem.AddHs(substrate)
    ranks = Chem.CanonicalRankAtoms(substrate, breakTies=False)
    sites = {}
    for atom in substrate.GetAtoms():
        if atom.GetIsAromatic() and atom.GetAtomicNum() == 6 and any(
            neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()
        ):
            sites.setdefault(ranks[atom.GetIdx()], atom.GetIdx())
    metal_index = complex_mol.GetNumAtoms() - 1
    for site in sites.values():
        editable = Chem.RWMol(Chem.CombineMols(complex_mol, substrate))
        editable.AddBond(complex_mol.GetNumAtoms() + site, metal_index, Chem.BondType.DATIVE)
        reactant = editable.GetMol()
        reactant.RemoveAllConformers()
        Chem.SanitizeMol(reactant)
        for atom in reactant.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        yield f"site_{site + 1:03d}", reactant


def create_products(reactant, metal):
    """Apply the notebook H-transfer SMARTS over resonance forms; preserve atom order."""
    mapped_metal = f"{metal[:-1]}:3]"
    reaction = AllChem.ReactionFromSmarts(
        f"([c:1](-[#1:4])->{mapped_metal}.[#7,#8,#16;-1:2])"
        f">>([c-:1]->{mapped_metal}.[*+0:2]-[#1:4])"
    )
    seen = set()
    for resonance, form in enumerate(Chem.ResonanceMolSupplier(reactant)):
        for (raw,) in reaction.RunReactants((form,)):
            try:
                source = [a.GetIntProp("react_atom_idx") for a in raw.GetAtoms()]
                if sorted(source) != list(range(reactant.GetNumAtoms())):
                    raise ValueError("Reaction changed the atom inventory")
                acceptor = next(a.GetIntProp("react_atom_idx") for a in raw.GetAtoms()
                                if a.HasProp("old_mapno") and a.GetIntProp("old_mapno") == 2)
                product = Chem.RenumberAtoms(raw, [source.index(i) for i in range(len(source))])
                Chem.SanitizeMol(product)
            except Chem.MolSanitizeException:
                continue
            for atom in product.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
            key_mol = Chem.Mol(product)
            for atom in key_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            key = Chem.MolToSmiles(key_mol)
            if key not in seen:
                seen.add(key)
                product.RemoveAllConformers()
                yield resonance, acceptor, product


def generate_reaction_pairs(row, output, *, executable=None, enumerate_only=False, seed=12345, prune_rmsd=0.5):
    """Generate mapped pairs for one ligand, filtering scored isomers at 20 kcal/mol."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "reaction_pairs.csv").exists():
        raise FileExistsError(f"Reaction pairs already exist in {output}; load reaction_pairs.csv to reuse them")
    for role in ("reactant", "product", "ts_guess"):
        (output / role).mkdir(exist_ok=True)
    modes, pairs, failures = [], [], []
    for mode_index, (donors, mol) in enumerate(create_permutations(
        row["ligand_smiles"], row["metal"], row["denticity"]
    )):
        if enumerate_only:
            modes.append((mode_index, donors, None, None, mol, None))
            continue
        import rxembed

        try:
            isomers = rxembed.metal(mol, geometry=row["metal_geom"])
            for isomer_index, isomer in enumerate(isomers):
                try:
                    ensemble = rxembed.embed(isomer, n=1, seed=seed, threads=1)
                    from conformer_sampling import embedded_conformers
                    conformers = [mol for batch in embedded_conformers(ensemble, prune_rmsd) for mol in batch]
                    if not conformers:
                        raise ValueError("Embedding produced no conformer")
                    geometry = conformers[0]
                    energy, _ = run_gxtb(geometry, output / "coordination", executable,
                                         name=f"mode_{mode_index + 1:02d}_isomer_{isomer_index + 1:02d}")
                    modes.append((mode_index, donors, isomer_index, energy, mol, rxembed.cxsmiles(isomer)))
                except (RuntimeError, ValueError, OSError) as exc:
                    failures.append({"mode": mode_index, "isomer": isomer_index, "error": str(exc)})
        except (RuntimeError, ValueError) as exc:
            failures.append({"mode": mode_index, "error": str(exc)})
    scored = [mode[3] for mode in modes if mode[3] is not None]
    minimum = min(scored) if scored else None
    mode_rows = []
    for mode_index, donors, isomer, energy, mol, coordination_smiles in modes:
        relative = (energy - minimum) * HARTREE_TO_KCAL if energy is not None else None
        retained = retain_coordination_mode(donors, Chem.AddHs(Chem.MolFromSmiles(row["ligand_smiles"])), relative)
        mode_rows.append({"mode": mode_index, "donors": donors, "isomer": isomer,
                          "energy_eh": energy, "relative_kcal_mol": relative, "retained": retained,
                          "coordination_smiles": coordination_smiles})
        if not retained:
            continue
        pairs.extend(reaction_pairs_for_coordination_mode(
            row, mode_index, donors, isomer, mol, coordination_smiles))
    (output / "coordination_modes.json").write_text(json.dumps(mode_rows, indent=2) + "\n")
    (output / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    columns = ["ligand_id", "mode", "donors", "isomer", "metal", "metal_geom", "coordination_smiles", "substrate", "site", "substrate_site",
               "resonance", "acceptor_map", "reactant_smiles", "product_smiles"]
    with (output / "reaction_pairs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(pairs)
    (output / "input.json").write_text(json.dumps({**row, "enumerate_only": enumerate_only,
                                                "seed": seed, "solvent": None}, indent=2) + "\n")
    LOG.info("Ligand %s: %d modes, %d pairs, %d failures", row["ID"], len(modes), len(pairs), len(failures))
    return pairs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligands-csv", type=Path, default=Path("ligands.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--enumerate-only", action="store_true", help="Skip embedding and energy filtering")
    parser.add_argument("--max-ligands", type=int)
    parser.add_argument("--xtb", default=os.environ.get("XTB_EXE", "xtb"))
    args = parser.parse_args(argv)
    if args.max_ligands is not None and args.max_ligands < 1:
        parser.error("--max-ligands must be positive")
    executable = None if args.enumerate_only else shutil.which(args.xtb)
    if not args.enumerate_only and executable is None:
        parser.error("Set XTB_EXE or --xtb to your g-xTB capable xTB executable")
    rows = read_ligands(args.ligands_csv)[:args.max_ligands]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    success = True
    for row in rows:
        pairs = generate_reaction_pairs(row, ligand_directory(args.output_dir, row["ID"]),
                            executable=executable, enumerate_only=args.enumerate_only)
        success = bool(pairs) and success
    return 0 if rows and success else 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
