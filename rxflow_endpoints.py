"""Sample endpoint minima and create TS guesses from mapped reaction pairs."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
from numbers import Integral
import os
from pathlib import Path
import shutil

# Set the native GFN-FF thread budget before importing its numerical dependencies.
os.environ.setdefault("OMP_NUM_THREADS", "1")

from rdkit import Chem
import rxembed

from rxflow import HARTREE_TO_KCAL, ligand_directory, run_gxtb


LOG = logging.getLogger(__name__)


def mapped_molecule(smiles):
    """Read explicit H atoms and put atoms in map order for XYZ export."""
    mol = rxembed.parse_smiles(smiles, remove_hs=False)
    if mol is None:
        raise ValueError("Invalid mapped endpoint SMILES")
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    if sorted(maps) != list(range(1, mol.GetNumAtoms() + 1)):
        raise ValueError("Expected unique contiguous atom maps starting at 1")
    return Chem.RenumberAtoms(mol, sorted(range(len(maps)), key=maps.__getitem__))


def create_endpoint_isomers(pair):
    """Seat the substrate in each vacancy while retaining the screened ligand arrangement."""
    if not pair.get("coordination_smiles"):
        raise ValueError("Endpoint sampling requires scored coordination isomers; rerun screening without --enumerate-only")
    selected = rxembed.metal(mapped_molecule(pair["coordination_smiles"]), geometry=pair["metal_geom"])
    if len(selected) != 1:
        raise ValueError("Stored ligand arrangement must specify exactly one isomer")
    source = selected[0]
    reactant, product = [mapped_molecule(pair[name]) for name in ("reactant_smiles", "product_smiles")]
    if [(a.GetAtomMapNum(), a.GetAtomicNum()) for a in reactant.GetAtoms()] != [
        (a.GetAtomMapNum(), a.GetAtomicNum()) for a in product.GetAtoms()
    ]:
        raise ValueError("Endpoint atom inventories differ")
    new_donors = [b.GetBeginAtomIdx() for b in reactant.GetBonds()
                  if b.GetBondType() == Chem.BondType.DATIVE and b.GetEndAtomIdx() == source.metal
                  and b.GetBeginAtomIdx() not in source.donors]
    if len(new_donors) != 1:
        raise ValueError(f"Expected one added substrate donor, found {new_donors}")
    vertices = source.centres[0].vertices
    # A haptic site can be named by any of its face atoms in the public Isomer API.
    sites = {slot: donor if isinstance(donor, int) else donor.atoms[0]
             for slot, donor in enumerate(vertices) if donor is not None}
    for vacant, donor in enumerate(vertices):
        if donor is None:
            seating = {**sites, vacant: new_donors[0]}
            yield vacant, rxembed.Isomer(reactant, pair["metal_geom"], seating), \
                rxembed.Isomer(product, pair["metal_geom"], seating)


def connectivity_changes(mol, reference, *, distance_change=2.0):
    """Report graph changes only when the associated distance moved by more than 2 angstrom."""
    from rxembed.pipeline.metrics import connectivity, coordination_changed

    if [a.GetAtomicNum() for a in mol.GetAtoms()] != [a.GetAtomicNum() for a in reference.GetAtoms()]:
        raise ValueError("Connectivity comparison requires the same atoms in the same order")
    start, end = reference.GetConformer(), mol.GetConformer()

    def large_distance_change(i, j):
        before = (start.GetAtomPosition(i) - start.GetAtomPosition(j)).Length()
        after = (end.GetAtomPosition(i) - end.GetAtomPosition(j)).Length()
        return abs(after - before) > distance_change

    formed, broken = connectivity(mol, mol.GetConformer().GetId(), charge=Chem.GetFormalCharge(mol))
    formed = [pair for pair in formed if large_distance_change(*pair)]
    broken = [pair for pair in broken if large_distance_change(*pair)]
    coordination = []
    metals = {b.GetEndAtomIdx() for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE}
    for metal in sorted(metals):
        donors = [b.GetBeginAtomIdx() for b in mol.GetBonds()
                  if b.GetBondType() == Chem.BondType.DATIVE and b.GetEndAtomIdx() == metal]
        left, joined = coordination_changed(mol, mol.GetConformer().GetId(), metal, donors)
        left = [donor for donor in left if large_distance_change(metal, donor)]
        joined = [donor for donor in joined if large_distance_change(metal, donor)]
        if left or joined:
            coordination.append({"metal": metal, "left": left, "joined": joined})
    return {"formed": formed, "broken": broken, "coordination": coordination}


def optimize_minima(isomer, output, executable, *, prefix="endpoint", n=3, seed=12345):
    """Optimize conformers, reject graph changes with >2 A distance shifts, retain within 3 kcal/mol."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if (output / f"{prefix}_minima.json").exists() or any(output.glob(f"{prefix}_conformer_*")):
        raise FileExistsError(f"Minimum sampling {prefix} already exists in {output}; choose another prefix")
    ensemble = rxembed.embed(isomer, n=n, seed=seed, threads=1)
    records, accepted = [], []
    for index, cid in enumerate(ensemble.ids):
        try:
            geometry = ensemble.restored_mol(cid)
            energy, optimized = run_gxtb(geometry, output, executable, optimize=True,
                                         name=f"{prefix}_conformer_{index + 1:03d}")
            changes = connectivity_changes(optimized, geometry)
            record = {"conformer": index, "energy_eh": energy, "connectivity": changes,
                      "retained": False}
            if any(changes.values()):
                record["error"] = "Connectivity changed with an atom-pair distance shift greater than 2 angstrom"
            else:
                accepted.append((record, optimized))
            records.append(record)
        except (RuntimeError, ValueError, OSError) as exc:
            records.append({"conformer": index, "retained": False, "error": str(exc)})
    if accepted:
        minimum = min(record["energy_eh"] for record, _ in accepted)
        for record, mol in accepted:
            relative = (record["energy_eh"] - minimum) * HARTREE_TO_KCAL
            record.update(relative_kcal_mol=relative, retained=relative <= 3.0 + 1e-8)
            if record["retained"]:
                path = output / f"{prefix}_minimum_{record['conformer'] + 1:03d}.xyz"
                path.write_text(Chem.MolToXYZBlock(mol))
                record["xyz"] = str(path.resolve())
    (output / f"{prefix}_minima.json").write_text(json.dumps({"requested": n, "embedded": len(ensemble.ids),
                                                   "connectivity_distance_change_angstrom": 2.0,
                                                   "conformers": records}, indent=2) + "\n")
    return records


def hydrogen_transfer(reactant, product):
    """Return donor, hydrogen and acceptor indices in the current atom order."""
    transfers = []
    for atom in reactant.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        before = [n.GetIdx() for n in atom.GetNeighbors()]
        after = [n.GetIdx() for n in product.GetAtomWithIdx(atom.GetIdx()).GetNeighbors()]
        if before != after:
            if len(before) != 1 or len(after) != 1:
                raise ValueError("Transfer requires one donor and one acceptor")
            transfers.append((before[0], atom.GetIdx(), after[0]))
    if len(transfers) != 1:
        raise ValueError(f"Expected one transferring hydrogen, found {transfers}")
    return transfers[0]


def check_rxts_bonds(path, expected_transfer, *, guess=None):
    """Require rxts's reported changes to include the mapped C-H transfer."""
    donor, hydrogen, acceptor = expected_transfer
    required_breaking = tuple(sorted((donor, hydrogen)))
    required_forming = tuple(sorted((acceptor, hydrogen)))
    problems = []
    atom_count = path.reactant.GetNumAtoms()

    def reported_pairs(name):
        values = getattr(guess if guess is not None else path, name, None)
        if values is None:
            problems.append(f"rxts did not report {name}")
            return []
        result = set()
        try:
            for pair in values:
                if (len(pair) != 2 or any(isinstance(i, bool) or not isinstance(i, Integral)
                                         or not 0 <= i < atom_count for i in pair) or pair[0] == pair[1]):
                    raise ValueError
                result.add(tuple(sorted(int(i) for i in pair)))
        except (TypeError, ValueError):
            problems.append(f"rxts reported invalid {name}; expected zero-based atom-index pairs")
        return sorted(result)

    breaking = reported_pairs("breaking_bonds")
    forming = reported_pairs("forming_bonds")
    if (any(type(i) is not int or not 0 <= i < atom_count for i in expected_transfer)
            or len(set(expected_transfer)) != 3 or path.product.GetNumAtoms() != atom_count):
        problems.append("The expected H-transfer indices do not match the rxts endpoints")
    else:
        if path.reactant.GetAtomWithIdx(donor).GetAtomicNum() != 6:
            problems.append("The expected breaking bond is not C-H: its donor is not carbon")
        if path.reactant.GetAtomWithIdx(hydrogen).GetAtomicNum() != 1:
            problems.append("The expected transferring atom is not hydrogen")
        try:
            if hydrogen_transfer(path.reactant, path.product) != tuple(expected_transfer):
                problems.append("The rxts endpoints identify a different H transfer")
        except ValueError as exc:
            problems.append(f"The rxts endpoints do not identify one H transfer: {exc}")
    if required_breaking not in breaking:
        problems.append(f"rxts did not identify expected C-H bond {list(required_breaking)} breaking")
    if required_forming not in forming:
        problems.append(f"rxts did not identify expected H-acceptor bond {list(required_forming)} forming")
    return dict(valid=not problems, index_base=0,
                expected_breaking_bonds=[list(required_breaking)], expected_forming_bonds=[list(required_forming)],
                breaking_bonds=[list(pair) for pair in breaking], forming_bonds=[list(pair) for pair in forming],
                reason="; ".join(problems) if problems else "rxts identifies the expected C-H breaking and H-acceptor formation")


def check_rxts_guess_motion(path, expected_transfer):
    """Use rxts's native bond-length series to check motion around its FF guess."""
    donor, hydrogen, acceptor = expected_transfer
    peak = path.peak_index("ff")
    report = dict(valid=False, index_base=0, analysis="rxts_bond_lengths",
                  peak_frame=peak, expected_breaking_bond=sorted((donor, hydrogen)),
                  expected_forming_bond=sorted((acceptor, hydrogen)), detected_atoms=[], bond_changes=[])
    if not 0 < peak < len(path.frames) - 1:
        report["reason"] = "The rxts FF peak has no scan frames on both sides to identify local bond motion"
        return report
    required_breaking = tuple(sorted((donor, hydrogen)))
    required_forming = tuple(sorted((acceptor, hydrogen)))
    series = {(kind, tuple(sorted(pair))): [float(value) for value in values]
              for (kind, pair), values in path.bond_lengths().items()}
    if any(len(values) != len(path.frames) or not all(math.isfinite(v) for v in values)
           for values in series.values()):
        raise ValueError("rxts returned incomplete or nonfinite bond-length data")
    if ("breaking", required_breaking) not in series or ("forming", required_forming) not in series:
        report["reason"] = "rxts bond-length data do not include both expected C-H and H-acceptor bonds"
        return report
    distances = {name: [series[key][i] for i in (peak - 1, peak, peak + 1)]
                 for name, key in (("C_H", ("breaking", required_breaking)),
                                   ("H_acceptor", ("forming", required_forming)))}
    breaking_delta = distances["C_H"][2] - distances["C_H"][0]
    forming_delta = distances["H_acceptor"][2] - distances["H_acceptor"][0]
    problems = []
    if breaking_delta <= 1e-8 or forming_delta >= -1e-8:
        problems.append("Local rxts motion must lengthen C-H while shortening the same H-acceptor bond")
    changes = [dict(kind=kind, atoms=list(pair), delta_angstrom=values[peak + 1] - values[peak - 1],
                    initial_angstrom=values[peak])
               for (kind, pair), values in sorted(series.items())
               if abs(values[peak + 1] - values[peak - 1]) > 1e-8]
    report.update(valid=not problems, frame_indices=[peak - 1, peak, peak + 1],
                  frame_lambdas=[float(path.frames[i].lam) for i in (peak - 1, peak, peak + 1)],
                  local_distances_angstrom=distances,
                  local_bond_deltas_angstrom=dict(C_H=breaking_delta, H_acceptor=forming_delta),
                  detected_atoms=sorted({i for change in changes for i in change["atoms"]}),
                  bond_changes=changes,
                  reason="; ".join(problems) if problems else "Local rxts motion breaks C-H and forms the expected H-acceptor bond")
    return report


def write_ts_input(guess, reactant, product, output, *, method="r2SCAN-3c", multiplicity=1):
    """Write an ORCA gas-phase OptTS input with a numerical Hessian and frequency calculation."""
    # ORCA syntax: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations_TS.html
    donor, hydrogen, acceptor = hydrogen_transfer(reactant, product)
    coordinates = "\n".join(Chem.MolToXYZBlock(guess).splitlines()[2:])
    Path(output).write_text(
        f"! {method} OptTS NumFreq\n"
        "%geom\n  Calc_Hess true\n  NumHess true\n"
        f"  TS_Mode {{B {donor} {hydrogen}}} end\n"
        f"  TS_Active_Atoms {{{donor} {hydrogen} {acceptor}}} end\n"
        "  Modify_Internal\n"
        f"    {{B {donor} {hydrogen} A}}\n    {{B {acceptor} {hydrogen} A}}\n"
        "  end\nend\n"
        f"* xyz {Chem.GetFormalCharge(reactant)} {multiplicity}\n{coordinates}\n*\n"
    )
    return {"donor": donor, "hydrogen": hydrogen, "acceptor": acceptor, "index_base": 0,
            "numhess": True, "calcfreq": True, "method": method}


def generate_ts_guess(reactant, product, output, *, prefix="reaction", n=10, seed=12345,
                      multiplicity=1, method="r2SCAN-3c", keep=1):
    """Make checked guesses, re-coembedding twice with new seeds if necessary."""
    import rxts

    if keep < 1:
        raise ValueError("At least one TS candidate is required")
    output = Path(output)
    for role in ("reactant", "product", "ts_guess"):
        directory = output / role
        directory.mkdir(parents=True, exist_ok=True)
        if (role == "ts_guess" and (
            any(directory.glob(f"{prefix}_*")) or (directory / f"{prefix}.xyz").exists()
            or (directory / f"{prefix}.inp").exists()
        )) or (role != "ts_guess" and (directory / f"{prefix}_reactive.xyz").exists()):
            raise FileExistsError(f"Reaction {prefix} already exists in {directory}; choose another prefix")
    paired = rxembed.coembed(reactant, product, n=n, seed=seed, threads=1)
    selection = {"requested": n, "coembedded": len(paired), "ts_rmsds": paired.ts_rmsds,
                 "coembed_errors": paired.errors, "seed": seed, "reembed_retries": 2}
    selection_path = output / "ts_guess" / f"{prefix}_selection.json"

    def save():
        selection_path.write_text(json.dumps(selection, indent=2, default=str) + "\n")

    if not paired.ts_rmsds:
        save()
        raise ValueError("No coembedded endpoint pairs survived")
    ranked = sorted(paired.ts_rmsds, key=paired.ts_rmsds.get)[:keep]
    ts_directory = output / "ts_guess"
    candidates = []
    selected = None
    for rank, initial_cid in enumerate(ranked, 1):
        name = prefix if keep == 1 else f"{prefix}_conformer_{rank:03d}"
        candidate = dict(rank=rank, attempts=[])
        candidates.append(candidate)
        if keep > 1:
            selection["candidates"] = candidates
        for retry in range(3):
            attempt_seed = seed if retry == 0 else (int(seed) + 2 * (rank - 1) + retry) % 2147483647
            attempt = dict(attempt=retry + 1, seed=attempt_seed)
            candidate["attempts"].append(attempt)
            attempt_name = f"{name}_attempt_{retry + 1:02d}"
            save()
            try:
                pool = paired if retry == 0 else rxembed.coembed(
                    reactant, product, n=n, seed=attempt_seed, threads=1)
                if not pool.ts_rmsds:
                    raise ValueError("No coembedded endpoint pairs survived")
                cid = initial_cid if retry == 0 else min(pool.ts_rmsds, key=pool.ts_rmsds.get)
                attempt.update(selected_conformer_id=cid, selected_ts_rmsd=pool.ts_rmsds[cid],
                               coembed_errors=pool.errors)
                left, right = pool[pool.ids.index(cid)]
                correspondence = dict(pool.correspondence)
                right = Chem.RenumberAtoms(right, [correspondence[i] for i in range(left.GetNumAtoms())])
                transfer = hydrogen_transfer(left, right)
                # Preserve each pair and profile, including rejected attempts.
                for role, molecule in (("reactant", left), ("product", right)):
                    endpoint_file = output / role / f"{attempt_name}_reactive.xyz"
                    endpoint_file.write_text(Chem.MolToXYZBlock(molecule))
                    attempt[f"{role}_xyz"] = str(endpoint_file.resolve())
                guesser = rxts.TransitionStateGuesser(
                    engine="gfnff", gfnff_topology="graph", charge=Chem.GetFormalCharge(left),
                    multiplicity=multiplicity, n_cores=1, n_conformers=1, n_embed=0,
                    mapping=list(range(left.GetNumAtoms())),
                )
                path = guesser.find_transition_state(left, right)
                guess = path.guess("ff")
                guess_file = ts_directory / f"{attempt_name}.xyz"
                profile_file = ts_directory / f"{attempt_name}_path.xyz"
                guess.dump(guess_file)
                path.dump(profile_file)
                attempt.update(guess_xyz=str(guess_file.resolve()), profile_xyz=str(profile_file.resolve()),
                               ff_barrier_kcal_mol=path.barrier("ff"), peak_lambda=path.peak_lambda("ff"))
                # Keep canonical inspection files even when an attempt is rejected.
                shutil.copy2(guess_file, ts_directory / f"{name}.xyz")
                shutil.copy2(profile_file, ts_directory / f"{name}_path.xyz")
                attempt["rxts_bond_check"] = check_rxts_bonds(path, transfer, guess=guess)
                if not attempt["rxts_bond_check"]["valid"]:
                    raise ValueError(attempt["rxts_bond_check"]["reason"])
                attempt["rxts_motion_check"] = check_rxts_guess_motion(path, transfer)
                if not attempt["rxts_motion_check"]["valid"]:
                    raise ValueError(attempt["rxts_motion_check"]["reason"])
                indices = write_ts_input(guess.mol, path.reactant, path.product,
                                         ts_directory / f"{name}.inp", method=method, multiplicity=multiplicity)
                attempt.update(valid=True, ts_mode=indices,
                               ts_input=str((ts_directory / f"{name}.inp").resolve()))
                candidate.update({key: value for key, value in attempt.items() if key != "attempt"})
                candidate.pop("error", None)
                if selected is None:
                    selected = candidate
                    for role in ("reactant", "product"):
                        target = output / role / f"{prefix}_reactive.xyz"
                        shutil.copy2(attempt[f"{role}_xyz"], target)
                    if name != prefix:
                        shutil.copy2(ts_directory / f"{name}.xyz", ts_directory / f"{prefix}.xyz")
                break
            except (RuntimeError, ValueError, OSError) as exc:
                attempt.update(valid=False, error=str(exc))
                candidate.update({key: value for key, value in attempt.items() if key != "attempt"})
            finally:
                if selected is None:
                    selection.update(candidate)
                save()
        if selected is not None:
            selection.update(selected)
            selection.pop("error", None)
        save()
    if selected is None:
        raise ValueError("No TS guesses survived after at most 3 attempts each: " +
                         "; ".join(c["error"] for c in candidates))
    return selection


def sample_reaction_endpoints(pair, output, executable, *, pair_index=1, n=3, coembed_conformers=10,
                              seed=12345, stage="all", method="r2SCAN-3c"):
    """Sample minima and the lowest-TS-RMSD reactive pair into one ligand directory."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    name = f"pair_{pair_index:03d}"
    input_path = output / f"{name}_{stage}_input.json"
    if input_path.exists():
        raise FileExistsError(f"Sampling input already exists: {input_path}")
    input_path.write_text(json.dumps({**pair, "n_conformers": n, "coembed_conformers": coembed_conformers,
                                     "seed": seed, "solvent": None, "stage": stage}, indent=2) + "\n")
    records = []
    for vacant, reactant, product in create_endpoint_isomers(pair):
        prefix = f"{name}_seating_{vacant + 1:02d}"
        (output / f"{prefix}_endpoints.json").write_text(json.dumps({"reactant": rxembed.cxsmiles(reactant),
            "product": rxembed.cxsmiles(product), "substrate": pair.get("substrate")}, indent=2) + "\n")
        if stage in ("all", "minima"):
            for role, isomer in (("reactant", reactant), ("product", product)):
                try:
                    minima = optimize_minima(isomer, output / role, executable, prefix=prefix, n=n, seed=seed)
                    records.append({"seating": vacant, "stage": role,
                                    "retained": sum(r["retained"] for r in minima)})
                except (RuntimeError, ValueError, OSError) as exc:
                    records.append({"seating": vacant, "stage": role, "error": str(exc)})
        if stage in ("all", "ts"):
            try:
                guess = generate_ts_guess(reactant, product, output, prefix=prefix,
                                          n=coembed_conformers, seed=seed, method=method)
                records.append({"seating": vacant, "stage": "ts", "guesses": 1,
                                "selected_ts_rmsd": guess["selected_ts_rmsd"]})
            except (RuntimeError, ValueError, OSError) as exc:
                records.append({"seating": vacant, "stage": "ts", "error": str(exc)})
    (output / f"{name}_{stage}_summary.json").write_text(json.dumps(records, indent=2) + "\n")
    LOG.info("Pair %s %s: %s", pair["ligand_id"], pair["site"], records)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pairs_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--conformers", type=int, default=3, help="Minimum conformers per endpoint")
    parser.add_argument("--coembed-conformers", type=int, default=10)
    parser.add_argument("--max-pairs", type=int)
    parser.add_argument("--stage", choices=("all", "minima", "ts"), default="all")
    parser.add_argument("--orca-method", default="r2SCAN-3c")
    parser.add_argument("--xtb", default=os.environ.get("XTB_EXE", "xtb"))
    args = parser.parse_args(argv)
    if min(args.conformers, args.coembed_conformers) < 1 or (args.max_pairs is not None and args.max_pairs < 1):
        parser.error("Counts must be positive")
    executable = shutil.which(args.xtb)
    if args.stage != "ts" and executable is None:
        parser.error("Set XTB_EXE or --xtb to your g-xTB capable xTB executable")
    with args.pairs_csv.open(newline="") as handle:
        pairs = list(csv.DictReader(handle))[:args.max_pairs]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    successful = bool(pairs)
    for index, pair in enumerate(pairs, 1):
        directory = ligand_directory(args.output_dir, pair["ligand_id"])
        directory.mkdir(parents=True, exist_ok=True)
        try:
            records = sample_reaction_endpoints(pair, directory, executable, pair_index=index,
                n=args.conformers, coembed_conformers=args.coembed_conformers, stage=args.stage, method=args.orca_method)
        except (RuntimeError, ValueError, OSError) as exc:
            records = [{"stage": "pair", "error": str(exc)}]
            (directory / f"pair_{index:03d}_{args.stage}_failure.json").write_text(json.dumps(records, indent=2) + "\n")
            LOG.error("Pair %s: %s", index, exc)
        successful = successful and bool(records) and all(
            not r.get("error") and r.get("retained", r.get("guesses", 0)) > 0 for r in records)
    return 0 if successful else 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    raise SystemExit(main())
