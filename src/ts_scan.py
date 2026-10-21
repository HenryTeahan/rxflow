"""Score saved rxts scans and select an interior H-transfer peak for recovery."""

import json
from pathlib import Path
import re
import shutil
import time

import numpy as np
from ase.io import read, write
import rxts


def scan_inputs(job, root, source=None):
    """Locate the scan and its own coembedded pair, including historical seeds."""
    root = Path(root)
    metadata = json.loads(job.get("ts_guess_json") or "{}")
    guess = root / job["input_xyz"]
    stem = guess.stem.removesuffix("_input")
    selection = None
    if metadata.get("selection_file"):
        selection = root / metadata["selection_file"]
    elif source is not None:
        folder = Path(source) / f"ligand_{int(job['ligand_id']):02d}"
        prefix = f"pair_{job['pair_index']:03d}_seating_{job['seating']:02d}"
        selection = folder / "ts_guess" / f"{prefix}_selection.json"
    report = json.loads(selection.read_text()) if selection and selection.is_file() else {}
    if not isinstance(report, dict):
        report = {}
    if job.get("conformer") and report.get("candidates"):
        report = next((c for c in report["candidates"] if c["rank"] == job["conformer"]), {})

    def resolve(value):
        path = Path(value)
        if path.is_absolute():
            return path
        choices = [root / path]
        if selection is not None:
            choices.extend([selection.parent / path, selection.parent.parent / path])
        return next((p for p in choices if p.is_file()), choices[0])

    profile = resolve(metadata["profile_xyz"]) if metadata.get("profile_xyz") else guess.with_name(f"{stem}_path.xyz")
    if not profile.is_file() and metadata.get("previous_input_xyz"):
        previous = resolve(metadata["previous_input_xyz"])
        profile = previous.with_name(previous.stem.removesuffix("_input") + "_path.xyz")
    if not profile.is_file() and report.get("profile_xyz"):
        profile = resolve(report["profile_xyz"])
    files = {"profile_xyz": profile}
    for role in ("reactant", "product"):
        key = f"{role}_xyz"
        if metadata.get(key):
            files[key] = resolve(metadata[key])
        elif report.get(key):
            files[key] = resolve(report[key])
        else:
            # The old canonical pair is valid only for its group fallback row.
            if job.get("conformer", 0) > 0:
                raise ValueError("Candidate's matching coembedded endpoints are unavailable; regenerate its rxts guess")
            prefix = re.sub(r"_conformer_\d+$", "", stem)
            files[key] = guess.parent.parent / role / f"{prefix}_reactive.xyz"
    for name, path in files.items():
        if not path.is_file():
            raise ValueError(f"Missing rxts {name}: {path}")
    return files


def interior_peak(frames, energies, donor, hydrogen, acceptor):
    """Choose the highest local maximum with the expected local H-transfer motion."""
    energies = np.asarray(energies, dtype=float)
    if energies.shape != (len(frames),) or not np.isfinite(energies).all():
        raise ValueError("g-xTB scan has incomplete or nonfinite energies")
    peaks = [i for i in range(1, len(frames) - 1)
             if energies[i] > energies[i - 1] + 1e-8 and energies[i] > energies[i + 1] + 1e-8]
    for i in sorted(peaks, key=lambda index: (-energies[index], index)):
        def distance(frame, a, b):
            return float(np.linalg.norm(frames[frame].positions[a] - frames[frame].positions[b]))
        ch_change = distance(i + 1, donor, hydrogen) - distance(i - 1, donor, hydrogen)
        ha_change = distance(i + 1, hydrogen, acceptor) - distance(i - 1, hydrogen, acceptor)
        if ch_change > 1e-8 and ha_change < -1e-8:
            return i, dict(C_H_change_angstrom=ch_change, H_acceptor_change_angstrom=ha_change)
    raise ValueError("No interior g-xTB local maximum with C-H breaking and H-acceptor formation")


def score_ts_scan(job, root, settings, source=None):
    """Preserve inputs and energy-only scores; return a checked recovery guess."""
    from orca_hpc import molecule_with_xyz

    root = Path(root)
    folder = root / Path(job["input_xyz"]).parent / f"job_{job['id']}_run_{job['attempt']:02d}_gxtb_scan"
    folder.mkdir(exist_ok=False)
    report_file = folder / "scoring.json"
    report = dict(valid=False, method="gxtb", solvent=None, job_id=job["id"],
                  report_file=str(report_file.relative_to(root)),
                  selection="highest_interior_local_maximum", original_input_xyz=job["input_xyz"])
    started = time.monotonic()
    try:
        inputs = scan_inputs(job, root, source)
        names = dict(profile_xyz="ff_path.xyz", reactant_xyz="reactant.xyz", product_xyz="product.xyz")
        for key, file in inputs.items():
            target = folder / names[key]
            shutil.copy2(file, target)
            report[key] = str(target.relative_to(root))
        frames = read(folder / "ff_path.xyz", index=":", format="xyz")
        if len(frames) < 3:
            raise ValueError("rxts scan needs at least three frames")
        # Check every frame against the recorded atom inventory/order before scoring.
        symbols = [a.GetSymbol() for a in molecule_with_xyz(job["smiles"],
                   (root / job["input_xyz"]).read_text()).GetAtoms()]
        if any(frame.get_chemical_symbols() != symbols or not np.isfinite(frame.positions).all()
               for frame in frames):
            raise ValueError("rxts scan atom inventory/order or coordinates differ from the TS job")
        indices = [job.get(key) for key in ("donor", "hydrogen", "acceptor")]
        if (len(set(indices)) != 3 or any(type(i) is not int or not 0 <= i < len(symbols) for i in indices)
                or symbols[indices[0]] != "C" or symbols[indices[1]] != "H"):
            raise ValueError("Invalid mapped C-H transfer indices for g-xTB scan scoring")
        if settings.get("xtb_exe"):
            rxts.set_xtb(settings["xtb_exe"])
        scores = rxts.scoring.score(frames, method="gxtb", solvent=None,
                                   charge=job["charge"], multiplicity=job["multiplicity"],
                                   n_cores=settings.get("cpus", 1))
        report.update(energies_ev=scores.energies_ev, failures=scores.failures,
                      n_evaluations=scores.n_evaluations, n_gradients=scores.n_gradients)
        if scores.failures:
            raise ValueError(f"g-xTB single points failed at scan frames: {scores.failures}")
        peak, motion = interior_peak(frames, scores.energies_ev, *indices)
        guess = folder / "guess.xyz"
        write(guess, frames[peak], format="xyz",
              comment=f"g-xTB interior local peak, frame={peak}; charge={job['charge']} multiplicity={job['multiplicity']}")
        report.update(valid=True, guess_xyz=str(guess.relative_to(root)), selected_frame=peak,
                      selected_energy_ev=scores.energies_ev[peak], local_motion=motion,
                      frames=len(frames))
        return report
    except Exception as exc:
        report["error"] = str(exc)
        raise
    finally:
        report["seconds"] = time.monotonic() - started
        report_file.write_text(json.dumps(report, indent=2) + "\n")


def best_scored_candidate(connection, job, root):
    """Use the lowest electronic-energy checked peak among this group's candidates."""
    choices = []
    for row in connection.execute("""SELECT id,ts_scan_json FROM jobs WHERE ligand_id=? AND pair_index=?
            AND seating=? AND endpoint='TS' AND ts_strategy='direct' AND status='failed'""",
            (job["ligand_id"], job["pair_index"], job["seating"])):
        report = json.loads(row["ts_scan_json"] or "{}")
        if report.get("valid") and all((Path(root) / report.get(key, "")).is_file()
                for key in ("guess_xyz", "reactant_xyz", "product_xyz")):
            choices.append((report["selected_energy_ev"], row["id"], report))
    if not choices:
        return None
    _, candidate_id, report = min(choices)
    return {**report, "source_job_id": candidate_id}
