"""ORCA external g-xTB inputs, execution and NumFreq result collection.

Follow inspiration/workflow_orca.py without its campaign/store machinery.
"""

from __future__ import annotations

import math
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import tempfile

from rdkit import Chem

from rxflow_endpoints import mapped_molecule

NUMBER = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)"
ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY(?:\s+\(From external program\))?\s+" + NUMBER)
FREQUENCY = re.compile(r"^\s*(\d+):\s*" + NUMBER + r"\s+cm\*\*-1", re.MULTILINE)
TS_MODE_POLICY_VERSION = 1


def molecule_with_xyz(smiles, xyz):
    """Restore coordinates onto the mapped graph without changing atom order."""
    mol = mapped_molecule(smiles)
    coordinates = Chem.MolFromXYZBlock(xyz)
    if coordinates is None or [a.GetAtomicNum() for a in coordinates.GetAtoms()] != [
        a.GetAtomicNum() for a in mol.GetAtoms()
    ]:
        raise ValueError("XYZ atom inventory/order differs from the mapped graph")
    if not all(math.isfinite(v) for position in coordinates.GetConformer().GetPositions() for v in position):
        raise ValueError("XYZ contains nonfinite coordinates")
    mol.AddConformer(Chem.Conformer(coordinates.GetConformer()), assignId=True)
    return mol


def ts_preoptimization_atoms(job):
    """Freeze the transfer atoms and the shortest acceptor-to-metal graph path."""
    mol = mapped_molecule(job["smiles"])
    reactive = [job.get(name) for name in ("donor", "hydrogen", "acceptor")]
    if (any(type(index) is not int or not 0 <= index < mol.GetNumAtoms() for index in reactive)
            or len(set(reactive)) != 3):
        raise ValueError("TS preoptimization requires distinct donor, hydrogen and acceptor indices")
    metals = {bond.GetEndAtomIdx() for bond in mol.GetBonds()
              if bond.GetBondType() == Chem.BondType.DATIVE}
    if len(metals) != 1:
        raise ValueError("TS preoptimization requires one central coordination metal")
    metal, = metals
    path = Chem.GetShortestPath(mol, job["acceptor"], metal)
    if not path:
        raise ValueError("The H acceptor has no graph path to the central metal")
    return tuple(sorted(set(reactive) | set(path)))


def orca_input(xyz, stage, charge, multiplicity, settings, atoms=()):
    """Write ExtOpt inputs; use a separate NumFreq job after optimization."""
    if stage == "singlepoint":
        count = int(xyz.splitlines()[0])
        return ("! PBE0 def2-TZVP D4 SP\n"
                f"%pal nprocs {settings['cpus']} end\n"
                f"%maxcore {int(settings['memory_mb'] * 0.75 / settings['cpus'])}\n"
                '%cpcm\n smd true\n SMDsolvent "water"\nend\n'
                f"* xyz {charge} {multiplicity}\n"
                + "\n".join(xyz.splitlines()[2:2 + count]) + "\n*\n")
    parameters = (settings["external_params"] if settings.get("method") == "uma" else
                  shlex.join(["-e", settings["xtb_exe"], "-s", settings["solvent"]]))
    external = settings["external_exe"]
    if any(c in external + parameters for c in ('"', '\n', '\r')):
        raise ValueError("ORCA external paths/parameters contain quote or newline delimiters")
    keyword = {"opt": "Opt", "verytightopt": "VeryTightOpt", "relax": "Opt", "ts_preopt": "Opt",
               "neb_reactant_opt": "Opt", "neb_product_opt": "Opt", "optts": "OptTS",
               "neb": "NEB-TS", "freq": "NumFreq", "external_sp": "SP"}[stage]
    header = (f"! ExtOpt {keyword}\n%pal nprocs {settings['cpus']} end\n"
              f"%maxcore {int(settings['memory_mb'] * 0.75 / settings['cpus'])}\n"
              f'%method\n ProgExt "{external}"\n Ext_Params "{parameters}"\nend\n')
    if stage in ("relax", "ts_preopt", "optts"):
        header += "%geom\n"
        if stage == "ts_preopt":
            header += " MaxIter 50\n"
        if stage in ("relax", "ts_preopt"):
            header += " Constraints\n" + "".join(f"  {{C {i} C}}\n" for i in atoms) + " end\n"
        if stage == "optts":
            donor, hydrogen, acceptor = atoms
            header += (" Calc_Hess true\n NumHess true\n"
                       f" TS_Mode {{B {donor} {hydrogen}}} end\n"
                       f" TS_Active_Atoms {{{donor} {hydrogen} {acceptor}}} end\n"
                       f" Modify_Internal\n  {{B {donor} {hydrogen} A}}\n"
                       f"  {{B {acceptor} {hydrogen} A}}\n end\n")
        header += "end\n"
    if stage == "neb":
        return (header + '%neb\n NEB_END_XYZFILE "product.xyz"\n NEB_TS_XYZFILE "guess.xyz"\nend\n'
                f"* xyzfile {charge} {multiplicity} reactant.xyz\n")
    count = int(xyz.splitlines()[0])
    return header + f"* xyz {charge} {multiplicity}\n" + "\n".join(xyz.splitlines()[2:2 + count]) + "\n*\n"


def collect_numfreq(text, atom_count):
    """Read reference electronic energy, Gibbs energy and the complete last frequency table."""
    if "ORCA TERMINATED NORMALLY" not in text:
        raise ValueError("ORCA NumFreq did not terminate normally")
    if "VIBRATIONAL FREQUENCIES" not in text:
        raise ValueError("ORCA output has no frequency calculation")
    section = text.rsplit("VIBRATIONAL FREQUENCIES", 1)[1].split("NORMAL MODES", 1)[0]
    modes = FREQUENCY.findall(section)
    if len(modes) != 3 * atom_count or {int(i) for i, value in modes} != set(range(3 * atom_count)):
        raise ValueError("ORCA NumFreq has an incomplete frequency table")
    thermal = text.rsplit("THERMOCHEMISTRY", 1)[-1]
    electronic = re.findall(r"Electronic energy\s*\.\.\.\s*" + NUMBER + r"\s+Eh", thermal)
    gibbs = re.findall(r"Final Gibbs free energy\s*\.\.\.\s*" + NUMBER + r"\s+Eh", thermal, re.IGNORECASE)
    if not electronic or not gibbs:
        raise ValueError("ORCA NumFreq is missing its thermochemical electronic/Gibbs energy")
    convert = lambda value: float(value.replace("D", "E").replace("d", "e"))
    frequencies = [convert(value) for i, value in sorted(modes, key=lambda item: int(item[0]))]
    energy, free_energy = convert(electronic[-1]), convert(gibbs[-1])
    if not all(math.isfinite(value) for value in [energy, free_energy, *frequencies]):
        raise ValueError("ORCA NumFreq returned a nonfinite result")
    return {"electronic_energy_eh": energy, "gibbs_free_energy_eh": free_energy,
            "frequencies": frequencies, "imaginary_count": sum(value < 0 for value in frequencies)}


def collect_singlepoint(text):
    """Accept only a normally terminated, converged SCF with a finite final energy."""
    if "ORCA TERMINATED NORMALLY" not in text:
        raise ValueError("ORCA single point did not terminate normally")
    upper = text.upper()
    if ("SCF CONVERGED AFTER" not in upper or "SCF NOT CONVERGED" in upper
            or "SCF DID NOT CONVERGE" in upper):
        raise ValueError("ORCA single-point SCF did not converge")
    energies = ENERGY.findall(text)
    if not energies:
        raise ValueError("ORCA single point has no final electronic energy")
    energy = float(energies[-1].replace("D", "E").replace("d", "e"))
    if not math.isfinite(energy):
        raise ValueError("ORCA single point returned a nonfinite electronic energy")
    return {"electronic_energy_eh": energy}


def ts_mode_fingerprint(output, job):
    """Invalidate a cached check when the output, atom mapping or analysis changes."""
    required = [job.get(name) for name in ("donor", "hydrogen", "acceptor")]
    if any(type(index) is not int or index < 0 for index in required) or len(set(required)) != 3:
        raise ValueError("TS mode check requires distinct zero-based donor, hydrogen and acceptor indices")
    payload = dict(policy_version=TS_MODE_POLICY_VERSION, graphrc_version=version("graphrc"),
                   xyzgraph_version=version("xyzgraph"), expected_atoms=sorted(required),
                   smiles=job["smiles"], charge=job["charge"], multiplicity=job["multiplicity"],
                   output_file=str(Path(output).resolve()),
                   output_sha256=hashlib.sha256(Path(output).read_bytes()).hexdigest())
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def check_ts_mode(output, job):
    """Require the imaginary mode's graphRC coordinates to contain every TS atom.

    Additional atoms are allowed. graphRC's default significance thresholds are
    used, and no trajectories or displaced geometries are written to disk.
    """
    from graphrc import run_vib_analysis

    fingerprint = ts_mode_fingerprint(output, job)
    required = sorted(job[name] for name in ("donor", "hydrogen", "acceptor"))
    molecule = mapped_molecule(job["smiles"])
    if required[-1] >= molecule.GetNumAtoms():
        raise ValueError("TS mode indices are outside the mapped molecular graph")
    text = Path(output).read_text(errors="replace")
    thermochemistry = collect_numfreq(text, molecule.GetNumAtoms())
    frequencies = thermochemistry["frequencies"]
    report = dict(policy_version=TS_MODE_POLICY_VERSION, graphrc_version=version("graphrc"),
                  fingerprint=fingerprint, output_file=str(Path(output).resolve()), index_base=0,
                  expected_atoms=required, detected_atoms=[], missing_atoms=required,
                  valid=False, imaginary_count=thermochemistry["imaginary_count"])
    if thermochemistry["imaginary_count"] != 1:
        report["reason"] = f"Expected 1 imaginary frequency, found {thermochemistry['imaginary_count']}"
        return report
    orca_mode = next(i for i, value in enumerate(frequencies) if value < 0)
    nonzero = [(i, value) for i, value in enumerate(frequencies) if abs(value) > 1e-5]
    mode = next((i for i, (original, value) in enumerate(nonzero) if original == orca_mode), None)
    if mode is None or mode + len(frequencies) - len(nonzero) != orca_mode:
        raise ValueError("The imaginary mode cannot be mapped to graphRC's non-zero mode numbering")
    analysis = run_vib_analysis(str(output), mode=mode, charge=job["charge"], multiplicity=job["multiplicity"],
                                save_trajectory=False, save_displacement=False, print_output=False)
    trajectory = analysis["trajectory"]
    if trajectory["frames"][0]["symbols"] != [atom.GetSymbol() for atom in molecule.GetAtoms()]:
        raise ValueError("Normal-mode atom inventory/order differs from the mapped TS graph")
    if abs(float(trajectory["frequencies"][mode]) - frequencies[orca_mode]) > 1e-4:
        raise ValueError("graphRC analyzed a different frequency than the imaginary TS mode")
    changes = {name: [dict(atoms=list(atoms), change=float(values[0]), initial=float(values[1]))
                      for atoms, values in sorted(analysis["vibrational"][name].items())]
               for name in ("bond_changes", "angle_changes", "dihedral_changes")}
    detected = sorted({index for coordinates in changes.values() for coordinate in coordinates
                       for index in coordinate["atoms"]})
    missing = sorted(set(required) - set(detected))
    report.update(valid=not missing, detected_atoms=detected, missing_atoms=missing,
                  orca_mode_index=orca_mode, graphrc_mode_index=mode,
                  frequency_cm1=frequencies[orca_mode], coordinate_changes=changes,
                  thresholds_relaxed=analysis["vibrational"].get("thresholds_relaxed", False))
    report["reason"] = ("All required TS atoms participate" if not missing else
                        f"TS imaginary mode is missing required atoms {missing}; expected {required}, detected {detected}")
    return report


def run_orca(job, root, settings, stage, xyz, *, attempt=1, neb_endpoints=None):
    """Run one ORCA stage in scratch and retain its flat input/output/geometry artifacts."""
    folder = root / Path(job["input_xyz"]).parent
    stem = f"{Path(job['input_xyz']).stem.removesuffix('_input')}_run_{job['attempt']:02d}_{stage}_{attempt}"
    atoms = ((job["donor"], job["hydrogen"], job["acceptor"])
             if job["endpoint"] == "TS" and stage != "singlepoint" else ())
    if stage == "ts_preopt":
        atoms = ts_preoptimization_atoms(job)
    inp = orca_input(xyz, stage, job["charge"], job["multiplicity"], settings, atoms)
    (folder / f"{stem}.inp").write_text(inp)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    orca_directory = str(Path(settings["orca_exe"]).parent)
    env["PATH"] = orca_directory + os.pathsep + env.get("PATH", "")
    if stage == "singlepoint":
        env["LD_LIBRARY_PATH"] = orca_directory + (
            os.pathsep + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    if settings.get("mpirun_exe"):
        mpi_bin = Path(settings["mpirun_exe"]).parent
        env["PATH"] = str(mpi_bin) + os.pathsep + env["PATH"]
        mpi_libraries = [str(mpi_bin.parent / "lib")]
        if (mpi_bin.parent / "lib64").is_dir():
            mpi_libraries.append(str(mpi_bin.parent / "lib64"))
        env["LD_LIBRARY_PATH"] = os.pathsep.join(mpi_libraries +
            ([env["LD_LIBRARY_PATH"]] if env.get("LD_LIBRARY_PATH") else []))
    result = {"out_file": str((folder / f"{stem}.out").relative_to(root)), "xyz": None}
    scratch_root = settings.get("scratch")
    if stage == "singlepoint" and not scratch_root:
        scratch_root = env.get("SCRATCH") or env.get("SLURM_TMPDIR")
    if scratch_root:
        Path(scratch_root).mkdir(parents=True, exist_ok=True)
    command = [settings["orca_exe"], f"{stem}.inp"]
    if stage == "singlepoint" and settings["cpus"] > 1:
        command.append("--use-hwthread-cpus")
    with tempfile.TemporaryDirectory(prefix="rxflow-orca-", dir=scratch_root or None) as scratch:
        work = Path(scratch)
        (work / f"{stem}.inp").write_text(inp)
        if stage == "neb":
            for name, contents in zip(("reactant", "product", "guess"), [*neb_endpoints, xyz]):
                (work / f"{name}.xyz").write_text(contents)
                (folder / f"{stem}_{name}.xyz").write_text(contents)
        try:
            with (work / f"{stem}.out").open("w") as stdout, (work / f"{stem}.err").open("w") as stderr:
                process = subprocess.Popen(command, cwd=work, env=env,
                    stdout=stdout, stderr=stderr, start_new_session=True)
                try:
                    code = process.wait(timeout=settings["timeout_seconds"])
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    code = -1  # Retain partial geometry/output for the next OptTS attempt.
                except BaseException:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    raise
            text = (work / f"{stem}.out").read_text(errors="replace")
            normal = code == 0 and "ORCA TERMINATED NORMALLY" in text
            result.update(normal=normal, converged=normal and (
                "THE OPTIMIZATION HAS CONVERGED" in text.upper() or "OPTIMIZATION CONVERGED" in text.upper()))
            result["max_cycles_reached"] = (
                "THE OPTIMIZATION DID NOT CONVERGE BUT REACHED THE MAXIMUM NUMBER OF" in text.upper()
                and "OPTIMIZATION CYCLES" in text.upper()
            )
            energies = ENERGY.findall(text)
            energy = float(energies[-1].replace("D", "E").replace("d", "e")) if energies else None
            result["electronic_energy_eh"] = energy if energy is not None and math.isfinite(energy) else None
            final = work / (f"{stem}_NEB-TS_converged.xyz" if stage == "neb" else f"{stem}.xyz")
            if final.exists():
                molecule_with_xyz(job["smiles"], final.read_text())
                result["xyz"] = final.read_text()
            if stage == "neb":
                result["converged"] = normal and final.exists()
            result["converged"] = result["converged"] and result["xyz"] is not None
            if stage == "freq":
                result.update(collect_numfreq(text, int(xyz.splitlines()[0])))
            if stage in ("singlepoint", "external_sp"):
                if not normal:
                    raise ValueError("ORCA single point failed or timed out")
                result.update(collect_singlepoint(text))
        finally:
            for path in work.glob(stem + "*"):
                if path.is_file() and path.suffix in (".inp", ".out", ".err", ".xyz", ".hess", ".trj"):
                    shutil.copy2(path, folder / path.name)
    return result
