"""ORCA external g-xTB inputs, execution and NumFreq result collection.

Follow inspiration/workflow_orca.py without its campaign/store machinery.
"""

from __future__ import annotations

import math
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


def orca_input(xyz, stage, charge, multiplicity, settings, atoms=()):
    """Write ExtOpt inputs; use a separate NumFreq job after optimization."""
    parameters = shlex.join(["-e", settings["xtb_exe"], "-s", settings["solvent"]])
    external = settings["external_exe"]
    if any(c in external + parameters for c in ('"', '\n', '\r')):
        raise ValueError("ORCA external paths/parameters contain quote or newline delimiters")
    keyword = {"opt": "Opt", "tightopt": "TightOpt", "relax": "Opt", "optts": "OptTS",
               "neb": "NEB-TS", "freq": "NumFreq"}[stage]
    header = (f"! ExtOpt {keyword}\n%pal nprocs {settings['cpus']} end\n"
              f"%maxcore {int(settings['memory_mb'] * 0.75 / settings['cpus'])}\n"
              f'%method\n ProgExt "{external}"\n Ext_Params "{parameters}"\nend\n')
    if stage in ("opt", "tightopt", "relax", "optts"):
        header += f"%geom\n MaxIter {settings['max_iterations']}\n"
        if stage == "relax":
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


def run_orca(job, root, settings, stage, xyz, *, attempt=1, neb_endpoints=None):
    """Run one ORCA stage in scratch and retain its flat input/output/geometry artifacts."""
    folder = root / Path(job["input_xyz"]).parent
    stem = f"{Path(job['input_xyz']).stem.removesuffix('_input')}_run_{job['attempt']:02d}_{stage}_{attempt}"
    atoms = (job["donor"], job["hydrogen"], job["acceptor"]) if job["endpoint"] == "TS" else ()
    inp = orca_input(xyz, stage, job["charge"], job["multiplicity"], settings, atoms)
    (folder / f"{stem}.inp").write_text(inp)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env["PATH"] = str(Path(settings["orca_exe"]).parent) + os.pathsep + env.get("PATH", "")
    result = {"out_file": str((folder / f"{stem}.out").relative_to(root)), "xyz": None}
    with tempfile.TemporaryDirectory(prefix="rxflow-orca-", dir=settings.get("scratch") or None) as scratch:
        work = Path(scratch)
        (work / f"{stem}.inp").write_text(inp)
        if stage == "neb":
            for name, contents in zip(("reactant", "product", "guess"), [*neb_endpoints, xyz]):
                (work / f"{name}.xyz").write_text(contents)
                (folder / f"{stem}_{name}.xyz").write_text(contents)
        try:
            with (work / f"{stem}.out").open("w") as stdout, (work / f"{stem}.err").open("w") as stderr:
                process = subprocess.Popen([settings["orca_exe"], f"{stem}.inp"], cwd=work, env=env,
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
        finally:
            for path in work.glob(stem + "*"):
                if path.is_file() and path.suffix in (".inp", ".out", ".err", ".xyz", ".hess", ".trj"):
                    shutil.copy2(path, folder / path.name)
    return result
