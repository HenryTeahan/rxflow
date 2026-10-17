"""Prepare a flat ORCA job database and drain it with independent HPC workers."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tomllib

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from rdkit import Chem

from orca_hpc import molecule_with_xyz, run_orca
from rxflow import HARTREE_TO_KCAL, SUBSTRATE, ligand_directory, read_ligands
from rxflow_endpoints import connectivity_changes, mapped_molecule


def connect_database(path):
    """Open the single jobs table using short transactions, as in Projects/workflow."""
    connection = sqlite3.connect(path, timeout=60, isolation_level=None)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA busy_timeout=60000")
    # Keep rollback journaling: WAL needs shared memory and is unsuitable across HPC nodes.
    connection.execute("PRAGMA journal_mode=DELETE")
    connection.execute("""CREATE TABLE IF NOT EXISTS jobs (
        id INTEGER PRIMARY KEY, ligand_id TEXT NOT NULL, pair_index INTEGER NOT NULL,
        seating INTEGER NOT NULL, endpoint TEXT NOT NULL CHECK(endpoint IN ('reactant','product','TS')),
        conformer INTEGER NOT NULL, status TEXT NOT NULL DEFAULT 'pending', stage TEXT NOT NULL DEFAULT 'opt',
        attempt INTEGER NOT NULL DEFAULT 0, input_xyz TEXT NOT NULL, xyz_file TEXT,
        smiles TEXT NOT NULL, charge INTEGER NOT NULL, multiplicity INTEGER NOT NULL DEFAULT 1,
        donor INTEGER, hydrogen INTEGER, acceptor INTEGER, settings_json TEXT NOT NULL,
        electronic_energy_eh REAL, gibbs_free_energy_eh REAL, frequencies_json TEXT,
        imaginary_count INTEGER, retained INTEGER, out_file TEXT, worker_id TEXT, error TEXT,
        updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
        UNIQUE(ligand_id,pair_index,seating,endpoint,conformer))""")
    # Upgrade already prepared databases without replacing rows or calculated results.
    columns = {row["name"] for row in connection.execute("PRAGMA table_info(jobs)")}
    for name, kind in {"substrate": "TEXT", "binding_permutation": "INTEGER",
                       "binding_mode": "INTEGER", "coordination_isomer": "INTEGER",
                       "donor_maps": "TEXT", "substrate_site": "INTEGER", "resonance": "INTEGER"}.items():
        if name not in columns:
            connection.execute(f"ALTER TABLE jobs ADD COLUMN {name} {kind}")
    return connection


def update_job(connection, job_id, **values):
    assignments = ",".join(f"{key}=?" for key in values)
    connection.execute(f"UPDATE jobs SET {assignments},updated_at=CURRENT_TIMESTAMP WHERE id=?",
                       [*values.values(), job_id])


def prepare_jobs(connection, source, root, settings):
    """Import each usable seed; retain missing/invalid endpoint and TS seeds as failed rows."""
    added, groups = 0, 0
    pair_tables = {}
    for descriptor in sorted(source.glob("ligand_*/pair_*_seating_*_endpoints.json")):
        prefix = descriptor.stem.removesuffix("_endpoints")
        selection_path = descriptor.parent / "ts_guess" / f"{prefix}_selection.json"
        selection = {}
        try:
            if selection_path.exists():
                selection = json.loads(selection_path.read_text())
        except json.JSONDecodeError:
            selection = {"error": "Incomplete TS selection JSON"}
        pair, seating = map(int, re.fullmatch(r"pair_(\d+)_seating_(\d+)", prefix).groups())
        ligand_input = json.loads((descriptor.parent / "input.json").read_text())
        ligand_id = str(ligand_input["ID"])
        if descriptor.parent not in pair_tables:
            with (descriptor.parent / "reaction_pairs.csv").open(newline="") as handle:
                pair_tables[descriptor.parent] = list(csv.DictReader(handle))
        reaction = pair_tables[descriptor.parent][pair - 1]
        graphs = json.loads(descriptor.read_text())
        # Preserve the selected coembedded endpoint pair for a possible NEB fallback.
        for role in ("reactant", "product"):
            reactive = descriptor.parent / role / f"{prefix}_reactive.xyz"
            if reactive.exists():
                xyz = reactive.read_text()
                molecule_with_xyz(graphs[role], xyz)
                target = root / descriptor.parent.name / role / reactive.name
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists() and target.read_text() != xyz:
                    raise ValueError(f"Existing coembedded input differs: {target}")
                if not target.exists():
                    target.write_text(xyz)
        summary_path = descriptor.parent / "sampling.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        errors = [e for e in summary.get("errors", []) if e.get("pair") == pair and e.get("seating") == seating]
        substrate = reaction.get("substrate") or ligand_input.get("substrate") or SUBSTRATE
        if graphs.get("substrate") and graphs["substrate"] != substrate:
            raise ValueError(f"Substrate metadata differs from sampled endpoints: {descriptor}")
        site = reaction.get("substrate_site")
        if not site:  # Historical alpha/beta labels: recover the actual bound substrate carbon.
            offset = mapped_molecule(reaction["coordination_smiles"]).GetNumAtoms()
            reactant = mapped_molecule(reaction["reactant_smiles"])
            site, = [b.GetBeginAtomIdx() - offset + 1 for b in reactant.GetBonds()
                     if b.GetBondType() == Chem.BondType.DATIVE and b.GetBeginAtomIdx() >= offset]
        metadata = dict(substrate=substrate, binding_mode=int(reaction["mode"]) + 1,
                        coordination_isomer=int(reaction["isomer"]) + 1,
                        donor_maps=json.dumps([i + 1 for i in json.loads(reaction["donors"])]),
                        substrate_site=int(site), resonance=int(reaction["resonance"]))
        binding = connection.execute("""SELECT binding_permutation FROM jobs WHERE ligand_id=?
            AND binding_mode=? AND coordination_isomer=? AND seating=? LIMIT 1""",
            (ligand_id, metadata["binding_mode"], metadata["coordination_isomer"], seating)).fetchone()
        metadata["binding_permutation"] = binding[0] if binding else connection.execute(
            "SELECT coalesce(max(binding_permutation),0)+1 FROM jobs WHERE ligand_id=?", (ligand_id,)).fetchone()[0]
        inputs = []
        for role in ("reactant", "product"):
            available = [int(re.search(r"conformer_(\d+)_input", path.name)[1])
                         for path in (descriptor.parent / role).glob(f"{prefix}_conformer_*_input.xyz")]
            count = max(summary.get("endpoint_conformers", 1), max(available, default=1))
            inputs.extend((role, conformer, descriptor.parent / role / f"{prefix}_conformer_{conformer:03d}_input.xyz")
                          for conformer in range(1, count + 1))
        inputs.append(("TS", 0, descriptor.parent / "ts_guess" / f"{prefix}.xyz"))
        for role, conformer, source_xyz in inputs:
            smiles = graphs["reactant" if role == "TS" else role]
            mol = mapped_molecule(smiles)
            error, xyz = None, None
            try:
                if role == "TS" and ("ts_input" not in selection or selection.get("error")):
                    raise ValueError(selection.get("error") or "No completed TS guess")
                xyz = source_xyz.read_text()
                molecule_with_xyz(smiles, xyz)
            except (ValueError, OSError) as exc:
                detail = next((e["error"] for e in errors if e.get("endpoint", role) == role), str(exc))
                error = f"{role} sampling failed: {detail}"
                xyz = None
            folder = root / descriptor.parent.name / ("ts_guess" if role == "TS" else role)
            folder.mkdir(parents=True, exist_ok=True)
            target = folder / (f"{prefix}_input.xyz" if role == "TS" else source_xyz.name)
            key = (ligand_id, pair, seating, role, conformer)
            previous = connection.execute("SELECT * FROM jobs WHERE ligand_id=? AND pair_index=? AND seating=? "
                                          "AND endpoint=? AND conformer=?", key).fetchone()
            protocol = json.dumps(settings, sort_keys=True)
            if previous and (previous["settings_json"] != protocol or previous["smiles"] != smiles):
                raise ValueError("Existing jobs have different chemistry/settings; use a new HPC results directory")
            if previous and any(previous[name] is not None and previous[name] != value
                                for name, value in metadata.items()):
                raise ValueError("Existing jobs have different substrate/binding metadata; use a new results directory")
            if xyz is not None and target.exists():
                if target.read_text() != xyz:
                    raise ValueError(f"Existing input differs: {target}")
            elif xyz is not None:
                target.write_text(xyz)
            atoms = selection.get("ts_mode", {})
            cursor = connection.execute("""INSERT OR IGNORE INTO jobs
                (ligand_id,pair_index,seating,endpoint,conformer,input_xyz,smiles,charge,multiplicity,
                 donor,hydrogen,acceptor,settings_json,substrate,binding_mode,coordination_isomer,
                 donor_maps,substrate_site,resonance,binding_permutation,status,stage,error,retained)
                 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (*key, str(target.relative_to(root)), smiles, Chem.GetFormalCharge(mol), 1,
                 atoms.get("donor"), atoms.get("hydrogen"), atoms.get("acceptor"), protocol, *metadata.values(),
                 "failed" if error else "pending", "sampling" if error else "opt", error, 0 if error else None))
            added += cursor.rowcount
            if previous and any(previous[name] != value for name, value in metadata.items()):
                update_job(connection, previous["id"], **metadata)
            if previous and previous["status"] == "failed" and previous["stage"] == "sampling":
                update_job(connection, previous["id"], status="failed" if error else "pending",
                           stage="sampling" if error else "opt", error=error, retained=0 if error else None,
                           donor=atoms.get("donor"), hydrogen=atoms.get("hydrogen"), acceptor=atoms.get("acceptor"))
        groups += 1
    return added, groups


def claim_job(connection):
    """Atomically claim a conformer; a TS waits until its endpoint jobs have finished."""
    connection.execute("BEGIN IMMEDIATE")
    try:
        job = connection.execute("""UPDATE jobs SET status='running',attempt=attempt+1,error=NULL,
            worker_id=?,updated_at=CURRENT_TIMESTAMP WHERE id=(
              SELECT j.id FROM jobs j WHERE j.status='pending' AND (j.endpoint!='TS' OR NOT EXISTS (
                SELECT 1 FROM jobs e WHERE e.ligand_id=j.ligand_id AND e.pair_index=j.pair_index
                  AND e.seating=j.seating AND e.endpoint!='TS' AND e.status IN ('pending','running')))
              ORDER BY (j.endpoint='TS'),j.id LIMIT 1) RETURNING *""",
            (os.environ.get("SLURM_JOB_ID", "local") + f":{os.getpid()}",)).fetchone()
        connection.commit()
        return dict(job) if job else None
    except BaseException:
        connection.rollback()
        raise


def select_endpoints(connection, job):
    """Apply the energy window inside the transaction that completes an endpoint."""
    group = (job["ligand_id"], job["pair_index"], job["seating"], job["endpoint"])
    connection.execute("""UPDATE jobs SET retained=(electronic_energy_eh <= (
        SELECT min(electronic_energy_eh) FROM jobs WHERE ligand_id=? AND pair_index=?
          AND seating=? AND endpoint=? AND status='succeeded') + ?),updated_at=CURRENT_TIMESTAMP
        WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint=? AND status='succeeded'""",
        (*group, 3.0 / HARTREE_TO_KCAL + 1e-12, *group))


def run_job(connection, job, root, executables):
    """Optimize endpoints or try direct OptTS before relaxation, then separate NumFreq."""
    settings = {**json.loads(job["settings_json"]), **executables}
    settings.setdefault("opt_attempts", settings.get("ts_attempts", 2))

    def calculate(stage, attempt=1, endpoints=None):
        update_job(connection, job["id"], stage=stage)
        result = run_orca(job, root, settings, stage, current, attempt=attempt, neb_endpoints=endpoints)
        update_job(connection, job["id"], out_file=result["out_file"])
        if result.get("electronic_energy_eh") is not None:
            update_job(connection, job["id"], electronic_energy_eh=result["electronic_energy_eh"])
        return result

    try:
        original = (root / job["input_xyz"]).read_text()
        current = original
        if job["endpoint"] != "TS":
            result = calculate("opt", 1)
            for attempt in range(2, settings["opt_attempts"] + 1):
                if result["converged"] or not result.get("max_cycles_reached") or not result.get("xyz"):
                    break
                current = result["xyz"]
                result = calculate("opt", attempt)
            if not result["converged"] or not result["xyz"]:
                raise RuntimeError("ORCA endpoint optimization did not converge")
            current = result["xyz"]
            changes = connectivity_changes(molecule_with_xyz(job["smiles"], current),
                                           molecule_with_xyz(job["smiles"], original))
            if any(changes.values()):
                raise ValueError(f"Graph change with a distance shift exceeding 2 angstrom: {changes}")
        else:
            result = calculate("optts", 1)
            if result["converged"] and result["xyz"]:
                current = result["xyz"]
            else:
                # Relax the original guess, not the failed direct OptTS geometry.
                current = original
                result = calculate("relax", 1)
                for attempt in range(2, settings["ts_attempts"] + 1):
                    if result["converged"] or not result.get("max_cycles_reached") or not result.get("xyz"):
                        break
                    current = result["xyz"]
                    result = calculate("relax", attempt)
                if not result["converged"] or not result["xyz"]:
                    raise RuntimeError("Constrained TS relaxation did not converge")
                current = result["xyz"]
                for attempt in range(2, settings["ts_attempts"] + 2):
                    result = calculate("optts", attempt)
                    current = result["xyz"] or current
                    if result["converged"] and result["xyz"]:
                        break
            if not result["converged"] or not result["xyz"]:
                if not settings["neb_fallback"]:
                    raise RuntimeError("OptTS attempts did not converge")
                endpoints = []
                ts_folder = Path(job["input_xyz"]).parent
                prefix = Path(job["input_xyz"]).name.removesuffix("_input.xyz")
                for role in ("reactant", "product"):
                    path = root / ts_folder.parent / role / f"{prefix}_reactive.xyz"
                    xyz = path.read_text()
                    endpoint = connection.execute("""SELECT smiles FROM jobs WHERE ligand_id=? AND pair_index=?
                        AND seating=? AND endpoint=? LIMIT 1""",
                        (job["ligand_id"], job["pair_index"], job["seating"], role)).fetchone()
                    if endpoint is None:
                        raise ValueError(f"NEB has no sampled {role} graph")
                    molecule_with_xyz(endpoint["smiles"], xyz)
                    endpoints.append(xyz)
                # NEB uses the original rxts TS guess and its coembedded endpoint pair.
                molecule_with_xyz(job["smiles"], original)
                current = original
                result = calculate("neb", endpoints=endpoints)
                if not result["converged"] or not result["xyz"]:
                    raise RuntimeError("NEB-TS fallback did not converge")
                current = result["xyz"]
        path = root / Path(job["input_xyz"]).parent / f"job_{job['id']}_run_{job['attempt']:02d}_optimized.xyz"
        path.write_text(current)
        update_job(connection, job["id"], xyz_file=str(path.relative_to(root)))
        result = calculate("freq")
        if job["endpoint"] != "TS" and result["imaginary_count"] > 0:
            tight = calculate("tightopt")
            if not tight["converged"] or not tight["xyz"]:
                raise RuntimeError("ORCA TightOpt retry did not converge")
            current = tight["xyz"]
            changes = connectivity_changes(molecule_with_xyz(job["smiles"], current),
                                           molecule_with_xyz(job["smiles"], original))
            if any(changes.values()):
                raise ValueError(f"Graph change with a distance shift exceeding 2 angstrom: {changes}")
            path.write_text(current)
            update_job(connection, job["id"], xyz_file=str(path.relative_to(root)))
            result = calculate("freq", 2)
        update_job(connection, job["id"], electronic_energy_eh=result["electronic_energy_eh"],
                   gibbs_free_energy_eh=result["gibbs_free_energy_eh"], imaginary_count=result["imaginary_count"],
                   frequencies_json=json.dumps(result["frequencies"]))
        expected = 1 if job["endpoint"] == "TS" else 0
        if result["imaginary_count"] != expected:
            raise ValueError(f"Expected {expected} imaginary frequencies, found {result['imaginary_count']}")
        connection.execute("BEGIN IMMEDIATE")
        try:
            update_job(connection, job["id"], status="succeeded", stage="complete", retained=1)
            if job["endpoint"] != "TS":
                select_endpoints(connection, job)
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
    except Exception as exc:
        update_job(connection, job["id"], status="failed", error=str(exc), retained=0)
        print(f"Job {job['id']} ({job['endpoint']}): {exc}", flush=True)


def submit_workflow(config_path, config, *, workers=100, account=None, partition=None, resume=False):
    """Submit the full workflow, or only ORCA workers for an existing database."""
    base = config_path.resolve().parent
    if workers < 1:
        raise ValueError("At least one worker is required")
    root = base / config["paths"]["results"]
    if resume:
        if not (root / "workflow.sqlite").is_file():
            raise ValueError("No prepared workflow database found; use submit for a new run")
    else:
        rows = read_ligands(base / config["paths"]["ligands_csv"])
        if not rows:
            raise ValueError("At least one CSV row is required")
        source = base / config["paths"]["source"]
        if (source.exists() and any(source.iterdir())) or (root / "workflow.sqlite").exists():
            raise ValueError("Choose fresh source/results directories for submit; use resume for prepared jobs")
    command = shlex.join([sys.executable, "-m", "hpc_workflow", "--config", str(config_path.resolve())])
    wrap = "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; exec " + command
    common = ["sbatch", "--parsable", f"--chdir={base}", "--export=ALL",
              f"--mem={config['calculation']['memory_mb']}", "--time=24:00:00"]
    for flag, value in (("account", account), ("partition", partition)):
        if value:
            common.append(f"--{flag}={value}")

    def submit(name, stage, *, cpus=1, dependency=None, dependency_type="afterok", array=None):
        args = [*common, f"--job-name=rxflow-{name}", f"--cpus-per-task={cpus}",
                f"--output={base / (name + '-%A_%a.out')}"]
        if dependency:
            args.append(f"--dependency={dependency_type}:{dependency}")
        if array:
            args.append(f"--array={array}")
        return subprocess.check_output([*args, "--wrap", wrap + " " + stage], text=True).strip().split(";")[0]

    if resume:
        return dict(orca=submit("orca", "worker", cpus=config["calculation"]["cpus"],
                                array=f"0-{workers - 1}%{workers}"))

    sampling = submit("sampling", 'sample --ligand-index "$SLURM_ARRAY_TASK_ID"', array=f"0-{len(rows) - 1}")
    prepared = submit("prepare", "prepare", dependency=sampling, dependency_type="afterany")
    orca = submit("orca", "worker", cpus=config["calculation"]["cpus"],
                  dependency=prepared, array=f"0-{workers - 1}%{workers}")
    return dict(sampling=sampling, prepare=prepared, orca=orca)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("hpc.toml"))
    parser.add_argument("command", choices=("submit", "resume", "sample", "prepare", "worker", "status", "retry"))
    parser.add_argument("--id", type=int, help="Specific stopped/failed row to retry")
    parser.add_argument("--ligand-index", type=int, help="Zero-based CSV row for a sampling array task")
    parser.add_argument("--workers", type=int, default=100, help="ORCA workers for submit/resume (default 100)")
    parser.add_argument("--account", help="Slurm account for submit/resume")
    parser.add_argument("--partition", help="Slurm partition for submit/resume")
    args = parser.parse_args(argv)
    config = tomllib.loads(args.config.read_text())
    base = args.config.resolve().parent
    root = (base / config["paths"]["results"]).resolve()
    root.mkdir(parents=True, exist_ok=True)
    settings = config["calculation"]
    settings.setdefault("opt_attempts", settings.get("ts_attempts", 2))
    for key in ("cpus", "memory_mb", "max_iterations", "opt_attempts", "ts_attempts", "timeout_seconds"):
        if settings[key] <= 0:
            parser.error(f"calculation.{key} must be positive")
    executables = {}
    if args.command in ("submit", "resume", "sample", "worker"):
        keys = ("xtb_exe",) if args.command == "sample" else ("orca_exe", "external_exe", "xtb_exe")
        for key in keys:
            value = os.path.expanduser(os.path.expandvars(config["executables"][key]))
            executable = shutil.which(value if not value.startswith(".") else str(base / value))
            if not executable:
                parser.error(f"Executable {key} not found: {value!r}; set it in {args.config} or its environment variable")
            executables[key] = str(Path(executable).resolve())
    if args.command in ("submit", "resume"):
        if not shutil.which("sbatch"):
            parser.error(f"{args.command} requires Slurm's sbatch command")
        try:
            print(submit_workflow(args.config, config, workers=args.workers,
                                  account=args.account, partition=args.partition, resume=args.command == "resume"))
        except (ValueError, subprocess.CalledProcessError) as exc:
            parser.error(str(exc))
        return
    if args.command == "sample":
        from hpc_sampling import sample_ligand

        rows = read_ligands(base / config["paths"]["ligands_csv"])
        if args.ligand_index is not None:
            if not 0 <= args.ligand_index < len(rows):
                parser.error("--ligand-index is outside the CSV row range")
            rows = [rows[args.ligand_index]]
        sampling = config.get("sampling", {})
        if any(sampling.get(key, default) < 1 for key, default in (
            ("endpoint_conformers", 3), ("coembed_conformers", 10))):
            parser.error("Sampling conformer counts must be positive")
        source = (base / config["paths"]["source"]).resolve()
        for row in rows:
            print(sample_ligand(row, ligand_directory(source, row["ID"]), executables["xtb_exe"], **sampling), flush=True)
        return
    connection = connect_database(root / "workflow.sqlite")
    try:
        if args.command == "prepare":
            print("New jobs, reaction/seating groups:", prepare_jobs(connection,
                (base / config["paths"]["source"]).resolve(), root, settings))
        elif args.command == "worker":
            while (job := claim_job(connection)) is not None:
                print(f"Job {job['id']}: ligand {job['ligand_id']} {job['endpoint']}", flush=True)
                run_job(connection, job, root, executables)
        elif args.command == "retry":
            if args.id is None:
                parser.error("retry requires --id; stop the previous worker before retrying a running row")
            row = connection.execute("SELECT stage FROM jobs WHERE id=?", (args.id,)).fetchone()
            if row and row["stage"] == "sampling":
                parser.error("Sampling failures need corrected seeds followed by prepare; they cannot be retried as ORCA jobs")
            connection.execute("UPDATE jobs SET status='pending',stage='opt',error=NULL,retained=NULL,"
                               "electronic_energy_eh=NULL,gibbs_free_energy_eh=NULL,frequencies_json=NULL,"
                               "imaginary_count=NULL,xyz_file=NULL,out_file=NULL,worker_id=NULL "
                               "WHERE id=? AND status IN ('failed','running')", (args.id,))
        else:
            for row in connection.execute("SELECT ligand_id,substrate,binding_permutation,substrate_site,endpoint,status,count(*) AS jobs FROM jobs "
                                          "GROUP BY ligand_id,substrate,binding_permutation,substrate_site,endpoint,status "
                                          "ORDER BY ligand_id,binding_permutation,substrate_site,endpoint,status"):
                print(dict(row))
    finally:
        connection.close()


if __name__ == "__main__":
    main()
