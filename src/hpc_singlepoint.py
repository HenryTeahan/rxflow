"""PBE0/def2-TZVP-D4/SMD(water) energies on selected gxTB geometries."""

from __future__ import annotations

from contextlib import closing
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tempfile

from hpc_precursors import reference_root
from orca_hpc import molecule_with_xyz, orca_input, run_orca
from rxflow import HARTREE_TO_KCAL

METHOD = "PBE0/def2-TZVP-D4/SMD(water)"
TS_ENERGY_WINDOW_KCAL = 3.0
TS_WINDOW_SKIP_REASON = "Outside current 3 kcal/mol TS Gibbs-energy window or source changed"


def read_jobs(path):
    path = Path(path)
    if not path.is_file():
        return []
    with closing(sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        return [dict(row) for row in db.execute("SELECT * FROM jobs")]


def state_key(row):
    return (*tuple(row[name] for name in ("ligand_id", "pair_index", "seating", "endpoint")),
            row.get("protonation_site", 0))


def lowest_jobs(rows):
    """One lowest-G successful NumFreq result per species/state, ignoring retained."""
    best = {}
    for row in rows:
        if (row.get("status") != "succeeded" or row.get("gibbs_free_energy_eh") is None
                or not math.isfinite(row["gibbs_free_energy_eh"])):
            continue
        key = state_key(row)
        if key not in best or (row["gibbs_free_energy_eh"], row["id"]) < (best[key]["gibbs_free_energy_eh"], best[key]["id"]):
            best[key] = row
    return list(best.values())


def selected_jobs(rows):
    """Lowest-G minima and every TS within 3 kcal/mol of its group's lowest gxTB G."""
    rows = list(rows)
    groups = {}
    for row in rows:
        if (row.get("status") == "succeeded" and row.get("gibbs_free_energy_eh") is not None
                and math.isfinite(row["gibbs_free_energy_eh"])):
            groups.setdefault(state_key(row), []).append(row)
    selected = []
    for group in groups.values():
        ranked = sorted(group, key=lambda row: (row["gibbs_free_energy_eh"], row["id"]))
        minimum = ranked[0]
        if minimum["endpoint"] != "TS":
            selected.append(minimum)
            continue
        cutoff = minimum["gibbs_free_energy_eh"] + TS_ENERGY_WINDOW_KCAL / HARTREE_TO_KCAL
        selected.extend(row for row in ranked if row["gibbs_free_energy_eh"] <= cutoff + 1e-12)
    return selected


def source_fingerprint(row, xyz):
    """Invalidate an SP if its source geometry, run, graph or thermal correction changes."""
    for name in ("electronic_energy_eh", "gibbs_free_energy_eh"):
        if row.get(name) is None or not math.isfinite(row[name]):
            raise ValueError(f"Source job {row['id']} lacks a finite NumFreq {name}")
    payload = {name: row[name] for name in ("id", "attempt", "smiles", "charge", "multiplicity",
                                          "electronic_energy_eh", "gibbs_free_energy_eh")}
    payload["geometry_sha256"] = hashlib.sha256(xyz.encode()).hexdigest()
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def source_geometry(row, database):
    if not row.get("xyz_file"):
        raise ValueError(f"Source job {row['id']} has no optimized geometry")
    xyz = (database.parent / row["xyz_file"]).read_text()
    molecule_with_xyz(row["smiles"], xyz)
    return xyz


def calculation_settings(config):
    settings = dict(cpus=4, memory_mb=8000, timeout_seconds=86400)
    overrides = config.get("singlepoint", {})
    unknown = set(overrides) - set(settings)
    if unknown:
        raise ValueError(f"Unknown singlepoint settings: {', '.join(sorted(unknown))}")
    settings.update(overrides)
    if any(not isinstance(value, int) or value <= 0 for value in settings.values()):
        raise ValueError("singlepoint cpus, memory_mb and timeout_seconds must be positive integers")
    if int(settings["memory_mb"] * 0.75 / settings["cpus"]) < 1:
        raise ValueError("singlepoint memory_mb is too small for the requested CPUs")
    return settings


def mpi_executable(config, base):
    """Find the node's MPI installation without loading environment modules."""
    for value in (config["executables"].get("mpirun_exe"),
                  config["executables"].get("mpirun_fallback_exe"), "mpirun"):
        if not value:
            continue
        value = os.path.expanduser(os.path.expandvars(value))
        if os.path.dirname(value):
            value = str((base / value).absolute())
        executable = shutil.which(value)
        if executable:
            # Preserve bin/mpirun rather than resolving its orterun symlink:
            # the sibling lib directory belongs to this MPI installation.
            return str(Path(executable).absolute())
    return None


def sources(config, base, root):
    references = reference_root(config, base)
    return [("reaction", root / "workflow.sqlite", root / "singlepoints.sqlite"),
            ("precursor", root / "precursors.sqlite", root / "singlepoints.sqlite"),
            ("bound_reactant", root / "acid_base.sqlite", root / "singlepoints.sqlite"),
            ("encounter", root / "encounters.sqlite", root / "singlepoints.sqlite"),
            ("reference", references / "references.sqlite", references / "singlepoints.sqlite")]


def connect_singlepoints(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=60, isolation_level=None)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA journal_mode=DELETE")
    db.execute("""CREATE TABLE IF NOT EXISTS jobs (
        id INTEGER PRIMARY KEY, source_database TEXT NOT NULL, source_job_id INTEGER NOT NULL,
        source_fingerprint TEXT NOT NULL, source_kind TEXT NOT NULL, source_json TEXT NOT NULL,
        method TEXT NOT NULL, settings_json TEXT NOT NULL, ligand_id TEXT NOT NULL,
        pair_index INTEGER NOT NULL, seating INTEGER NOT NULL, endpoint TEXT NOT NULL,
        conformer INTEGER NOT NULL, charge INTEGER NOT NULL, multiplicity INTEGER NOT NULL,
        smiles TEXT NOT NULL, input_xyz TEXT NOT NULL, geometry_sha256 TEXT NOT NULL,
        gxtb_electronic_energy_eh REAL NOT NULL, gxtb_gibbs_free_energy_eh REAL NOT NULL,
        status TEXT NOT NULL DEFAULT 'pending', attempt INTEGER NOT NULL DEFAULT 0,
        electronic_energy_eh REAL, corrected_gibbs_free_energy_eh REAL,
        out_file TEXT, error TEXT, worker_id TEXT, updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
        UNIQUE(source_database,source_job_id,source_fingerprint,method))""")
    if 'reused_from_id' not in {row[1] for row in db.execute('PRAGMA table_info(jobs)')}:
        try:
            db.execute('ALTER TABLE jobs ADD COLUMN reused_from_id INTEGER')
        except sqlite3.OperationalError:
            if 'reused_from_id' not in {row[1] for row in db.execute('PRAGMA table_info(jobs)')}:
                raise
    return db


def plan_singlepoints(config, base, root):
    if not (root / "workflow.sqlite").is_file():
        raise ValueError("No reaction database exists; prepare/run the gxTB workflow first")
    calculation_settings(config)
    summary = []
    for kind, source, destination in sources(config, base, root):
        rows = selected_jobs(read_jobs(source))
        summary.append(dict(source_kind=kind, source_database=str(source),
                            singlepoint_database=str(destination), selected_states=len(rows),
                            source_available=source.is_file()))
    return dict(method=METHOD, ts_energy_window_kcal=TS_ENERGY_WINDOW_KCAL, sources=summary,
                selected_states=sum(item["selected_states"] for item in summary))


def freeze_geometry(target, xyz):
    """Publish a complete immutable snapshot, including during concurrent run submissions."""
    if target.exists():
        if target.read_text() != xyz:
            raise ValueError(f"Existing single-point geometry differs: {target}")
        return
    with tempfile.NamedTemporaryFile(mode="w", dir=target.parent, prefix=".rxflow-sp-", delete=False) as handle:
        handle.write(xyz)
        temporary = Path(handle.name)
    try:
        try:
            os.link(temporary, target)
        except FileExistsError:
            if target.read_text() != xyz:
                raise ValueError(f"Existing single-point geometry differs: {target}")
    finally:
        temporary.unlink()


def prepare_singlepoints(config, base, root):
    plan = plan_singlepoints(config, base, root)
    settings = calculation_settings(config)
    added = existing = reused = skipped = reactivated = 0
    errors = []
    for kind, source, destination in sources(config, base, root):
        rows = selected_jobs(read_jobs(source))
        if not rows and not destination.is_file():
            continue
        with closing(connect_singlepoints(destination)) as db:
            selected_ts = set()
            for row in rows:
                try:
                    xyz = source_geometry(row, source)
                    fingerprint = source_fingerprint(row, xyz)
                except (ValueError, OSError) as exc:
                    errors.append(dict(source_kind=kind, source_job_id=row["id"], error=str(exc)))
                    continue
                identity = (str(source.resolve()), row["id"], fingerprint, METHOD)
                if row["endpoint"] == "TS":
                    selected_ts.add((row["id"], fingerprint))
                if db.execute("SELECT id FROM jobs WHERE source_database=? AND source_job_id=? AND source_fingerprint=? AND method=?", identity).fetchone():
                    existing += 1
                    continue
                folder = Path(row["xyz_file"]).parent / "singlepoint"
                method_tag = hashlib.sha256(METHOD.encode()).hexdigest()[:8]
                stem = f"job_{row['id']}_{fingerprint[:12]}_{method_tag}"
                relative = folder / f"{stem}_input.xyz"
                target = destination.parent / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                freeze_geometry(target, xyz)
                target.with_suffix(".inp").write_text(orca_input(xyz, "singlepoint", row["charge"], row["multiplicity"], settings))
                inserted = db.execute("""INSERT OR IGNORE INTO jobs (source_database,source_job_id,source_fingerprint,method,
                    source_kind,source_json,settings_json,ligand_id,pair_index,seating,endpoint,conformer,
                    charge,multiplicity,smiles,input_xyz,geometry_sha256,gxtb_electronic_energy_eh,gxtb_gibbs_free_energy_eh)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (*identity, kind, json.dumps(row, sort_keys=True), json.dumps(settings, sort_keys=True),
                     *(row[name] for name in ("ligand_id", "pair_index", "seating", "endpoint", "conformer", "charge", "multiplicity", "smiles")),
                     str(relative), hashlib.sha256(xyz.encode()).hexdigest(), row["electronic_energy_eh"], row["gibbs_free_energy_eh"]))
                added += inserted.rowcount
                existing += int(inserted.rowcount == 0)
                if inserted.rowcount:
                    # A graph/state label does not affect an electronic single point.
                    # Reuse only the identical frozen XYZ, charge, spin and method.
                    for previous in db.execute('''SELECT * FROM jobs WHERE status='succeeded'
                            AND method=? AND geometry_sha256=? AND charge=? AND multiplicity=?
                            AND electronic_energy_eh IS NOT NULL ORDER BY id''',
                            (METHOD, hashlib.sha256(xyz.encode()).hexdigest(), row['charge'], row['multiplicity'])):
                        frozen = destination.parent / previous['input_xyz']
                        if (not math.isfinite(previous['electronic_energy_eh']) or not frozen.is_file()
                                or frozen.read_text() != xyz):
                            continue
                        energy = previous['electronic_energy_eh']
                        corrected = row['gibbs_free_energy_eh'] - row['electronic_energy_eh'] + energy
                        db.execute('''UPDATE jobs SET status='succeeded',electronic_energy_eh=?,
                            corrected_gibbs_free_energy_eh=?,out_file=?,reused_from_id=? WHERE id=?''',
                            (energy, corrected, previous['out_file'], previous['id'], inserted.lastrowid))
                        reused += 1
                        break
            # A newly found lower TS can move an old pending correction outside
            # the window. Preserve completed work, but do not calculate stale seeds.
            previous = db.execute("""SELECT id,source_job_id,source_fingerprint,status,error FROM jobs
                WHERE source_database=? AND endpoint='TS' AND method=?
                AND status IN ('pending','skipped')""", (str(source.resolve()), METHOD)).fetchall()
            for queued in previous:
                selected = (queued["source_job_id"], queued["source_fingerprint"]) in selected_ts
                if queued["status"] == "pending" and not selected:
                    db.execute("UPDATE jobs SET status='skipped',error=?,updated_at=CURRENT_TIMESTAMP WHERE id=? AND status='pending'",
                               (TS_WINDOW_SKIP_REASON, queued["id"]))
                    skipped += 1
                elif queued["status"] == "skipped" and selected and queued["error"] == TS_WINDOW_SKIP_REASON:
                    db.execute("UPDATE jobs SET status='pending',error=NULL,updated_at=CURRENT_TIMESTAMP WHERE id=? AND status='skipped'",
                               (queued["id"],))
                    reactivated += 1
    result = dict(**plan, added=added, existing=existing, reused=reused, skipped=skipped,
                  reactivated=reactivated, errors=errors, settings=settings)
    root.mkdir(parents=True, exist_ok=True)
    (root / "singlepoints.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def claim_singlepoint(db):
    db.execute("BEGIN IMMEDIATE")
    try:
        row = db.execute("""UPDATE jobs SET status='running',attempt=attempt+1,error=NULL,
            worker_id=?,updated_at=CURRENT_TIMESTAMP WHERE id=(SELECT id FROM jobs
            WHERE status='pending' AND method=? ORDER BY id LIMIT 1) RETURNING *""",
            (os.environ.get("SLURM_JOB_ID", "local") + f":{os.getpid()}", METHOD)).fetchone()
        db.commit()
        return dict(row) if row else None
    except BaseException:
        db.rollback()
        raise


def run_singlepoint(db, job, root, orca_exe, mpirun_exe=None):
    stem = f"{Path(job['input_xyz']).stem.removesuffix('_input')}_run_{job['attempt']:02d}_singlepoint_1"
    output = Path(job["input_xyz"]).parent / f"{stem}.out"
    try:
        xyz = (root / job["input_xyz"]).read_text()
        if hashlib.sha256(xyz.encode()).hexdigest() != job["geometry_sha256"]:
            raise ValueError("Frozen single-point geometry has changed")
        molecule_with_xyz(job["smiles"], xyz)
        settings = {**json.loads(job["settings_json"]), "orca_exe": orca_exe}
        if mpirun_exe:
            settings["mpirun_exe"] = mpirun_exe
        result = run_orca(job, root, settings, "singlepoint", xyz)
        energy = result["electronic_energy_eh"]
        corrected = job["gxtb_gibbs_free_energy_eh"] - job["gxtb_electronic_energy_eh"] + energy
        if not math.isfinite(corrected):
            raise ValueError("Nonfinite corrected Gibbs energy")
        db.execute("""UPDATE jobs SET status='succeeded',electronic_energy_eh=?,
            corrected_gibbs_free_energy_eh=?,out_file=?,error=NULL,updated_at=CURRENT_TIMESTAMP WHERE id=?""",
            (energy, corrected, result["out_file"], job["id"]))
    except (ValueError, OSError, RuntimeError) as exc:
        db.execute("""UPDATE jobs SET status='failed',electronic_energy_eh=NULL,
            corrected_gibbs_free_energy_eh=NULL,out_file=?,error=?,updated_at=CURRENT_TIMESTAMP WHERE id=?""",
            (str(output), str(exc), job["id"]))
        print(f"Single point {job['id']} ({job['source_kind']}): {exc}", flush=True)


def database_paths(config, base, root):
    return [reference_root(config, base) / "singlepoints.sqlite", root / "singlepoints.sqlite"]


def apply_singlepoint_energies(data, reaction_database, precursor_database, reference_database, mode="auto"):
    """Use a complete corrected energy model, or explicitly requested partial DFT results."""
    if mode not in ("auto", "gxtb", "dft"):
        raise ValueError("Energy mode must be auto, gxtb or dft")
    data.setdefault("bound_reactant_jobs", [])
    data.setdefault("encounter_jobs", [])
    selections = {}
    counts = {}
    for key, source in (("jobs", reaction_database), ("precursor_jobs", precursor_database),
                        ("bound_reactant_jobs", reaction_database.with_name("acid_base.sqlite")),
                        ("encounter_jobs", reaction_database.with_name("encounters.sqlite")),
                        ("reference_jobs", reference_database)):
        records = read_jobs(source.parent / "singlepoints.sqlite")
        index = {(row["source_database"], row["source_job_id"], row["source_fingerprint"]): row
                 for row in records if row["method"] == METHOD}
        matches = {}
        for row in selected_jobs(data[key]):
            state, sp, detail = "not_queued", None, None
            try:
                xyz = source_geometry(row, source)
                fingerprint = source_fingerprint(row, xyz)
                sp = index.get((str(source.resolve()), row["id"], fingerprint))
                if sp:
                    state = sp["status"]
                    if state == "succeeded":
                        energy = sp["electronic_energy_eh"]
                        corrected = sp["corrected_gibbs_free_energy_eh"]
                        if (energy is None or corrected is None or not math.isfinite(energy)
                                or not math.isfinite(corrected)
                                or abs(corrected - (row["gibbs_free_energy_eh"] - row["electronic_energy_eh"] + energy)) > 1e-9):
                            state, detail = "invalid", "Stored single-point correction is invalid"
            except (ValueError, OSError) as exc:
                state, detail = "invalid_source", str(exc)
            matches[row["id"]] = dict(status=state, job=sp, error=detail)
            counts[state] = counts.get(state, 0) + 1
        selections[key] = matches
    total = sum(counts.values())
    ready = counts.get("succeeded", 0)
    use_dft = mode == "dft" or (mode == "auto" and total > 0 and ready == total)
    data["singlepoint_summary"] = dict(method=METHOD, requested_mode=mode,
        active_mode="dft" if use_dft else "gxtb", selected_states=total, ready=ready, counts=counts)
    data["energy_model"] = METHOD + " + gxTB thermal correction" if use_dft else "gxTB NumFreq"
    if not use_dft:
        return data
    for key in selections:
        updated = []
        for row in data[key]:
            row = dict(row)
            row["gxtb_gibbs_free_energy_eh"] = row.get("gibbs_free_energy_eh")
            row["gxtb_electronic_energy_eh"] = row.get("electronic_energy_eh")
            match = selections[key].get(row["id"])
            row["singlepoint_selected"] = match is not None
            row["singlepoint_status"] = match["status"] if match else "not_selected"
            row["gibbs_free_energy_eh"] = None
            if match and match["status"] == "succeeded":
                sp = match["job"]
                row.update(gibbs_free_energy_eh=sp["corrected_gibbs_free_energy_eh"],
                           dft_electronic_energy_eh=sp["electronic_energy_eh"],
                           singlepoint_job_id=sp["id"], singlepoint_out_file=sp["out_file"])
            updated.append(row)
        data[key] = updated
    precursors = {tuple(row[name] for name in ("ligand_id", "pair_index", "seating", "endpoint")): row
                  for row in lowest_jobs(data["precursor_jobs"])}
    references = {row["ligand_id"]: row for row in lowest_jobs(data["reference_jobs"])}
    bound = {tuple(row[name] for name in ("ligand_id", "pair_index", "seating", "endpoint")): row
             for row in lowest_jobs(data["jobs"])}

    def energy(row):
        return row["gibbs_free_energy_eh"] if row else None

    def total_energy(*values):
        return None if any(value is None for value in values) else sum(values)

    acid, acetate = energy(references.get("acetic_acid")), energy(references.get("acetate"))
    links = []
    for original in data["diagram_references"]:
        link = dict(original)
        d = precursors.get((link["ligand_id"], link["binding_permutation"], 0, "deprotonated"))
        p = precursors.get((link["ligand_id"], link["binding_permutation"], link["protonation_site"], "protonated"))
        s = precursors.get((link["substrate_id"], 0, 0, "substrate"))
        link.update(deprotonated_job_id=d["id"] if d else None, protonated_job_id=p["id"] if p else None,
                    substrate_job_id=s["id"] if s else None, deprotonated_complex_gibbs_eh=energy(d),
                    protonated_complex_gibbs_eh=energy(p), substrate_gibbs_eh=energy(s),
                    acetic_acid_gibbs_eh=acid, acetate_gibbs_eh=acetate,
                    deprotonated_complex_plus_substrate_eh=total_energy(energy(d), energy(s)),
                    protonated_complex_plus_substrate_eh=total_energy(energy(p), energy(s)),
                    deprotonated_start_gibbs_eh=total_energy(energy(d), energy(s), acid),
                    protonated_start_gibbs_eh=total_energy(energy(p), energy(s), acetate))
        for role in ("reactant", "TS", "product"):
            row = bound.get((link["ligand_id"], link["pair_index"], link["seating"], role))
            link[f"{role}_with_acid_gibbs_eh"] = total_energy(energy(row), acid)
        links.append(link)
    from hpc_precursors import add_bound_reactant_references
    from hpc_encounters import add_encounter_references
    data["diagram_references"] = add_encounter_references(
        add_bound_reactant_references(links, data["bound_reactant_jobs"]), data["encounter_jobs"])
    data["reference_minima"] = references
    return data


def singlepoint_status(config, base, root):
    result = plan_singlepoints(config, base, root)
    result["queues"] = []
    for path in database_paths(config, base, root):
        rows = read_jobs(path)
        counts = {}
        for row in rows:
            counts[row["status"]] = counts.get(row["status"], 0) + 1
        failures = [{name: row[name] for name in ("id", "source_kind", "source_job_id", "ligand_id",
                                                  "pair_index", "seating", "endpoint", "error", "out_file")}
                    for row in rows if row["status"] == "failed"]
        result["queues"].append(dict(database=str(path), counts=counts, failed=failures))
    from hpc_precursors import load_diagram_inputs
    result["selection"] = load_diagram_inputs(root / "workflow.sqlite", root / "precursors.sqlite",
                                             reference_root(config, base) / "references.sqlite")["singlepoint_summary"]
    return result


def submit_singlepoints(config_path, config, *, workers=10, account=None, partition=None, resume=False):
    if workers < 1:
        raise ValueError("Single-point worker count must be positive")
    base = config_path.resolve().parent
    root = (base / config["paths"]["results"]).resolve()
    if not resume:
        prepared = prepare_singlepoints(config, base, root)
        if prepared["errors"]:
            raise ValueError(f"Cannot prepare all selected states: {prepared['errors']}")
    paths = database_paths(config, base, root)
    if not any(path.is_file() for path in paths):
        raise ValueError("No single-point database exists; use singlepoint prepare or submit first")
    settings = calculation_settings(config)
    pending = [row for path in paths for row in read_jobs(path) if row["status"] == "pending" and row["method"] == METHOD]
    if not pending:
        return dict(workers=None, pending=0, method=METHOD)
    for row in pending:
        frozen = json.loads(row["settings_json"])
        if any(frozen[name] > settings[name] for name in ("cpus", "memory_mb", "timeout_seconds")):
            raise ValueError("Pending single points require larger resources than the configured singlepoint allocation")
    workers = min(workers, len(pending))
    command = shlex.join([sys.executable, "-m", "hpc_workflow", "--config", str(config_path.resolve()), "singlepoint", "worker"])
    logs = base / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    args = ["sbatch", "--parsable", f"--chdir={base}", "--export=ALL", "--job-name=rxflow-singlepoint",
            f"--array=0-{workers - 1}%{workers}", "--nodes=1",
            f"--ntasks={settings['cpus']}", "--cpus-per-task=1", "--no-requeue",
            f"--mem={settings['memory_mb']}", f"--time={math.ceil(settings['timeout_seconds'] / 60)}",
            f"--output={logs / 'singlepoint-%A_%a.out'}",
            f"--error={logs / 'singlepoint-%A_%a.err'}"]
    for flag, value in (("account", account), ("partition", partition)):
        if value:
            args.append(f"--{flag}={value}")
    wrap = "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; exec " + command
    job_id = subprocess.check_output([*args, "--wrap", wrap], text=True).strip().split(";")[0]
    return dict(workers=job_id, pending=len(pending), array_workers=workers, method=METHOD, settings=settings)


def handle_singlepoints(args, config, base, root, parser):
    action = args.precursor_command
    if action is None:
        parser.error("singlepoint requires plan, prepare, submit, worker, status, retry, or resume")
    try:
        if action in ("plan", "prepare", "status"):
            function = {"plan": plan_singlepoints, "prepare": prepare_singlepoints, "status": singlepoint_status}[action]
            print(json.dumps(function(config, base, root), indent=2))
            return
        if action == "retry":
            if args.id is None:
                raise ValueError("singlepoint retry requires --id; stop the old worker first")
            path = (reference_root(config, base) if args.references else root) / "singlepoints.sqlite"
            if not path.is_file():
                raise ValueError(f"No single-point database: {path}")
            with closing(connect_singlepoints(path)) as db:
                db.execute("""UPDATE jobs SET status='pending',error=NULL,worker_id=NULL,
                    electronic_energy_eh=NULL,corrected_gibbs_free_energy_eh=NULL,out_file=NULL,
                    updated_at=CURRENT_TIMESTAMP WHERE id=? AND status IN ('failed','running')""", (args.id,))
            return
        value = config["executables"].get("singlepoint_orca_exe", config["executables"]["orca_exe"])
        value = os.path.expanduser(os.path.expandvars(value))
        executable = shutil.which(str(base / value) if value.startswith(".") else value)
        if not executable:
            raise ValueError(f"ORCA executable not found: {value!r}")
        executable = str(Path(executable).resolve())
        if action in ("submit", "resume"):
            if not shutil.which("sbatch"):
                raise ValueError("singlepoint submit/resume requires Slurm's sbatch command")
            print(json.dumps(submit_singlepoints(args.config, config, workers=args.workers,
                  account=args.account, partition=args.partition, resume=action == "resume"), indent=2))
            return
        mpi = mpi_executable(config, base) if calculation_settings(config)["cpus"] > 1 else None
        if calculation_settings(config)["cpus"] > 1 and mpi is None:
            raise ValueError("No usable mpirun on this node: checked configured MPI paths and PATH; jobs were left pending")
        if mpi:
            print(f"MPI launcher: {mpi}", flush=True)
        paths = [path for path in database_paths(config, base, root) if path.is_file()]
        if not paths:
            raise ValueError("No single-point database exists; prepare or submit first")
        allocation = calculation_settings(config)
        connections = []
        try:
            for path in paths:
                db = connect_singlepoints(path)
                connections.append((path, db))
                for row in db.execute("SELECT settings_json FROM jobs WHERE status='pending' AND method=?", (METHOD,)):
                    frozen = json.loads(row[0])
                    if any(frozen[name] > allocation[name] for name in ("cpus", "memory_mb", "timeout_seconds")):
                        raise ValueError("Pending single points exceed this worker's configured allocation")
            while True:
                for path, db in connections:
                    job = claim_singlepoint(db)
                    if job is not None:
                        print(f"Single point {job['id']}: {job['source_kind']} source {job['source_job_id']}", flush=True)
                        run_singlepoint(db, job, path.parent, executable, mpi)
                        break
                else:
                    break
        finally:
            for _, db in connections:
                db.close()
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.error(str(exc))
