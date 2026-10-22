"""Portable ORCA/UMA reoptimizations seeded from the plotted state minima.

The source databases are read only. A prepared run contains its own geometries,
queues, recovery inputs and protocol; workers never need the source cluster.
"""

from __future__ import annotations

from collections import Counter
from contextlib import closing, ExitStack
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import sqlite3
import subprocess
import sys


QUEUES = {
    "reaction": ("jobs", "workflow.sqlite"),
    "encounter": ("encounter_jobs", "encounters.sqlite"),
    "precursor": ("precursor_jobs", "precursors.sqlite"),
    "bound": ("bound_reactant_jobs", "acid_base.sqlite"),
    "reference": ("reference_jobs", "references/references.sqlite"),
}
MANIFEST = "uma_manifest.json"


def results_root(config, base, source_root):
    value = config.get("uma", {}).get("results")
    return (base / os.path.expandvars(os.path.expanduser(value))).resolve() if value else source_root / "uma"


def calculation_settings(config):
    options = config.get("uma", {})
    settings = dict(cpus=1, memory_mb=8000, opt_attempts=2, ts_attempts=2,
                    neb_fallback=True, timeout_seconds=86400)
    settings.update(config.get("calculation", {}))
    settings.update(options.get("calculation", {}))
    # The UMA wrapper owns the model/solvation. Never send gxTB flags to it.
    settings.pop("solvent", None)
    settings.pop("max_iterations", None)
    parameters = options.get("external_params", [])
    if isinstance(parameters, str):
        parameters = shlex.split(parameters)
    if not isinstance(parameters, list) or any(not isinstance(p, str) for p in parameters):
        raise ValueError("uma.external_params must be a string or a list of strings")
    for key in ("cpus", "memory_mb", "opt_attempts", "ts_attempts", "timeout_seconds"):
        if type(settings[key]) is not int or settings[key] <= 0:
            raise ValueError(f"uma.calculation.{key} must be a positive integer")
    settings.update(method="uma", external_params=shlex.join(parameters),
                    model_label=options.get("model_label", "UMA"))
    if not isinstance(settings["model_label"], str) or not settings["model_label"].strip():
        raise ValueError("uma.model_label must be a nonempty string identifying the model/protocol")
    return settings


def method_protocol(settings):
    # Executable locations and allocations can change when moving clusters.
    # Chemistry arguments are frozen; device is provided by the wrapper/env.
    if settings.get("method") != "uma" or any(key not in settings for key in ("external_params", "model_label")):
        raise ValueError("Result is not identified as a complete UMA protocol")
    return {key: settings[key] for key in ("method", "external_params", "model_label")}


def open_queue(root, name):
    from hpc_workflow import connect_database
    from hpc_precursors import connect_precursors
    path = root / QUEUES[name][1]
    path.parent.mkdir(parents=True, exist_ok=True)
    if name == "reaction":
        return connect_database(path)
    return connect_precursors(path, bound=name == "bound", encounter=name == "encounter")


def load_manifest(root, settings=None):
    path = root / MANIFEST
    if not path.is_file():
        raise ValueError(f"No prepared UMA run at {root}; run 'rxflow ... uma prepare' first")
    manifest = json.loads(path.read_text())
    if manifest.get("version") != 1 or manifest.get("method") != "uma":
        raise ValueError("Unrecognized UMA run manifest")
    if settings is not None and manifest["protocol"] != method_protocol(settings):
        raise ValueError("UMA model/arguments differ from this prepared run; use its frozen configuration or a new results directory")
    return manifest


def validate_seed_files(root, manifest):
    from hpc_singlepoint import geometry_hash
    for record in manifest["sources"]:
        queue_root = root / "references" if record["queue"] == "reference" else root
        if geometry_hash((queue_root / record["input_xyz"]).read_text()) != record["geometry_sha256"]:
            raise ValueError(f"Prepared UMA seed changed: {record['input_xyz']}")


def selected_minima(rows):
    """One lowest corrected-G conformer per plotted state, without sink filtering."""
    from hpc_precursors import canonical_substrate
    best = {}
    for row in rows:
        energy = row.get("gibbs_free_energy_eh")
        if (row.get("status") != "succeeded" or energy is None or not math.isfinite(energy)
                or row.get("singlepoint_status") != "succeeded" or row.get("singlepoint_job_id") is None):
            continue
        substrate = row.get("substrate")
        comparable_substrate = (canonical_substrate(substrate) if substrate and row["endpoint"] in
                                ("reactant", "TS", "product", "encounter") else substrate)
        key = (row["ligand_id"], comparable_substrate,
               row["charge"], row["multiplicity"], row.get("substrate_site"),
               row["pair_index"], row["seating"], row.get("binding_permutation"),
               row["endpoint"], row.get("protonation_site", 0))
        previous = best.get(key)
        if previous is None or (energy, row["id"]) < (previous["gibbs_free_energy_eh"], previous["id"]):
            best[key] = row
    return sorted(best.values(), key=lambda row: row["id"])


def plan_uma(config, base, source_root):
    from hpc_precursors import load_diagram_inputs, reference_root
    from hpc_singlepoint import source_geometry, geometry_hash
    from rxflow import ligand_directory
    from ts_scan import scan_inputs

    source_reference = reference_root(config, base) / "references.sqlite"
    data = load_diagram_inputs(source_root / "workflow.sqlite", source_root / "precursors.sqlite",
                               source_reference, energy_mode="dft")
    settings = calculation_settings(config)
    root = results_root(config, base, source_root)
    if root == source_root or source_root.is_relative_to(root):
        raise ValueError("UMA results must be separate from the source results directory")
    include_balanced = config.get("uma", {}).get("include_balanced_references", True)
    records = []
    reaction_index = {row["id"]: row for row in data["jobs"]}
    for name, (key, database_name) in QUEUES.items():
        if not include_balanced and name not in ("reaction", "encounter"):
            continue
        source_db = source_reference if name == "reference" else source_root / database_name
        for row in selected_minima(data[key]):
            xyz = source_geometry(row, source_db)
            folder = ligand_directory(Path(), row["ligand_id"]) / ("ts_guess" if row["endpoint"] == "TS" else row["endpoint"])
            target = folder / f"source_{row['id']}_input.xyz"
            recovery = {}
            if row["endpoint"] == "TS":
                candidate = row
                seen = set()
                while candidate.get("ts_parent_id") in reaction_index and candidate["ts_parent_id"] not in seen:
                    seen.add(candidate["ts_parent_id"])
                    candidate = reaction_index[candidate["ts_parent_id"]]
                try:
                    inputs = scan_inputs(candidate, source_root, base / config["paths"]["source"])
                    recovery = {key: path.read_text() for key, path in inputs.items()}
                except (ValueError, OSError, KeyError):
                    # Selected TS minima may lack their archived scan. OptTS and
                    # the normal group NEB fallback still use portable endpoints.
                    pass
            records.append(dict(queue=name, source_database=str(source_db.resolve()),
                                source_job_id=row["id"], source_gibbs_eh=row["gibbs_free_energy_eh"],
                                geometry_sha256=geometry_hash(xyz), input_xyz=str(target),
                                row=row, xyz=xyz, recovery=recovery))
    if not any(record["queue"] == "reaction" for record in records):
        raise ValueError("No plotted states with successful matching PBE0 single points are available to seed UMA")
    return dict(root=str(root), method="uma", protocol=method_protocol(settings),
                source_energy_model=data["energy_model"],
                counts=dict(Counter(record["queue"] for record in records)), records=records)


def portable_config(config, settings):
    """A destination config needs no source CSV, source database or original path."""
    def value(v):
        return str(v).lower() if isinstance(v, bool) else json.dumps(v)
    lines = ['[paths]', 'results = "."', '', '[uma]', 'results = "."',
             'orca_exe = "${ORCA_EXE}"', 'external_exe = "${ORCA_UMA_EXTERNAL_TOOLS}"']
    for key in ("model_label",):
        lines.append(f"{key} = {value(settings[key])}")
    lines.append(f'external_params = {value(shlex.split(settings["external_params"]))}')
    options = config.get("uma", {})
    for key, default in (("gpus_per_task", 1), ("time", "24:00:00")):
        lines.append(f"{key} = {value(options.get(key, default))}")
    lines.extend(['', '[uma.calculation]'])
    for key in ("cpus", "memory_mb", "opt_attempts", "ts_attempts", "neb_fallback", "timeout_seconds"):
        lines.append(f"{key} = {value(settings[key])}")
    lines.extend(['', '[sampling]'])
    for key, default in (("ts_conformers", 20), ("prune_rmsd", 0.5), ("seed", 12345)):
        lines.append(f"{key} = {value(config.get('sampling', {}).get(key, default))}")
    return "\n".join(lines) + "\n"


def prepare_uma(config, base, source_root):
    settings = calculation_settings(config)
    root = results_root(config, base, source_root)
    if (root / MANIFEST).is_file():
        manifest = load_manifest(root, settings)
        validate_seed_files(root, manifest)
        return dict(root=str(root), existing=True, counts=manifest["counts"], method="uma")
    if root.exists() and any(root.iterdir()):
        raise ValueError("UMA results directory is nonempty without a complete manifest; use a fresh directory")
    plan = plan_uma(config, base, source_root)
    root.mkdir(parents=True, exist_ok=True)
    records = plan["records"]
    with ExitStack() as stack:
        queues = {name: stack.enter_context(closing(open_queue(root, name))) for name in QUEUES}
        for record in records:
            db, row = queues[record["queue"]], record["row"]
            queue_root = root / "references" if record["queue"] == "reference" else root
            target = queue_root / record["input_xyz"]
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(record["xyz"])
            columns = {item["name"] for item in db.execute("PRAGMA table_info(jobs)")}
            keep = ("id", "ligand_id", "pair_index", "seating", "endpoint", "conformer", "smiles", "charge",
                    "multiplicity", "donor", "hydrogen", "acceptor", "substrate", "binding_permutation",
                    "binding_mode", "coordination_isomer", "donor_maps", "substrate_site", "resonance", "protonation_site")
            values = {key: row[key] for key in keep if key in row and key in columns}
            job_settings = dict(settings)
            if row["endpoint"] == "TS":
                # Reserve conformer 0 for the existing group NEB fallback.
                values.update(conformer=1, ts_strategy="direct")
                metadata = {}
                for role, name in (("profile_xyz", "scan.xyz"), ("reactant_xyz", "reactant.xyz"), ("product_xyz", "product.xyz")):
                    if role in record["recovery"]:
                        path = target.parent / f"source_{row['id']}_recovery" / name
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(record["recovery"][role])
                        metadata[role] = str(path.relative_to(root))
                job_settings["ts_scan_recovery"] = "profile_xyz" in metadata
                values["ts_guess_json"] = json.dumps(metadata)
            values.update(input_xyz=record["input_xyz"], settings_json=json.dumps(job_settings, sort_keys=True))
            db.execute(f"INSERT INTO jobs ({','.join(values)}) VALUES ({','.join('?' for _ in values)})", tuple(values.values()))

        reaction = queues["reaction"]
        for record in records:
            row = record["row"]
            if record["queue"] != "reaction" or row["endpoint"] != "TS":
                continue
            metadata = json.loads(reaction.execute("SELECT ts_guess_json FROM jobs WHERE id=?", (row["id"],)).fetchone()[0])
            for role in ("reactant", "product"):
                endpoint = reaction.execute("SELECT input_xyz FROM jobs WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint=?",
                                            (row["ligand_id"], row["pair_index"], row["seating"], role)).fetchone()
                if endpoint and f"{role}_xyz" not in metadata:
                    metadata[f"{role}_xyz"] = endpoint[0]
            reaction.execute("UPDATE jobs SET ts_guess_json=? WHERE id=?", (json.dumps(metadata), row["id"]))
            if settings["neb_fallback"] and all(f"{role}_xyz" in metadata for role in ("reactant", "product")):
                seed = dict(reaction.execute("SELECT * FROM jobs WHERE id=?", (row["id"],)).fetchone())
                keep = ("ligand_id", "pair_index", "seating", "endpoint", "input_xyz", "smiles", "charge", "multiplicity",
                        "donor", "hydrogen", "acceptor", "settings_json", "substrate", "binding_permutation", "binding_mode",
                        "coordination_isomer", "donor_maps", "substrate_site", "resonance", "ts_guess_json")
                values = {key: seed[key] for key in keep}
                values.update(conformer=0, stage="neb", ts_strategy="neb")
                reaction.execute(f"INSERT INTO jobs ({','.join(values)}) VALUES ({','.join('?' for _ in values)})", tuple(values.values()))

        # Preserve balanced link identities, never the source energy values.
        for name, table in (("precursor", "precursor_links"), ("bound", "bound_reactant_links")):
            original = source_root / QUEUES[name][1]
            if not original.is_file():
                continue
            with closing(sqlite3.connect(f"file:{original}?mode=ro", uri=True)) as source:
                source.row_factory = sqlite3.Row
                if not source.execute("SELECT 1 FROM sqlite_master WHERE name=?", (table,)).fetchone():
                    continue
                for row in source.execute(f"SELECT * FROM {table}"):
                    values = dict(row)
                    queues[name].execute(f"INSERT OR IGNORE INTO {table} ({','.join(values)}) VALUES ({','.join('?' for _ in values)})", tuple(values.values()))
    manifest = {key: value for key, value in plan.items() if key not in ("records", "root")}
    manifest.update(version=1, sources=[{key: value for key, value in record.items()
                                       if key not in ("row", "xyz", "recovery")} for record in records])
    (root / "run_uma.toml").write_text(portable_config(config, settings))
    (root / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    return dict(root=str(root), existing=False, counts=plan["counts"], method="uma",
                portable_config=str(root / "run_uma.toml"))


def status_uma(root):
    manifest = load_manifest(root)
    counts = []
    for name, (_, filename) in QUEUES.items():
        path = root / filename
        if not path.is_file():
            continue
        with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as db:
            db.row_factory = sqlite3.Row
            counts.extend(dict(queue=name, **dict(row)) for row in db.execute(
                "SELECT endpoint,status,stage,count(*) jobs FROM jobs GROUP BY endpoint,status,stage"))
    return dict(root=str(root), method=manifest["protocol"]["model_label"], counts=counts)


def resolve_executables(config, base):
    options = config.get("uma", {})
    shared = config.get("executables", {})
    result = {}
    for key, default in (("orca_exe", shared.get("orca_exe", "${ORCA_EXE}")),
                         ("external_exe", "${ORCA_UMA_EXTERNAL_TOOLS}")):
        value = os.path.expanduser(os.path.expandvars(options.get(key, default)))
        found = shutil.which(str(base / value) if not Path(value).is_absolute() and "/" in value else value)
        if not found:
            raise ValueError(f"UMA executable {key} not found: {value!r}; set uma.{key} or its environment variable on the compute cluster")
        result[key] = str(Path(found).resolve())
    for key in ("mpirun_exe",):
        value = options.get(key)
        if value:
            found = shutil.which(os.path.expandvars(os.path.expanduser(value)))
            if not found:
                raise ValueError(f"UMA executable {key} not found: {value!r}")
            result[key] = found
    return result


def worker_uma(config, base, root):
    from hpc_workflow import claim_job, run_job
    settings = calculation_settings(config)
    manifest = load_manifest(root, settings)
    validate_seed_files(root, manifest)
    executables = resolve_executables(config, base)
    # Allocation/device arguments can be resolved on the destination cluster.
    executables.update({key: settings[key] for key in ("cpus", "memory_mb", "timeout_seconds")})
    executables["external_params"] = os.path.expandvars(settings["external_params"])
    with ExitStack() as stack:
        queues = {name: stack.enter_context(closing(open_queue(root, name))) for name in QUEUES}
        while True:
            for name, db in queues.items():
                job = claim_job(db, config.get("sampling", {}) if name == "reaction" else None)
                if job is not None:
                    print(f"UMA {name} job {job['id']}: {job['ligand_id']} {job['endpoint']}", flush=True)
                    run_job(db, job, root / "references" if name == "reference" else root, executables)
                    break
            else:
                break


def submit_uma(config_path, config, base, root, *, workers=10, partition=None, account=None, dry_run=False):
    settings = calculation_settings(config)
    load_manifest(root, settings)
    options = config.get("uma", {})
    gpus = options.get("gpus_per_task", 1)
    if type(gpus) is not int or gpus < 0 or type(workers) is not int or workers < 1:
        raise ValueError("UMA workers must be positive and uma.gpus_per_task must be a nonnegative integer (0 for CPU)")
    pending = sum(row["jobs"] for row in status_uma(root)["counts"] if row["status"] == "pending")
    if not pending:
        return dict(root=str(root), pending=0, workers=None)
    workers = min(workers, pending)
    logs = root / "logs"
    command = shlex.join([sys.executable, "-m", "hpc_workflow", "--config", str(config_path.resolve()), "uma", "worker"])
    args = ["sbatch", "--parsable", f"--chdir={base}", "--export=ALL", "--job-name=rxflow-uma",
            "--nodes=1", "--ntasks=1", f"--cpus-per-task={settings['cpus']}",
            f"--array=0-{workers - 1}%{workers}", f"--mem={settings['memory_mb']}",
            f"--time={options.get('time', '24:00:00')}",
            f"--output={logs / 'uma-%A_%a.out'}", f"--error={logs / 'uma-%A_%a.err'}"]
    if gpus:
        args.append(f"--gpus-per-task={gpus}")
    for key, value in (("partition", partition or options.get("partition")), ("account", account or options.get("account"))):
        if value:
            args.append(f"--{key}={value}")
    args += ["--wrap", "exec " + command]
    if dry_run:
        return dict(root=str(root), pending=pending, command=args, dry_run=True)
    logs.mkdir(parents=True, exist_ok=True)
    job_id = subprocess.check_output(args, text=True).strip().split(";")[0]
    return dict(root=str(root), pending=pending, workers=job_id, array_workers=workers)


def load_uma_diagram_inputs(reaction_database, uma_root=None):
    from hpc_precursors import load_diagram_inputs
    root = Path(uma_root).resolve() if uma_root is not None else reaction_database.parent / "uma"
    manifest = load_manifest(root)
    data = load_diagram_inputs(root / "workflow.sqlite", root / "precursors.sqlite",
                               root / "references/references.sqlite", energy_mode="gxtb")
    counts = Counter()
    for key, _ in QUEUES.values():
        for row in data[key]:
            settings = json.loads(row["settings_json"])
            if method_protocol(settings) != manifest["protocol"]:
                raise ValueError(f"UMA database contains a different energy protocol: {key} job {row['id']}")
            row.update(energy_method="uma", method_status=row["status"])
            if row["status"] == "succeeded" and (row.get("gibbs_free_energy_eh") is None
                    or not math.isfinite(row["gibbs_free_energy_eh"])):
                row["method_status"] = "invalid"
            if row.get("ts_strategy") != "sampling" and row["status"] != "skipped":
                counts[row["method_status"]] += 1
    data["singlepoint_summary"] = dict(active_mode="uma", requested_mode="uma", method=manifest["protocol"]["model_label"],
        selected_states=sum(counts.values()), ready=counts["succeeded"], counts=dict(counts))
    data["energy_model"] = manifest["protocol"]["model_label"] + " + UMA NumFreq thermal correction"
    data["results_database"] = str(root / "workflow.sqlite")
    data["results_root"] = str(root)
    return data


def handle_uma(args, config, base, source_root, parser):
    action = args.precursor_command or "plan"
    root = results_root(config, base, source_root)
    try:
        if action == "plan":
            if (root / MANIFEST).is_file():
                result = status_uma(root)
            else:
                result = {key: value for key, value in plan_uma(config, base, source_root).items() if key != "records"}
        elif action == "prepare":
            if args.dry_run:
                result = {key: value for key, value in plan_uma(config, base, source_root).items() if key != "records"}
            else:
                result = prepare_uma(config, base, source_root)
        elif action == "status":
            result = status_uma(root)
        elif action in ("submit", "resume"):
            if not (root / MANIFEST).is_file():
                if action == "resume" or args.dry_run:
                    raise ValueError("Prepare a portable UMA run before resume or submission dry-run")
                prepare_uma(config, base, source_root)
            result = submit_uma(args.config, config, base, root, workers=args.workers,
                                partition=args.partition, account=args.account, dry_run=args.dry_run)
        elif action == "worker":
            worker_uma(config, base, root)
            return
        elif action == "verify-ts":
            from hpc_workflow import audit_transition_states
            load_manifest(root)
            result = audit_transition_states(root, dry_run=args.dry_run, job_id=args.id)
        elif action == "retry":
            from hpc_workflow import update_job
            if args.id is None:
                raise ValueError("uma retry requires --id and optionally --queue")
            load_manifest(root)
            with closing(open_queue(root, args.queue)) as db:
                row = db.execute("SELECT * FROM jobs WHERE id=?", (args.id,)).fetchone()
                if row is None or row["status"] not in ("failed", "running"):
                    raise ValueError("Only failed/interrupted UMA rows can be retried; stop active workers before recovery")
                stage = "ts_sampling" if row["ts_strategy"] == "sampling" else "neb" if row["ts_strategy"] == "neb" else "opt"
                update_job(db, args.id, status="pending", stage=stage, error=None, worker_id=None,
                           retained=None, xyz_file=None, out_file=None, electronic_energy_eh=None,
                           gibbs_free_energy_eh=None, frequencies_json=None, imaginary_count=None,
                           ts_mode_status=None, ts_mode_fingerprint=None, ts_mode_json=None, ts_scan_json=None)
                result = dict(queue=args.queue, id=args.id, status="pending")
        else:
            raise ValueError(f"Unsupported UMA action: {action}")
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, KeyError, sqlite3.Error, subprocess.CalledProcessError) as exc:
        parser.error(str(exc))
