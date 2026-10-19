"""Prepare a flat ORCA job database and drain it with independent HPC workers."""

from __future__ import annotations

import argparse
from contextlib import closing
import csv
from datetime import datetime, timezone
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

from orca_hpc import check_ts_mode, molecule_with_xyz, run_orca, ts_mode_fingerprint
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
                       "donor_maps": "TEXT", "substrate_site": "INTEGER", "resonance": "INTEGER",
                       "ts_rmsd": "REAL", "ts_strategy": "TEXT", "ts_mode_status": "TEXT",
                       "ts_mode_fingerprint": "TEXT", "ts_mode_json": "TEXT",
                       "ts_mode_recovery": "TEXT"}.items():
        if name not in columns:
            connection.execute(f"ALTER TABLE jobs ADD COLUMN {name} {kind}")
    return connection


def update_job(connection, job_id, **values):
    assignments = ",".join(f"{key}=?" for key in values)
    connection.execute(f"UPDATE jobs SET {assignments},updated_at=CURRENT_TIMESTAMP WHERE id=?",
                       [*values.values(), job_id])


def prepare_jobs(connection, source, root, settings):
    """Import each usable seed; retain missing/invalid endpoint and TS seeds as failed rows."""
    settings = {key: value for key, value in settings.items() if key != "max_iterations"}
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
        candidates = selection.get("candidates", [])
        for candidate in candidates:
            rank = candidate["rank"]
            inputs.append(("TS", rank, descriptor.parent / "ts_guess" / f"{prefix}_conformer_{rank:03d}.xyz"))
        if not candidates or settings["neb_fallback"]:
            inputs.append(("TS", 0, descriptor.parent / "ts_guess" / f"{prefix}.xyz"))
        for role, conformer, source_xyz in inputs:
            ts_selection = next((c for c in candidates if c["rank"] == conformer), selection)
            ts_strategy = ("direct" if conformer else "neb") if role == "TS" and candidates else None
            ts_rmsd = ts_selection.get("selected_ts_rmsd") if role == "TS" else None
            smiles = graphs["reactant" if role == "TS" else role]
            mol = mapped_molecule(smiles)
            error, xyz = None, None
            try:
                if role == "TS" and ("ts_input" not in ts_selection or ts_selection.get("error")):
                    raise ValueError(ts_selection.get("error") or "No completed TS guess")
                xyz = source_xyz.read_text()
                molecule_with_xyz(smiles, xyz)
            except (ValueError, OSError) as exc:
                detail = next((e["error"] for e in errors if e.get("endpoint", role) == role), str(exc))
                error = f"{role} sampling failed: {detail}"
                xyz = None
            folder = root / descriptor.parent.name / ("ts_guess" if role == "TS" else role)
            folder.mkdir(parents=True, exist_ok=True)
            target = folder / (f"{source_xyz.stem}_input.xyz" if role == "TS" else source_xyz.name)
            key = (ligand_id, pair, seating, role, conformer)
            previous = connection.execute("SELECT * FROM jobs WHERE ligand_id=? AND pair_index=? AND seating=? "
                                          "AND endpoint=? AND conformer=?", key).fetchone()
            protocol = json.dumps(settings, sort_keys=True)
            if previous and ({key: value for key, value in json.loads(previous["settings_json"]).items()
                              if key != "max_iterations"} != settings or previous["smiles"] != smiles):
                raise ValueError("Existing jobs have different chemistry/settings; use a new HPC results directory")
            if previous and any(previous[name] is not None and previous[name] != value
                                for name, value in metadata.items()):
                raise ValueError("Existing jobs have different substrate/binding metadata; use a new results directory")
            if xyz is not None and target.exists():
                if target.read_text() != xyz:
                    raise ValueError(f"Existing input differs: {target}")
            elif xyz is not None:
                target.write_text(xyz)
            atoms = ts_selection.get("ts_mode", {})
            cursor = connection.execute("""INSERT OR IGNORE INTO jobs
                (ligand_id,pair_index,seating,endpoint,conformer,input_xyz,smiles,charge,multiplicity,
                 donor,hydrogen,acceptor,settings_json,substrate,binding_mode,coordination_isomer,
                 donor_maps,substrate_site,resonance,binding_permutation,status,stage,error,retained,ts_rmsd,ts_strategy)
                 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (*key, str(target.relative_to(root)), smiles, Chem.GetFormalCharge(mol), 1,
                 atoms.get("donor"), atoms.get("hydrogen"), atoms.get("acceptor"), protocol, *metadata.values(),
                 "failed" if error else "pending", "sampling" if error else ("neb" if ts_strategy == "neb" else "opt"),
                 error, 0 if error else None, ts_rmsd, ts_strategy))
            added += cursor.rowcount
            if previous and any(previous[name] != value for name, value in metadata.items()):
                update_job(connection, previous["id"], **metadata)
            if previous and previous["status"] == "failed" and previous["stage"] == "sampling":
                update_job(connection, previous["id"], status="failed" if error else "pending",
                           stage="sampling" if error else "opt", error=error, retained=0 if error else None,
                           donor=atoms.get("donor"), hydrogen=atoms.get("hydrogen"), acceptor=atoms.get("acceptor"),
                           ts_rmsd=ts_rmsd, ts_strategy=ts_strategy)
        groups += 1
    return added, groups


def claim_job(connection):
    """Claim ready jobs; NEB waits for all TS candidates and runs only if none succeeded."""
    connection.execute("BEGIN IMMEDIATE")
    try:
        connection.execute("""UPDATE jobs AS j SET status='skipped',stage='complete',retained=0,
            updated_at=CURRENT_TIMESTAMP WHERE j.status='pending' AND j.ts_strategy='neb'
            AND EXISTS (SELECT 1 FROM jobs t WHERE t.ligand_id=j.ligand_id AND t.pair_index=j.pair_index
              AND t.seating=j.seating AND t.endpoint='TS' AND t.status='succeeded')""")
        job = connection.execute("""UPDATE jobs SET status='running',attempt=attempt+1,error=NULL,
            worker_id=?,updated_at=CURRENT_TIMESTAMP WHERE id=(
              SELECT j.id FROM jobs j WHERE j.status='pending' AND (j.endpoint!='TS' OR NOT EXISTS (
                SELECT 1 FROM jobs e WHERE e.ligand_id=j.ligand_id AND e.pair_index=j.pair_index
                  AND e.seating=j.seating AND e.endpoint!='TS' AND e.status IN ('pending','running')))
              AND (coalesce(j.ts_strategy,'')!='neb' OR NOT EXISTS (
                SELECT 1 FROM jobs t WHERE t.ligand_id=j.ligand_id AND t.pair_index=j.pair_index
                  AND t.seating=j.seating AND t.endpoint='TS' AND t.ts_strategy='direct'
                  AND t.status IN ('pending','running')))
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


def select_transition_state(connection, job):
    group = (job["ligand_id"], job["pair_index"], job["seating"])
    connection.execute("""UPDATE jobs SET retained=(id=(SELECT id FROM jobs
        WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint='TS' AND status='succeeded'
        ORDER BY gibbs_free_energy_eh,id LIMIT 1)),updated_at=CURRENT_TIMESTAMP
        WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint='TS' AND status='succeeded'""",
        (*group, *group))


def run_job(connection, job, root, executables):
    """Preoptimize TS guesses, try OptTS, and prepare coembedded endpoints for NEB."""
    settings = {**json.loads(job["settings_json"]), **executables}
    settings.setdefault("opt_attempts", settings.get("ts_attempts", 2))

    def calculate(stage, attempt=1, endpoints=None, calculation_job=None):
        update_job(connection, job["id"], stage=stage)
        result = run_orca(calculation_job or job, root, settings, stage, current,
                          attempt=attempt, neb_endpoints=endpoints)
        update_job(connection, job["id"], out_file=result["out_file"])
        if result.get("electronic_energy_eh") is not None:
            update_job(connection, job["id"], electronic_energy_eh=result["electronic_energy_eh"])
        return result

    def calculate_frequency(attempt=1):
        path.write_text(current)
        update_job(connection, job["id"], xyz_file=str(path.relative_to(root)))
        result = calculate("freq", attempt)
        if job["endpoint"] == "TS":
            update_job(connection, job["id"], electronic_energy_eh=result["electronic_energy_eh"],
                       gibbs_free_energy_eh=result["gibbs_free_energy_eh"], imaginary_count=result["imaginary_count"],
                       frequencies_json=json.dumps(result["frequencies"]),
                       ts_mode_status=None, ts_mode_fingerprint=None, ts_mode_json=None)
            if result["imaginary_count"] != 1:
                raise ValueError(f"Expected 1 imaginary frequencies, found {result['imaginary_count']}")
            try:
                report = check_ts_mode(root / result["out_file"], job)
            except ImportError:
                raise
            except Exception as exc:
                # An unreadable mode is never accepted; use the same recovery
                # as a successfully analyzed but incorrect imaginary mode.
                reason = f"TS mode analysis failed: {exc}"
                update_job(connection, job["id"], ts_mode_status="error",
                           ts_mode_json=json.dumps(dict(valid=False, reason=reason)))
                raise ValueError(reason) from exc
            update_job(connection, job["id"], ts_mode_status="passed" if report["valid"] else "failed",
                       ts_mode_fingerprint=report["fingerprint"], ts_mode_json=json.dumps(report, sort_keys=True))
            if not report["valid"]:
                raise ValueError(report["reason"])
        return result

    def optimize_with_restarts(stage, first_attempt=1):
        nonlocal current
        for attempt in range(first_attempt, first_attempt + settings.get("ts_attempts", 2)):
            result = calculate(stage, attempt)
            if result.get("xyz"):
                current = result["xyz"]
            if result["converged"] and result.get("xyz"):
                return True
        return False

    def optimize_endpoint(stage, endpoint_job):
        """Use the regular endpoint optimization/restart and connectivity checks."""
        nonlocal current
        reference = current
        result = calculate(stage, calculation_job=endpoint_job)
        for attempt in range(2, settings["opt_attempts"] + 1):
            if result["converged"] or not result.get("max_cycles_reached") or not result.get("xyz"):
                break
            current = result["xyz"]
            result = calculate(stage, attempt, calculation_job=endpoint_job)
        if not result["converged"] or not result["xyz"]:
            raise RuntimeError("ORCA endpoint optimization did not converge")
        current = result["xyz"]
        changes = connectivity_changes(molecule_with_xyz(endpoint_job["smiles"], current),
                                       molecule_with_xyz(endpoint_job["smiles"], reference))
        if any(changes.values()):
            raise ValueError(f"Graph change with a distance shift exceeding 2 angstrom: {changes}")

    try:
        original = (root / job["input_xyz"]).read_text()
        current = original
        path = root / Path(job["input_xyz"]).parent / f"job_{job['id']}_run_{job['attempt']:02d}_optimized.xyz"
        frequency_result = None
        if job["endpoint"] != "TS":
            optimize_endpoint("opt", job)
        else:
            strategy = job.get("ts_strategy")
            recovery_only = job.get("ts_mode_recovery") == "relax"
            failure = None
            if strategy != "neb" or recovery_only:
                # One capped preoptimization, then the normal OptTS retry budget.
                # Reaching 50 cycles still supplies the last usable geometry.
                result = calculate("ts_preopt")
                if not result.get("xyz") or not (result["converged"] or result.get("max_cycles_reached")):
                    failure = RuntimeError("Constrained TS preoptimization failed before its 50-cycle limit")
                else:
                    current = result["xyz"]
                    if not optimize_with_restarts("optts"):
                        failure = RuntimeError("OptTS attempts did not converge")
                    else:
                        try:
                            frequency_result = calculate_frequency()
                        except ValueError as exc:
                            failure = exc
            if frequency_result is None:
                # Modern candidates share one group-level NEB job, claimed only
                # after every candidate fails. Historical jobs use inline NEB.
                if strategy == "direct" or not settings["neb_fallback"]:
                    raise failure or RuntimeError("OptTS attempts did not converge")
                update_job(connection, job["id"], xyz_file=None, electronic_energy_eh=None,
                           gibbs_free_energy_eh=None, frequencies_json=None, imaginary_count=None)
                endpoints = []
                ts_folder = Path(job["input_xyz"]).parent
                prefix = Path(job["input_xyz"]).name.removesuffix("_input.xyz")
                for role in ("reactant", "product"):
                    reactive_path = root / ts_folder.parent / role / f"{prefix}_reactive.xyz"
                    current = reactive_path.read_text()
                    endpoint = connection.execute("""SELECT smiles,charge,multiplicity FROM jobs
                        WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint=? LIMIT 1""",
                        (job["ligand_id"], job["pair_index"], job["seating"], role)).fetchone()
                    if endpoint is None:
                        raise ValueError(f"NEB has no sampled {role} graph")
                    endpoint_job = {**job, **dict(endpoint), "endpoint": role}
                    molecule_with_xyz(endpoint_job["smiles"], current)
                    optimize_endpoint(f"neb_{role}_opt", endpoint_job)
                    endpoints.append(current)
                # Retain the original aligned TS guess; endpoints are separately
                # optimized from its original coembedded pair, without overwriting it.
                molecule_with_xyz(job["smiles"], original)
                current = original
                result = calculate("neb", endpoints=endpoints)
                if not result["converged"] or not result["xyz"]:
                    raise RuntimeError("NEB-TS fallback did not converge")
                current = result["xyz"]
                frequency_result = calculate_frequency(2 if recovery_only or strategy != "neb" else 1)
        result = frequency_result if frequency_result is not None else calculate_frequency()
        if job["endpoint"] != "TS" and result["imaginary_count"] > 0:
            tight = calculate("verytightopt")
            if not tight["converged"] or not tight["xyz"]:
                raise RuntimeError("ORCA VeryTightOpt retry did not converge")
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
            update_job(connection, job["id"], status="succeeded", stage="complete", retained=1,
                       ts_mode_recovery=None)
            if job["endpoint"] != "TS":
                select_endpoints(connection, job)
            else:
                select_transition_state(connection, job)
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
    except Exception as exc:
        update_job(connection, job["id"], status="failed", error=str(exc), retained=0)
        print(f"Job {job['id']} ({job['endpoint']}): {exc}", flush=True)


def audit_transition_states(root, *, dry_run=False, job_id=None):
    """Idempotently validate completed TS modes and queue rejected results for recovery."""
    # Dependency errors must stop the audit before touching any calculation rows.
    try:
        import graphrc
    except ImportError as exc:
        raise ValueError("TS mode verification requires graphrc; install the project dependencies") from exc
    root = Path(root).resolve()
    database = root / "workflow.sqlite"
    if not database.is_file():
        raise ValueError("No reaction database exists; prepare/run the workflow first")
    connection = (sqlite3.connect(f"file:{database}?mode=ro", uri=True, timeout=60, isolation_level=None)
                  if dry_run else connect_database(database))
    connection.row_factory = sqlite3.Row
    summary = dict(dry_run=dry_run, graphrc_version=graphrc.__version__, examined=0, checked=0,
                   passed=0, invalid=0, errors=0, cached=0, queued_for_recovery=0,
                   unrecoverable=0, changed_during_audit=0, neb_reactivated=0,
                   backup_file=None, report_file=None, rejected_jobs=[], jobs=[])
    try:
        query = "SELECT * FROM jobs WHERE endpoint='TS' AND status='succeeded'"
        rows = connection.execute(query + (" AND id=?" if job_id is not None else "") + " ORDER BY id",
                                  (job_id,) if job_id is not None else ()).fetchall()
        for original_row in rows:
            job = dict(original_row)
            summary["examined"] += 1
            fingerprint = None
            error = False
            try:
                if not job["out_file"]:
                    raise ValueError("The TS has no NumFreq output file")
                output = root / job["out_file"]
                fingerprint = ts_mode_fingerprint(output, job)
                if (job.get("ts_mode_status") == "passed" and job.get("ts_mode_fingerprint") == fingerprint
                        and job.get("ts_mode_json")):
                    report = json.loads(job["ts_mode_json"])
                    if report.get("valid") is True:
                        summary["cached"] += 1
                        summary["passed"] += 1
                        summary["jobs"].append({**report, "job_id": job["id"], "cached": True})
                        continue
                report = check_ts_mode(output, job)
            except Exception as exc:
                error = True
                report = dict(valid=False, fingerprint=fingerprint, reason=f"TS mode analysis failed: {exc}",
                              output_file=job["out_file"], index_base=0,
                              expected_atoms=[job.get(name) for name in ("donor", "hydrogen", "acceptor")],
                              detected_atoms=[], missing_atoms=[])
            summary["checked"] += 1
            summary["errors" if error else ("passed" if report["valid"] else "invalid")] += 1
            record = dict(job_id=job["id"], cached=False, **report)
            if not report["valid"]:
                record["previous_result"] = {key: job.get(key) for key in (
                    "status", "stage", "attempt", "out_file", "xyz_file", "electronic_energy_eh",
                    "gibbs_free_energy_eh", "imaginary_count", "frequencies_json", "ts_strategy")}
                required = [job.get(name) for name in ("donor", "hydrogen", "acceptor")]
                try:
                    seed = molecule_with_xyz(job["smiles"], (root / job["input_xyz"]).read_text())
                    recoverable = (all(type(index) is int and 0 <= index < seed.GetNumAtoms() for index in required)
                                   and len(set(required)) == 3)
                except (ValueError, OSError):
                    recoverable = False
                record["recovery_available"] = recoverable
                summary["rejected_jobs"].append({key: record[key] for key in (
                    "job_id", "reason", "expected_atoms", "detected_atoms", "missing_atoms", "recovery_available")})
            summary["jobs"].append(record)
            if dry_run:
                continue
            if summary["backup_file"] is None:
                backup_dir = root / "ts_mode_audits" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                backup_dir.mkdir(parents=True, exist_ok=False)
                backup = backup_dir / "workflow.sqlite"
                with closing(sqlite3.connect(backup)) as saved:
                    connection.backup(saved)
                summary["backup_file"] = str(backup)
            connection.execute("BEGIN IMMEDIATE")
            try:
                current = connection.execute("SELECT status,attempt,out_file FROM jobs WHERE id=?", (job["id"],)).fetchone()
                if tuple(current) != (job["status"], job["attempt"], job["out_file"]):
                    summary["changed_during_audit"] += 1
                    connection.rollback()
                    continue
                values = dict(ts_mode_status="error" if error else ("passed" if report["valid"] else "failed"),
                              ts_mode_fingerprint=fingerprint, ts_mode_json=json.dumps(record, sort_keys=True))
                if not report["valid"]:
                    values.update(status="pending" if recoverable else "failed", stage="relax" if recoverable else "ts_mode_check",
                                  ts_mode_recovery="relax" if recoverable else None,
                                  error=None if recoverable else report["reason"] + "; original seed or TS indices unavailable",
                                  retained=None if recoverable else 0, worker_id=None, xyz_file=None,
                                  electronic_energy_eh=None, gibbs_free_energy_eh=None,
                                  imaginary_count=None, frequencies_json=None)
                update_job(connection, job["id"], **values)
                if not report["valid"]:
                    summary["queued_for_recovery" if recoverable else "unrecoverable"] += 1
                    select_transition_state(connection, job)
                    group = (job["ligand_id"], job["pair_index"], job["seating"])
                    # Keep reactivation atomic with rejection: even an interrupted
                    # audit leaves the group ready for its normal fallback.
                    cursor = connection.execute("""UPDATE jobs SET status='pending',stage='neb',retained=NULL,
                        error=NULL,worker_id=NULL,updated_at=CURRENT_TIMESTAMP
                        WHERE ligand_id=? AND pair_index=? AND seating=? AND ts_strategy='neb' AND status='skipped'
                        AND NOT EXISTS (SELECT 1 FROM jobs t WHERE t.ligand_id=? AND t.pair_index=? AND t.seating=?
                                        AND t.endpoint='TS' AND t.status='succeeded')""", (*group, *group))
                    summary["neb_reactivated"] += cursor.rowcount
                connection.commit()
            except BaseException:
                connection.rollback()
                raise
        if not dry_run:
            report_file = root / "ts_mode_audit.json"
            summary["report_file"] = str(report_file)
            report_file.write_text(json.dumps(summary, indent=2) + "\n")
        return summary
    finally:
        connection.close()


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
    parser.add_argument("command", choices=("submit", "resume", "sample", "prepare", "worker", "status", "retry", "precursors", "acid-base", "singlepoint", "verify-ts"))
    parser.add_argument("precursor_command", nargs="?", choices=("submit", "prepare", "worker", "status", "retry", "resume", "plan"),
                        help="Action for the precursors, acid-base or singlepoint workflow")
    parser.add_argument("--conformers", type=int, default=10, help="Fresh conformers per precursor species (default 10)")
    parser.add_argument("--id", type=int, help="Specific row to retry or verify with verify-ts")
    parser.add_argument("--ligand-index", type=int, help="Zero-based CSV row for a sampling array task")
    parser.add_argument("--workers", type=int, help="Workers for submit/resume (default 10 for singlepoint, otherwise 100)")
    parser.add_argument("--references", action="store_true", help="Retry an ID in the shared single-point reference database")
    parser.add_argument("--bound", action="store_true", help="Target the per-run protonated bound complexes for acid-base retry")
    parser.add_argument("--account", help="Slurm account for submit/resume")
    parser.add_argument("--partition", help="Slurm partition for submit/resume")
    parser.add_argument("--dry-run", action="store_true", help="For verify-ts, report checks without changing jobs or writing files")
    args = parser.parse_args(argv)
    if args.workers is None:
        args.workers = 10 if args.command == "singlepoint" else 100
    config = tomllib.loads(args.config.read_text())
    base = args.config.resolve().parent
    root = (base / config["paths"]["results"]).resolve()
    if args.command != "verify-ts":
        root.mkdir(parents=True, exist_ok=True)
    settings = config["calculation"]
    settings.pop("max_iterations", None)  # Ignore the removed override in older config files.
    settings.setdefault("opt_attempts", settings.get("ts_attempts", 2))
    for key in ("cpus", "memory_mb", "opt_attempts", "ts_attempts", "timeout_seconds"):
        if settings[key] <= 0:
            parser.error(f"calculation.{key} must be positive")
    if args.command == "singlepoint":
        from hpc_singlepoint import handle_singlepoints

        return handle_singlepoints(args, config, base, root, parser)
    if args.command in ("precursors", "acid-base"):
        from hpc_precursors import handle_precursors

        return handle_precursors(args, config, base, root, parser)
    if args.precursor_command is not None:
        parser.error("The extra action argument is only supported for precursors, acid-base or singlepoint")
    if args.command == "verify-ts":
        try:
            result = audit_transition_states(root, dry_run=args.dry_run, job_id=args.id)
            print(json.dumps({key: value for key, value in result.items() if key != "jobs"}, indent=2))
        except (ValueError, OSError, sqlite3.Error) as exc:
            parser.error(str(exc))
        return
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
            ("endpoint_conformers", 3), ("coembed_conformers", 10), ("ts_candidates", 5))):
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
