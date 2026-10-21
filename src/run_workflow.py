"""Run ligand workflows in parallel, with persistent progress in results/."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import shutil
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import rxembed

from rxflow import generate_reaction_pairs, ligand_directory, read_ligands, rxembed_revision
from rxflow_endpoints import sample_reaction_endpoints


def save_status(path, status):
    """Replace a progress file atomically so it can be read during calculations."""
    status["updated_at"] = datetime.now(timezone.utc).isoformat()
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(status, indent=2) + "\n")
    temporary.replace(path)


def run_ligand(row, results, executable, conformers, coembed_conformers):
    """Screen one ligand and sample every retained reaction pair."""
    output = ligand_directory(results, row["ID"])
    output.mkdir(parents=True, exist_ok=True)
    status = {"ligand_id": row["ID"], "status": "running", "stage": "screening",
              "pairs_total": 0, "pairs_completed": 0, "ts_guesses": 0, "errors": []}
    save_status(output / "status.json", status)
    try:
        pairs = generate_reaction_pairs(row, output, executable=executable)
        status.update(stage="endpoints", pairs_total=len(pairs))
        save_status(output / "status.json", status)
        for index, pair in enumerate(pairs, 1):
            status["current_pair"] = index
            save_status(output / "status.json", status)
            try:
                records = sample_reaction_endpoints(pair, output, executable, pair_index=index,
                    n=conformers, coembed_conformers=coembed_conformers)
                status["ts_guesses"] += sum(r.get("guesses", 0) for r in records)
                for record in records:
                    if record.get("error") or record.get("retained") == 0:
                        status["errors"].append({"pair": index, **record})
                if not records:
                    status["errors"].append({"pair": index, "error": "No compatible endpoint seating"})
            except (RuntimeError, ValueError, OSError) as exc:
                failure = {"pair": index, "error": str(exc)}
                status["errors"].append(failure)
                save_status(output / f"pair_{index:03d}_failure.json", failure)
            status["pairs_completed"] = index
            save_status(output / "status.json", status)
        status.update(stage="finished", status="completed_with_errors" if status["errors"] else
                      "complete" if pairs else "no_candidates")
    except Exception as exc:
        logging.exception("Ligand %s failed", row["ID"])
        status.update(status="failed", error=str(exc))
    save_status(output / "status.json", status)
    return status


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligands-csv", type=Path, default=Path("ligands.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--xtb", default=os.environ.get("XTB_EXE", "xtb"))
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--coembed-conformers", type=int, default=10)
    parser.add_argument("--workers", type=int, help="Default: one worker per ligand")
    args = parser.parse_args(argv)
    rows = read_ligands(args.ligands_csv)
    workers = args.workers if args.workers is not None else len(rows)
    if not rows or min(workers, args.conformers, args.coembed_conformers) < 1:
        parser.error("At least one ligand and positive counts are required")
    executable = shutil.which(args.xtb)
    if executable is None:
        parser.error("Set XTB_EXE or --xtb to your g-xTB capable xTB executable")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "run.json").exists() or any(output.glob("ligand_*")):
        parser.error(f"A workflow already exists in {output}; choose an empty output directory")
    revision = rxembed_revision()
    run = {"status": "running", "pid": os.getpid(), "python": sys.executable,
           "rxembed_source": rxembed.__file__, "rxembed_commit": revision,
           "ligands": len(rows), "workers": workers, "threads_per_worker": 1,
           "endpoint_conformers": args.conformers, "coembed_conformers": args.coembed_conformers,
           "solvent": None, "xtb": executable, "ligands_completed": 0, "results": []}
    save_status(output / "run.json", run)
    logging.info("Starting %d ligands with %d workers; %d endpoint / %d coembedded conformers",
                 len(rows), workers, args.conformers, args.coembed_conformers)
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_ligand, row, output, executable, args.conformers,
                                   args.coembed_conformers): row["ID"] for row in rows}
            for future in as_completed(futures):
                ligand = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {"ligand_id": ligand, "status": "failed", "error": str(exc)}
                run["results"].append(result)
                run["ligands_completed"] += 1
                save_status(output / "run.json", run)
                logging.info("Ligand %s: %s (%d/%d finished)", ligand, result["status"],
                             run["ligands_completed"], len(rows))
        run["status"] = "complete" if all(r["status"] == "complete" for r in run["results"]) else "completed_with_errors"
    except BaseException as exc:
        run.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed", error=str(exc))
        raise
    finally:
        save_status(output / "run.json", run)
    return 0 if run["status"] == "complete" else 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(processName)s %(levelname)s %(message)s")
    raise SystemExit(main())
