"""Retest the six original stereo sampling failures without changing workflow rows."""
import argparse
import csv
from importlib.metadata import distribution
import json
import os
from pathlib import Path
import sqlite3
import sys
import time
import traceback

os.environ.setdefault("OMP_NUM_THREADS", "1")
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from rxflow_endpoints import create_endpoint_isomers, generate_ts_guess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='experiments/rxembed_730eae5_stereo_check')
    parser.add_argument('--keep', type=int, default=1)
    args = parser.parse_args()
    out = BASE / args.output
    out.mkdir(exist_ok=False)
    summary = dict(rxembed=json.loads(distribution('rxembed').read_text('direct_url.json')), cases=[])
    for root, seeds in [('hpc_tm_results', 'hpc_tm_seeds'), ('ED_dF_hpc_results', 'ED_dF_hpc_seeds'),
                        ('ED_pyr_hpc_results', 'ED_pyr_hpc_seeds')]:
        with sqlite3.connect(f'file:{BASE/root/"workflow.sqlite"}?mode=ro', uri=True) as conn:
            conn.row_factory = sqlite3.Row
            rows = [dict(row) for row in conn.execute(
                "SELECT * FROM jobs WHERE endpoint='TS' AND stage='sampling' AND status='failed' ORDER BY id")]
        for row in rows:
            started = time.monotonic()
            case = dict(system=root, job_id=row['id'], pair_index=row['pair_index'],
                        seating=row['seating'], previous_error=row['error'])
            print(f"Testing {root} job {row['id']}", flush=True)
            try:
                source = BASE / seeds / f"ligand_{int(row['ligand_id']):02d}"
                with (source/'reaction_pairs.csv').open() as handle:
                    pairs = list(csv.DictReader(handle))
                pair = pairs[row['pair_index']-1]
                left, right = next((left,right) for seating,left,right in create_endpoint_isomers(pair)
                                   if seating+1 == row['seating'])
                folder = out / root / f"job_{row['id']}"
                result = generate_ts_guess(left, right, folder, n=10, seed=12345,
                                           multiplicity=row['multiplicity'], keep=args.keep)
                case.update(succeeded=True, guess=str(folder/'ts_guess/reaction.xyz'),
                            selected_ts_rmsd=result['selected_ts_rmsd'], ts_mode=result['ts_mode'],
                            endpoint_stereo=[left.stereo_label, right.stereo_label],
                            candidate_count=len(result.get('candidates', [])))
            except Exception as exc:
                case.update(succeeded=False, error=str(exc), traceback=traceback.format_exc())
            case['seconds'] = time.monotonic()-started
            summary['cases'].append(case)
            (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
            print(json.dumps(case,indent=2),flush=True)
    print('Succeeded:',sum(case['succeeded'] for case in summary['cases']), '/',len(summary['cases']),flush=True)


if __name__ == '__main__':
    main()
