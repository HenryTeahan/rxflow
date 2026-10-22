"""Import the six checked stereo repairs, preserving all other workflow rows."""
import json
from datetime import datetime, timezone
from pathlib import Path
import shutil
import sqlite3
import sys
import tomllib

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from hpc_workflow import connect_database, prepare_jobs


def main():
    experiment = BASE / 'experiments/rxembed_730eae5_stereo_fixed'
    summary = json.loads((experiment / 'summary.json').read_text())
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_stereo_repair'
    reports = []
    for config_name in ('hpc.toml', 'hpc_ED_dF.toml', 'hpc_ED_pyr.toml'):
        config = tomllib.loads((BASE / config_name).read_text())
        source = BASE / config['paths']['source']
        root = BASE / config['paths']['results']
        cases = [case for case in summary['cases'] if case['system'] == root.name]
        assert len(cases) == 2 and all(case['succeeded'] for case in cases)
        expected = {case['job_id'] for case in cases}
        connection = connect_database(root / 'workflow.sqlite')
        before = {row['id']: dict(row) for row in connection.execute('SELECT * FROM jobs')}
        assert not any(row['status'] in ('running', 'pending') for row in before.values())
        assert all(before[job]['status'] == 'failed' and before[job]['stage'] == 'sampling'
                   for job in expected)
        archive = root / 'ts_workflow_resets' / stamp
        archive.mkdir(parents=True, exist_ok=False)
        with sqlite3.connect(archive / 'workflow.sqlite') as backup:
            connection.backup(backup)
        ligand = source / 'ligand_03'
        shutil.copy2(ligand / 'sampling.json', archive / 'sampling.json')
        for case in cases:
            folder = experiment / root.name / f"job_{case['job_id']}"
            prefix = f"pair_{case['pair_index']:03d}_seating_04"
            selection = json.loads((folder / 'ts_guess/reaction_selection.json').read_text())
            assert len(selection['candidates']) == 5
            assert all(candidate['valid'] for candidate in selection['candidates'])
            paths = {}
            for role in ('reactant', 'product', 'ts_guess'):
                for path in (folder / role).iterdir():
                    assert path.is_file() and path.name.startswith('reaction')
                    target = ligand / role / (prefix + path.name[len('reaction'):])
                    assert not target.exists(), target
                    paths[str(path.resolve())] = str(target.resolve())

            def relocate(value):
                if isinstance(value, dict):
                    return {key: relocate(item) for key, item in value.items()}
                if isinstance(value, list):
                    return [relocate(item) for item in value]
                return paths.get(value, value) if isinstance(value, str) else value

            for original, target in paths.items():
                original, target = Path(original), Path(target)
                if original.suffix == '.json':
                    data = relocate(json.loads(original.read_text()))
                    if original.name == 'reaction_selection.json':
                        data['endpoint_stereo'] = dict(zip(('reactant', 'product'), case['endpoint_stereo']))
                    target.write_text(json.dumps(data, indent=2) + '\n')
                else:
                    shutil.copy2(original, target)
        sampling = json.loads((ligand / 'sampling.json').read_text())
        repaired = {(case['pair_index'], case['seating']) for case in cases}
        old_errors = sampling['errors']
        sampling['errors'] = [error for error in old_errors if
                              (error.get('pair'), error.get('seating')) not in repaired]
        sampling['groups'] += len(cases)
        sampling['status'] = 'completed_with_errors' if sampling['errors'] else 'complete'
        (ligand / 'sampling.json').write_text(json.dumps(sampling, indent=2) + '\n')
        settings = config['calculation'].copy()
        settings.pop('max_iterations', None)
        settings.setdefault('opt_attempts', settings.get('ts_attempts', 2))
        connection.execute('BEGIN IMMEDIATE')
        try:
            added, groups = prepare_jobs(connection, source, root, settings)
            after = {row['id']: dict(row) for row in connection.execute('SELECT * FROM jobs')}
            changed = {job for job in before if before[job] != after[job]}
            new = [row for job, row in after.items() if job not in before]
            assert changed == expected, changed
            assert added == len(new) == 10, (added, len(new))
            assert all(row['endpoint'] == 'TS' and row['ts_strategy'] == 'direct'
                       and row['status'] == 'pending' for row in new)
            assert all(after[job]['status'] == 'pending' and after[job]['ts_strategy'] == 'neb'
                       for job in expected)
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
        finally:
            connection.close()
        report = dict(config=config_name, archive=str(archive), groups=groups,
                      added_direct_ids=[row['id'] for row in new], repaired_fallback_ids=sorted(expected),
                      prior_jobs=[before[job] for job in sorted(expected)],
                      other_existing_rows_unchanged=True)
        (archive / 'repair.json').write_text(json.dumps(report, indent=2) + '\n')
        reports.append({key: value for key, value in report.items() if key != 'prior_jobs'})
        print(json.dumps(reports[-1], indent=2), flush=True)
    (experiment / 'import.json').write_text(json.dumps(reports, indent=2) + '\n')


if __name__ == '__main__':
    main()
