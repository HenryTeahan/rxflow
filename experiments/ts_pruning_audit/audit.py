"""Read-only audit of embedded and optimized TS diversity; never modifies queues."""
from collections import Counter, defaultdict
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import sqlite3
import sys
import tomllib

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

BASE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BASE / 'src'))
from conformer_sampling import prune_conformers
from orca_hpc import molecule_with_xyz


def read_rows(path, sql):
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True, timeout=5)) as db:
        db.row_factory = sqlite3.Row
        return [dict(row) for row in db.execute(sql)]


def pair_rmsd(left, right):
    atoms = [(atom.GetIdx(), atom.GetIdx()) for atom in left.GetAtoms() if atom.GetAtomicNum() != 1]
    return rdMolAlign.GetAlignmentTransform(left, right, atomMap=atoms, reflect=False)[0]


def main():
    summary = dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                   metric='Aligned heavy-atom RMSD; original atom mapping; separate molecular identities',
                   energy_order='Ascending g-xTB NumFreq Gibbs free energy', runs=[])
    details = []
    for config_name in ('hpc.toml', 'hpc_ED_dF.toml', 'hpc_ED_pyr.toml', 'hpc_ni.toml'):
        config = tomllib.loads((BASE / config_name).read_text())
        root = (BASE / config['paths']['results']).resolve()
        rows = read_rows(root / 'workflow.sqlite', "SELECT * FROM jobs WHERE endpoint='TS'")
        # Audit the energy window before the production optimized-geometry prune.
        energy_groups = defaultdict(list)
        for row in rows:
            if row['status'] == 'succeeded' and row['gibbs_free_energy_eh'] is not None:
                energy_groups[(row['ligand_id'], row['pair_index'], row['seating'])].append(row)
        selected = []
        for group in energy_groups.values():
            minimum = min(row['gibbs_free_energy_eh'] for row in group)
            selected.extend(row for row in group if (row['gibbs_free_energy_eh'] - minimum) * 627.509474 <= 3. + 1e-9)
        groups = defaultdict(list)
        for row in selected:
            groups[(row['ligand_id'], row['pair_index'], row['seating'])].append(row)
        totals = Counter()
        mismatched_seed_groups = []
        sample_reports = [json.loads(row['ts_sampling_json']) for row in rows
                          if row['ts_strategy'] == 'sampling' and row['status'] == 'skipped'
                          and row['ts_sampling_json']]
        embedded = sum(report['embedded'] for report in sample_reports)
        retained = sum(report['retained'] for report in sample_reports)
        seed_groups = defaultdict(list)
        for row in rows:
            if row['ts_strategy'] == 'sampled':
                seed_groups[row['ts_parent_id']].append(row)
        seed_counts = Counter()
        for children in seed_groups.values():
            children.sort(key=lambda row: -row['conformer'])
            mols = [molecule_with_xyz(row['smiles'], (root / row['input_xyz']).read_text()) for row in children]
            count = len(prune_conformers(mols, 0.5))
            seed_counts.update(original=len(mols), repruned=count)
            if count != len(mols):
                mismatched_seed_groups.append(dict(parent_id=children[0]['ts_parent_id'], retained=len(mols), repruned=count))
        removed_ids = set()
        representatives = []
        for key, members in sorted(groups.items()):
            members.sort(key=lambda row: (row['gibbs_free_energy_eh'], row['id']))
            mols = [molecule_with_xyz(row['smiles'], (root / row['xyz_file']).read_text()) for row in members]
            lookup = {id(mol): row for mol, row in zip(mols, members)}
            counts = {}
            kept_half = []
            for threshold in (0.1, 0.25, 0.5):
                kept = prune_conformers(mols, threshold)
                counts[str(threshold)] = len(kept)
                totals[str(threshold)] += len(kept)
                if threshold == 0.5:
                    kept_half = kept
            keep_ids = {lookup[id(mol)]['id'] for mol in kept_half}
            removed_ids.update(row['id'] for row in members if row['id'] not in keep_ids)
            closest = None
            for i, mol in enumerate(mols):
                for j in range(i):
                    # Exact atom order is required; different labels/graphs are preserved by prism's wrapper.
                    rmsd = pair_rmsd(mol, mols[j])
                    if closest is None or rmsd < closest['rmsd_angstrom']:
                        closest = dict(job_ids=[members[j]['id'], members[i]['id']], rmsd_angstrom=rmsd,
                                       delta_G_kcal=abs(members[i]['gibbs_free_energy_eh'] - members[j]['gibbs_free_energy_eh']) * 627.509474)
            record = dict(run=config_name, ligand=key[0], pair=key[1], seating=key[2],
                          selected=len(members), unique_0p1=counts['0.1'], unique_0p25=counts['0.25'],
                          unique_0p5=counts['0.5'], kept_job_ids=sorted(keep_ids),
                          unique_smiles=len({Chem.MolToSmiles(mol) for mol in mols}), closest_pair=closest)
            details.append(record)
            if closest and len(members) >= 10:
                representatives.append(record)
        dft = read_rows(root / 'singlepoints.sqlite', "SELECT source_job_id,status FROM jobs WHERE endpoint='TS'")
        duplicate_status = Counter(row['status'] for row in dft if row['source_job_id'] in removed_ids)
        representatives.sort(key=lambda row: (-row['selected'], row['unique_0p5']))
        report = dict(run=config_name, groups=len(groups), published_embedding_batches=len(sample_reports),
                      embedded=embedded, retained_embedded=retained,
                      retained_seeds_rechecked=dict(seed_counts), seed_groups_with_remaining_duplicates=mismatched_seed_groups,
                      selected_in_3kcal=len(selected), optimized_unique_by_RMSD=dict(totals),
                      potentially_redundant_DFT_status=dict(duplicate_status),
                      examples=representatives[:3])
        summary['runs'].append(report)
        print(json.dumps(report), flush=True)
    out = Path(__file__).resolve().parent
    (out / 'report.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out / 'groups.json').write_text(json.dumps(details, indent=2) + '\n')
    with (out / 'groups.csv').open('w', newline='') as handle:
        columns = ['run', 'ligand', 'pair', 'seating', 'selected', 'unique_0p1', 'unique_0p25', 'unique_0p5', 'unique_smiles']
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(details)


if __name__ == '__main__':
    main()
