"""Audit optimized substrate coordination without modifying workflow databases."""
import argparse
from collections import Counter, defaultdict
import csv
from importlib.metadata import version
import json
from pathlib import Path
import sqlite3
import sys

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
import networkx as nx
import numpy as np
from rdkit import Chem
import xyzgraph
from hpc_singlepoint import selected_jobs
from orca_hpc import molecule_with_xyz
from rxembed.pipeline.metrics import coordination_changed
from rxflow_endpoints import connectivity_changes, mapped_molecule


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=BASE / 'hpc_tm_results')
    parser.add_argument('--ligand', default='11')
    parser.add_argument('--output', type=Path, default=BASE / 'experiments/ligand_11_binding_audit')
    args = parser.parse_args()
    root = args.results.resolve()
    with sqlite3.connect(f'file:{root / "workflow.sqlite"}?mode=ro', uri=True) as db:
        db.row_factory = sqlite3.Row
        rows = [dict(row) for row in db.execute('SELECT * FROM jobs WHERE ligand_id=? ORDER BY id',
                                                (args.ligand,))]
    selected = {row['id'] for row in selected_jobs(rows)}
    groups = defaultdict(list)
    for row in rows:
        groups[(row['pair_index'], row['seating'], row['binding_permutation'])].append(row)
    records = []
    periodic = Chem.GetPeriodicTable()
    for key, group in groups.items():
        reactant = next(row for row in group if row['endpoint'] == 'reactant')
        reference = mapped_molecule(reactant['smiles'])
        dative = [bond for bond in reference.GetBonds() if bond.GetBondType() == Chem.BondType.DATIVE]
        metal, = {bond.GetEndAtomIdx() for bond in dative}
        ligand_donors = {number - 1 for number in json.loads(reactant['donor_maps'])}
        substrate_donor, = {bond.GetBeginAtomIdx() for bond in dative} - ligand_donors
        disconnected = Chem.RWMol(reference)
        for bond in dative:
            disconnected.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        substrate, = [set(fragment) for fragment in Chem.GetMolFrags(disconnected.GetMol())
                       if substrate_donor in fragment]
        heavy = sorted(i for i in substrate if reference.GetAtomWithIdx(i).GetAtomicNum() > 1)
        for row in group:
            if row['status'] != 'succeeded' or not row['xyz_file']:
                continue
            molecule = molecule_with_xyz(row['smiles'], (root / row['xyz_file']).read_text())
            original = molecule_with_xyz(row['smiles'], (root / row['input_xyz']).read_text())
            positions = molecule.GetConformer().GetPositions()
            initial = original.GetConformer().GetPositions()
            distance = lambda i: float(np.linalg.norm(positions[metal] - positions[i]))
            before = float(np.linalg.norm(initial[metal] - initial[substrate_donor]))
            atoms = [(atom.GetSymbol(), tuple(positions[atom.GetIdx()])) for atom in molecule.GetAtoms()]
            graph = xyzgraph.build_graph(atoms, charge=row['charge'], multiplicity=row['multiplicity'], quick=True)
            contacts = [i for i in heavy if graph.has_edge(metal, i)]
            radii = {i: periodic.GetRcovalent(molecule.GetAtomWithIdx(i).GetAtomicNum()) for i in [metal, *heavy]}
            near = [i for i in heavy if distance(i) <= 1.3 * (radii[metal] + radii[i])]
            donors = [bond.GetBeginAtomIdx() for bond in molecule.GetBonds()
                      if bond.GetBondType() == Chem.BondType.DATIVE and bond.GetEndAtomIdx() == metal]
            left, joined = coordination_changed(molecule, molecule.GetConformer().GetId(), metal, donors)
            if not contacts and not near:
                state = 'uncoordinated_substrate'
            elif substrate_donor in contacts:
                state = 'intended_coordination'
            elif contacts and substrate_donor in left:
                state = 'changed_binding_site'
            else:
                state = 'borderline_or_perception_disagreement'
            nearest = min(heavy, key=distance)
            current = connectivity_changes(molecule, original) if row['endpoint'] != 'TS' else None
            metal_component = nx.node_connected_component(graph, metal)
            record = {name: row[name] for name in ('id', 'pair_index', 'seating', 'binding_permutation',
                                                   'substrate_site', 'resonance', 'endpoint', 'conformer',
                                                   'retained', 'gibbs_free_energy_eh', 'input_xyz', 'xyz_file')}
            record.update(selected_for_notebook=row['id'] in selected, classification=state,
                          metal_index=metal, substrate_donor_index=substrate_donor,
                          initial_metal_site_distance_angstrom=before,
                          final_metal_site_distance_angstrom=distance(substrate_donor),
                          metal_site_distance_change_angstrom=distance(substrate_donor) - before,
                          nearest_substrate_heavy_atom_index=nearest,
                          nearest_substrate_heavy_distance_angstrom=distance(nearest),
                          xyzgraph_metal_substrate_contacts=contacts,
                          radius_supported_substrate_contacts=near,
                          xyzgraph_substrate_connected_to_metal=bool(set(heavy) & metal_component),
                          native_coordination_lost=left, native_coordination_joined=joined,
                          new_ligand_contacts=[dict(index=i, symbol=molecule.GetAtomWithIdx(i).GetSymbol(),
                                                    distance_angstrom=distance(i))
                                               for i in joined if i not in substrate],
                          current_workflow_connectivity_check=current)
            records.append(record)
    counts = Counter((r['endpoint'], r['classification']) for r in records)
    selected_records = [r for r in records if r['selected_for_notebook']]
    report = dict(results=str(root), ligand_id=args.ligand, index_base=0,
                  xyzgraph_version=version('xyzgraph'), rxembed_version=version('rxembed'),
                  method=dict(graph='xyzgraph.build_graph, quick=True, default geometry parameters',
                              distance_check='Any substrate heavy atom within 1.3 times the metal/atom covalent-radius sum',
                              uncoordinated='No graph metal-substrate heavy-atom edge AND no radius-supported substrate contact',
                              borderline='Retained separately; never classified as dissociated from graph absence alone'),
                  successful_structures=len(records),
                  counts=[dict(endpoint=role, classification=state, count=count)
                          for (role, state), count in sorted(counts.items())],
                  selected_states=selected_records, structures=records)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    for name, data in [('conformers.csv', records), ('selected_states.csv', selected_records)]:
        with (args.output / name).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows({key: json.dumps(value) if isinstance(value, (dict, list)) else value
                             for key, value in row.items()} for row in data)
    print(json.dumps({key: value for key, value in report.items() if key not in ('structures', 'selected_states')}, indent=2))
    for record in selected_records:
        print(record['id'], record['pair_index'], record['endpoint'], record['classification'],
              f"Pd-site={record['final_metal_site_distance_angstrom']:.3f}",
              f"nearest={record['nearest_substrate_heavy_distance_angstrom']:.3f}", flush=True)


if __name__ == '__main__':
    main()
