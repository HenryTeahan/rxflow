"""Uncoordinated substrate encounters with each deprotonated ligand–metal arrangement."""
from contextlib import closing
from datetime import datetime, timezone
import json
import math
import sqlite3

from rdkit import Chem
import rxembed
from conformer_sampling import embedded_conformers, prune_conformers

from orca_hpc import molecule_with_xyz
from rxflow import read_ligands
from rxflow_endpoints import mapped_molecule


def encounter_plan(config, base, root):
    from hpc_precursors import canonical_substrate, precursor_plan
    rows = {row['ID']: row for row in read_ligands(base / config['paths']['ligands_csv'])}
    return [dict(spec, endpoint='encounter', substrate=canonical_substrate(rows[spec['ligand_id']]['substrate']))
            for spec in precursor_plan(config, base, root) if spec['endpoint'] == 'deprotonated']


def encounter_graph(spec):
    from hpc_precursors import precursor_graph
    from rxembed.stereo import defined_stereo_label
    _, native = precursor_graph(dict(spec, endpoint='deprotonated'))
    core = mapped_molecule(rxembed.cxsmiles(native))
    molecule = Chem.CombineMols(core, Chem.AddHs(Chem.MolFromSmiles(spec['substrate'])))
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    Chem.SanitizeMol(molecule)
    sites = {slot: donor if isinstance(donor, int) else donor.atoms[0]
             for slot, donor in enumerate(native.centres[0].vertices) if donor is not None}
    source = rxembed.Isomer(molecule, spec['metal_geom'], sites)
    source.stereo_label = defined_stereo_label(molecule)
    if Chem.GetFormalCharge(molecule) != 1 or len(Chem.GetMolFrags(molecule)) != 2:
        raise ValueError('An encounter must contain a +1 ligand–metal complex and a neutral substrate')
    return molecule, source


def sample_encounter(spec, *, n=10, seed=12345, prune_rmsd=0.5):
    """Sample native noncovalent contact modes without adding a substrate–metal bond."""
    molecule, source = encounter_graph(spec)
    modes = rxembed.nci_modes(molecule, kinds=('HB', 'CHPI', 'HBPI', 'ANPI', 'CATLP'),
                             inter_fragment=True, max_modes=min(n, 4), seed=seed)
    if not modes:
        raise ValueError('No interfragment encounter contact modes were found')
    batches, failures = [], []
    for index, contact in enumerate(modes.values()):
        try:
            # Several modes can fail: retain the others and request enough seeds
            # for ten total final jobs even when only one mode survives.
            embedded = rxembed.embed(source, contacts=contact, n=n, seed=seed + index, threads=1)
            batches.append([mol for batch in embedded_conformers(embedded, prune_rmsd) for mol in batch])
        except (ValueError, RuntimeError, TimeoutError) as exc:
            failures.append(str(exc))
            batches.append([])
    conformers = []
    for rank in range(max(map(len, batches), default=0)):
        for batch in batches:
            if rank < len(batch):
                geometry = batch[rank]
                if ([a.GetAtomicNum() for a in geometry.GetAtoms()] != [a.GetAtomicNum() for a in molecule.GetAtoms()]
                        or [a.GetAtomMapNum() for a in geometry.GetAtoms()] != list(range(1, molecule.GetNumAtoms() + 1))
                        or Chem.GetFormalCharge(geometry) != spec['charge']):
                    raise ValueError('Encounter embedding changed atom inventory, order or charge')
                conformers.append(geometry)
    if not conformers and failures:
        raise ValueError('All encounter contact modes failed: ' + '; '.join(failures))
    return molecule, prune_conformers(conformers, prune_rmsd)[:n]


def uncoordinated_substrate(molecule):
    """Qualify old detached reactants for reuse; never filter encounter optimizations."""
    import xyzgraph
    metals = {b.GetEndAtomIdx() for b in molecule.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE}
    if len(metals) != 1:
        return False
    metal, = metals
    core = next(set(f) for f in Chem.GetMolFrags(molecule) if metal in f)
    substrate = set(range(molecule.GetNumAtoms())) - core
    if not substrate:
        return False
    positions = molecule.GetConformer().GetPositions()
    graph = xyzgraph.build_graph([(a.GetSymbol(), tuple(positions[a.GetIdx()])) for a in molecule.GetAtoms()],
                                charge=Chem.GetFormalCharge(molecule), quick=True)
    periodic = Chem.GetPeriodicTable()
    for i in substrate:
        atom = molecule.GetAtomWithIdx(i)
        if atom.GetAtomicNum() <= 1:
            continue
        distance = (molecule.GetConformer().GetAtomPosition(metal) - molecule.GetConformer().GetAtomPosition(i)).Length()
        if graph.has_edge(metal, i) or distance <= 1.3 * (
                periodic.GetRcovalent(molecule.GetAtomWithIdx(metal).GetAtomicNum())
                + periodic.GetRcovalent(atom.GetAtomicNum())):
            return False
    return True


def import_encounter_minima(config, root, plan, *, conformers=10, ligand_id=None):
    """Reuse verified detached optimizations, preserving their original files and provenance."""
    from hpc_precursors import canonical_substrate, connect_precursors, precursor_location
    database = root / 'workflow.sqlite'
    with closing(sqlite3.connect(f'file:{database}?mode=ro', uri=True)) as db:
        db.row_factory = sqlite3.Row
        columns = {r[1] for r in db.execute('PRAGMA table_info(jobs)')}
        if not {'gibbs_free_energy_eh', 'xyz_file', 'imaginary_count'} <= columns:
            return 0
        candidates = [dict(r) for r in db.execute("""SELECT * FROM jobs WHERE endpoint='reactant'
            AND (status='succeeded' OR (status='failed' AND stage='connectivity_filter'))
            AND imaginary_count=0 AND gibbs_free_energy_eh IS NOT NULL
            ORDER BY gibbs_free_energy_eh,id""")]
    settings = {k: v for k, v in config['calculation'].items() if k != 'max_iterations'}
    settings.setdefault('opt_attempts', settings.get('ts_attempts', 2))
    imported = 0
    with closing(connect_precursors(root / 'encounters.sqlite', encounter=True)) as db:
        db.execute('''CREATE TABLE IF NOT EXISTS encounter_imports (
            job_id INTEGER PRIMARY KEY, source_database TEXT NOT NULL, source_job_id INTEGER NOT NULL,
            UNIQUE(source_database,source_job_id))''')
        for spec in plan:
            if ligand_id is not None and str(spec['ligand_id']) != str(ligand_id):
                continue
            if db.execute('SELECT 1 FROM jobs WHERE ligand_id=? AND pair_index=? LIMIT 1',
                          (spec['ligand_id'], spec['binding_permutation'])).fetchone():
                continue
            accepted = []
            native, _ = encounter_graph(spec)
            for row in candidates:
                if (row['ligand_id'], row['binding_permutation']) != (spec['ligand_id'], spec['binding_permutation']):
                    continue
                if not row['xyz_file'] or not all(row.get(k) is not None and math.isfinite(row[k])
                                               for k in ('electronic_energy_eh', 'gibbs_free_energy_eh')):
                    continue
                try:
                    original = mapped_molecule(row['smiles'])
                    metal = next(b.GetEndAtomIdx() for b in native.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE)
                    core_size = len(next(f for f in Chem.GetMolFrags(native) if metal in f))
                    if ([a.GetAtomicNum() for a in list(original.GetAtoms())[:core_size]] !=
                            [a.GetAtomicNum() for a in list(native.GetAtoms())[:core_size]]
                            or row['charge'] != spec['charge'] or row['multiplicity'] != spec['multiplicity']
                            or json.loads(row['settings_json']) != settings):
                        continue
                    def core_connectivity(mol):
                        return {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                                 max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                                 b.GetBondType() == Chem.BondType.DATIVE)
                                for b in mol.GetBonds()
                                if b.GetBeginAtomIdx() < core_size and b.GetEndAtomIdx() < core_size}
                    if core_connectivity(original) != core_connectivity(native):
                        continue
                    editable = Chem.RWMol(original)
                    for bond in original.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DATIVE and bond.GetBeginAtomIdx() >= core_size:
                            editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    graph = editable.GetMol()
                    fragments = Chem.GetMolFrags(graph, asMols=True)
                    if (len(fragments) != 2 or canonical_substrate(Chem.MolToSmiles(fragments[1])) != spec['substrate']):
                        continue
                    from rxembed.stereo import defined_stereo_label
                    # Legacy graphs missing the screened donor hand cannot be imported.
                    native_core = Chem.GetMolFrags(native, asMols=True)[0]
                    if defined_stereo_label(fragments[0]) != defined_stereo_label(native_core):
                        continue
                    smiles = Chem.MolToSmiles(graph)
                    xyz = (root / row['xyz_file']).read_text()
                    current = molecule_with_xyz(smiles, xyz)
                    if not uncoordinated_substrate(current):
                        continue
                except (ValueError, RuntimeError, OSError):
                    continue
                accepted.append((row, smiles, xyz))
                if len(accepted) == conformers:
                    break
            for number, (row, smiles, xyz) in enumerate(accepted, 1):
                relative = precursor_location(spec, number)
                target = root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists() and target.read_text() != xyz:
                    raise ValueError(f'Existing encounter geometry differs: {target}')
                target.write_text(xyz)
                optimized = relative.with_name(relative.stem.removesuffix('_input') + '_imported.xyz')
                (root / optimized).write_text(xyz)
                cursor = db.execute('''INSERT INTO jobs (ligand_id,pair_index,seating,endpoint,conformer,
                    input_xyz,xyz_file,smiles,charge,multiplicity,settings_json,binding_permutation,
                    binding_mode,coordination_isomer,donor_maps,substrate,status,stage,
                    electronic_energy_eh,gibbs_free_energy_eh,frequencies_json,imaginary_count,out_file)
                    VALUES (?,?,0,'encounter',?,?,?,?,?,?,?,?,?,?,?,?,'succeeded','complete',?,?,?,?,?)''',
                    (spec['ligand_id'], spec['binding_permutation'], number, str(relative), str(optimized),
                     smiles, row['charge'], row['multiplicity'], json.dumps(settings, sort_keys=True),
                     spec['binding_permutation'], spec['binding_mode'], spec['coordination_isomer'],
                     spec['donor_maps'], spec['substrate'], row['electronic_energy_eh'], row['gibbs_free_energy_eh'],
                     row['frequencies_json'], 0, row['out_file']))
                db.execute('INSERT INTO encounter_imports VALUES (?,?,?)',
                           (cursor.lastrowid, str(database.resolve()), row['id']))
                imported += 1
    return imported


def prepare_encounters(config, base, root, *, conformers=10, ligand_id=None):
    from hpc_precursors import encounter_protocol_extension, prepare_species, species_protocol
    if conformers < 1:
        raise ValueError('Encounter conformer count must be positive')
    plan = encounter_plan(config, base, root)
    if ligand_id is not None and not any(str(s['ligand_id']) == str(ligand_id) for s in plan):
        raise ValueError(f'No encounter binding arrangements for ligand {ligand_id}')
    # Validate before importing or changing any queue rows. Only an upward
    # conformer-count change is allowed for an existing encounter database.
    database = root / 'encounters.sqlite'
    if database.is_file():
        with closing(sqlite3.connect(f'file:{database}?mode=ro', uri=True)) as db:
            frozen = db.execute("SELECT value FROM precursor_metadata WHERE key='protocol'").fetchone()
        requested = species_protocol(config, plan, [], conformers)[2]
        if frozen and frozen[0] != requested and not encounter_protocol_extension(frozen[0], requested):
            raise ValueError('Encounter chemistry/settings changed; use a fresh results directory')
    requeued = requeue_filtered_encounters(root, ligand_id=ligand_id)
    imported = import_encounter_minima(config, root, plan, conformers=conformers, ligand_id=ligand_id)
    summary = prepare_species(config, root, plan, [], conformers=conformers, encounter=True, ligand_id=ligand_id)
    summary['imported'] = imported
    summary['requeued_geometry_failures'] = requeued
    (root / 'encounters.json').write_text(json.dumps(summary, indent=2) + '\n')
    return {k: v for k, v in summary.items() if k not in ('plan', 'links')}


def requeue_filtered_encounters(root, *, ligand_id=None):
    """Retry old geometry rejections; reuse converged optimizations for NumFreq."""
    from hpc_precursors import connect_precursors
    from hpc_workflow import update_job
    database = root / 'encounters.sqlite'
    if not database.is_file():
        return 0
    with closing(connect_precursors(database, encounter=True)) as db:
        sql = "SELECT * FROM jobs WHERE status='failed' AND error LIKE 'Graph change or bond-length change%'"
        rows = [dict(r) for r in db.execute(sql + (' AND ligand_id=?' if ligand_id is not None else ''),
                                          (str(ligand_id),) if ligand_id is not None else ())]
        if not rows:
            return 0
        archive = root / 'encounter_filter_resets' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive.mkdir(parents=True)
        with closing(sqlite3.connect(archive / 'encounters.sqlite')) as snapshot:
            db.backup(snapshot)
        (archive / 'previous_rows.json').write_text(json.dumps(rows, indent=2) + '\n')
        retries = []
        for row in rows:
            restart = None
            if row['out_file']:
                output = root / row['out_file']
                optimized = output.with_suffix('.xyz')
                try:
                    text = output.read_text()
                    if ('ORCA TERMINATED NORMALLY' in text and
                            ('THE OPTIMIZATION HAS CONVERGED' in text.upper() or 'OPTIMIZATION CONVERGED' in text.upper())):
                        molecule_with_xyz(row['smiles'], optimized.read_text())
                        restart = str(optimized.relative_to(root))
                except (OSError, ValueError):
                    pass
            retries.append((row, restart))
        db.execute('BEGIN IMMEDIATE')
        try:
            changed = []
            for row, restart in retries:
                current = db.execute('SELECT status,error FROM jobs WHERE id=?', (row['id'],)).fetchone()
                if current['status'] != 'failed' or current['error'] != row['error']:
                    continue
                update_job(db, row['id'], status='pending', stage='freq' if restart else 'opt',
                           error=None, retained=None, worker_id=None, xyz_file=restart,
                           electronic_energy_eh=None, gibbs_free_energy_eh=None,
                           frequencies_json=None, imaginary_count=None, out_file=None)
                changed.append(dict(id=row['id'], restart_xyz=restart, stage='freq' if restart else 'opt'))
            db.commit()
        except BaseException:
            db.rollback()
            raise
        (archive / 'retry.json').write_text(json.dumps(changed, indent=2) + '\n')
    return len(changed)


def add_encounter_references(links, jobs):
    from hpc_singlepoint import lowest_jobs
    minima = {(r['ligand_id'], r['binding_permutation']): r for r in lowest_jobs(jobs)}
    requested = {(r['ligand_id'], r['binding_permutation']) for r in jobs}
    result = []
    for link in links:
        row = dict(link)
        key = (row['ligand_id'], row['binding_permutation'])
        best = minima.get(key)
        energy = best['gibbs_free_energy_eh'] if best else None
        acid = row.get('acetic_acid_gibbs_eh')
        row.update(encounter_requested=key in requested, encounter_job_id=best['id'] if best else None,
                   encounter_gibbs_eh=energy, encounter_with_acid_gibbs_eh=(energy + acid
                       if energy is not None and acid is not None else None))
        result.append(row)
    return result
