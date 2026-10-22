"""Fresh ligand–metal, substrate and shared acid/base solution reference energies."""

from __future__ import annotations

from contextlib import closing
import json
import os
from pathlib import Path
import shlex
import shutil
import sqlite3
import subprocess
import sys

from rdkit import Chem
import rxembed

from conformer_sampling import embedded_conformers
from hpc_workflow import claim_job, connect_database, run_job, update_job
from rxflow import ligand_directory, read_ligands
from rxflow_endpoints import mapped_molecule


def reference_root(config, base):
    return (base / config["paths"].get("references", "solution_references")).resolve()


def acid_base_plan():
    """One shared neutral acid and its closed-shell -1 conjugate base."""
    return [dict(ligand_id=name, binding_permutation=0, binding_mode=None,
                 coordination_isomer=None, seating=0, donor_maps="[]", multiplicity=1,
                 endpoint="substrate", protonation_site=0, charge=charge,
                 substrate=smiles, reference_name=name)
            for name, smiles, charge in (("acetic_acid", "CC(=O)O", 0),
                                         ("acetate", "CC(=O)[O-]", -1))]


def connect_precursors(path, *, bound=False, encounter=False):
    """Keep reference calculations separate from the reaction jobs and their role constraint."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=60, isolation_level=None)
    try:
        roles = "'encounter'" if encounter else ("'protonated_bound'" if bound else "'deprotonated','protonated','substrate'")
        site_column = "protonation_site INTEGER NOT NULL DEFAULT 0," if bound else ""
        site_key = ",protonation_site" if bound else ""
        connection.execute(f"""CREATE TABLE IF NOT EXISTS jobs (
            id INTEGER PRIMARY KEY, ligand_id TEXT NOT NULL, pair_index INTEGER NOT NULL,
            seating INTEGER NOT NULL, endpoint TEXT NOT NULL
                CHECK(endpoint IN ({roles})),
            conformer INTEGER NOT NULL, status TEXT NOT NULL DEFAULT 'pending',
            stage TEXT NOT NULL DEFAULT 'opt', attempt INTEGER NOT NULL DEFAULT 0,
            input_xyz TEXT NOT NULL, xyz_file TEXT, smiles TEXT NOT NULL,
            charge INTEGER NOT NULL, multiplicity INTEGER NOT NULL DEFAULT 1,
            donor INTEGER, hydrogen INTEGER, acceptor INTEGER, settings_json TEXT NOT NULL,
            electronic_energy_eh REAL, gibbs_free_energy_eh REAL, frequencies_json TEXT,
            imaginary_count INTEGER, retained INTEGER, out_file TEXT, worker_id TEXT, error TEXT,
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP, {site_column}
            UNIQUE(ligand_id,pair_index,seating,endpoint,conformer{site_key}))""")
        connection.execute("CREATE TABLE IF NOT EXISTS precursor_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        connection.execute("""CREATE TABLE IF NOT EXISTS precursor_links (
            ligand_id TEXT NOT NULL, pair_index INTEGER NOT NULL, seating INTEGER NOT NULL,
            binding_permutation INTEGER NOT NULL, substrate_site INTEGER, resonance INTEGER,
            protonation_site INTEGER NOT NULL, substrate_id TEXT NOT NULL,
            PRIMARY KEY(ligand_id,pair_index,seating,protonation_site))""")
    finally:
        connection.close()
    connection = connect_database(path)
    if bound:
        connection.execute("""CREATE TABLE IF NOT EXISTS bound_reactant_links (
            ligand_id TEXT NOT NULL, pair_index INTEGER NOT NULL, seating INTEGER NOT NULL,
            binding_permutation INTEGER NOT NULL, substrate_site INTEGER NOT NULL, resonance INTEGER,
            protonation_site INTEGER NOT NULL, bound_pair_index INTEGER NOT NULL,
            PRIMARY KEY(ligand_id,pair_index,seating,protonation_site))""")
        return connection
    connection.execute("""CREATE VIEW IF NOT EXISTS precursor_minima AS SELECT * FROM (
        SELECT id,ligand_id,pair_index,seating,endpoint,binding_permutation,conformer,
               charge,gibbs_free_energy_eh,xyz_file,
               row_number() OVER (PARTITION BY ligand_id,pair_index,seating,endpoint
                                  ORDER BY gibbs_free_energy_eh,id) energy_rank
        FROM jobs WHERE status='succeeded' AND gibbs_free_energy_eh IS NOT NULL)
        WHERE energy_rank=1""")
    connection.execute("""CREATE VIEW IF NOT EXISTS precursor_diagram_references AS
        SELECT l.*, d.id AS deprotonated_job_id, p.id AS protonated_job_id, s.id AS substrate_job_id,
               d.gibbs_free_energy_eh AS deprotonated_complex_gibbs_eh,
               p.gibbs_free_energy_eh AS protonated_complex_gibbs_eh,
               s.gibbs_free_energy_eh AS substrate_gibbs_eh,
               d.gibbs_free_energy_eh + s.gibbs_free_energy_eh AS deprotonated_complex_plus_substrate_eh,
               p.gibbs_free_energy_eh + s.gibbs_free_energy_eh AS protonated_complex_plus_substrate_eh
        FROM precursor_links l
        LEFT JOIN precursor_minima d ON d.ligand_id=l.ligand_id AND d.pair_index=l.binding_permutation
            AND d.seating=0 AND d.endpoint='deprotonated'
        LEFT JOIN precursor_minima p ON p.ligand_id=l.ligand_id AND p.pair_index=l.binding_permutation
            AND p.seating=l.protonation_site AND p.endpoint='protonated'
        LEFT JOIN precursor_minima s ON s.ligand_id=l.substrate_id AND s.endpoint='substrate'""")
    return connection


def canonical_substrate(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError("Invalid free-substrate SMILES")
    if Chem.GetFormalCharge(molecule) != 0 or any(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()):
        raise ValueError("Free substrates must be neutral, fully protonated closed-shell molecules")
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule)


def precursor_plan(config, base, root):
    """Reuse binding identities only; build every new molecular graph from its CSV ligand."""
    rows = {row["ID"]: row for row in read_ligands(base / config["paths"]["ligands_csv"])}
    reaction_db = root / "workflow.sqlite"
    if not reaction_db.is_file():
        raise ValueError("Prepare the reaction workflow first; precursor binding IDs come from workflow.sqlite")
    with closing(sqlite3.connect(f"file:{reaction_db}?mode=ro", uri=True)) as connection:
        connection.row_factory = sqlite3.Row
        bindings = [dict(row) for row in connection.execute("""SELECT DISTINCT ligand_id,
            binding_permutation,binding_mode,coordination_isomer,seating,donor_maps,multiplicity
            FROM jobs ORDER BY ligand_id,binding_permutation""")]
    plan = []
    seen = set()
    for binding in bindings:
        ligand_id = binding["ligand_id"]
        if ligand_id not in rows:
            raise ValueError(f"Ligand {ligand_id} is absent from the configured CSV")
        row = rows[ligand_id]
        if any(binding[key] is None for key in ("binding_permutation", "coordination_isomer", "donor_maps")):
            raise ValueError("Missing binding metadata; backfill the reaction database with prepare first")
        identity = (ligand_id, binding["binding_permutation"])
        if identity in seen:
            raise ValueError(f"Conflicting metadata for binding {identity}")
        seen.add(identity)
        ligand = Chem.MolFromSmiles(row["ligand_smiles"])
        if ligand is None or not ligand.GetNumAtoms():
            raise ValueError(f"Ligand {ligand_id}: invalid ligand SMILES")
        ligand = Chem.AddHs(ligand)
        if Chem.GetFormalCharge(ligand) != -1:
            raise ValueError(f"Ligand {ligand_id} must have charge -1 in its deprotonated input state")
        metal = Chem.MolFromSmiles(row["metal"])
        if metal.GetNumAtoms() != 1 or Chem.GetFormalCharge(metal) != 2:
            raise ValueError(f"Ligand {ligand_id}: precursor complexes require a +2 central metal")
        sites = [atom.GetIdx() + 1 for atom in ligand.GetAtoms() if atom.GetFormalCharge() < 0]
        if not sites:
            raise ValueError(f"Ligand {ligand_id} has no explicitly deprotonated atom")
        common = dict(binding, ligand_smiles=row["ligand_smiles"], metal=row["metal"],
                      metal_geom=row["metal_geom"])
        plan.append(dict(common, endpoint="deprotonated", protonation_site=0, charge=1))
        for site in sites:
            plan.append(dict(common, endpoint="protonated", protonation_site=site, charge=2))
    substrates = {}
    for row in rows.values():
        substrates.setdefault(canonical_substrate(row["substrate"]), []).append(row["ID"])
    for index, (smiles, ligand_ids) in enumerate(sorted(substrates.items()), 1):
        plan.append(dict(ligand_id=f"substrate_{index:03d}", binding_permutation=0,
                         binding_mode=None, coordination_isomer=None, seating=0, donor_maps="[]",
                         multiplicity=1, endpoint="substrate", protonation_site=0, charge=0,
                         substrate=smiles, ligand_ids=ligand_ids))
    if not bindings:
        raise ValueError("The reaction database contains no prepared binding arrangements")
    return plan


def precursor_links(plan, reaction_db):
    """Keep every reaction pair/resonance, while sharing identical reference species."""
    substrates = {spec["substrate"]: spec["ligand_id"] for spec in plan if spec["endpoint"] == "substrate"}
    sites = {}
    for spec in plan:
        if spec["endpoint"] == "protonated":
            sites.setdefault((spec["ligand_id"], spec["binding_permutation"]), []).append(spec["protonation_site"])
    with closing(sqlite3.connect(f"file:{reaction_db}?mode=ro", uri=True)) as connection:
        connection.row_factory = sqlite3.Row
        groups = connection.execute("""SELECT DISTINCT ligand_id,pair_index,seating,
            binding_permutation,substrate_site,resonance,substrate FROM jobs""").fetchall()
    links = []
    for row in groups:
        substrate_id = substrates.get(canonical_substrate(row["substrate"]))
        if substrate_id is None:
            raise ValueError(f"Reaction substrate for ligand {row['ligand_id']} is absent from the configured CSV")
        for site in sites[(row["ligand_id"], row["binding_permutation"])]:
            links.append(dict(ligand_id=row["ligand_id"], pair_index=row["pair_index"], seating=row["seating"],
                              binding_permutation=row["binding_permutation"], substrate_site=row["substrate_site"],
                              resonance=row["resonance"], protonation_site=site, substrate_id=substrate_id))
    return sorted(links, key=lambda row: (row["ligand_id"], row["pair_index"], row["seating"], row["protonation_site"]))


def bound_reactant_plan(config, base, root):
    """Every recorded bound reactant, protonated at each native ligand anion site."""
    protonated = {}
    for spec in precursor_plan(config, base, root):
        if spec["endpoint"] == "protonated":
            protonated.setdefault((spec["ligand_id"], spec["binding_permutation"]), []).append(spec)
    with closing(sqlite3.connect(f"file:{root / 'workflow.sqlite'}?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        reactants = db.execute("""SELECT DISTINCT ligand_id,pair_index,seating,binding_permutation,
            substrate_site,resonance,substrate,smiles FROM jobs WHERE endpoint='reactant'
            ORDER BY ligand_id,pair_index,seating""").fetchall()
    plan, seen = [], {}
    for row in reactants:
        for spec in protonated[(row["ligand_id"], row["binding_permutation"])]:
            identity = (row["ligand_id"], row["binding_permutation"], row["seating"],
                        row["substrate_site"], spec["protonation_site"])
            if identity in seen:
                shared = seen[identity]
            else:
                shared = dict(spec, **{key: row[key] for key in (
                    "pair_index", "seating", "substrate_site", "substrate")},
                    resonance=None, endpoint="protonated_bound", reactant_smiles=row["smiles"], reaction_pairs=[])
                seen[identity] = shared
                plan.append(shared)
            pathway = dict(pair_index=row["pair_index"], resonance=row["resonance"])
            if pathway in shared["reaction_pairs"]:
                raise ValueError(f"Conflicting bound-reactant metadata for {identity}")
            shared["reaction_pairs"].append(pathway)
    if not plan:
        raise ValueError("No bound reactant graphs are available in workflow.sqlite")
    return plan


def protonate_native(molecule, site):
    molecule = Chem.Mol(molecule)
    target = molecule.GetAtomWithIdx(site - 1)
    if target.GetFormalCharge() >= 0:
        raise ValueError("Protonation must target an explicitly deprotonated ligand atom")
    target.SetFormalCharge(target.GetFormalCharge() + 1)
    target.SetNumExplicitHs(target.GetNumExplicitHs() + 1)
    target.SetNoImplicit(True)
    Chem.SanitizeMol(molecule)
    molecule = Chem.AddHs(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return molecule


def protonated_bound_graph(spec, native_complex):
    """Use the native protonated ligand and the recorded substrate and donor slots.

    Some stored reactants use a different ligand resonance/tautomer form. Native
    protonation sites therefore refer to the CSV graph, not its shifted anion.
    """
    original = mapped_molecule(spec["reactant_smiles"])
    original.RemoveAllConformers()
    n_core = native_complex.GetNumAtoms()
    if Chem.GetFormalCharge(original) != 1:
        raise ValueError("The original bound reactant must have charge +1")
    if [a.GetAtomicNum() for a in list(original.GetAtoms())[:n_core]] != [
            a.GetAtomicNum() for a in native_complex.GetAtoms()]:
        raise ValueError("Native ligand/metal atom order differs from the bound reactant")
    substrate_donors = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in original.GetBonds()
                        if b.GetBondType() == Chem.BondType.DATIVE and b.GetBeginAtomIdx() >= n_core]
    if len(substrate_donors) != 1:
        raise ValueError("Expected one recorded substrate-to-metal bond")
    substrate = Chem.RWMol(original)
    for index in reversed(range(n_core)):
        substrate.RemoveAtom(index)
    substrate = substrate.GetMol()
    Chem.SanitizeMol(substrate)
    if canonical_substrate(Chem.MolToSmiles(substrate)) != canonical_substrate(spec["substrate"]):
        raise ValueError("The bound substrate must remain fully protonated and match the configured substrate")
    protonated = protonate_native(native_complex, spec["protonation_site"])
    combined = Chem.CombineMols(protonated, substrate)
    # Retain the original bound atom-map order and append only the new ligand H.
    order = list(range(n_core)) + list(range(n_core + 1, combined.GetNumAtoms())) + [n_core]
    editable = Chem.RWMol(Chem.RenumberAtoms(combined, order))
    for donor, metal in substrate_donors:
        editable.AddBond(donor, metal, Chem.BondType.DATIVE)
    molecule = editable.GetMol()
    Chem.SanitizeMol(molecule)
    sites = {}
    for bond in original.GetBonds():
        if bond.GetBondType() != Chem.BondType.DATIVE:
            continue
        donor = bond.GetBeginAtom()
        note = donor.GetProp("atomNote") if donor.HasProp("atomNote") else ""
        if not note.startswith("s") or not note[1:].isdigit():
            raise ValueError("The recorded bound reactant lacks donor-slot annotations")
        sites.setdefault(int(note[1:]), donor.GetIdx())
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    if Chem.GetFormalCharge(molecule) != spec["charge"] or spec["charge"] != 2:
        raise ValueError("Protonated bound reactants must have charge +2")
    return molecule, rxembed.Isomer(molecule, spec["metal_geom"], sites)


def precursor_graph(spec):
    """Return a real-metal mapped graph and fresh-embedding source, with no input coordinates."""
    if spec["endpoint"] == "substrate":
        molecule = Chem.AddHs(Chem.MolFromSmiles(spec["substrate"]))
        if Chem.GetFormalCharge(molecule) != spec["charge"]:
            raise ValueError("Free reference charge mismatch")
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        return molecule, molecule
    if spec["endpoint"] == "encounter":
        from hpc_encounters import encounter_graph
        return encounter_graph(spec)
    ligand = Chem.AddHs(Chem.MolFromSmiles(spec["ligand_smiles"]))
    metal = Chem.MolFromSmiles(spec["metal"])
    editable = Chem.RWMol(Chem.CombineMols(ligand, metal))
    donor_maps = json.loads(spec["donor_maps"])
    for donor_map in donor_maps:
        if not 1 <= donor_map <= ligand.GetNumAtoms():
            raise ValueError("Donor map is outside the ligand")
        editable.AddBond(donor_map - 1, ligand.GetNumAtoms(), Chem.BondType.DATIVE)
    deprotonated = editable.GetMol()
    Chem.SanitizeMol(deprotonated)
    for atom in deprotonated.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    if spec["endpoint"] == "protonated_bound":
        return protonated_bound_graph(spec, deprotonated)
    isomers = rxembed.metal(deprotonated, geometry=spec["metal_geom"])
    index = spec["coordination_isomer"] - 1
    if not 0 <= index < len(isomers):
        raise ValueError("Recorded coordination isomer is unavailable for this ligand graph")
    isomer = isomers[index]
    if spec["endpoint"] == "deprotonated":
        molecule, source = deprotonated, isomer
    else:
        molecule = protonate_native(deprotonated, spec["protonation_site"])
        sites = {slot: donor if isinstance(donor, int) else donor.atoms[0]
                 for slot, donor in enumerate(isomer.centres[0].vertices) if donor is not None}
        source = rxembed.Isomer(molecule, spec["metal_geom"], sites)
    if Chem.GetFormalCharge(molecule) != spec["charge"]:
        raise ValueError(f"Precursor charge mismatch: expected {spec['charge']}")
    if molecule.GetNumConformers():
        raise ValueError("Precursors must start from a graph without inherited coordinates")
    return molecule, source


def sample_precursor(spec, *, n=10, seed=12345, prune_rmsd=0.5):
    """Freshly embed n total conformers; retain each successful stereo ensemble fairly."""
    if spec["endpoint"] == "encounter":
        from hpc_encounters import sample_encounter
        return sample_encounter(spec, n=n, seed=seed, prune_rmsd=prune_rmsd)
    molecule, source = precursor_graph(spec)
    if spec.get("reference_name"):
        # Request the full seed budget; the explicit prism step removes duplicates.
        embedded = rxembed.embed(source, n=n, params=rxembed.EmbedParams(seed=seed, threads=1, prune_rms=-1))
    else:
        embedded = rxembed.embed(source, n=n, seed=seed, threads=1)
    batches = embedded_conformers(embedded, prune_rmsd)
    # Round-robin across undefined ligand stereo alternatives, without discarding one hand.
    conformers = []
    for rank in range(max(map(len, batches), default=0)):
        for batch in batches:
            if rank < len(batch):
                geometry = batch[rank]
                if Chem.GetFormalCharge(geometry) != spec["charge"]:
                    raise ValueError("Embedding returned a precursor with the wrong total charge")
                if [a.GetAtomicNum() for a in geometry.GetAtoms()] != [a.GetAtomicNum() for a in molecule.GetAtoms()]:
                    raise ValueError("Embedding changed the precursor atom inventory/order")
                if [a.GetAtomMapNum() for a in geometry.GetAtoms()] != [a.GetAtomMapNum() for a in molecule.GetAtoms()]:
                    raise ValueError("Embedding changed the precursor atom-map order")
                conformers.append(geometry)
                if len(conformers) == n:
                    return molecule, conformers
    return molecule, conformers


def precursor_location(spec, conformer):
    if spec["endpoint"] == "encounter":
        return ligand_directory(Path(), spec["ligand_id"]) / "complex" / "encounter" / (
            f"binding_{spec['binding_permutation']:03d}_conformer_{conformer:03d}_input.xyz")
    if spec["endpoint"] == "protonated_bound":
        return ligand_directory(Path(), spec["ligand_id"]) / "complex" / "protonated_bound" / (
            f"pair_{spec['pair_index']:03d}_seating_{spec['seating']:02d}_"
            f"site_{spec['protonation_site']:03d}_conformer_{conformer:03d}_input.xyz")
    if spec.get("reference_name"):
        name = spec["reference_name"]
        return Path(name) / f"{name}_conformer_{conformer:03d}_input.xyz"
    if spec["endpoint"] == "substrate":
        return Path("substrate") / f"{spec['ligand_id']}_neutral_conformer_{conformer:03d}_input.xyz"
    state = spec["endpoint"]
    if state == "protonated":
        state += f"_site_{spec['protonation_site']:03d}"
    return ligand_directory(Path(), spec["ligand_id"]) / "complex" / (
        f"binding_{spec['binding_permutation']:03d}_{state}_conformer_{conformer:03d}_input.xyz")


def prepare_precursors(config, base, root, *, conformers=10, reference=False, ligand_id=None):
    if conformers < 1:
        raise ValueError("Precursor conformer count must be positive")
    plan = acid_base_plan() if reference else precursor_plan(config, base, root)
    links = [] if reference else precursor_links(plan, root / "workflow.sqlite")
    summary = prepare_species(config, root, plan, links, conformers=conformers, reference=reference)
    if not reference:
        from hpc_encounters import prepare_encounters
        encounter_conformers = config.get("sampling", {}).get("encounter_conformers", conformers)
        summary["encounters"] = prepare_encounters(config, base, root, conformers=encounter_conformers,
                                                    ligand_id=ligand_id)
    return summary


def prepare_bound_reactants(config, base, root, *, conformers=10):
    if conformers < 1:
        raise ValueError("Precursor conformer count must be positive")
    plan = bound_reactant_plan(config, base, root)
    links = [dict(ligand_id=spec["ligand_id"], pair_index=pathway["pair_index"], seating=spec["seating"],
                  binding_permutation=spec["binding_permutation"], substrate_site=spec["substrate_site"],
                  resonance=pathway["resonance"], protonation_site=spec["protonation_site"],
                  bound_pair_index=spec["pair_index"])
             for spec in plan for pathway in spec["reaction_pairs"]]
    return prepare_species(config, root, plan, links,
                           conformers=conformers, bound=True)


def species_protocol(config, plan, links, conformers):
    settings = {key: value for key, value in config["calculation"].items() if key != "max_iterations"}
    settings.setdefault("opt_attempts", settings.get("ts_attempts", 2))
    seed = config.get("sampling", {}).get("seed", 12345)
    sampling = config.get("sampling", {})
    protocol = json.dumps(dict(version=1, conformers=conformers, seed=seed, settings=settings, plan=plan,
                               links=links, **({"prune_rmsd": sampling["prune_rmsd"]}
                                              if sampling.get("prune_rmsd", 0.5) != 0.5 else {})), sort_keys=True)
    return settings, seed, protocol


def encounter_protocol_extension(frozen_protocol, requested_protocol):
    """Allow append-only encounter-plan growth and upward conformer extensions."""
    if frozen_protocol == requested_protocol:
        return False
    return _species_protocol_extension(frozen_protocol, requested_protocol, allow_count_increase=True)


def precursor_protocol_extension(frozen_protocol, requested_protocol):
    """Allow adding new precursor species/links without changing existing chemistry."""
    return _species_protocol_extension(frozen_protocol, requested_protocol, allow_count_increase=False)


def _species_protocol_extension(frozen_protocol, requested_protocol, *, allow_count_increase):
    old, new = json.loads(frozen_protocol), json.loads(requested_protocol)
    old_count, new_count = old.pop("conformers", None), new.pop("conformers", None)
    old_plan, new_plan = old.pop("plan", None), new.pop("plan", None)
    old_links, new_links = old.pop("links", None), new.pop("links", None)
    if old != new or not isinstance(old_count, int) or not isinstance(new_count, int):
        return False
    if new_count < old_count or (new_count != old_count and not allow_count_increase):
        return False
    if not isinstance(old_plan, list) or not isinstance(new_plan, list):
        return False
    if not isinstance(old_links, list) or not isinstance(new_links, list):
        return False
    encode = lambda values: {json.dumps(value, sort_keys=True, separators=(",", ":")) for value in values}
    old_plan_set, new_plan_set = encode(old_plan), encode(new_plan)
    old_links_set, new_links_set = encode(old_links), encode(new_links)
    if not old_plan_set <= new_plan_set or not old_links_set <= new_links_set:
        return False
    return (new_count > old_count or old_plan_set != new_plan_set or old_links_set != new_links_set)


def prepare_species(config, root, plan, links, *, conformers=10, reference=False, bound=False,
                    encounter=False, ligand_id=None):
    settings, seed, protocol = species_protocol(config, plan, links, conformers)
    database_name = "encounters.sqlite" if encounter else ("acid_base.sqlite" if bound else ("references.sqlite" if reference else "precursors.sqlite"))
    connection = connect_precursors(root / database_name, bound=bound, encounter=encounter)
    added = repaired = 0
    errors = []
    sampling_records = []
    try:
        frozen = connection.execute("SELECT value FROM precursor_metadata WHERE key='protocol'").fetchone()
        extending_encounters = False
        if frozen and frozen[0] != protocol:
            extending_encounters = encounter and encounter_protocol_extension(frozen[0], protocol)
            extending_precursors = (not encounter and not reference and not bound
                                    and precursor_protocol_extension(frozen[0], protocol))
            if not (extending_encounters or extending_precursors):
                raise ValueError("Precursor chemistry/settings changed; use a fresh results directory")
        if not frozen:
            connection.execute("INSERT INTO precursor_metadata VALUES ('protocol',?)", (protocol,))
        for link in links:
            table = "bound_reactant_links" if bound else "precursor_links"
            connection.execute(f"INSERT OR IGNORE INTO {table} VALUES (?,?,?,?,?,?,?,?)", tuple(link.values()))
        for spec in plan:
            if ligand_id is not None and str(spec["ligand_id"]) != str(ligand_id):
                continue
            identity = (spec["ligand_id"], spec["pair_index"] if bound else spec["binding_permutation"],
                        spec["seating"] if bound else spec["protonation_site"], spec["endpoint"])
            site_filter = " AND protonation_site=?" if bound else ""
            query_identity = (*identity, spec["protonation_site"]) if bound else identity
            previous = {r["conformer"]: dict(r) for r in connection.execute(
                "SELECT * FROM jobs WHERE ligand_id=? AND pair_index=? AND seating=? AND endpoint=?" + site_filter,
                query_identity)}
            needed = [i for i in range(1, conformers + 1)
                      if i not in previous or (previous[i]["status"] == "failed" and previous[i]["stage"] == "sampling")]
            if not needed:
                continue
            molecule, _ = precursor_graph(spec)
            geometries, error = [], None
            try:
                pruning = config.get("sampling", {}).get("prune_rmsd", 0.5)
                batch_seed = seed + (max(previous, default=0) * 1_000_003 if extending_encounters else 0)
                molecule, geometries = sample_precursor(spec, n=len(needed), seed=batch_seed,
                                                        prune_rmsd=pruning)
                if not geometries:
                    raise ValueError("No precursor conformers survived embedding and pruning")
            except (ValueError, RuntimeError, TimeoutError) as exc:
                error = str(exc)
            sampling_records.append(dict(ligand_id=spec["ligand_id"], endpoint=spec["endpoint"],
                pair_index=spec["pair_index"] if bound else spec["binding_permutation"],
                seating=spec["seating"] if bound else spec["protonation_site"],
                requested=len(needed), retained=len(geometries), error=error))
            smiles = Chem.MolToSmiles(molecule)
            for batch_index, conformer in enumerate(needed):
                relative = precursor_location(spec, conformer)
                target = root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                valid = batch_index < len(geometries)
                detail = None if valid or error is None else error
                if valid:
                    text = Chem.MolToXYZBlock(geometries[batch_index])
                    if target.exists() and target.read_text() != text:
                        raise ValueError(f"Existing precursor seed differs: {target}")
                    if not target.exists():
                        target.write_text(text)
                values = dict(status="pending" if valid else ("failed" if error else "skipped"),
                              stage="opt" if valid else ("sampling" if error else "pruned"),
                              error=f"Precursor embedding failed: {detail}" if detail else None)
                if conformer in previous:
                    update_job(connection, previous[conformer]["id"], **values)
                    repaired += int(valid)
                else:
                    extra_columns = ",protonation_site,substrate_site,resonance" if bound else ""
                    extra_parameters = ",?,?,?" if bound else ""
                    extra_values = (spec["protonation_site"], spec["substrate_site"], spec["resonance"]) if bound else ()
                    connection.execute(f"""INSERT INTO jobs
                        (ligand_id,pair_index,seating,endpoint,conformer,input_xyz,smiles,charge,multiplicity,
                         settings_json,binding_permutation,binding_mode,coordination_isomer,donor_maps,substrate,
                         status,stage,error{extra_columns}) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?{extra_parameters})""",
                        (*identity, conformer, str(relative), smiles, spec["charge"], spec["multiplicity"],
                         json.dumps(settings, sort_keys=True), spec["binding_permutation"], spec["binding_mode"],
                         spec["coordination_isomer"], spec["donor_maps"], spec.get("substrate"),
                         values["status"], values["stage"], values["error"], *extra_values))
                    added += 1
                if detail:
                    errors.append(dict(ligand_id=spec["ligand_id"], binding_permutation=spec["binding_permutation"],
                                       state=spec["endpoint"], conformer=conformer, error=detail))
        if frozen and frozen[0] != protocol:
            connection.execute("UPDATE precursor_metadata SET value=? WHERE key='protocol'", (protocol,))
        summary = dict(species=len(plan), conformers_per_species=conformers, added=added, repaired=repaired,
                       errors=errors, sampling=sampling_records, plan=plan, links=links)
        summary_name = "encounters.json" if encounter else ("acid_base.json" if bound else ("references.json" if reference else "precursors.json"))
        (root / summary_name).write_text(json.dumps(summary, indent=2) + "\n")
        return summary
    finally:
        connection.close()


def bound_reactant_minima(rows):
    best = {}
    for row in rows:
        if row["status"] != "succeeded" or row["gibbs_free_energy_eh"] is None:
            continue
        key = tuple(row[name] for name in ("ligand_id", "binding_permutation", "seating", "substrate_site", "protonation_site"))
        if key not in best or (row["gibbs_free_energy_eh"], row["id"]) < (best[key]["gibbs_free_energy_eh"], best[key]["id"]):
            best[key] = row
    return best


def bound_reactant_status(database):
    if not database.is_file():
        return dict(database=str(database), available=False, counts=[], lowest_gibbs=[])
    with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        rows = [dict(row) for row in db.execute("SELECT * FROM jobs")]
        counts = [dict(row) for row in db.execute(
            "SELECT endpoint,status,stage,count(*) jobs FROM jobs GROUP BY endpoint,status,stage")]
        links = [dict(row) for row in db.execute("SELECT * FROM bound_reactant_links")]
    return dict(database=str(database), available=True, counts=counts,
                lowest_gibbs=list(bound_reactant_minima(rows).values()), pathway_links=links)


def add_bound_reactant_references(links, jobs):
    minima = bound_reactant_minima(jobs)
    requested = {tuple(row[name] for name in ("ligand_id", "binding_permutation", "seating", "substrate_site", "protonation_site"))
                 for row in jobs}
    result = []
    for original in links:
        row = dict(original)
        key = tuple(row[name] for name in ("ligand_id", "binding_permutation", "seating", "substrate_site", "protonation_site"))
        bound = minima.get(key)
        gibbs = bound["gibbs_free_energy_eh"] if bound else None
        acetate = row["acetate_gibbs_eh"]
        row.update(protonated_bound_requested=key in requested,
                   protonated_bound_job_id=bound["id"] if bound else None,
                   protonated_bound_gibbs_eh=gibbs,
                   protonated_bound_with_acetate_gibbs_eh=(gibbs + acetate
                       if gibbs is not None and acetate is not None else None))
        result.append(row)
    return result


def balanced_diagram_references(connection, reference_database, reaction_database):
    """Join global acid/base minima and bound reaction energies without changing either DB."""
    references = {"acetic_acid": None, "acetate": None}
    if reference_database.is_file():
        with closing(sqlite3.connect(f"file:{reference_database}?mode=ro", uri=True)) as db:
            reference_protocol = db.execute("SELECT value FROM precursor_metadata WHERE key='protocol'").fetchone()
            precursor_protocol = connection.execute("SELECT value FROM precursor_metadata WHERE key='protocol'").fetchone()
            if reference_protocol and precursor_protocol:
                if json.loads(reference_protocol[0])["settings"] != json.loads(precursor_protocol[0])["settings"]:
                    raise ValueError("Acid/base and precursor calculation settings differ; use matching references")
            references.update(db.execute("""SELECT ligand_id,min(gibbs_free_energy_eh)
                FROM jobs WHERE status='succeeded' AND gibbs_free_energy_eh IS NOT NULL
                GROUP BY ligand_id""").fetchall())
    bound = {}
    if reaction_database.is_file():
        with closing(sqlite3.connect(f"file:{reaction_database}?mode=ro", uri=True)) as db:
            for ligand, pair, seating, role, energy in db.execute("""SELECT ligand_id,pair_index,seating,
                endpoint,min(gibbs_free_energy_eh) FROM jobs WHERE status='succeeded'
                AND gibbs_free_energy_eh IS NOT NULL GROUP BY ligand_id,pair_index,seating,endpoint"""):
                bound[(ligand, pair, seating, role)] = energy
    def total(*values):
        return None if any(value is None for value in values) else sum(values)
    result = []
    for row in connection.execute("SELECT * FROM precursor_diagram_references"):
        row = dict(row)
        row.update(acetic_acid_gibbs_eh=references["acetic_acid"], acetate_gibbs_eh=references["acetate"])
        row["protonated_start_gibbs_eh"] = total(row["protonated_complex_plus_substrate_eh"], references["acetate"])
        row["deprotonated_start_gibbs_eh"] = total(row["deprotonated_complex_plus_substrate_eh"], references["acetic_acid"])
        identity = (row["ligand_id"], row["pair_index"], row["seating"])
        for role in ("reactant", "TS", "product"):
            row[f"{role}_with_acid_gibbs_eh"] = total(bound.get((*identity, role)), references["acetic_acid"])
        result.append(row)
    bound_database = reaction_database.with_name("acid_base.sqlite")
    bound_jobs = []
    if bound_database.is_file():
        with closing(sqlite3.connect(f"file:{bound_database}?mode=ro", uri=True)) as db:
            db.row_factory = sqlite3.Row
            bound_jobs = [dict(row) for row in db.execute("SELECT * FROM jobs")]
    from hpc_encounters import add_encounter_references
    from hpc_singlepoint import read_jobs
    return add_encounter_references(add_bound_reactant_references(result, bound_jobs),
                                    read_jobs(reaction_database.with_name("encounters.sqlite")))


def has_plotted_gibbs(row):
    """Admit only a validated result from the notebook's explicitly loaded method."""
    import math
    energy = row.get("gibbs_free_energy_eh")
    if row.get("status") != "succeeded" or energy is None or not math.isfinite(energy):
        return False
    if row.get("energy_method") == "uma":
        return row.get("method_status") == "succeeded"
    return row.get("singlepoint_status") == "succeeded" and row.get("singlepoint_job_id") is not None


def load_diagram_inputs(reaction_database, precursor_database=None, reference_database=None, *, energy_mode="auto", uma_root=None):
    """Read notebook inputs without creating, upgrading or writing any calculation database."""
    reaction_database = Path(reaction_database).resolve()
    if energy_mode.lower() == "uma":
        from hpc_uma import load_uma_diagram_inputs
        return load_uma_diagram_inputs(reaction_database, uma_root)
    precursor_database = (Path(precursor_database).resolve() if precursor_database is not None
                          else reaction_database.with_name("precursors.sqlite"))
    reference_database = (Path(reference_database).resolve() if reference_database is not None
                          else reaction_database.parent.parent / "solution_references" / "references.sqlite")
    if not reaction_database.is_file():
        raise FileNotFoundError(f"Database not found: {reaction_database}")

    def read_jobs(path):
        if not path.is_file():
            return []
        with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as db:
            db.row_factory = sqlite3.Row
            return [dict(row) for row in db.execute("SELECT * FROM jobs")]

    jobs = read_jobs(reaction_database)
    precursor_jobs = read_jobs(precursor_database)
    reference_jobs = read_jobs(reference_database)
    bound_reactant_jobs = read_jobs(reaction_database.with_name("acid_base.sqlite"))
    encounter_jobs = read_jobs(reaction_database.with_name("encounters.sqlite"))
    links = []
    if precursor_database.is_file():
        with closing(sqlite3.connect(f"file:{precursor_database}?mode=ro", uri=True)) as db:
            db.row_factory = sqlite3.Row
            links = balanced_diagram_references(db, reference_database, reaction_database)
    minima = {}
    for row in reference_jobs:
        if row["status"] == "succeeded" and row["gibbs_free_energy_eh"] is not None:
            previous = minima.get(row["ligand_id"])
            if previous is None or (row["gibbs_free_energy_eh"], row["id"]) < (previous["gibbs_free_energy_eh"], previous["id"]):
                minima[row["ligand_id"]] = row
    data = dict(jobs=jobs, precursor_jobs=precursor_jobs, reference_jobs=reference_jobs,
                bound_reactant_jobs=bound_reactant_jobs, encounter_jobs=encounter_jobs,
                diagram_references=links, reference_minima=minima)
    from hpc_singlepoint import apply_singlepoint_energies
    return apply_singlepoint_energies(data, reaction_database, precursor_database, reference_database, energy_mode)


def precursor_status(connection, *, reference_database=None, reaction_database=None):
    counts = [dict(row) for row in connection.execute(
        "SELECT endpoint,status,stage,count(*) jobs FROM jobs GROUP BY endpoint,status,stage")]
    best = [dict(row) for row in connection.execute("""SELECT id,ligand_id,binding_permutation,
        seating AS protonation_site,endpoint,conformer,charge,gibbs_free_energy_eh,xyz_file
        FROM precursor_minima ORDER BY ligand_id,binding_permutation,endpoint,protonation_site""")]
    links = (balanced_diagram_references(connection, reference_database, reaction_database)
             if reference_database is not None and reaction_database is not None else
             [dict(row) for row in connection.execute("SELECT * FROM precursor_diagram_references")])
    return dict(counts=counts, lowest_gibbs=best, diagram_references=links)


def submit_precursors(config_path, config, *, workers=100, conformers=10, account=None, partition=None, resume=False, reference=False):
    if workers < 1 or conformers < 1:
        raise ValueError("Worker and conformer counts must be positive")
    base = config_path.resolve().parent
    root = reference_root(config, base) if reference else base / config["paths"]["results"]
    database = root / ("references.sqlite" if reference else "precursors.sqlite")
    workflow = "acid-base" if reference else "precursors"
    if resume and not database.is_file():
        raise ValueError(f"No database exists; use {workflow} submit or prepare first")
    if not resume and database.exists() and not reference:
        raise ValueError(f"Database already exists; use {workflow} prepare to finish sampling, or resume")
    plan = acid_base_plan() if reference else precursor_plan(config, base, root)
    if reference:
        plan += bound_reactant_plan(config, base, base / config["paths"]["results"])
    else:
        from hpc_encounters import encounter_plan
        plan += encounter_plan(config, base, root)
    command = shlex.join([sys.executable, "-m", "hpc_workflow", "--config", str(config_path.resolve()), workflow])
    wrap = "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; exec " + command
    logs = base / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    common = ["sbatch", "--parsable", f"--chdir={base}", "--export=ALL",
              f"--mem={config['calculation']['memory_mb']}", "--time=24:00:00"]
    for flag, value in (("account", account), ("partition", partition)):
        if value:
            common.append(f"--{flag}={value}")

    def submit(stage, *, cpus=1, dependency=None, array=None):
        args = [*common, f"--job-name=rxflow-{workflow}-{stage}", f"--cpus-per-task={cpus}",
                f"--output={logs / (workflow + '-' + stage + '-%A_%a.out')}",
                f"--error={logs / (workflow + '-' + stage + '-%A_%a.err')}"]
        if dependency:
            args.append(f"--dependency=afterok:{dependency}")
        if array:
            args.append(f"--array={array}")
        suffix = f" --conformers {conformers}" if stage == "prepare" else ""
        return subprocess.check_output([*args, "--wrap", f"{wrap} {stage}{suffix}"], text=True).strip().split(";")[0]

    prepared = submit("prepare")
    workers_id = submit("worker", cpus=config["calculation"]["cpus"], dependency=prepared,
                        array=f"0-{workers - 1}%{workers}")
    return dict(prepare=prepared, workers=workers_id, species=len(plan), conformers_per_species=conformers)


def handle_precursors(args, config, base, root, parser):
    reference = args.command == "acid-base"
    workflow = args.command
    run_root = root
    if reference:
        root = reference_root(config, base)
    action = args.precursor_command
    if getattr(args, "bound", False) and (not reference or action != "retry"):
        parser.error("--bound is supported only for acid-base retry")
    if action is None:
        parser.error(f"{workflow} requires submit, prepare, worker, status, retry, or resume")
    if action == "plan":
        parser.error("plan is supported by singlepoint; use status for this workflow")
    try:
        if action == "prepare":
            summary = prepare_precursors(config, base, root, conformers=args.conformers, reference=reference,
                                         ligand_id=getattr(args, "ligand_id", None))
            if reference:
                bound_summary = prepare_bound_reactants(config, base, run_root, conformers=args.conformers)
                summary = dict(references=summary, bound_reactants=bound_summary)
                summary = {key: {k: v for k, v in value.items() if k not in ("plan", "links")}
                           for key, value in summary.items()}
            print(json.dumps({key: value for key, value in summary.items() if key not in ("plan", "links")}, indent=2))
            return
        executables = {}
        if action in ("submit", "resume", "worker"):
            for key in ("orca_exe", "external_exe", "xtb_exe"):
                value = os.path.expanduser(os.path.expandvars(config["executables"][key]))
                executable = shutil.which(value if not value.startswith(".") else str(base / value))
                if not executable:
                    raise ValueError(f"Executable {key} not found: {value!r}")
                executables[key] = str(Path(executable).resolve())
        if action in ("submit", "resume"):
            if not shutil.which("sbatch"):
                raise ValueError(f"{workflow} submit/resume requires Slurm's sbatch command")
            print(json.dumps(submit_precursors(args.config, config, workers=args.workers,
                  conformers=args.conformers, account=args.account, partition=args.partition,
                  resume=action == "resume", reference=reference), indent=2))
            return
        database = root / ("references.sqlite" if reference else "precursors.sqlite")
        if reference and getattr(args, "bound", False):
            database = run_root / "acid_base.sqlite"
        encounter = getattr(args, "encounters", False)
        if encounter:
            if reference or action != "retry":
                raise ValueError("--encounters is supported only for precursors retry")
            database = root / "encounters.sqlite"
        if not database.is_file():
            raise ValueError(f"No database exists; use {workflow} submit or prepare first")
        connection = connect_precursors(database, bound=reference and getattr(args, "bound", False), encounter=encounter)
        try:
            if action == "worker":
                queues = [(connection, root)]
                bound_connection = None
                if reference and (run_root / "acid_base.sqlite").is_file():
                    bound_connection = connect_precursors(run_root / "acid_base.sqlite", bound=True)
                    queues.append((bound_connection, run_root))
                if not reference and (root / "encounters.sqlite").is_file():
                    bound_connection = connect_precursors(root / "encounters.sqlite", encounter=True)
                    queues.append((bound_connection, root))
                try:
                    while True:
                        for queue, queue_root in queues:
                            job = claim_job(queue)
                            if job is not None:
                                print(f"Precursor {job['id']}: {job['ligand_id']} {job['endpoint']}", flush=True)
                                run_job(queue, job, queue_root, executables)
                                break
                        else:
                            break
                finally:
                    if bound_connection is not None:
                        bound_connection.close()
            elif action == "retry":
                if args.id is None:
                    raise ValueError(f"{workflow} retry requires --id; stop any previous worker first")
                row = connection.execute("SELECT stage FROM jobs WHERE id=?", (args.id,)).fetchone()
                if row and row["stage"] == "sampling":
                    raise ValueError(f"Retry embedding failures with {workflow} prepare")
                connection.execute("""UPDATE jobs SET status='pending',stage='opt',error=NULL,retained=NULL,
                    electronic_energy_eh=NULL,gibbs_free_energy_eh=NULL,frequencies_json=NULL,
                    imaginary_count=NULL,xyz_file=NULL,out_file=NULL,worker_id=NULL
                    WHERE id=? AND status IN ('failed','running')""", (args.id,))
            else:
                if reference:
                    summary = dict(references=precursor_status(connection),
                                   bound_reactants=bound_reactant_status(run_root / "acid_base.sqlite"))
                    print(json.dumps(summary, indent=2))
                    return
                paths = {} if reference else dict(reference_database=reference_root(config, base) / "references.sqlite",
                                                  reaction_database=root / "workflow.sqlite")
                summary = precursor_status(connection, **paths)
                from hpc_singlepoint import read_jobs
                encounter_rows = read_jobs(root / "encounters.sqlite")
                from collections import Counter
                summary["encounters"] = dict(database=str(root / "encounters.sqlite"),
                    counts=[dict(status=status, jobs=count) for status, count in Counter(
                        row["status"] for row in encounter_rows).items()])
                print(json.dumps(summary, indent=2))
        finally:
            connection.close()
    except (ValueError, subprocess.CalledProcessError) as exc:
        parser.error(str(exc))
