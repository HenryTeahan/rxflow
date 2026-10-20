"""Energy accounting and pathway preservation in the reaction-diagram notebook."""

from collections import Counter, defaultdict
from contextlib import closing
import hashlib
from html import escape
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import rdDepictor

from hpc_precursors import load_diagram_inputs
from rxflow import add_substrate


class ReactionDiagramTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.notebook = json.loads(Path(__file__).with_name('reaction_diagrams.ipynb').read_text())
        cls.keys = ('ligand_id', 'pair_index', 'seating', 'binding_permutation')
        namespace = dict(Counter=Counter, defaultdict=defaultdict, escape=escape, jobs=[], profiles=[], keys=cls.keys,
                         HARTREE_TO_KCAL=627.509474, BARRIER_LIMIT_KCAL=25.0,
                         Chem=Chem, rdDepictor=rdDepictor)
        site_key = next(cell for cell in cls.notebook['cells'] if cell.get('id') == 'substrate-site-key')
        exec(compile(''.join(site_key['source']), 'substrate site key', 'exec'), namespace)
        helper = next(cell for cell in cls.notebook['cells'] if cell.get('id') == 'permutation-coverage')
        exec(compile(''.join(helper['source']), 'notebook helpers', 'exec'), namespace)
        coverage = next(cell for cell in cls.notebook['cells'] if cell.get('id') == 'reactant-ts-coverage-helpers')
        exec(compile(''.join(coverage['source']), 'reactant TS coverage helpers', 'exec'), namespace)
        cls.namespace = namespace

    def setUp(self):
        self.profiles = [dict(ligand_id='11', pair_index=pair, seating=4, binding_permutation=1,
                              substrate_site=site, resonance=resonance, reactant=-100., TS=-99.98,
                              product=-100.01)
                         for pair, site, resonance in [(1, 6, 0), (2, 6, 1), (3, 7, 0), (4, 7, 1)]]
        self.links = [dict(**{key: profile[key] for key in self.keys},
                          substrate_site=profile['substrate_site'], resonance=profile['resonance'],
                          protonation_site=1, protonated_start_gibbs_eh=-130.,
                          deprotonated_start_gibbs_eh=-130., protonated_complex_gibbs_eh=-101.,
                          deprotonated_complex_gibbs_eh=-100., substrate_gibbs_eh=-20.,
                          acetate_gibbs_eh=-9., encounter_gibbs_eh=-100.1,
                          encounter_with_acid_gibbs_eh=-110.1) for profile in self.profiles]

    def diagram(self, profiles=None, links=None, acid=-10.):
        return self.namespace['energy_diagram']('11', self.profiles if profiles is None else profiles,
                                               self.links if links is None else links, acid)

    def test_encounter_state_balances_spectators_without_counting_substrate_twice(self):
        diagram = self.diagram()
        self.assertEqual(diagram['reference'], -130.)
        self.assertEqual(diagram['state_keys'], ['protonated_start', 'deprotonated_start', 'encounter', 'reactant', 'TS', 'product'])
        route = diagram['routes'][0]
        self.assertEqual(route['energies']['encounter'], -110.1)
        self.assertIn(('deprotonated_start', 'encounter'), diagram['edges'])
        self.assertIn(('encounter', 'reactant'), diagram['edges'])
        self.assertEqual(route['energies']['reactant'], -110.)
        self.assertEqual(route['energies']['TS'], -109.98)
        self.assertAlmostEqual(route['relative']['TS'] - route['relative']['reactant'], .02 * 627.509474)
        self.assertEqual(route['relative']['protonated_start'], 0.)

    def test_all_sites_use_reference_from_all_recorded_precursors(self):
        links = self.links + [dict(self.links[0], pair_index=9, binding_permutation=2,
                                  protonated_start_gibbs_eh=-131.)]
        diagram = self.diagram(links=links)
        self.assertEqual(diagram['reference'], -131.)
        self.assertTrue(all(route['relative']['protonated_start'] == 627.509474 for route in diagram['routes']))

    def test_protonated_bound_state_uses_alternative_formation_edges(self):
        links = [dict(row, protonated_bound_requested=True,
                      protonated_bound_with_acetate_gibbs_eh=-129.) for row in self.links]
        diagram = self.diagram(links=links)
        self.assertIn('protonated_bound', diagram['state_keys'])
        self.assertIn(('protonated_start', 'protonated_bound'), diagram['edges'])
        self.assertIn(('protonated_bound', 'reactant'), diagram['edges'])
        self.assertIn(('deprotonated_start', 'encounter'), diagram['edges'])
        self.assertIn(('encounter', 'reactant'), diagram['edges'])
        self.assertEqual(diagram['state_keys'].index('encounter') + 1, diagram['state_keys'].index('reactant'))
        self.assertNotIn(('deprotonated_start', 'protonated_bound'), diagram['edges'])
        self.assertNotIn(('protonated_bound', 'deprotonated_start'), diagram['edges'])
        self.assertEqual(len(diagram['routes']), 4)
        for route in diagram['routes']:
            self.assertEqual(route['energies']['protonated_bound'], -129.)
            self.assertAlmostEqual(route['relative']['protonated_bound'], 627.509474)

    def test_missing_protonated_bound_energy_remains_blank(self):
        links = [dict(row, protonated_bound_requested=True,
                      protonated_bound_with_acetate_gibbs_eh=None) for row in self.links]
        diagram = self.diagram(links=links)
        self.assertIn('protonated_bound', diagram['state_keys'])
        for route in diagram['routes']:
            self.assertIsNone(route['relative']['protonated_bound'])
            self.assertIn('protonated bound reactant + acetate', route['missing'])

    def test_missing_protonated_input_remains_blank_with_explicit_fallback(self):
        links = [dict(row, protonated_start_gibbs_eh=None, acetate_gibbs_eh=None) for row in self.links]
        diagram = self.diagram(links=links)
        self.assertEqual(diagram['reference'], -110.)
        self.assertIn('reference unavailable', diagram['reference_label'])
        for route in diagram['routes']:
            self.assertIsNone(route['relative']['protonated_start'])
            self.assertIn('acetate', route['missing'])
            self.assertIsNotNone(route['relative']['deprotonated_start'])

    def test_missing_acid_uses_original_diagram_without_mixing_energy_conventions(self):
        diagram = self.diagram(acid=None)
        self.assertFalse(diagram['extended'])
        self.assertEqual(diagram['state_keys'], ['encounter', 'reactant', 'TS', 'product'])
        self.assertEqual(diagram['routes'][0]['energies']['encounter'], -100.1)
        self.assertEqual(diagram['reference'], -100.)
        self.assertEqual(diagram['routes'][0]['energies']['reactant'], -100.)
        self.assertEqual(diagram['routes'][0]['relative']['reactant'], 0.)

    def test_missing_link_does_not_hide_completed_reaction_pathway(self):
        diagram = self.diagram(links=[])
        self.assertEqual(len(diagram['routes']), 4)
        self.assertIsNone(diagram['routes'][0]['relative']['deprotonated_start'])
        self.assertEqual(diagram['routes'][0]['relative']['reactant'], 0.)
        self.assertTrue(any('No precursor link' in message for message in diagram['routes'][0]['missing']))

    def test_missing_encounter_stays_blank_and_does_not_change_local_barriers(self):
        links = [dict(row, encounter_gibbs_eh=None, encounter_with_acid_gibbs_eh=None) for row in self.links]
        diagram = self.diagram(links=links)
        for route in diagram['routes']:
            self.assertIsNone(route['relative']['encounter'])
            self.assertIn('relaxed encounter complex', route['missing'])
            self.assertAlmostEqual(route['forward_barrier_kcal_mol'], .02 * 627.509474)

    def test_resonance_and_native_protonation_routes_are_not_collapsed(self):
        links = self.links + [dict(self.links[0], protonation_site=3, protonated_start_gibbs_eh=-129.)]
        diagram = self.diagram(links=links)
        self.assertEqual(len(diagram['routes']), 5)
        self.assertEqual({route['profile']['pair_index'] for route in diagram['routes']}, {1, 2, 3, 4})
        self.assertEqual({route['profile']['resonance'] for route in diagram['routes']}, {0, 1})
        self.assertEqual({route['protonation_site'] for route in diagram['routes']}, {1, 3})

    def test_loader_handles_absent_optional_databases_without_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'workflow.sqlite'
            with closing(sqlite3.connect(path)) as db, db:
                db.execute('CREATE TABLE jobs (id INTEGER PRIMARY KEY)')
                db.execute('INSERT INTO jobs VALUES (1)')
            before = hashlib.sha256(path.read_bytes()).digest()
            result = load_diagram_inputs(path)
            self.assertEqual(result['jobs'], [{'id': 1}])
            self.assertEqual(result['diagram_references'], [])
            self.assertEqual(result['reference_minima'], {})
            self.assertEqual(hashlib.sha256(path.read_bytes()).digest(), before)
            self.assertFalse(path.with_name('precursors.sqlite').exists())
            with self.assertRaises(FileNotFoundError):
                load_diagram_inputs(path.with_name('absent.sqlite'))

    def test_missing_energy_display_and_all_cells_compile(self):
        self.assertEqual(self.namespace['format_energy'](None), 'missing')
        for cell in self.notebook['cells']:
            if cell['cell_type'] == 'code':
                compile(''.join(cell['source']), cell['id'], 'exec')

    def test_substrate_key_matches_sampling_atom_order_and_symmetry_for_all_systems(self):
        for smiles in ('[H]C1=CC(OC)=C(OC)C(OC)=C1', 'COC(=O)C1=C(F)C=CC=C1F',
                       'COC(=O)C1=CC=NC=C1'):
            with self.subTest(substrate=smiles):
                recorded = {int(site.removeprefix('site_'))
                            for site, _ in add_substrate(Chem.MolFromSmiles('[Pd+2]'), smiles)}
                key = self.namespace['substrate_site_key'](smiles, recorded)
                self.assertEqual({site['substrate_site'] for site in key['sites']}, recorded)
                self.assertEqual([site['name'] for site in key['sites']], ['A', 'B'])
                for site in key['sites']:
                    self.assertIn(site['substrate_site'] - 1, site['atoms'])
                    for index in site['atoms']:
                        atom = key['molecule'].GetAtomWithIdx(index)
                        self.assertFalse(atom.HasProp('atomNote'))
                        self.assertTrue(atom.GetIsAromatic())
                        self.assertEqual(atom.GetAtomicNum(), 6)
                        self.assertGreater(atom.GetTotalNumHs(), 0)
                # Missing profiles must not rename B to A.
                partial = self.namespace['substrate_site_key'](smiles, {max(recorded)})
                self.assertEqual(key['sites'], partial['sites'])

    def test_substrate_key_rejects_wrong_atom_indices(self):
        with self.assertRaisesRegex(ValueError, 'do not match substrate atom order'):
            self.namespace['substrate_site_key']('COC(=O)C1=CC=NC=C1', {1})

    def test_site_key_precedes_profiles_and_is_not_repeated_per_ligand(self):
        cells = self.notebook['cells']
        indices = {cell['id']: index for index, cell in enumerate(cells)}
        self.assertLess(indices['substrate-site-key'], indices['560d679e'])
        self.assertEqual(sum(cell['id'] == 'substrate-site-key' for cell in cells), 1)

    def test_both_forward_and_reverse_barriers_must_pass_inclusive_limit(self):
        for forward, reverse, passes in [(24., 24., True), (25., 25., True),
                                         (24., 26., False), (26., 24., False),
                                         (25.004, 24., False)]:
            with self.subTest(forward=forward, reverse=reverse):
                profile = dict(TS=0., reactant=-forward / 627.509474, product=-reverse / 627.509474)
                result = self.namespace['proton_transfer_barriers'](profile)
                self.assertAlmostEqual(result['forward_barrier_kcal_mol'], forward)
                self.assertAlmostEqual(result['reverse_barrier_kcal_mol'], reverse)
                self.assertEqual(result['both_barriers_within_limit'], passes)

    def test_barriers_include_encounter_and_ignore_acid_and_plot_energy_zero(self):
        original = self.diagram()
        shifted_links = [dict(row, protonated_start_gibbs_eh=-200.) for row in self.links]
        shifted = self.diagram(links=shifted_links, acid=-99.)
        for before, after in zip(original['routes'], shifted['routes']):
            for key in ('forward_barrier_kcal_mol', 'reverse_barrier_kcal_mol', 'both_barriers_within_limit'):
                self.assertEqual(before[key], after[key])
        self.assertEqual(len(shifted['routes']), len(self.profiles))

    def test_forward_reference_is_lower_encounter_or_bound_reactant(self):
        scale = 627.509474
        for encounter, expected, reference in [(None, 24., 'reactant'),
                                                (-20. / scale, 24., 'reactant'),
                                                (-24. / scale, 24., 'reactant'),
                                                (-26. / scale, 26., 'encounter')]:
            with self.subTest(encounter=encounter):
                profile = dict(TS=0., reactant=-24. / scale, product=-23. / scale, encounter=encounter)
                result = self.namespace['proton_transfer_barriers'](profile)
                self.assertAlmostEqual(result['forward_barrier_kcal_mol'], expected)
                self.assertAlmostEqual(result['bound_forward_barrier_kcal_mol'], 24.)
                self.assertAlmostEqual(result['reverse_barrier_kcal_mol'], 23.)
                self.assertEqual(result['forward_reference_state'], reference)
                self.assertEqual(result['both_barriers_within_limit'], expected <= 25.)

    def test_route_and_grid_caption_use_same_encounter_adjusted_barrier(self):
        profile = dict(self.profiles[0], encounter=-100.1)
        expected = self.namespace['proton_transfer_barriers'](profile)
        route = self.diagram(profiles=[profile])['routes'][0]
        self.assertEqual(route['forward_reference_state'], 'encounter')
        self.assertEqual(route['forward_reference_gibbs_eh'], -100.1)
        self.assertAlmostEqual(route['forward_barrier_kcal_mol'], .12 * 627.509474)
        self.assertAlmostEqual(route['bound_forward_barrier_kcal_mol'], .02 * 627.509474)
        self.assertAlmostEqual(route['reverse_barrier_kcal_mol'], .03 * 627.509474)
        self.assertEqual(route['forward_barrier_kcal_mol'], expected['forward_barrier_kcal_mol'])
        legend = self.namespace['pathway_grid_legend'](profile, 'N', 17, 0.)
        self.assertIn(f"Forward {expected['forward_barrier_kcal_mol']:.2f}", legend)

    def test_profile_selection_matches_encounter_to_ligand_and_binding_in_active_model(self):
        jobs = []
        for pair, binding in ((1, 1), (2, 2)):
            group = self.coverage_group(pair, ts_status='succeeded', binding_permutation=binding)
            group.append(dict(group[0], id=pair * 10 + 2, endpoint='product'))
            for row in group:
                row.update(binding_mode=1, coordination_isomer=1)
            jobs.extend(group)
        encounters = [dict(id=1, ligand_id='11', binding_permutation=1, status='succeeded', gibbs_free_energy_eh=-100.2),
                      dict(id=2, ligand_id='11', binding_permutation=2, status='succeeded', gibbs_free_energy_eh=-99.5),
                      dict(id=3, ligand_id='12', binding_permutation=1, status='succeeded', gibbs_free_energy_eh=-999.),
                      dict(id=4, ligand_id='11', binding_permutation=1, status='failed', gibbs_free_energy_eh=-999.),
                      dict(id=5, ligand_id='11', binding_permutation=1, status='succeeded', gibbs_free_energy_eh=None)]
        namespace = dict(jobs=jobs, encounter_jobs=encounters, defaultdict=defaultdict)
        cell = next(c for c in self.notebook['cells'] if c.get('id') == '560d679e')
        exec(compile(''.join(cell['source']), 'profile selection', 'exec'), namespace)
        first, second = namespace['profiles']
        self.assertEqual((first['encounter'], first['encounter_job_id']), (-100.2, 1))
        self.assertEqual((second['encounter'], second['encounter_job_id']), (-99.5, 2))
        barrier = self.namespace['proton_transfer_barriers']
        self.assertEqual(barrier(first)['forward_reference_state'], 'encounter')
        self.assertEqual(barrier(second)['forward_reference_state'], 'reactant')

    def coverage_group(self, pair, delta=0., ts_status='failed', ligand='11', **metadata):
        shared = dict(ligand_id=ligand, pair_index=pair, seating=4, binding_permutation=1,
                      substrate_site=6, resonance=0, substrate='CC', charge=1, multiplicity=1,
                      conformer=1, retained=1, updated_at='2026-10-19 12:00:00')
        shared.update(metadata)
        reactant = dict(shared, id=pair * 10, endpoint='reactant', status='succeeded', stage='complete',
                        gibbs_free_energy_eh=-100. + delta / 627.509474)
        ts = dict(shared, id=pair * 10 + 1, endpoint='TS', status=ts_status, ts_strategy='direct',
                  stage='complete' if ts_status == 'succeeded' else 'optts',
                  gibbs_free_energy_eh=-99.9 if ts_status == 'succeeded' else None,
                  error=None if ts_status == 'succeeded' else 'OptTS attempts did not converge')
        return [reactant, ts]

    def test_reactant_reference_includes_failed_paths_and_unretained_conformers(self):
        failed = self.coverage_group(1)
        failed[0].update(retained=0, conformer=7)
        higher_conformer = dict(failed[0], id=12, conformer=1, retained=1, gibbs_free_energy_eh=-99.9)
        failed_energy = dict(failed[0], id=13, status='failed', gibbs_free_energy_eh=-101.)
        rows = failed + [higher_conformer, failed_energy] + self.coverage_group(2, delta=4., ts_status='succeeded')
        original = json.dumps(rows, sort_keys=True)
        report = self.namespace['reactant_ts_coverage'](rows)
        low, high = report['paths']
        self.assertEqual(low['reactant_job_id'], failed[0]['id'])
        self.assertEqual(low['reactant_conformer'], 7)
        self.assertEqual(low['delta_reactant_kcal_mol'], 0.)
        self.assertAlmostEqual(high['delta_reactant_kcal_mol'], 4.)
        self.assertEqual(len(report['low_energy_failed']), 1)
        self.assertEqual(report['ligands'][0]['minimum_failed'], 1)
        self.assertEqual(json.dumps(rows, sort_keys=True), original)

    def test_grid_caption_uses_failed_group_minimum_and_keeps_both_local_barriers(self):
        rows = self.coverage_group(1) + self.coverage_group(2, delta=4., ts_status='succeeded')
        report = self.namespace['reactant_ts_coverage'](rows)
        relative = {tuple(row[key] for key in self.keys): row['delta_reactant_kcal_mol']
                    for row in report['paths']}
        reactant = rows[2]['gibbs_free_energy_eh']
        profile = dict(ligand_id='11', pair_index=2, seating=4, binding_permutation=1,
                       substrate_site=6, resonance=0, reactant=reactant,
                       TS=reactant + 12. / 627.509474, product=reactant - 3. / 627.509474)
        legend = self.namespace['pathway_grid_legend'](
            profile, 'N', 17, relative[tuple(profile[key] for key in self.keys)])
        self.assertIn('Forward 12.00; reverse 15.00 kcal/mol', legend)
        self.assertIn('Relative reactant G: +4.00 kcal/mol', legend)

    def test_TS_recovery_is_counted_once_and_uses_terminal_NEB_reason(self):
        rows = self.coverage_group(1)
        rows.append(dict(rows[-1], id=12, ts_strategy='neb', status='succeeded', error=None))
        second = self.coverage_group(2, delta=1.)
        second.append(dict(second[-1], id=22, ts_strategy='neb', stage='freq', error='Wrong final NEB mode'))
        report = self.namespace['reactant_ts_coverage'](rows + second)
        self.assertEqual(report['totals']['complete'], 1)
        self.assertEqual(report['totals']['failed'], 1)
        self.assertEqual(report['low_energy_failed'][0]['terminal_job_id'], 22)
        self.assertEqual(report['low_energy_failed'][0]['terminal_reason'], 'Wrong final NEB mode')

    def test_pending_and_running_NEB_recovery_are_not_terminal_failures(self):
        rows = []
        for pair, status in [(1, 'pending'), (2, 'running')]:
            group = self.coverage_group(pair)
            group.append(dict(group[-1], id=pair * 10 + 2, ts_strategy='neb', status=status))
            rows.extend(group)
        report = self.namespace['reactant_ts_coverage'](rows)
        self.assertEqual(report['totals']['active'], 2)
        self.assertEqual(report['totals']['failed'], 0)
        self.assertFalse(report['low_energy_failed'])

    def test_missing_DFT_singlepoint_does_not_turn_a_completed_TS_into_a_failure(self):
        rows = self.coverage_group(1, ts_status='succeeded')
        rows[-1].update(gibbs_free_energy_eh=None, gxtb_gibbs_free_energy_eh=-99.9,
                        singlepoint_selected=True, singlepoint_status='not_queued')
        second = self.coverage_group(2)
        second[0].update(gibbs_free_energy_eh=None, gxtb_gibbs_free_energy_eh=-100.)
        report = self.namespace['reactant_ts_coverage'](rows + second)
        self.assertEqual(report['totals']['complete'], 1)
        self.assertEqual(report['totals']['failed'], 1)
        self.assertFalse(report['low_energy_failed'])
        self.assertEqual(len(report['unranked']), 1)
        self.assertIn('DFT single point unavailable', report['unranked'][0]['reactant_energy_reason'])
        self.assertEqual(report['windows'][-1]['total'], 1)

    def test_resonance_paths_are_preserved_and_sites_share_the_ligand_reference(self):
        rows = (self.coverage_group(1) + self.coverage_group(2, delta=2., resonance=1)
                + self.coverage_group(3, delta=5., substrate_site=7, binding_permutation=2))
        report = self.namespace['reactant_ts_coverage'](rows)
        self.assertEqual(len(report['paths']), 3)
        self.assertEqual([row['pair_index'] for row in report['low_energy_failed']], [1, 2])
        self.assertAlmostEqual(report['paths'][-1]['delta_reactant_kcal_mol'], 5.)

    def test_energy_references_keep_ligands_and_different_compositions_separate(self):
        rows = (self.coverage_group(1) + self.coverage_group(2, delta=9., ligand='12')
                + self.coverage_group(3, delta=20., substrate='CCC')
                + self.coverage_group(4, delta=30., charge=2)
                + self.coverage_group(5, delta=40., multiplicity=3))
        report = self.namespace['reactant_ts_coverage'](rows)
        self.assertTrue(all(row['delta_reactant_kcal_mol'] == 0 for row in report['paths']))

    def test_energy_windows_are_inclusive_and_missing_energies_are_separate(self):
        rows = self.coverage_group(1) + self.coverage_group(2, delta=3.) + self.coverage_group(3, delta=3.0001)
        missing = self.coverage_group(4)
        missing[0].update(status='failed', gibbs_free_energy_eh=-1000.)
        rows.extend(missing)
        report = self.namespace['reactant_ts_coverage'](rows, window_kcal=3.)
        self.assertEqual([row['pair_index'] for row in report['low_energy_failed']], [1, 2])
        self.assertEqual(report['totals']['failed'], 4)
        self.assertEqual(len(report['unranked']), 1)
        self.assertEqual(next(row for row in report['windows'] if row['limit_kcal_mol'] == 0.)['failed'], 1)
        self.assertEqual(report['windows'][-1]['total'], 3)
        for invalid in [-1., float('nan'), float('inf')]:
            with self.assertRaises(ValueError):
                self.namespace['reactant_ts_coverage'](rows, invalid)

    def test_unprepared_TS_and_nonfinite_reactant_energies_are_reported(self):
        rows = self.coverage_group(1)[:1]
        rows[0]['gibbs_free_energy_eh'] = float('nan')
        rows.append(dict(rows[0], id=12, gibbs_free_energy_eh=float('inf')))
        report = self.namespace['reactant_ts_coverage'](rows)
        self.assertEqual(report['totals']['not_prepared'], 1)
        self.assertEqual(report['totals']['failed'], 0)
        self.assertEqual(len(report['unranked']), 1)


if __name__ == '__main__':
    unittest.main()
