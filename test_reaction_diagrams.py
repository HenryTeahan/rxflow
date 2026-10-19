"""Energy accounting and pathway preservation in the reaction-diagram notebook."""

from collections import Counter, defaultdict
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from hpc_precursors import load_diagram_inputs


class ReactionDiagramTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.notebook = json.loads(Path(__file__).with_name('reaction_diagrams.ipynb').read_text())
        cls.keys = ('ligand_id', 'pair_index', 'seating', 'binding_permutation')
        namespace = dict(Counter=Counter, defaultdict=defaultdict, jobs=[], profiles=[], keys=cls.keys,
                         HARTREE_TO_KCAL=627.509474, BARRIER_LIMIT_KCAL=25.0)
        helper = next(cell for cell in cls.notebook['cells'] if cell.get('id') == 'permutation-coverage')
        exec(compile(''.join(helper['source']), 'notebook helpers', 'exec'), namespace)
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
                          acetate_gibbs_eh=-9.) for profile in self.profiles]

    def diagram(self, profiles=None, links=None, acid=-10.):
        return self.namespace['energy_diagram']('11', self.profiles if profiles is None else profiles,
                                               self.links if links is None else links, acid)

    def test_five_states_balance_spectators_without_counting_substrate_twice(self):
        diagram = self.diagram()
        self.assertEqual(diagram['reference'], -130.)
        self.assertEqual(diagram['state_keys'], ['protonated_start', 'deprotonated_start', 'reactant', 'TS', 'product'])
        route = diagram['routes'][0]
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
        self.assertIn(('deprotonated_start', 'reactant'), diagram['edges'])
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
        self.assertEqual(diagram['state_keys'], ['reactant', 'TS', 'product'])
        self.assertEqual(diagram['reference'], -100.)
        self.assertEqual(diagram['routes'][0]['energies']['reactant'], -100.)
        self.assertEqual(diagram['routes'][0]['relative']['reactant'], 0.)

    def test_missing_link_does_not_hide_completed_reaction_pathway(self):
        diagram = self.diagram(links=[])
        self.assertEqual(len(diagram['routes']), 4)
        self.assertIsNone(diagram['routes'][0]['relative']['deprotonated_start'])
        self.assertEqual(diagram['routes'][0]['relative']['reactant'], 0.)
        self.assertIn('No precursor link', diagram['routes'][0]['missing'][0])

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

    def test_barriers_use_bound_reaction_step_and_ignore_precursor_energy_zero(self):
        original = self.diagram()
        shifted_links = [dict(row, protonated_start_gibbs_eh=-200.) for row in self.links]
        shifted = self.diagram(links=shifted_links, acid=-99.)
        for before, after in zip(original['routes'], shifted['routes']):
            for key in ('forward_barrier_kcal_mol', 'reverse_barrier_kcal_mol', 'both_barriers_within_limit'):
                self.assertEqual(before[key], after[key])
        self.assertEqual(len(shifted['routes']), len(self.profiles))


if __name__ == '__main__':
    unittest.main()
