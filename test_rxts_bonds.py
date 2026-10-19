"""Verify rxts's returned bond changes before accepting a TS guess."""

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
import rxts
from rdkit import Chem
from rxts.structure import Frame, ReactionPath

from rxflow_endpoints import check_rxts_bonds, check_rxts_guess_motion, generate_ts_guess, mapped_molecule


def endpoints():
    reactant = mapped_molecule('[C:1]([H:2])#[N:3]')
    product = mapped_molecule('[C-:1]#[N+:3][H:2]')
    for molecule, hydrogen_x in ((reactant, 1.1), (product, 1.6)):
        conformer = Chem.Conformer(3)
        for i, position in enumerate(((0., 0., 0.), (hydrogen_x, 0., 0.), (2.7, 0., 0.))):
            conformer.SetAtomPosition(i, position)
        molecule.AddConformer(conformer)
    return reactant, product


def scan_path(reactant, product):
    start = reactant.GetConformer().GetPositions()
    end = product.GetConformer().GetPositions()
    frames = [Frame(lam=lam, positions=(1 - lam) * start + lam * end,
                    energy=energy, e_reactant=energy, e_product=energy)
              for lam, energy in ((0., 0.), (.5, 1.), (1., 0.))]
    return ReactionPath(frames=frames, reactant=reactant, product=product,
                        forming_bonds=[(1, 2)], breaking_bonds=[(0, 1)], mapping=[0, 1, 2],
                        engine='gfnff', charge=0, multiplicity=1, solvent=None,
                        n_cores=1, mol=reactant)


class ReturnedBondTests(unittest.TestCase):
    def setUp(self):
        self.reactant, self.product = endpoints()
        self.path = SimpleNamespace(reactant=self.reactant, product=self.product,
                                    breaking_bonds=[(0, 1)], forming_bonds=[(1, 2)])

    def test_rxts_really_infers_expected_bonds_from_mapped_endpoints(self):
        guesser = rxts.TransitionStateGuesser(engine='gfnff', gfnff_topology='graph',
                                              charge=0, multiplicity=1, n_embed=0, mapping=[0, 1, 2])
        prepared = guesser.prepare_endpoints(self.reactant, self.product)
        report = check_rxts_bonds(prepared, (0, 1, 2))
        self.assertTrue(report['valid'], report)
        self.assertEqual(report['breaking_bonds'], [[0, 1]])
        self.assertEqual(report['forming_bonds'], [[1, 2]])

    def test_additional_changes_and_reversed_pair_order_are_allowed(self):
        self.path.breaking_bonds = [(1, 0), (2, 0)]
        self.path.forming_bonds = [(2, 1), (0, 2)]
        report = check_rxts_bonds(self.path, (0, 1, 2))
        self.assertTrue(report['valid'], report)

    def test_native_guess_bond_lists_are_checked(self):
        path = scan_path(self.reactant, self.product)
        guess = path.guess('ff')
        self.assertTrue(check_rxts_bonds(path, (0, 1, 2), guess=guess)['valid'])
        guess.breaking_bonds.clear()
        self.assertTrue(path.breaking_bonds)
        self.assertFalse(check_rxts_bonds(path, (0, 1, 2), guess=guess)['valid'])

    def test_missing_expected_ch_break_is_rejected(self):
        self.path.breaking_bonds = [(0, 2)]
        report = check_rxts_bonds(self.path, (0, 1, 2))
        self.assertFalse(report['valid'])
        self.assertIn('C-H bond [0, 1] breaking', report['reason'])

    def test_forming_bond_must_use_the_same_h_and_expected_acceptor(self):
        self.path.forming_bonds = [(0, 2)]
        report = check_rxts_bonds(self.path, (0, 1, 2))
        self.assertFalse(report['valid'])
        self.assertIn('H-acceptor bond [1, 2] forming', report['reason'])

    def test_missing_malformed_or_out_of_range_report_is_rejected(self):
        for value in (None, [(0,)], [(0, 9)], [(True, 1)], ['bad']):
            with self.subTest(value=value):
                self.path.breaking_bonds = value
                self.assertFalse(check_rxts_bonds(self.path, (0, 1, 2))['valid'])

    def test_non_carbon_h_donor_is_rejected_even_if_reported_pairs_match(self):
        reactant = mapped_molecule('[O:1]([H:2])[H:3].[N:4]([H:5])([H:6])[H:7]')
        product = mapped_molecule('[O-:1][H:3].[N+:4]([H:2])([H:5])([H:6])[H:7]')
        path = SimpleNamespace(reactant=reactant, product=product,
                               breaking_bonds=[(0, 1)], forming_bonds=[(1, 3)])
        report = check_rxts_bonds(path, (0, 1, 3))
        self.assertFalse(report['valid'])
        self.assertIn('donor is not carbon', report['reason'])


class GuessMotionTests(unittest.TestCase):
    def setUp(self):
        self.path = scan_path(*endpoints())

    def test_native_bond_length_series_identifies_both_changing_bonds(self):
        report = check_rxts_guess_motion(self.path, (0, 1, 2))
        self.assertTrue(report['valid'], report)
        self.assertEqual(report['detected_atoms'], [0, 1, 2])
        self.assertEqual({tuple(b['atoms']) for b in report['bond_changes']}, {(0, 1), (1, 2)})
        self.assertGreater(report['local_bond_deltas_angstrom']['C_H'], 0)
        self.assertLess(report['local_bond_deltas_angstrom']['H_acceptor'], 0)

    def test_correct_metadata_does_not_pass_when_ch_is_not_breaking(self):
        self.assertTrue(check_rxts_bonds(self.path, (0, 1, 2))['valid'])
        for frame in self.path.frames:
            frame.positions[1, 0] = 1.35
            frame.positions[2, 0] = 2.6 + .2 * frame.lam
        report = check_rxts_guess_motion(self.path, (0, 1, 2))
        self.assertFalse(report['valid'], report)
        self.assertIn('lengthen C-H', report['reason'])

    def test_both_expected_bond_length_series_are_required(self):
        with patch.object(self.path, 'bond_lengths', return_value={('breaking', (0, 1)): [1., 1.3, 1.6]}):
            report = check_rxts_guess_motion(self.path, (0, 1, 2))
        self.assertFalse(report['valid'])
        self.assertIn('both expected', report['reason'])

    def test_reverse_transfer_direction_is_rejected(self):
        for frame in self.path.frames:
            frame.positions[1, 0] = 1.6 - .5 * frame.lam
        report = check_rxts_guess_motion(self.path, (0, 1, 2))
        self.assertFalse(report['valid'])
        self.assertIn('lengthen C-H', report['reason'])

    def test_endpoint_peak_or_no_local_motion_is_rejected(self):
        self.path.frames[0].energy = 10.
        self.assertIn('both sides', check_rxts_guess_motion(self.path, (0, 1, 2))['reason'])
        self.path.frames[0].energy = 0.
        for frame in self.path.frames:
            frame.positions = np.array(self.path.frames[1].positions)
        self.assertFalse(check_rxts_guess_motion(self.path, (0, 1, 2))['valid'])

    def test_rigid_translation_and_rotation_do_not_count_as_internal_motion(self):
        reference = np.array(self.path.frames[1].positions)
        for frame in self.path.frames:
            angle = frame.lam
            rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                                 [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
            frame.positions = reference @ rotation + frame.lam * np.array([2., 3., 4.])
        self.assertFalse(check_rxts_guess_motion(self.path, (0, 1, 2))['valid'])


class GuessCreationTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.reactant, self.product = endpoints()
        self.paired = MagicMock()
        self.paired.ts_rmsds = {7: .1}
        self.paired.errors = []
        self.paired.correspondence = {0: 0, 1: 1, 2: 2}
        self.paired.ids = [7]
        self.paired.__len__.return_value = 1
        self.paired.__getitem__.return_value = (self.reactant, self.product)
        self.path = scan_path(self.reactant, self.product)

    def generate(self):
        with patch('rxflow_endpoints.rxembed.coembed', return_value=self.paired), \
                patch('rxts.TransitionStateGuesser') as guesser:
            guesser.return_value.find_transition_state.return_value = self.path
            result = generate_ts_guess(self.reactant, self.product, self.root)
            self.assertNotIn('forming_bonds', guesser.call_args.kwargs)
            self.assertNotIn('breaking_bonds', guesser.call_args.kwargs)
            self.assertEqual(guesser.call_args.kwargs['mapping'], [0, 1, 2])
        return result

    def test_accepted_guess_records_returned_bond_report_with_profile(self):
        selection = self.generate()
        self.assertTrue(selection['rxts_bond_check']['valid'])
        saved = json.loads((self.root/'ts_guess/reaction_selection.json').read_text())
        self.assertEqual(saved['rxts_bond_check']['breaking_bonds'], [[0, 1]])
        self.assertEqual(saved['rxts_bond_check']['forming_bonds'], [[1, 2]])
        self.assertTrue(saved['rxts_motion_check']['valid'])
        for name in ('reaction.xyz', 'reaction_path.xyz', 'reaction.inp'):
            self.assertTrue((self.root/'ts_guess'/name).is_file())

    def test_two_recoembedding_retries_use_distinct_seeds_and_accept_third_guess(self):
        bad = scan_path(self.reactant, self.product)
        bad.breaking_bonds = [(0, 2)]
        with patch('rxflow_endpoints.rxembed.coembed', return_value=self.paired) as coembed, \
                patch('rxts.TransitionStateGuesser') as guesser:
            guesser.return_value.find_transition_state.side_effect = [bad, bad, self.path]
            selection = generate_ts_guess(self.reactant, self.product, self.root, seed=42)
        self.assertEqual([c.kwargs['seed'] for c in coembed.call_args_list], [42, 43, 44])
        self.assertTrue(all(c.args == (self.reactant, self.product) for c in coembed.call_args_list))
        self.assertEqual([a['valid'] for a in selection['attempts']], [False, False, True])
        self.assertEqual(selection['seed'], 44)
        self.assertNotIn('error', selection)
        for attempt in (1, 2, 3):
            self.assertTrue((self.root/'ts_guess'/f'reaction_attempt_{attempt:02d}_path.xyz').exists())
        self.assertEqual((self.root/'reactant/reaction_reactive.xyz').read_text(),
                         Path(selection['reactant_xyz']).read_text())

    def test_three_bad_guesses_fail_without_a_fourth_attempt(self):
        self.path.forming_bonds = [(0, 2)]
        with patch('rxflow_endpoints.rxembed.coembed', return_value=self.paired) as coembed, \
                patch('rxts.TransitionStateGuesser') as guesser:
            guesser.return_value.find_transition_state.return_value = self.path
            with self.assertRaisesRegex(ValueError, 'at most 3 attempts'):
                generate_ts_guess(self.reactant, self.product, self.root, seed=42)
        self.assertEqual(coembed.call_count, 3)
        self.assertEqual(guesser.return_value.find_transition_state.call_count, 3)
        saved = json.loads((self.root/'ts_guess/reaction_selection.json').read_text())
        self.assertEqual(len(saved['attempts']), 3)
        self.assertTrue(all(not a['valid'] for a in saved['attempts']))
        self.assertFalse((self.root/'ts_guess/reaction.inp').exists())

    def test_successful_first_guess_does_not_recoembed(self):
        with patch('rxflow_endpoints.rxembed.coembed', return_value=self.paired) as coembed, \
                patch('rxts.TransitionStateGuesser') as guesser:
            guesser.return_value.find_transition_state.return_value = self.path
            selection = generate_ts_guess(self.reactant, self.product, self.root)
        self.assertEqual(coembed.call_count, 1)
        self.assertEqual(len(selection['attempts']), 1)

    def test_failed_first_candidate_does_not_poison_a_valid_second_candidate(self):
        self.paired.ts_rmsds = {7: .1, 8: .2}
        self.paired.ids = [7, 8]
        bad = scan_path(self.reactant, self.product)
        bad.breaking_bonds = [(0, 2)]
        with patch('rxflow_endpoints.rxembed.coembed', return_value=self.paired), \
                patch('rxts.TransitionStateGuesser') as guesser:
            guesser.return_value.find_transition_state.side_effect = [bad, bad, bad, self.path]
            selection = generate_ts_guess(self.reactant, self.product, self.root, keep=2)
        self.assertNotIn('error', selection)
        self.assertEqual(selection['rank'], 2)
        self.assertIn('error', selection['candidates'][0])
        self.assertTrue(selection['candidates'][1]['valid'])
        self.assertEqual((self.root/'ts_guess/reaction.xyz').read_text(),
                         (self.root/'ts_guess/reaction_conformer_002.xyz').read_text())

    def test_rejected_guess_preserves_profile_and_records_failure_before_ts_input(self):
        self.path.breaking_bonds = [(0, 2)]
        with self.assertRaisesRegex(ValueError, 'No TS guesses survived'):
            self.generate()
        saved = json.loads((self.root/'ts_guess/reaction_selection.json').read_text())
        self.assertFalse(saved['rxts_bond_check']['valid'])
        self.assertIn('C-H bond', saved['error'])
        self.assertNotIn('ts_input', saved)
        self.assertFalse((self.root/'ts_guess/reaction.inp').exists())
        self.assertTrue((self.root/'ts_guess/reaction_path.xyz').exists())
        self.assertTrue((self.root/'ts_guess/reaction.xyz').exists())

    def test_correct_reported_bonds_but_wrong_actual_motion_cannot_create_ts_input(self):
        for frame in self.path.frames:
            frame.positions[1, 0] = 1.35
            frame.positions[2, 0] = 2.6 + .2 * frame.lam
        with self.assertRaisesRegex(ValueError, 'No TS guesses survived'):
            self.generate()
        saved = json.loads((self.root/'ts_guess/reaction_selection.json').read_text())
        self.assertTrue(saved['rxts_bond_check']['valid'])
        self.assertFalse(saved['rxts_motion_check']['valid'])
        self.assertNotIn('ts_input', saved)
        self.assertTrue((self.root/'ts_guess/reaction_path.xyz').exists())
        self.assertFalse((self.root/'ts_guess/reaction.inp').exists())


if __name__ == '__main__':
    unittest.main()
