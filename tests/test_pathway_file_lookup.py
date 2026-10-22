"""Check that notebook pathway lookup preserves selected geometries and TS modes."""
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock


class PathwayFileLookupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name).resolve()
        profiles = []
        for ligand, ts_id, seating in [('5', 2120, 1), ('6', 2121, 1), ('5', 2122, 2)]:
            profile = dict(ligand_id=ligand, substrate_site=2, binding_permutation=1,
                           pair_index=2, resonance=0, seating=seating)
            for i, name in enumerate(('reactant', 'ts', 'product')):
                filename = f'{ts_id}_{name}.xyz'
                (root / filename).write_text('1\nselected geometry\nH 0 0 0\n')
                profile[f'{name}_row'] = dict(id=ts_id if name == 'ts' else ts_id + 100 + i,
                                               xyz_file=filename)
            profile['ts_row'].update(frequencies_json='[0, -1327.21, 20.0]',
                                     donor=0, hydrogen=1, acceptor=2)
            profiles.append(profile)
        renderer = Mock(return_value='rendered')
        self.ns = dict(DB_PATH=root / 'workflow.sqlite', profiles=profiles, math=math,
                       xyzrender=SimpleNamespace(render=renderer), display=Mock(), Markdown=str)
        notebook = json.loads((Path(__file__).resolve().parents[1] / 'reaction_diagrams.ipynb').read_text())
        cell = next(c for c in notebook['cells'] if c.get('id') == 'pathway-file-lookup')
        exec(compile(''.join(cell['source']), 'pathway-file-lookup', 'exec'), self.ns)
        self.root = root

    def test_numeric_and_readable_ids_return_exact_selected_files(self):
        files = self.ns['get_files'](2120)
        self.assertEqual(list(files), ['reactant', 'TS', 'product'])
        self.assertEqual(files['TS'], self.root / '2120_ts.xyz')
        self.assertEqual(files, self.ns['get_files']('L5/2-1-2-0/s1'))
        self.assertEqual(files, self.ns['get_files']('2120'))

    def test_ligand_and_seating_keep_otherwise_identical_keys_distinct(self):
        self.assertEqual(self.ns['get_files']('L6/2-1-2-0/s1')['TS'], self.root / '2121_ts.xyz')
        self.assertEqual(self.ns['get_files']('L5/2-1-2-0/s2')['TS'], self.root / '2122_ts.xyz')
        with self.assertRaises(KeyError):
            self.ns['get_files']('2-1-2-0')
        with self.assertRaises(KeyError):
            self.ns['get_files'](9999)

    def test_render_uses_selected_structures_and_ts_reacting_bonds(self):
        files = self.ns['get_files'](2120, render=True)
        calls = self.ns['xyzrender'].render.call_args_list
        self.assertEqual([call.args[0] for call in calls], list(files.values()))
        self.assertEqual(calls[1].kwargs['ts_bonds'], [(1, 2), (2, 3)])
        self.assertNotIn('ts_bonds', calls[0].kwargs)
        self.assertNotIn('ts_bonds', calls[2].kwargs)

    def test_missing_xyz_does_not_fall_back_to_another_conformer(self):
        (self.root / '2120_product.xyz').unlink()
        with self.assertRaisesRegex(FileNotFoundError, 'product'):
            self.ns['get_files'](2120, render=True)
        self.ns['xyzrender'].render.assert_not_called()

    def test_frequency_report_distinguishes_multiple_missing_and_zero_modes(self):
        profile = self.ns['profiles'][0]
        self.assertEqual(self.ns['format_ts_imaginary_frequencies'](profile), '-1327.21')
        profile['ts_row']['frequencies_json'] = '[-1327.21, -50.1, 0, 20]'
        self.assertEqual(self.ns['ts_imaginary_frequencies'](profile), (-1327.21, -50.1))
        profile['ts_row']['frequencies_json'] = None
        self.assertEqual(self.ns['format_ts_imaginary_frequencies'](profile), 'missing')
        profile['ts_row']['frequencies_json'] = '[0, 20]'
        self.assertEqual(self.ns['format_ts_imaginary_frequencies'](profile), 'none')


if __name__ == '__main__':
    unittest.main()
