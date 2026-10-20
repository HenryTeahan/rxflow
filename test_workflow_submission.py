"""Slurm dependencies and off-node precursor preparation for submit/resume."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hpc_workflow import submit_workflow


class WorkflowSubmissionTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.config_path = self.base / 'run config.toml'
        self.config = dict(paths=dict(results='results', source='seeds', ligands_csv='ligands.csv'),
                           calculation=dict(cpus=1, memory_mb=2000))
        self.root = self.base / 'results'
        self.root.mkdir()
        self.database = self.root / 'workflow.sqlite'
        self.database.write_bytes(b'unchanged reaction results')

    def submit(self, **kwargs):
        with patch('hpc_workflow.subprocess.check_output',
                   side_effect=[f'{i};cluster\n' for i in range(100, 110)]) as sbatch, \
                patch('hpc_precursors.prepare_precursors', side_effect=AssertionError('No local preparation')), \
                patch('hpc_encounters.sample_encounter', side_effect=AssertionError('No local embedding')), \
                patch('hpc_singlepoint.prepare_singlepoints', side_effect=AssertionError('No early SP selection')):
            result = submit_workflow(self.config_path, self.config, **kwargs)
        return result, [call.args[0] for call in sbatch.call_args_list]

    def test_resume_prepares_all_precursors_on_slurm_and_preserves_results(self):
        result, calls = self.submit(resume=True, workers=7, partition='kemi1', account='chemistry')
        self.assertEqual(result, dict(orca='100', precursor_prepare='101', precursors='102',
                                      singlepoint_submit='103'))
        reaction, prepare, precursors, singlepoints = calls
        self.assertIn(' worker', reaction[-1])
        self.assertIn('--array=0-6%7', reaction)
        self.assertIn('precursors prepare --conformers 10', prepare[-1])
        self.assertNotIn('--ligand-id', prepare[-1])
        self.assertFalse(any(arg.startswith('--array=') for arg in prepare))
        self.assertIn('--dependency=afterok:101', precursors)
        self.assertIn('precursors worker', precursors[-1])
        self.assertIn('--array=0-6%7', precursors)
        self.assertIn('--dependency=afterany:100:102', singlepoints)
        self.assertIn('singlepoint submit --workers 10 --account chemistry --partition kemi1', singlepoints[-1])
        for args in calls:
            self.assertIn('--partition=kemi1', args)
            self.assertIn('--account=chemistry', args)
            self.assertIn('--export=ALL', args)
            self.assertIn(str(self.config_path), args[-1])
        self.assertEqual(self.database.read_bytes(), b'unchanged reaction results')
        self.assertFalse((self.root / 'precursors.sqlite').exists())
        self.assertFalse((self.root / 'encounters.sqlite').exists())

    def test_new_run_prepares_precursors_only_after_reaction_database_import(self):
        self.database.unlink()
        with patch('hpc_workflow.read_ligands', return_value=[{'ID': '1'}, {'ID': '2'}]):
            result, calls = self.submit(workers=3, conformers=6)
        self.assertEqual(result, dict(sampling='100', prepare='101', orca='102',
            precursor_prepare='103', precursors='104', singlepoint_submit='105'))
        sampling, prepare, reaction, precursor_prepare, precursors, singlepoints = calls
        self.assertIn('--array=0-1', sampling)
        self.assertIn('--dependency=afterany:100', prepare)
        self.assertIn('--dependency=afterok:101', reaction)
        self.assertIn('--dependency=afterok:101', precursor_prepare)
        self.assertIn('precursors prepare --conformers 6', precursor_prepare[-1])
        self.assertIn('--dependency=afterok:103', precursors)
        self.assertIn('--dependency=afterany:102:104', singlepoints)

    def test_resume_does_not_require_an_existing_precursor_database(self):
        _, calls = self.submit(resume=True)
        self.assertEqual(len(calls), 4)
        self.assertIn('precursors prepare', calls[1][-1])

    def test_resume_requires_reaction_database_and_positive_resource_counts(self):
        with patch('hpc_workflow.subprocess.check_output') as sbatch:
            for invalid in (dict(workers=0), dict(conformers=0)):
                with self.assertRaises(ValueError):
                    submit_workflow(self.config_path, self.config, resume=True, **invalid)
            self.database.unlink()
            with self.assertRaisesRegex(ValueError, 'No prepared workflow database'):
                submit_workflow(self.config_path, self.config, resume=True)
            sbatch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
