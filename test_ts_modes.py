"""Normal-mode atom coverage, audit idempotency, and recovery scheduling."""

from contextlib import closing
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hpc_workflow import audit_transition_states, claim_job, connect_database
from orca_hpc import check_ts_mode
from hpc_singlepoint import source_fingerprint


SMILES = '[C:1]([H:2])#[N:3]'
XYZ = '3\nTS fixture\nC 0 0 0\nH 1.3 0 0\nN 2.6 0 0\n'


def orca_output(displacement=0.5, imaginary_column=6, spectator=False):
    """Small ORCA-style NumFreq output with a controlled H-transfer eigenvector."""
    atoms = [('C', 0., 0., 0.), ('H', 1.3, 0., 0.), ('N', 2.6, 0., 0.)]
    if spectator:
        atoms.append(('He', 0., 4., 0.))
    dof = 3 * len(atoms)
    lines = ['O   R   C   A', 'CARTESIAN COORDINATES (ANGSTROEM)', '--------------------']
    lines += [f'{symbol} {x:.6f} {y:.6f} {z:.6f}' for symbol, x, y, z in atoms]
    lines += ['', 'VIBRATIONAL FREQUENCIES', '--------------------']
    for mode in range(dof):
        frequency = 0. if mode < 6 else (-500. if mode == imaginary_column else 200. + mode)
        lines.append(f'{mode}: {frequency:.4f} cm**-1')
    lines += ['NORMAL MODES', '------------', 'Cartesian displacements', '']
    for start in range(0, dof, 6):
        columns = list(range(start, min(start + 6, dof)))
        lines.append(' '.join(str(column) for column in columns))
        for row in range(dof):
            values = [displacement if row == 3 and column == imaginary_column else 0. for column in columns]
            lines.append(f'{row} ' + ' '.join(f'{value:.6f}' for value in values))
        lines.append('')
    lines += ['THERMOCHEMISTRY', 'Electronic energy ... -10.000000 Eh',
              'Final Gibbs free energy ... -9.900000 Eh', 'ORCA TERMINATED NORMALLY']
    return '\n'.join(lines) + '\n'


class ModeCoverageTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.output = self.root / 'frequency.out'
        self.job = dict(donor=0, hydrogen=1, acceptor=2, smiles=SMILES, charge=0, multiplicity=1)

    def test_real_graphrc_reads_out_and_detects_all_expected_transfer_atoms(self):
        self.output.write_text(orca_output())
        result = check_ts_mode(self.output, self.job)
        self.assertTrue(result['valid'], result)
        self.assertEqual(result['missing_atoms'], [])
        self.assertEqual(result['detected_atoms'], [0, 1, 2])
        self.assertEqual(result['orca_mode_index'], 6)
        self.assertEqual(result['graphrc_mode_index'], 0)
        self.assertEqual(list(self.root.iterdir()), [self.output])

    def test_imaginary_mode_without_expected_motion_is_rejected(self):
        self.output.write_text(orca_output(displacement=0.001))
        result = check_ts_mode(self.output, self.job)
        self.assertFalse(result['valid'])
        self.assertEqual(result['missing_atoms'], [0, 1, 2])

    def test_additional_participating_atoms_are_allowed(self):
        self.output.write_text(orca_output(spectator=True))
        job = dict(self.job, smiles=SMILES + '.[He:4]')
        analysis = dict(trajectory=dict(frames=[dict(symbols=['C', 'H', 'N', 'He'])],
                                        frequencies=[-500., 207., 208., 209., 210., 211.]),
                        vibrational=dict(bond_changes={(0, 1): (1., 1.3), (1, 2): (1., 1.3), (0, 3): (.5, 4.)},
                                         angle_changes={}, dihedral_changes={}, thresholds_relaxed=False))
        with patch('graphrc.run_vib_analysis', return_value=analysis):
            result = check_ts_mode(self.output, job)
        self.assertTrue(result['valid'])
        self.assertEqual(result['detected_atoms'], [0, 1, 2, 3])

    def test_analyzes_the_actual_imaginary_mode_not_an_unrelated_first_mode(self):
        self.output.write_text(orca_output(imaginary_column=7))
        result = check_ts_mode(self.output, self.job)
        self.assertTrue(result['valid'])
        self.assertEqual(result['orca_mode_index'], 7)
        self.assertEqual(result['graphrc_mode_index'], 1)

    def test_missing_ts_indices_are_not_silently_accepted(self):
        self.output.write_text(orca_output())
        with self.assertRaisesRegex(ValueError, 'zero-based donor'):
            check_ts_mode(self.output, dict(self.job, hydrogen=None))


class ModeAuditTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.database = self.root / 'workflow.sqlite'
        with closing(connect_database(self.database)):
            pass

    def add(self, identifier, *, valid=True, strategy='direct', status='succeeded', gibbs=-9.9):
        seed = self.root / f'seed_{identifier}.xyz'
        output = self.root / f'frequency_{identifier}.out'
        seed.write_text(XYZ)
        output.write_text(orca_output(displacement=.5 if valid else .001))
        with closing(connect_database(self.database)) as db:
            db.execute('''INSERT INTO jobs (id,ligand_id,pair_index,seating,endpoint,conformer,status,stage,
                attempt,input_xyz,xyz_file,smiles,charge,multiplicity,settings_json,donor,hydrogen,acceptor,
                electronic_energy_eh,gibbs_free_energy_eh,imaginary_count,out_file,ts_strategy,retained)
                VALUES (?, '1',1,1,'TS',?,?, 'complete',1,?,?,?,0,1,'{}',0,1,2,-10,?,1,?,?,1)''',
                (identifier, 0 if strategy == 'neb' else identifier, status, seed.name, seed.name, SMILES,
                 gibbs, output.name, strategy))
        return output

    def row(self, identifier):
        with closing(connect_database(self.database)) as db:
            return dict(db.execute('SELECT * FROM jobs WHERE id=?', (identifier,)).fetchone())

    def test_dry_run_leaves_database_and_outputs_untouched(self):
        self.add(1)
        self.add(2, valid=False)
        before = hashlib.sha256(self.database.read_bytes()).digest()
        report = audit_transition_states(self.root, dry_run=True)
        self.assertEqual((report['passed'], report['invalid']), (1, 1))
        self.assertEqual(report['queued_for_recovery'], 0)
        self.assertEqual(hashlib.sha256(self.database.read_bytes()).digest(), before)
        self.assertFalse((self.root / 'ts_mode_audits').exists())
        self.assertFalse((self.root / 'ts_mode_audit.json').exists())

    def test_audit_queues_only_wrong_mode_and_is_idempotent(self):
        self.add(1)
        bad_output = self.add(2, valid=False, gibbs=-11.)
        self.add(3, strategy='neb', status='skipped')
        original = self.row(1)
        fingerprint = source_fingerprint(original, XYZ)
        report = audit_transition_states(self.root)
        self.assertEqual((report['passed'], report['invalid'], report['queued_for_recovery']), (1, 1, 1))
        self.assertTrue(Path(report['backup_file']).is_file())
        self.assertEqual(source_fingerprint(self.row(1), XYZ), fingerprint)
        rejected = self.row(2)
        self.assertEqual((rejected['status'], rejected['stage'], rejected['ts_mode_recovery']), ('pending', 'relax', 'relax'))
        self.assertEqual(rejected['attempt'], 1)
        self.assertIsNone(rejected['gibbs_free_energy_eh'])
        self.assertEqual(json.loads(rejected['ts_mode_json'])['previous_result']['gibbs_free_energy_eh'], -11.)
        self.assertTrue(bad_output.is_file())
        self.assertEqual(self.row(3)['status'], 'skipped')
        again = audit_transition_states(self.root)
        self.assertEqual((again['cached'], again['checked'], again['queued_for_recovery']), (1, 0, 0))
        self.assertIsNone(again['backup_file'])
        self.assertEqual(self.row(2), rejected)

    def test_skipped_neb_reactivates_only_when_group_has_no_valid_ts(self):
        self.add(1, valid=False)
        self.add(2, strategy='neb', status='skipped')
        report = audit_transition_states(self.root)
        self.assertEqual(report['neb_reactivated'], 1)
        self.assertEqual(self.row(2)['status'], 'pending')
        with closing(connect_database(self.database)) as db:
            candidate = claim_job(db)
            self.assertEqual(candidate['id'], 1)
            self.assertEqual(candidate['ts_mode_recovery'], 'relax')
            self.assertIsNone(claim_job(db))  # NEB waits for the recovering candidate.
            db.execute("UPDATE jobs SET status='failed' WHERE id=1")
            self.assertEqual(claim_job(db)['id'], 2)

    def test_changed_output_invalidates_previously_passed_cache(self):
        output = self.add(1)
        self.assertEqual(audit_transition_states(self.root)['passed'], 1)
        output.write_text(orca_output(displacement=0.001))
        report = audit_transition_states(self.root)
        self.assertEqual((report['cached'], report['invalid'], report['queued_for_recovery']), (0, 1, 1))

    def test_missing_normal_mode_output_is_reported_and_recoverable_from_seed(self):
        output = self.add(1)
        output.unlink()
        report = audit_transition_states(self.root)
        self.assertEqual((report['errors'], report['queued_for_recovery']), (1, 1))
        self.assertEqual(self.row(1)['ts_mode_status'], 'error')
        self.assertEqual(self.row(1)['status'], 'pending')

    def test_unrecoverable_wrong_mode_is_excluded_instead_of_accepted(self):
        self.add(1, valid=False)
        (self.root / 'seed_1.xyz').unlink()
        report = audit_transition_states(self.root)
        self.assertEqual((report['invalid'], report['queued_for_recovery'], report['unrecoverable']), (1, 0, 1))
        self.assertEqual(self.row(1)['status'], 'failed')


if __name__ == '__main__':
    unittest.main()
