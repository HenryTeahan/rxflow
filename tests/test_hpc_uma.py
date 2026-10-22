"""Portable method runs, exact plotted seeds, GPU submission and isolated plots."""
from contextlib import closing
import copy
import json
import math
from pathlib import Path
import shutil
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from hpc_precursors import has_plotted_gibbs, load_diagram_inputs
from hpc_uma import (calculation_settings, load_manifest, load_uma_diagram_inputs, plan_uma,
                     prepare_uma, selected_minima, status_uma, submit_uma, worker_uma)
from hpc_workflow import claim_job, connect_database, run_job
from orca_hpc import orca_input

XYZ = '3\nselected conformer\nO 0 0 0\nH 0.9 0 0\nH 0 0.9 0\n'
SMILES = '[O:1]([H:2])[H:3]'


class UmaWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name).resolve()
        self.source = self.base / 'source_results'
        self.source.mkdir()
        with closing(connect_database(self.source / 'workflow.sqlite')):
            pass
        self.before = (self.source / 'workflow.sqlite').read_bytes()
        self.config = dict(paths=dict(results='source_results', source='seeds', references='source_references'),
                           calculation=dict(cpus=1, memory_mb=2000),
                           sampling=dict(ts_conformers=2),
                           uma=dict(model_label='UMA test model', external_params=['--device', '${UMA_DEVICE}'],
                                    gpus_per_task=1))
        self.config_path = self.base / 'run.toml'
        self.config_path.write_text('[paths]\nresults="source_results"\n')
        self.data = dict(jobs=[], encounter_jobs=[], precursor_jobs=[], bound_reactant_jobs=[],
                         reference_jobs=[], diagram_references=[], energy_model='PBE0 corrected')
        for job_id, endpoint, gibbs in ((1, 'reactant', -10.), (2, 'TS', -9.9), (3, 'product', -10.01),
                                         (4, 'reactant', -9.), (5, 'encounter', -10.05)):
            filename = f'selected_{job_id}.xyz'
            (self.source / filename).write_text(XYZ)
            row = dict(id=job_id, ligand_id='5', pair_index=1, seating=1, endpoint=endpoint,
                       conformer=job_id, status='succeeded', stage='complete', attempt=1,
                       input_xyz=filename, xyz_file=filename, smiles=SMILES, charge=0, multiplicity=1,
                       substrate='O', binding_permutation=1, binding_mode=1, coordination_isomer=1,
                       donor_maps='[1]', substrate_site=1, resonance=0, donor=0, hydrogen=1, acceptor=2,
                       gibbs_free_energy_eh=gibbs, electronic_energy_eh=gibbs-.01,
                       singlepoint_status='succeeded', singlepoint_job_id=100+job_id,
                       settings_json='{}', ts_guess_json=None, ts_parent_id=None)
            self.data['encounter_jobs' if endpoint=='encounter' else 'jobs'].append(row)

    def prepare(self):
        with patch('hpc_precursors.load_diagram_inputs', return_value=copy.deepcopy(self.data)), \
             patch('ts_scan.scan_inputs', side_effect=ValueError('no archived scan')):
            return prepare_uma(self.config, self.base, self.source)

    def test_plan_selects_one_corrected_minimum_per_state_and_never_writes_sources(self):
        with patch('hpc_precursors.load_diagram_inputs', return_value=self.data), \
             patch('ts_scan.scan_inputs', side_effect=ValueError('no scan')):
            result = plan_uma(self.config, self.base, self.source)
        self.assertEqual(result['counts'], dict(reaction=3, encounter=1))
        self.assertEqual({r['source_job_id'] for r in result['records']}, {1,2,3,5})
        self.assertFalse((self.source / 'uma').exists())
        self.assertEqual((self.source / 'workflow.sqlite').read_bytes(), self.before)
        bad = dict(self.data['jobs'][0], id=88, gibbs_free_energy_eh=-100., singlepoint_status='pending')
        self.assertNotIn(bad, selected_minima([*self.data['jobs'],bad]))

    def test_prepared_run_survives_transfer_without_source_files_or_databases(self):
        result = self.prepare()
        root = Path(result['root'])
        self.assertEqual((self.source / 'workflow.sqlite').read_bytes(), self.before)
        destination = self.base / 'other_cluster' / 'run'
        shutil.copytree(root, destination)
        shutil.rmtree(self.source)
        portable = __import__('tomllib').loads((destination/'run_uma.toml').read_text())
        self.assertEqual(portable['paths']['results'], '.')
        status = status_uma(destination)
        self.assertEqual(sum(row['jobs'] for row in status['counts']),5)  # 4 seeds + NEB
        self.assertTrue(all(not Path(r['input_xyz']).is_absolute() for r in load_manifest(destination)['sources']))
        imported = load_uma_diagram_inputs(self.base/'absent/workflow.sqlite',destination)
        self.assertEqual(imported['singlepoint_summary']['active_mode'],'uma')
        self.assertEqual(imported['singlepoint_summary']['ready'],0)
        self.assertTrue(all(not has_plotted_gibbs(row) for row in imported['jobs']))
        with closing(sqlite3.connect(destination/'workflow.sqlite')) as db, db:
            db.row_factory=sqlite3.Row
            ts=dict(db.execute("SELECT * FROM jobs WHERE ts_strategy='direct'").fetchone())
            meta=json.loads(ts['ts_guess_json'])
            self.assertTrue((destination/meta['reactant_xyz']).is_file())
            self.assertTrue((destination/meta['product_xyz']).is_file())
            self.assertFalse(json.loads(ts['settings_json'])['ts_scan_recovery'])
        again=prepare_uma(portable,destination,destination)
        self.assertTrue(again['existing'])

    def test_preparation_is_idempotent_and_protocol_or_seed_changes_are_rejected(self):
        root=Path(self.prepare()['root'])
        self.assertTrue(prepare_uma(self.config,self.base,self.source)['existing'])
        other=copy.deepcopy(self.config)
        other['uma']['model_label']='another model'
        with self.assertRaisesRegex(ValueError,'differ'):
            prepare_uma(other,self.base,self.source)
        record=load_manifest(root)['sources'][0]
        (root/record['input_xyz']).write_text(XYZ.replace('0.9','0.8'))
        with self.assertRaisesRegex(ValueError,'seed changed'):
            prepare_uma(self.config,self.base,self.source)

    def test_uma_loader_uses_own_gibbs_geometries_and_frequencies_with_no_dft_fallback(self):
        root=Path(self.prepare()['root'])
        with closing(sqlite3.connect(root/'workflow.sqlite')) as db, db:
            db.execute("UPDATE jobs SET status='succeeded',stage='complete',gibbs_free_energy_eh=-20.,electronic_energy_eh=-20.1,xyz_file=input_xyz,frequencies_json='[-500,1,2]' WHERE id=2")
        inputs=load_diagram_inputs(self.source/'workflow.sqlite',energy_mode='uma')
        ts=next(row for row in inputs['jobs'] if row['id']==2)
        self.assertTrue(has_plotted_gibbs(ts))
        self.assertEqual(ts['gibbs_free_energy_eh'],-20.)
        self.assertEqual(ts['frequencies_json'],'[-500,1,2]')
        self.assertEqual(inputs['results_database'],str(root/'workflow.sqlite'))
        self.assertFalse(any(has_plotted_gibbs(r) for r in inputs['jobs'] if r['endpoint']=='product'))
        with closing(sqlite3.connect(root/'workflow.sqlite')) as db, db:
            db.execute("UPDATE jobs SET settings_json='{}' WHERE id=2")
        with self.assertRaises((ValueError,KeyError)):
            load_diagram_inputs(self.source/'workflow.sqlite',energy_mode='uma')

    def test_portable_cli_needs_no_source_config_sections(self):
        import contextlib
        import io
        from hpc_workflow import main
        root=Path(self.prepare()['root'])
        with contextlib.redirect_stdout(io.StringIO()) as output:
            main(['--config',str(root/'run_uma.toml'),'uma','status'])
        self.assertIn('UMA test model',output.getvalue())
        with contextlib.redirect_stdout(io.StringIO()) as output:
            main(['--config',str(root/'run_uma.toml'),'uma','submit','--dry-run'])
        self.assertIn('--gpus-per-task=1',output.getvalue())

    def test_submission_is_portable_and_requests_cpu_or_gpu_explicitly(self):
        root=Path(self.prepare()['root'])
        for gpus in (0,1):
            config=copy.deepcopy(self.config); config['uma']['gpus_per_task']=gpus
            with patch('hpc_uma.subprocess.check_output') as submit:
                report=submit_uma(self.config_path,config,self.base,root,workers=3,
                                   partition='gpuq',account='mygroup',dry_run=True)
                submit.assert_not_called()
            command=report['command']
            self.assertEqual(any(a.startswith('--gpus-per-task=') for a in command),bool(gpus))
            self.assertIn('--partition=gpuq',command)
            self.assertIn('--array=0-2%3',command)
            self.assertIn('uma worker',command[-1])
        with patch('hpc_uma.subprocess.check_output',return_value='123;cluster\n') as submit:
            report=submit_uma(self.config_path,self.config,self.base,root,workers=2)
        self.assertEqual(report['workers'],'123')
        self.assertTrue((root/'logs').is_dir())

    def test_both_notebook_selections_admit_only_returned_uma_states(self):
        from collections import defaultdict
        from hpc_precursors import canonical_substrate
        root=Path(self.prepare()['root'])
        for filename in ('workflow.sqlite','encounters.sqlite'):
            with closing(sqlite3.connect(root/filename)) as db, db:
                db.execute("UPDATE jobs SET status='succeeded',stage='complete',gibbs_free_energy_eh=CASE endpoint WHEN 'TS' THEN -9.9 WHEN 'product' THEN -10.01 WHEN 'encounter' THEN -10.05 ELSE -10 END,electronic_energy_eh=-10.1,xyz_file=input_xyz WHERE coalesce(ts_strategy,'')!='neb'")
        data=load_diagram_inputs(self.source/'workflow.sqlite',energy_mode='uma')
        base=Path(__file__).resolve().parents[1]
        reaction=json.loads((base/'reaction_diagrams.ipynb').read_text())
        cell=next(c for c in reaction['cells'] if c.get('id')=='560d679e')
        namespace=dict(jobs=data['jobs'],encounter_jobs=data['encounter_jobs'],
                       math=math,defaultdict=defaultdict,canonical_substrate=canonical_substrate)
        exec(''.join(cell['source']),namespace)
        self.assertEqual(len(namespace['profiles']),1)
        self.assertEqual(namespace['profiles'][0]['lowest_start_gibbs_eh'],-10.05)
        comparison=json.loads((base/'system_energy_comparison.ipynb').read_text())
        scope=dict(inputs_by_system={'UMA':data},site_label_by_identity={},HARTREE_TO_KCAL=627.509474,
                   MINIMUM_TOLERANCE_KCAL=1e-6,defaultdict=defaultdict,canonical_substrate=canonical_substrate,math=math)
        exec(''.join(comparison['cells'][2]['source']),scope)
        self.assertEqual(len(scope['pathways']),1)
        self.assertAlmostEqual(scope['pathways'][0]['limiting_barrier'],.15*627.509474)

    def test_generic_uma_inputs_keep_existing_ts_stages_without_gxtb_flags(self):
        settings={**calculation_settings(self.config),'external_exe':'/opt/uma_wrapper'}
        for stage in ('opt','verytightopt','ts_preopt','optts','neb','freq','external_sp'):
            text=orca_input(XYZ,stage,0,1,settings,atoms=(0,1,2))
            self.assertIn('ExtOpt',text)
            self.assertIn('/opt/uma_wrapper',text)
            self.assertNotIn('PBE0',text)
            self.assertNotIn('-s h2o',text)
            self.assertNotIn('xtb',text.lower())
        self.assertIn('NumFreq',orca_input(XYZ,'freq',0,1,settings))

    def test_encounter_connectivity_change_is_rejected_by_shared_worker(self):
        root=Path(self.prepare()['root'])
        from hpc_uma import open_queue
        with closing(open_queue(root,'encounter')) as db:
            job=claim_job(db)
            optimized=dict(out_file='opt.out',xyz=XYZ,converged=True,max_cycles_reached=False)
            with patch('hpc_workflow.run_orca',return_value=optimized) as calculate, \
                 patch('hpc_workflow.connectivity_changes',return_value={'broken':[(0,1)]}):
                run_job(db,job,root,{})
            row=db.execute('SELECT status,error FROM jobs WHERE id=?',(job['id'],)).fetchone()
            self.assertEqual(row['status'],'failed')
            self.assertIn('Graph change',row['error'])
            self.assertEqual(calculate.call_count,1)

    def test_ts_uses_same_preopt_optts_numfreq_and_mode_validation(self):
        root=Path(self.prepare()['root'])
        with closing(connect_database(root/'workflow.sqlite')) as db:
            db.execute("UPDATE jobs SET status='succeeded' WHERE endpoint!='TS'")
            job=claim_job(db)
            good=dict(out_file='out.out',xyz=XYZ,converged=True,max_cycles_reached=False)
            freq=dict(out_file='freq.out',electronic_energy_eh=-20.,gibbs_free_energy_eh=-19.9,
                      imaginary_count=1,frequencies=[-500.,1.,2.])
            mode=dict(valid=True,fingerprint='tested',reason='valid transfer')
            with patch('hpc_workflow.run_orca',side_effect=[good,good,freq]) as calculate, \
                 patch('hpc_workflow.connectivity_changes',return_value={'broken':[(0,1)],'formed':[(1,2)]}), \
                 patch('hpc_workflow.check_ts_mode',return_value=mode):
                run_job(db,job,root,{})
            self.assertEqual([call.args[3] for call in calculate.call_args_list],['ts_preopt','optts','freq'])
            saved=db.execute('SELECT * FROM jobs WHERE id=?',(job['id'],)).fetchone()
            self.assertEqual(saved['status'],'succeeded')
            self.assertEqual(saved['ts_mode_status'],'passed')
            self.assertEqual(saved['gibbs_free_energy_eh'],-19.9)

    def test_ts_unrelated_connectivity_change_is_rejected_before_frequency(self):
        root=Path(self.prepare()['root'])
        with closing(connect_database(root/'workflow.sqlite')) as db:
            db.execute("UPDATE jobs SET status='succeeded' WHERE endpoint!='TS'")
            job=claim_job(db)
            good=dict(out_file='out.out',xyz=XYZ,converged=True,max_cycles_reached=False)
            with patch('hpc_workflow.run_orca',return_value=good) as calculate, \
                 patch('hpc_workflow.connectivity_changes',return_value={'broken':[(0,2)]}):
                run_job(db,job,root,{})
            self.assertEqual([call.args[3] for call in calculate.call_args_list],['ts_preopt','optts'])
            saved=db.execute('SELECT * FROM jobs WHERE id=?',(job['id'],)).fetchone()
            self.assertEqual(saved['status'],'failed')
            self.assertIn('outside the intended H transfer',saved['error'])


if __name__=='__main__':
    unittest.main()
