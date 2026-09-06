import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('controller', Path(__file__).resolve().parents[2] / 'hpc/policy/control_policy_rerun.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

def task(host='midway', state='SUBMITTED', job='1'):
    return dict(id=job, state=state, deps=[], preferred=host,
                attempts=[dict(host=host, job=job, token='owned-'+job, submitted=0)])

class ControllerTests(unittest.TestCase):
    def test_missing_observation_is_not_failure(self):
        t=task(); ledger=dict(tasks=[t],events=[])
        with patch.object(c,'ssh',return_value=''): c.poll(ledger)
        self.assertEqual(t['state'],'SUBMITTED')

    def test_failed_exit_is_not_success(self):
        t=task(); ledger=dict(tasks=[t],events=[])
        with patch.object(c,'ssh',return_value='1|COMPLETED|1:0|12'): c.poll(ledger)
        self.assertEqual(t['state'],'FAILED')

    def test_migration_requires_confirmed_cancellation(self):
        t=task(state='CANCELLING'); ledger=dict(tasks=[t],events=[])
        with patch.object(c,'ssh',return_value='1|RUNNING|0:0|12'): c.poll(ledger)
        self.assertEqual(t['state'],'CANCELLING')
        with patch.object(c,'ssh',return_value='1|CANCELLED|0:0|12'): c.poll(ledger)
        self.assertEqual(t['state'],'READY')
        self.assertEqual(t['preferred'],'midway3')
        self.assertEqual(t['migration_destination'],'midway3')

    def test_migrated_job_cannot_return_to_emptier_original_queue(self):
        first=task(state='READY');first['migration_destination']='midway3'
        other=task('midway3',job='2');other['attempts'][0]['slurm_state']='RUNNING'
        ledger=dict(tasks=[first,other],events=[],max_active=2)
        with tempfile.TemporaryDirectory() as d, patch.object(c,'poll'), patch.object(c,'submit') as submit:
            c.step(ledger,Path(d)/'ledger.json',SimpleNamespace(migrate_after=1))
        self.assertEqual(submit.call_args.args[2],'midway3')

    def test_job_starting_during_cancel_check_is_not_cancelled(self):
        first=task(); first['attempts'][0]['slurm_state']='PENDING'
        other=task('midway3',job='2'); other['attempts'][0]['slurm_state']='RUNNING'
        ledger=dict(tasks=[first,other],events=[],max_active=2)
        with tempfile.TemporaryDirectory() as d, patch.object(c,'poll'), patch.object(c,'ssh',return_value='RUNNING') as ssh:
            c.step(ledger,Path(d)/'ledger.json',SimpleNamespace(migrate_after=1))
            self.assertEqual(ssh.call_count,1)
            self.assertTrue(ssh.call_args.args[1].startswith('squeue'))
        self.assertEqual(first['state'],'SUBMITTED')

    def test_lost_submission_is_reconciled_by_owned_name(self):
        t=task(state='SUBMITTING'); ledger=dict(tasks=[t],events=[])
        with patch.object(c,'ssh',return_value=''): c.reconcile_submissions(ledger)
        self.assertEqual(t['state'],'SUBMITTING')
        with patch.object(c,'ssh',return_value='77|owned-1'): c.reconcile_submissions(ledger)
        self.assertEqual(t['attempts'][-1]['job'],'77')
        self.assertEqual(t['state'],'SUBMITTED')

if __name__=='__main__': unittest.main()
