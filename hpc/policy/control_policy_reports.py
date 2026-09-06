#!/usr/bin/env python3
"""Schedule reporting after recorded production prerequisites complete."""
import argparse, fcntl, importlib.util, json, shlex, time
from pathlib import Path

spec=importlib.util.spec_from_file_location('main_controller',Path(__file__).with_name('control_policy_rerun.py'))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--production-ledger',required=True);p.add_argument('--ledger',required=True)
    p.add_argument('--code-root',required=True);p.add_argument('--interval',type=int,default=30)
    args=p.parse_args();path=Path(args.ledger)
    with open(str(path)+'.lock','w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if path.exists():ledger=json.loads(path.read_text())
        else:
            ledger=dict(tasks=[dict(id=s,stage=s,state='READY',attempts=[]) for s in ('median','cost-baseline','cost-bootstrap','assemble')],events=[])
            c.save(path,ledger)
        while True:
            production=json.loads(Path(args.production_ledger).read_text())
            if production.get('invalidated_reason'):raise RuntimeError('Production invalidated')
            c.reconcile_submissions(ledger);c.poll(ledger);c.save(path,ledger)
            states={t['id']:t['state'] for t in production['tasks']}
            own={t['id']:t['state'] for t in ledger['tasks']}
            for t in ledger['tasks']:
                if t['state']!='READY':continue
                required={'median':['cap-3500/benchmark/prepare'],
                  'cost-baseline':['cap-3500/benchmark/collect'],
                  'cost-bootstrap':['cap-3500/cluster-weighted/collect'],
                  'assemble':[k for k in states if k.endswith('/summarize')]}[t['stage']]
                if not all(states.get(k)=='COMPLETE' for k in required):continue
                if t['stage']=='assemble' and not all(own[s]=='COMPLETE' for s in ('median','cost-baseline','cost-bootstrap')):continue
                host='midway3' if t['stage']=='cost-bootstrap' else 'midway'
                root=production['run_root'];token='adult-report-'+t['stage']+'-'+str(int(time.time()))
                env=dict(REPO_ROOT=args.code_root,RUN_ROOT=root,STAGE=t['stage'],
                  DISTANCE_DATA=c.INPUTS+'/full-many-pots-experiment-1451.rds',
                  POLICY_CENSUS=c.INPUTS+'/takeup_census.RData',
                  POLICY_GUROBI_ROOT=c.BASE+'/scratch/policy-adult-population-20260906/runtime/gurobi952',
                  POLICY_CODE_REVISION=production['revision'])
                script='#!/usr/bin/env bash\nset -euo pipefail\n'+''.join('export '+k+'='+shlex.quote(v)+'\n' for k,v in env.items())
                script+='exec bash '+shlex.quote(args.code_root+'/hpc/policy/slurm_policy_adult_reports.sh')+'\n'
                command='sbatch --parsable --partition='+c.HOSTS[host]+' --account=pi-akaring --cpus-per-task=8 --mem=28G --time=02:00:00 --job-name='+token+' --output='+shlex.quote(root+'/_logs/'+token+'-%j.log')
                t['state']='SUBMITTING';a=dict(host=host,token=token,submitted=time.time(),command=command,script=script)
                t['attempts'].append(a);c.save(path,ledger)
                job=c.ssh(host,command,script).splitlines()[-1].split(';')[0]
                if not job.isdigit():raise RuntimeError('Unknown sbatch response: '+job)
                a['job']=job;t['state']='SUBMITTED';c.save(path,ledger)
            print(json.dumps({t['id']:t['state'] for t in ledger['tasks']}),flush=True)
            if all(t['state']=='COMPLETE' for t in ledger['tasks']):break
            if any(t['state']=='FAILED' for t in ledger['tasks']):raise RuntimeError('Reporting job failed; inspect ledger and log')
            time.sleep(args.interval)

if __name__=='__main__':main()
