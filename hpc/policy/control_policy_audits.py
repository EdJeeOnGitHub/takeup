#!/usr/bin/env python3
"""Audit each collected panel independently, with explicit full draw counts."""
import argparse, fcntl, importlib.util, json, shlex, time
from pathlib import Path
spec=importlib.util.spec_from_file_location('c',Path(__file__).with_name('control_policy_rerun.py'))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
p=argparse.ArgumentParser();p.add_argument('--production-ledger',required=True);p.add_argument('--ledger',required=True)
args=p.parse_args();path=Path(args.ledger)
with open(str(path)+'.lock','w') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    production=json.loads(Path(args.production_ledger).read_text())
    ledger=json.loads(path.read_text()) if path.exists() else dict(tasks=[dict(id=t['key'],prerequisite=t['id'],model=t['model'],cap=t['cap'],state='READY',attempts=[]) for t in production['tasks'] if t['stage']=='collect'],events=[])
    c.save(path,ledger)
    while True:
        production=json.loads(Path(args.production_ledger).read_text())
        if production.get('invalidated_reason'):raise RuntimeError('Production invalidated')
        c.reconcile_submissions(ledger);c.poll(ledger);c.save(path,ledger)
        states={t['id']:t['state'] for t in production['tasks']}
        for i,t in enumerate(ledger['tasks']):
            if t['state']!='READY' or states.get(t['prerequisite'])!='COMPLETE':continue
            host=list(c.HOSTS)[i%2];root=production['run_root'];token='adult-audit-'+t['model']+'-'+str(t['cap'])+'-'+str(int(time.time()))
            argv=['Rscript',c.BASE+'/scratch/policy-adult-population-20260906/check-policy-adult-outputs-final.R',
              '--run-root='+root,'--model='+t['model'],'--cap='+str(t['cap']),'--expected-combinations=1',
              '--output-path='+root+'/reports/validation/'+t['id'].replace('/','-'),
              '--distance-data='+c.INPUTS+'/full-many-pots-experiment-1451.rds',
              '--census-data='+c.INPUTS+'/takeup_census.RData',
              '--comparison-200-csv='+c.INPUTS+'/policy-original-200-parameters.csv']
            script='#!/usr/bin/env bash\nset -euo pipefail\nmodule load R/4.2.0\ncd '+shlex.quote(production['code_root'])+'\n'+ ' '.join(shlex.quote(a) for a in argv)+'\n'
            command='sbatch --parsable --partition='+c.HOSTS[host]+' --account=pi-akaring --cpus-per-task=1 --mem=28G --time=02:00:00 --job-name='+token+' --output='+shlex.quote(root+'/_logs/'+token+'-%j.log')
            t['state']='SUBMITTING';a=dict(host=host,token=token,submitted=time.time(),script=script,command=command);t['attempts'].append(a);c.save(path,ledger)
            job=c.ssh(host,command,script).splitlines()[-1].split(';')[0]
            if not job.isdigit():raise RuntimeError('Unknown sbatch response')
            a['job']=job;t['state']='SUBMITTED';c.save(path,ledger)
        print(json.dumps({s:sum(t['state']==s for t in ledger['tasks']) for s in sorted({t['state'] for t in ledger['tasks']})}),flush=True)
        if all(t['state']=='COMPLETE' for t in ledger['tasks']):break
        if any(t['state']=='FAILED' for t in ledger['tasks']):raise RuntimeError('Independent audit failed; inspect ledger/log')
        time.sleep(20)
