#!/usr/bin/env python3
"""Own a two-cluster policy DAG; migrate only terminally cancelled pending work.

A ledger is local to the controlling checkout. flock prevents concurrent drivers.
Every retry has a separate output root. Production requires an explicit smoke gate.
"""
import argparse, fcntl, hashlib, json, os, shlex, subprocess, time
from pathlib import Path

MODELS = {'benchmark':1600, 'private-distance-community-image':800,
          'full-information':800, 'exclude-dispersed':4000, 'cluster-shock':4000,
          'tight-multinomial':1600, 'second-order-observability':4000,
          'grouped-lambda':1600, 'arm-lambda':1600, 'student-t5':1600, 'finite-mixture':3200}
HOSTS = {'midway': 'broadwl', 'midway3': 'caslake'}
TERMINAL = {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE'}
BASE = '/project/akaring/takeup-data'
INPUTS = BASE + '/scratch/policy-adult-population-20260906/inputs'

def ssh(host, command, stdin=None):
    return subprocess.run(['ssh','-F',os.path.expanduser('~/.ssh/config'),'-o','BatchMode=yes',
                           '-o','ConnectTimeout=15',host,command], input=stdin, text=True,
                          capture_output=True, timeout=90, check=True).stdout.strip()

def save(path, ledger):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(ledger, indent=2)+'\n'); temp.replace(path)

def create(args):
    tasks=[]
    combinations=[(m,3500) for m in MODELS]+[('benchmark',c) for c in ([10000] if args.smoke else [4500,5500,10000])]+[('cluster-weighted',3500)]
    for j,(model,cap) in enumerate(combinations):
        key=f'cap-{cap}/{model}'; count=2 if args.smoke else MODELS.get(model,999)
        prep=f'{key}/prepare'; pred=f'{key}/predict'
        preferred=list(HOSTS)[j%2]
        def add(id,stage,deps,**kw):
            tasks.append(dict(id=id,stage=stage,deps=deps,key=key,model=model,cap=cap,
                              state='READY',preferred=preferred,attempts=[],**kw))
        add(prep,'prepare',[]);add(pred,'predict',[prep]); opt=[]
        for scenario in range(1,6):
            for first in range(1,count+1,args.shard_size):
                end=min(count,first+args.shard_size-1)
                id=f'{key}/scenario-{scenario}/{first:05d}-{end:05d}'
                add(id,'optimize',[pred],scenario=scenario,first=first,end=end)
                opt.append(id)
        collect=f'{key}/collect';add(collect,'collect',opt)
        add(f'{key}/summarize','summarize',[collect])
    return dict(version=1,created=time.time(),run_root=args.run_root,code_root=args.code_root,
                smoke=args.smoke,revision=args.revision,tasks=tasks,events=[],max_active=args.max_active)

def event(ledger, **kw):
    ledger['events'].append(dict(time=time.time(),**kw))

def poll(ledger):
    # sacct supplies terminal evidence; a missing row or SSH failure never means failed.
    for host in HOSTS:
        active=[t for t in ledger['tasks'] if t['state'] in ('SUBMITTED','CANCELLING') and t['attempts'][-1]['host']==host]
        if not active:continue
        ids=','.join(t['attempts'][-1]['job'] for t in active)
        try:
            rows=ssh(host, 'sacct -X -n -P -j '+shlex.quote(ids)+' --format=JobIDRaw,State,ExitCode,ElapsedRaw').splitlines()
            states={}
            for row in rows:
                fields=row.split('|')
                if len(fields)>=4:states[fields[0]]=(fields[1].split()[0].rstrip('+'),fields[2],fields[3])
            for task in active:
                attempt=task['attempts'][-1];record=states.get(attempt['job'])
                if not record:continue
                state,exit_code,elapsed=record
                attempt.update(slurm_state=state,exit_code=exit_code,elapsed=elapsed,last_seen=time.time())
                if state in TERMINAL:
                    if state=='COMPLETED' and exit_code=='0:0':task['state']='COMPLETE'
                    elif task['state']=='CANCELLING' and state=='CANCELLED':
                        task['state']='READY'; task['preferred']=next(h for h in HOSTS if h!=host)
                        task['migration_destination']=task['preferred']
                    else:task['state']='FAILED'
                    event(ledger,task=task['id'],job=attempt['job'],state=task['state'])
        except (subprocess.SubprocessError,ValueError) as error:
            event(ledger,host=host,observation_error=str(error))

def submit(ledger, task, host, path):
    attempt_id=len(task['attempts'])+1
    code_root=task.get('code_root',ledger['code_root'])
    allocation=f"{ledger['run_root']}/_shards/{task['key']}/scenario-{task.get('scenario',0)}/{task.get('first',0):05d}-{task.get('end',0):05d}/attempt-{attempt_id}"
    env=dict(REPO_ROOT=code_root,PROJECT_ROOT=code_root,
      MODEL_ID=task['model'],STAGE=task['stage'],OUTPUT_PATH=f"{ledger['run_root']}/{task['key']}",
      DISTANCE_DATA=INPUTS+'/full-many-pots-experiment-1451.rds',DISTANCE_CAP=str(task['cap']),
      CENSUS_DATA=INPUTS+'/takeup_census.RData',CLUSTER_WORKSPACE=INPUTS+'/dist_fit106.RData',
      BENCHMARK_FIT_PATH=INPUTS+'/assigned',BENCHMARK_DRAWS_PER_CHAIN='400',
      WEIGHTED_PATH=BASE+'/candidate-hpc-cd5f295-assigned/work/cluster-weight/modes',
      POLICY_GUROBI_ROOT=BASE+'/scratch/policy-adult-population-20260906/runtime/gurobi952',
      POPULATION_WEIGHTING='adult-census',TARGET_MODE='draw-specific-experimental-control',
      NUM_CORES='8',OPTIMIZE_CORES='8',DRAW_BATCH_SIZE='4',SOLVER_THREADS='1',SOLVER_SEED='0',
      POLICY_SOLVER='gurobi',MAX_DRAWS='2' if ledger['smoke'] else '0',
      POLICY_CODE_REVISION=task.get('revision',ledger['revision']),SOLVER_TIME_LIMIT='300',
      SHARD_ROOT=f"{ledger['run_root']}/_shards/{task['key']}")
    if task['stage']=='optimize':env.update(ALLOCATION_ROOT=allocation,SLURM_ARRAY_TASK_ID=str(task['scenario']),DRAW_START=str(task['first']),DRAW_END=str(task['end']))
    script='#!/usr/bin/env bash\nset -euo pipefail\n'+''.join(f'export {k}={shlex.quote(v)}\n' for k,v in env.items())
    script+='exec bash '+shlex.quote(code_root+'/hpc/policy/slurm_policy_model_robustness.sh')+'\n'
    token="adult-"+hashlib.sha256((ledger["run_root"]+task["id"]).encode()).hexdigest()[:10]+f"-{attempt_id}"
    log=f"{ledger['run_root']}/_logs/{token}-%j.log"
    command='mkdir -p '+shlex.quote(ledger['run_root']+'/_logs')+'; sbatch --parsable --partition='+HOSTS[host]+' --account=pi-akaring --cpus-per-task=8 --mem=28G --time=02:00:00 --job-name='+token+' --output='+shlex.quote(log)
    # Persist submission intent first. If SSH response is lost, reconcile by job name;
    # never blindly resubmit a task whose job ID is unknown.
    task['state']='SUBMITTING'
    attempt=dict(host=host,token=token,submitted=time.time(),command=command,script=script,allocation_root=allocation)
    task['attempts'].append(attempt);save(path,ledger)
    try:
        job=ssh(host,command,script).splitlines()[-1].split(';')[0]
        if not job.isdigit():raise ValueError('Invalid sbatch response: '+job)
        attempt['job']=job;task['state']='SUBMITTED';task.pop('migration_destination',None);event(ledger,task=task['id'],job=job,host=host,state='SUBMITTED')
    finally:save(path,ledger)

def reconcile_submissions(ledger):
    for task in ledger['tasks']:
        if task['state']!='SUBMITTING':continue
        a=task['attempts'][-1]
        # Search live and historical IDs by the unique recorded job name.
        rows=ssh(a['host'],'sacct -X -n -P -S '+time.strftime('%Y-%m-%d',time.localtime(a['submitted']))+' --name='+shlex.quote(a['token'])+' --format=JobIDRaw,JobName%100').splitlines()
        ids={r.split('|')[0] for r in rows if len(r.split('|'))>=2 and r.split('|')[1]==a['token']}
        if len(ids)==1:a['job']=ids.pop();task['state']='SUBMITTED';task.pop('migration_destination',None)
        elif len(ids)>1:raise RuntimeError('Ambiguous submission ownership: '+task['id'])
        # No observation is not proof that submission failed. Leave it pending.

def step(ledger,path,args):
    control=path.with_suffix('.control.json')
    if control.exists():
        maximum=int(json.loads(control.read_text())['max_active'])
        if maximum<1:raise ValueError('max_active must be positive')
        ledger['max_active']=maximum
    reconcile_submissions(ledger);poll(ledger);save(path,ledger)
    byid={t['id']:t for t in ledger['tasks']}
    active={h:sum(t['state'] in ('SUBMITTED','SUBMITTING','CANCELLING') and t['attempts'][-1]['host']==h for t in ledger['tasks']) for h in HOSTS}
    # Only migrate a job after rechecking it is PENDING; --state prevents cancellation
    # of a job that starts between observation and cancellation.
    for t in ledger['tasks']:
        if t['state']!='SUBMITTED':continue
        a=t['attempts'][-1];other=next(h for h in HOSTS if h!=a['host'])
        if (a.get('slurm_state')=='PENDING' and len(t['attempts']) < 3 and time.time()-a['submitted']>args.migrate_after and active[other]<ledger['max_active'] and any(x['attempts'] and x['attempts'][-1]['host']==other and x['attempts'][-1].get('slurm_state')=='RUNNING' for x in ledger['tasks'])):
            current=ssh(a['host'],'squeue -h -j '+a['job']+' -o "%T"').strip()
            if current=='PENDING':
                t['state']='CANCELLING';save(path,ledger)
                ssh(a['host'],'scancel --state=PENDING '+a['job']);event(ledger,task=t['id'],state='CANCELLING',destination=other)
    for t in ledger['tasks']:
        if t['state']!='READY' or not all(byid[d]['state']=='COMPLETE' for d in t['deps']):continue
        choices=[t['migration_destination']] if t.get('migration_destination') else sorted(HOSTS,key=lambda h:(active[h],h!=t['preferred']))
        host=choices[0]
        if active[host]>=ledger['max_active']:continue
        submit(ledger,t,host,path);active[host]+=1
    save(path,ledger)

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','step','watch','status','tune'])
    p.add_argument('--ledger',required=True);p.add_argument('--run-root');p.add_argument('--code-root')
    p.add_argument('--revision',default='uncommitted-snapshot');p.add_argument('--smoke',action='store_true')
    p.add_argument('--smoke-gate');p.add_argument('--shard-size',type=int,default=400)
    p.add_argument('--max-active',type=int,default=8);p.add_argument('--migrate-after',type=int,default=300)
    p.add_argument('--interval',type=int,default=60)
    args=p.parse_args();path=Path(args.ledger);path.parent.mkdir(parents=True,exist_ok=True)
    if args.action=='tune':
        if not path.exists() or args.max_active<1:raise ValueError('Existing ledger and positive maximum required')
        save(path.with_suffix('.control.json'),dict(max_active=args.max_active))
        return
    if args.action=='status':
        ledger=json.loads(path.read_text())
        counts={s:sum(t['state']==s for t in ledger['tasks']) for s in sorted({t['state'] for t in ledger['tasks']})}
        print(json.dumps(dict(time=time.time(),counts=counts,invalidated_reason=ledger.get('invalidated_reason'))))
        return
    with open(str(path)+'.lock','w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.action=='init':
            if path.exists():raise RuntimeError('Ledger already exists')
            if not args.run_root or not args.code_root or args.shard_size<1 or args.max_active<1:raise ValueError('Missing/invalid run configuration')
            if not args.smoke:
                if not args.smoke_gate:raise ValueError('Production requires --smoke-gate')
                gate=json.loads(Path(args.smoke_gate).read_text())
                if not gate.get('passed') or gate.get('revision')!=args.revision:raise ValueError('Smoke gate missing or does not match code revision')
            save(path,create(args));return
        ledger=json.loads(path.read_text())
        if ledger.get('invalidated_reason'):raise RuntimeError('Run invalidated: '+ledger['invalidated_reason'])
        while True:
            if args.action!='status':step(ledger,path,args)
            counts={s:sum(t['state']==s for t in ledger['tasks']) for s in sorted({t['state'] for t in ledger['tasks']})}
            print(json.dumps(dict(time=time.time(),counts=counts)),flush=True)
            if args.action!='watch' or all(t['state']=='COMPLETE' for t in ledger['tasks']):break
            time.sleep(args.interval)
if __name__=='__main__':main()
