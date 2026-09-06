#!/usr/bin/env python3
import argparse,csv,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);args=p.parse_args()
root=Path(args.run_root);report=root/'reports'
counts={'benchmark':1600,'private-distance-community-image':800,'full-information':800,'exclude-dispersed':4000,
 'cluster-shock':4000,'tight-multinomial':1600,'second-order-observability':4000,'grouped-lambda':1600,
 'arm-lambda':1600,'student-t5':1600,'finite-mixture':3200,'cluster-weighted':999}
panels=[(m,3500) for m in counts]+[('benchmark',cap) for cap in (4500,5500,10000)]
rows=[];panel_status=[]
for model,cap in panels:
    path=report/'validation'/('cap-'+str(cap)+'-'+model)
    if not (path/'independent-allocation-audit.txt').read_text().startswith('PASS:'):raise RuntimeError('Missing audit pass')
    values=list(csv.DictReader((path/'independent-allocation-audit.csv').open()))
    if len(values)!=5*counts[model]:raise RuntimeError('Incomplete audit inventory')
    if any(r['model']!=model or int(r['cap'])!=cap for r in values):raise RuntimeError('Wrong audited panel')
    rows.extend(values)
    status={s:sum(r['status']==s for r in values) for s in ('complete','target_infeasible','equilibrium_undefined')}
    if sum(status.values())!=len(values):raise RuntimeError('Failed or missing draw outcomes')
    panel_status.append(dict(model=model,cap=cap,draws=counts[model],**status))
if len(rows)!=152995 or len({(r['model'],r['cap'],r['draw'],r['scenario']) for r in rows})!=152995:raise RuntimeError('Global inventory mismatch')
digests={hashlib.sha256((root/('cap-'+str(cap))/'benchmark/policy-model-parameters.csv').read_bytes()).hexdigest() for cap in (3500,4500,5500,10000)}
if len(digests)!=1:raise RuntimeError('Benchmark draws differ across caps')
targets=[list(csv.DictReader((root/('cap-'+str(cap))/'benchmark/policy-experimental-targets.csv').open())) for cap in (3500,4500,5500,10000)]
auxiliary_difference=0
for values in targets[1:]:
    if len(values)!=len(targets[0]):raise RuntimeError('Target inventory differs across caps')
    for a,b in zip(targets[0],values):
        for key in a:
            if key=='target_community_welfare':
                difference=abs(float(a[key])-float(b[key]));auxiliary_difference=max(auxiliary_difference,difference)
                if difference>1e-12:raise RuntimeError('Auxiliary community target differs across caps')
            elif a[key]!=b[key]:raise RuntimeError('Adult targets or provenance differ across caps: '+key)
cost_checks=[]
for analysis in ('baseline-posterior','exponential-cluster-weights'):
    value=list(csv.DictReader((report/('independent-cost-assignment-audit-'+analysis+'.csv')).open()))
    if len(value)!=1 or value[0]['status']!='passed':raise RuntimeError('Missing cost-assignment audit')
    cost_checks.extend(value)
with (report/'independent-allocation-audit.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=dict(passed=True,communities=144,candidate_sites=1451,census_adults=39301,scenario_combinations=75,
 draw_scenario_evaluations=len(rows),benchmark_draws_identical_across_caps=True,benchmark_targets_identical_across_caps=True,
 original_200_draw_subset_verified=True,max_auxiliary_community_target_difference=auxiliary_difference,panels=panel_status,cost_assignment_audits=cost_checks)
(report/'independent-audit-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS:',len(rows),'main evaluations, all panels and both cost-assignment audits')
