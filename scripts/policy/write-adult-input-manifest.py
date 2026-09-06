#!/usr/bin/env python3
"""Record actual source paths/content hashes and selected-draw inventories."""
import argparse, csv, hashlib
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);p.add_argument('--input-root',required=True)
p.add_argument('--code-root',action='append',default=[]);p.add_argument('--audit-script',action='append',default=[]);args=p.parse_args()
root=Path(args.run_root);inputs=Path(args.input_root);out=root/'reports';out.mkdir(exist_ok=True)
cache={};rows=[];draws=[]
def record(path,role,model='',cap=''):
    path=Path(path)
    if not path.is_file():raise FileNotFoundError(path)
    name=str(path.resolve())
    if name not in cache:
        h=hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
        cache[name]=(h.hexdigest(),path.stat().st_size)
    digest,size=cache[name]
    rows.append(dict(role=role,model=model,cap_m=cap,path=name,sha256=digest,bytes=size))

for path in sorted(root.glob('cap-*/*/policy-model-parameters.csv')):
    model=path.parent.name;cap=path.parent.parent.name.split('-')[1]
    records=list(csv.DictReader(path.open()))
    for name in ('policy-model-parameters.csv','policy-model-parameters.rds','policy-population.csv',
                 'policy-experimental-targets.csv','policy-cache-manifest.rds','policy-edge-demand-draw-map.csv'):
        record(path.parent/name,'derived-input',model,cap)
    sources=set()
    for row in records:
        source=row.get('source_csv') or row.get('mode_csv')
        if not source or source=='NA':source=row.get('mode_csv')
        if not source or source=='NA':raise ValueError('Missing source provenance for '+str(path))
        sources.add(source)
        draws.append(dict(model=model,cap_m=cap,draw=row['draw'],replicate=row['replicate'],
          chain=row.get('chain',''),iteration=row.get('iteration',''),comparison_200=row.get('comparison_200',''),source_path=source))
    for source in sorted(sources):record(source,'structural-fit-or-mode',model,cap)
for name in ('takeup_census.RData','dist_fit106.RData','full-many-pots-experiment-1451.rds','policy-original-200-parameters.csv'):
    record(inputs/name,'shared-input')
record('/project/akaring/takeup-data/data/stan_analysis_data/dist_fit104.RData','household-geography')
for code in args.code_root:record(Path(code)/'SHA256SUMS','code-content-manifest')
for script in args.audit_script:record(script,'validation-code')
record(inputs.parent/'runtime/gurobi952/bin/gurobi_cl','solver-executable')
for name,data in (('input-source-manifest.csv',rows),('selected-draw-manifest.csv',draws)):
    with (out/name).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
print('Hashed',len(cache),'unique inputs; recorded',len(draws),'model/cap draw selections')
