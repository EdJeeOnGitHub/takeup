from pathlib import Path
import csv,re,json,math,shutil,hashlib
root=Path('/home/ed/projects/takeup'); out=root/'ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun'
new=root/'build/expected-distance-full-audit-20260929/work/realized';old=root/'build/expected-distance-historical-audit-20260929/work/realized';paper=Path('/home/ed/projects/overleaf/overleaf-takeup'); saved=root/'build/work/realized/temp-data/tidy-rf-tes'; prejoin=root/'build/expected-distance-prejoin-audit-20260929/work/realized'
rows=[];summary=[]
for nf in sorted((new/'temp-data/tidy-rf-tes').glob('*.csv')) + [new/'temp-data/temp.csv']:
 of=old/nf.relative_to(new)
 if not of.exists():continue
 nr=list(csv.DictReader(nf.open()));orr=list(csv.DictReader(of.open()))
 if not nr or not {'assigned_treatment','assigned_dist_group','estimate','std_error'}.issubset(nr[0]):continue
 keyed={(r['assigned_treatment'],r['assigned_dist_group']):r for r in orr}
 saved_file=saved/nf.name if nf.parent.name=='tidy-rf-tes' else root/'build/work/realized/temp-data'/nf.name
 saved_rows={(r['assigned_treatment'],r['assigned_dist_group']):r for r in csv.DictReader(saved_file.open())} if saved_file.exists() else {}
 model_rows=[]
 for r in nr:
  k=r['assigned_treatment'],r['assigned_dist_group'];o=keyed.get(k)
  if o is None:continue
  try:e,se,oe,ose=[float(v) for v in [r['estimate'],r['std_error'],o['estimate'],o['std_error']]]
  except (ValueError,TypeError):continue
  p=math.erfc(abs(e/se)/math.sqrt(2));op=math.erfc(abs(oe/ose)/math.sqrt(2))
  changes='; '.join(f'{t:g}: '+('gains' if p<t else 'loses') for t in [.01,.05,.1] if (p<t)!=(op<t))
  rec=dict(file=nf.name,arm=k[0],group=k[1],historical_estimate=oe,corrected_estimate=e,change_pp=100*(e-oe),historical_se=ose,corrected_se=se,historical_p=op,corrected_p=p,threshold_changes=changes)
  sr=saved_rows.get(k)
  rec.update(saved_500_estimate='',saved_500_se='',saved_500_p='',simulation_change_pp='',join_change_pp='')
  if sr:
   try:
    se500=float(sr['estimate']); sd500=float(sr['std_error'])
    rec.update(saved_500_estimate=se500,saved_500_se=sd500,saved_500_p=math.erfc(abs(se500/sd500)/math.sqrt(2)),simulation_change_pp=100*(se500-oe),join_change_pp=100*(e-se500))
   except (ValueError,TypeError):pass
  model_rows.append(rec);rows.append(rec)
 if model_rows:summary.append(dict(file=nf.name,contrasts=len(model_rows),max_abs_change_pp=max(abs(r['change_pp']) for r in model_rows),threshold_crossings=sum(bool(r['threshold_changes']) for r in model_rows)))
def write(name,records):
 if records:
  with (out/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
write('all-bootstrap-contrasts.csv',rows);write('model-summary.csv',summary)
write('significance-crossings.csv',[r for r in rows if r['threshold_changes']])
# Compare numeric cells while ignoring stars, captions, and purely cosmetic LaTeX.
def numeric_rows(path):
 result=[]
 for line in path.read_text().splitlines():
  if '&' not in line or any(x in line for x in ['multicolumn','cmidrule','Dependent variable','\\caption','\\label']):continue
  parts=line.split('&');label=re.sub(r'[^a-z0-9]', '', parts[0].strip().replace('Control mean','Control').lower())
  cells=[]
  for cell in parts[1:]:
   cell=cell.replace('{,}',',').replace('$<$','<').replace('\\textless{}','<')
   vals=re.findall(r'(?<![A-Za-z])[-+]?(?:\d[\d,]*\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',cell)
   cells.append(tuple(float(x.replace(',','')) for x in vals))
  result.append((label,cells))
 return result
included=json.loads((out/'included-rf-tables.json').read_text())
for fn in ['prior-deworming-robustness.tex','fob-baseline-imbalance-robustness.tex']:
 included.append(dict(source='online-appendix.tex',table='tables/'+fn,line=0))
coverage=[]
for item in included:
 name=Path(item['table']).name;nf=new/'presentations/rf-tables/main-specs'/name;of=old/'presentations/rf-tables/main-specs'/name;pf=paper/item['table']
 rec=dict(table=item['table'],corrected_generated=nf.exists(),historical_generated=of.exists(),historical_matches_paper_numeric_rows='',corrected_matches_paper_numeric_rows='',prejoin_500_matches_paper_numeric_rows='')
 if nf.exists() and of.exists():
  rec['historical_matches_paper_numeric_rows']=numeric_rows(pf)==numeric_rows(of)
  rec['corrected_matches_paper_numeric_rows']=numeric_rows(pf)==numeric_rows(nf)
  for version,src in [('corrected',nf),('historical',of),('manuscript',pf)]:
   dest=out/version/item['table'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
 pre=prejoin/'presentations/rf-tables/main-specs'/name
 if not pre.exists():pre=root/'build/realized/presentations/rf-tables/main-specs'/name
 if pre.exists():
  rec['prejoin_500_matches_paper_numeric_rows']=numeric_rows(pf)==numeric_rows(pre)
  dest=out/'prejoin-500'/item['table'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(pre,dest)
 coverage.append(rec)
write('included-table-coverage.csv',coverage)
hashes=json.loads((out/'manuscript-before-sha256.json').read_text());changed=[name for name,h in hashes.items() if hashlib.sha256((paper/name).read_bytes()).hexdigest()!=h]
(out/'manuscript-unchanged.json').write_text(json.dumps(dict(checked=len(hashes),changed=changed),indent=2)+'\n'); audited_changed=[name for name in changed if name in ['ECM ReStud.tex','experimental-design.tex','online-appendix.tex'] or name.startswith('rf-tables/')]; assert not audited_changed,audited_changed
print(f'{len(summary)} paired bootstrap outputs; {len(rows)} contrasts; {sum(bool(r["threshold_changes"]) for r in rows)} threshold crossings. Audited manuscript files unchanged; unrelated changes: {changed}.')
