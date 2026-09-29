from pathlib import Path
import csv,json,difflib,math
out=Path('ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun')
coverage={Path(r['table']).name:r for r in csv.DictReader(open(out/'included-table-coverage.csv'))}
mapping={
'discrete-dist-covs-tidy-tes.csv':'rf_discrete_dist_covs_tbl.tex',
'reducedform-discrete-fob-tidy-tes.csv':'rf_discrete_fob_spec_tbl.tex',
'reducedform-dist-cts-tidy-tes.csv':'rf_dist_cts_spec_tbl.tex',
'reducedform-cts-fob-tidy-tes.csv':'rf_cts_fob_spec_tbl.tex',
'reducedform-discrete-sob-tidy-tes.csv':'rf_discrete_sob_spec_tbl.tex',
'reducedform-discrete-pct-yesno-tidy-tes.csv':'rf_discrete_belief_decomposition_tbl.tex',
'reducedform-discrete-pct-yesnodk-tidy-tes.csv':'rf_discrete_belief_decomposition_tbl.tex',
'externality-knowledge-tidy-tes.csv':'rf_externality_knowledge_tbl.tex',
'predicted-endline-deworm-takeup-tidy-tes.csv':'predicted_endline_deworm_takeup_spec_tbl.tex',
'preference-for-bracelet-tidy-tes.csv':'preference_for_bracelet_abs_diff_spec_tbl.tex',
'reducedform-discrete-fob-no-covs-no-mu-d-tidy-tes.csv':'rf_discrete_fob_no_covs_no_mu_d_spec_tbl.tex',
'discrete-dist-no-covs-no-mu-d-tidy-tes.csv':'rf_discrete_dist_no_covs_no_mu_d_tbl.tex'}
records=[]
for r in csv.DictReader(open(out/'all-bootstrap-contrasts.csv')):
 fn=r['file']
 if fn not in mapping:continue
 arms=['control','bracelet','calendar','ink','signal','bracelet - calendar']
 if fn.startswith('predicted'): arms=arms[:4]
 if fn.startswith('preference'):arms=arms[:4]+['abs(calendar) - abs(bracelet)']
 if r['arm'] not in arms:continue
 c=coverage[mapping[fn]]
 baseline='100-draw reconstruction' if c['historical_matches_paper_numeric_rows']=='True' else '500-draw prejoin' if c['prejoin_500_matches_paper_numeric_rows']=='True' else 'unresolved manuscript source'
 prefix='saved_500' if baseline=='500-draw prejoin' else 'historical'
 try:e,p=float(r[prefix+'_estimate']),float(r[prefix+'_p'])
 except ValueError:continue
 ne,np=float(r['corrected_estimate']),float(r['corrected_p'])
 changes='; '.join(f'{int(100*t)}%: '+('gains' if np<t else 'loses') for t in [.01,.05,.1] if (np<t)!=(p<t))
 records.append(dict(table=mapping[fn],model=fn,arm=r['arm'],group=r['group'],baseline=baseline,before_estimate=e,after_estimate=ne,before_p=p,after_p=np,threshold_changes=changes))
def write(name,rs):
 if rs:
  with (out/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
write('manuscript-bootstrap-contrasts.csv',records);write('manuscript-significance-crossings.csv',[r for r in records if r['threshold_changes']])
lines=['# Significance changes in included bootstrap tables','', 'Comparisons use the reconstruction that matches each manuscript table’s displayed numerical cells. These are numerical significance thresholds, separate from correcting the star formatter. Effects are percentage points.','', '| Table / contrast | Before: effect; p | Corrected: effect; p | Threshold |','|---|---:|---:|---|']
for r in records:
 if not r['threshold_changes']:continue
 lines.append(f"| {r['table']} ({'Panel B' if 'pct-yesno-' in r['model'] else 'Panel C' if 'pct-yesnodk-' in r['model'] else 'main panel'}): {r['arm']}, {r['group']} | {100*r['before_estimate']:.2f}; {r['before_p']:.3f} | {100*r['after_estimate']:.2f}; {r['after_p']:.3f} | {r['threshold_changes']} |")
(out/'significance-review.md').write_text('\n'.join(lines)+'\n')
# Full, review-only diffs, including analytical attrition and robustness tables.
diffs=out/'table-diffs';diffs.mkdir(exist_ok=True)
for c in coverage.values():
 rel=Path(c['table']);n=out/'corrected'/rel;p=out/'manuscript'/rel
 if n.exists() and p.exists():
  diff=difflib.unified_diff(p.read_text().splitlines(True),n.read_text().splitlines(True),fromfile='manuscript/'+str(rel),tofile='corrected/'+str(rel))
  (diffs/(rel.stem+'.diff')).write_text(''.join(diff))
print(len(records),'included contrasts;',sum(bool(r['threshold_changes']) for r in records),'threshold crossings')
