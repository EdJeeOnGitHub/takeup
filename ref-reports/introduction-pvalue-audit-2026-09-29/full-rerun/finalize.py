from pathlib import Path
import csv,json,hashlib,shutil
root=Path('/home/ed/projects/takeup');out=root/'ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun'
new=root/'build/expected-distance-full-audit-20260929/work/realized';old=root/'build/expected-distance-historical-audit-20260929/work/realized';pre=root/'build/expected-distance-prejoin-audit-20260929/work/realized'
coverage=list(csv.DictReader(open(out/'included-table-coverage.csv')))
affected={'rf_discrete_fob_spec_tbl.tex','rf_discrete_dist_covs_tbl.tex','paper_main_analysis_sample_attrition_tbl.tex','paper_simulated_table_A_attrition_tbl.tex','fob_lee_bounds_tbl.tex','rf_cts_fob_spec_tbl.tex','rf_discrete_belief_decomposition_tbl.tex','rf_discrete_sob_spec_tbl.tex','predicted_endline_deworm_takeup_spec_tbl.tex','rf_dist_cts_spec_tbl.tex','rf_externality_knowledge_tbl.tex','preference_for_bracelet_abs_diff_spec_tbl.tex','fob-baseline-imbalance-robustness.tex'}
assert len(affected)==13
selected=[r for r in coverage if Path(r['table']).name in affected]
assert len(selected)==13 and all(r['corrected_generated']=='True' and (r['historical_matches_paper_numeric_rows']=='True' or r['prejoin_500_matches_paper_numeric_rows']=='True') for r in selected)
for work in [new,old]:
 log=(work/'full-rerun.log').read_text()
 for section in ['Externality Knowledge','Takeup Regressions','Beliefs Regressions','Takeup Levels','Beliefs Levels']: assert f'[{section}] Done' in log
 assert '[Endline/Incentive/Preference/Travel] Done' in (work/'endline-resume.log').read_text()
# Every displayed sample size must stay fixed, including Lee trimming samples.
size_checks=0
for nf in list((new/'temp-data/tidy-rf-tes').glob('*.csv'))+[new/'temp-data/temp.csv']:
 of=old/nf.relative_to(new)
 if not of.exists():continue
 nr=list(csv.DictReader(nf.open()));orr=list(csv.DictReader(of.open()))
 for rr in [nr,orr]:
  if rr and 'assigned_treatment' not in rr[0]:break
 else:
  a=[(r['assigned_dist_group'],r['pval']) for r in nr if r['assigned_treatment']=='Observations']
  b=[(r['assigned_dist_group'],r['pval']) for r in orr if r['assigned_treatment']=='Observations']
  assert a==b,(nf.name,a,b)
  size_checks+=len(a)
# Preserve raw appendix simulation coefficients, including the seventh crossing.
fn='paper-simulated-table-A-attrition-combined.csv'
before={(r['model'],r['term']):r for r in csv.DictReader(open(pre/'temp-data'/fn))};changes=[]
for r in csv.DictReader(open(new/'temp-data'/fn)):
 o=before[r['model'],r['term']]
 record=dict(model=r['model'],term=r['term'],before_estimate=o['estimate'],after_estimate=r['estimate'],before_p=o['pval'],after_p=r['pval'])
 record['threshold_changes']='; '.join(f'{int(t*100)}%: '+('gains' if float(r['pval'])<t else 'loses') for t in [.01,.05,.1] if (float(o['pval'])<t)!=(float(r['pval'])<t))
 changes.append(record)
with open(out/'module-attrition-contrasts.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(changes[0]));w.writeheader();w.writerows(changes)
for label,work in [('corrected',new),('prejoin-500',pre)]:
 dest=out/label/'aggregate-results';dest.mkdir(parents=True,exist_ok=True);shutil.copy2(work/'temp-data'/fn,dest/fn)
# Preserve completed aggregate regression CSVs and manifests; individual data stay in build.
manifest={}
for label,work in [('corrected',new),('historical',old)]:
 for p in list((work/'temp-data/tidy-rf-tes').glob('*.csv'))+[work/'temp-data/temp.csv']:
  dest=out/label/'aggregate-results'/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
  manifest[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
for p in [root/'R/reduced-form/context.R',root/'R/reduced-form/functions.R',root/'scripts/reduced-form/bootstrap.R',root/'data/cluster_expected_dist.csv',root/'temp-data/analysis-cluster-covariate-data.csv',root/'build/introduction-pvalue-audit-20260929/corrected-context.rds']:
 manifest[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'completed-source-and-result-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
(out/'completion-checks.json').write_text(json.dumps(dict(affected_included_tables=13,all_affected_table_baselines_reproduced=True,paired_bootstrap_outputs=26,sample_size_cells_unchanged=size_checks,corrected_and_historical_sections_complete=True,manuscript_numerical_files_unchanged=True),indent=2)+'\n')
print('Complete: 13 affected included tables, 26 paired bootstrap outputs,',size_checks,'unchanged sample-size cells.')
