from pathlib import Path
import re,json,collections,difflib,hashlib
paper=Path('/home/ed/projects/overleaf/overleaf-takeup'); audit=Path('ref-reports/introduction-pvalue-audit-2026-09-29');out=audit/'promotion';stage=out/'staged';backup=out/'before'
files=['rf_discrete_fob_spec_tbl.tex','rf_discrete_dist_covs_tbl.tex','paper_main_analysis_sample_attrition_tbl.tex','paper_simulated_table_A_attrition_tbl.tex','fob_lee_bounds_tbl.tex','rf_cts_fob_spec_tbl.tex','rf_discrete_belief_decomposition_tbl.tex','rf_discrete_sob_spec_tbl.tex','predicted_endline_deworm_takeup_spec_tbl.tex','rf_dist_cts_spec_tbl.tex','rf_externality_knowledge_tbl.tex','preference_for_bracelet_abs_diff_spec_tbl.tex']
paths=['rf-tables/main-specs/'+f for f in files]+['tables/fob-baseline-imbalance-robustness.tex']
manifest={}
def save(rel,old,new):
 for base,text in [(stage,new),(backup,old)]:
  dest=base/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
 manifest[rel]={'before_sha256':hashlib.sha256(old.encode()).hexdigest(),'after_sha256':hashlib.sha256(new.encode()).hexdigest()}
 (out/(Path(rel).name+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=rel,tofile=rel)))
def data_row(line):
 if '&' not in line or 'multicolumn' in line or 'Dependent variable' in line:return None
 label,body=line.split('&',1)
 if not re.search(r'\d',body) or not re.search(r'\\\\\s*$',body.rstrip()):return None
 key=re.sub('[^a-z0-9]','',label.replace('Control mean','Control').lower())
 return key,label,body
for rel in paths:
 old=(paper/rel).read_text();new=(audit/'full-rerun/corrected'/rel).read_text()
 rows=collections.defaultdict(list)
 for line in new.splitlines(True):
  item=data_row(line)
  if item:rows[item[0]].append(item[2])
 result=[];used=0
 for line in old.splitlines(True):
  item=data_row(line)
  if item:
   key,label,body=item
   assert rows[key],(rel,key)
   result.append(label+'&'+rows[key].pop(0));used+=1
  else:result.append(line)
 assert not any(rows.values()),(rel,{k:len(v) for k,v in rows.items() if v})
 save(rel,old,''.join(result));print(rel,used,'data rows')
# Active-text replacements are anchored to their paragraphs; comments are retained.
def edit_lines(rel, rules):
 old=(paper/rel).read_text();lines=old.splitlines(True)
 for anchor,pairs in rules:
  hits=[i for i,l in enumerate(lines) if not l.lstrip().startswith('%') and anchor in l]
  assert len(hits)==1,(rel,anchor,hits)
  i=hits[0]
  for before,after in pairs:
   assert lines[i].count(before)==1,(rel,anchor,before)
   lines[i]=lines[i].replace(before,after)
 save(rel,old,''.join(lines))
edit_lines('ECM ReStud.tex',[
 ('We first document how distance', [('75 percent','74 percent'),('p=0.056','p=0.020'),('6 percentage points in Far','7 percentage points in Far'),('p=\\)','p=0.665\\)'),('25 percentage points','26 percentage points'),('p=0.007','p=0.005'),('5--6 percentage points','5--7 percentage points'),('p=0.357','p=0.185'),('0.644','0.324'),('p=0.009','p=0.008')]),
 ('We then examine take-up', [('16 percentage points','16.5 percentage points'),('\\(-0.3\\) percentage points (\\(p=\\))','\\(0.4\\) percentage points (\\(p=0.940\\))'),('12 percentage point decline','13 percentage point decline'),('p=0.002','p=0.003'),('7.5 percentage points','7.8 percentage points'),('p=0.008','p=0.007'),('2.7 percentage points relative to Control (\\(p=\\))','2.6 percentage points relative to Control (\\(p=0.359\\))'),('p=0.048','p=0.034'),('p=0.35\\)','p=0.36\\)'),('p=0.083','p=0.088')]),
 ('Two alternative explanations', [('p=0.959','p=0.940')])])
edit_lines('experimental-design.tex',[
 ('(\\(p = 0.081\\)) but not',[('p = 0.081','p = 0.063'),('p = 0.195','p = 0.181')]),
 ('order of magnitude smaller',[('16.2','16.5')]),
 ('distance gap in take-up. Finally',[('p=0.511','p=0.488'),('p=0.754','p=0.672')]),
 ('arms: joint tests do not reject',[('p = 0.431','p = 0.270'),('p = 0.220','p = 0.158'),('The only significant difference is 8.8 percentage points lower missingness for Bracelet respondents in Close communities','At the 5 percent level, the only arm difference is 9.1 percentage points lower missingness for Bracelet respondents in Close communities'),('To account for this imbalance','The pooled Bracelet difference is 6.5 percentage points (\\(p=0.067\\)). To account for this imbalance')]),
 ('observability is bounded below',[('12.2','12.7')]),
 ('qualitative pattern of the main results is unchanged. Again',[('Again, we find no significant difference in omission between the bracelet and calendar arm, either overall (\\(p=0.316\\)) or across distance cells (\\(0.111\\)).','We do not reject equality of omission rates between the bracelet and calendar arms overall (\\(p=0.205\\)) or jointly across distance cells (\\(p=0.119\\)).')]),
 ('In the no-item Control arm, observability',[('68.1','67.9'),('74.5','74.4'),('59.9','59.6'),('-14.6','-14.8')]),
 ('Public signals substantially increase observability',[('14.4','14.9'),('6.4','6.6'),('24.8','25.6'),('5.2','5.3'),('16.9','16.7'),('4.3','4.4'),('p=0.009','p=0.008'),('p=0.001','p<0.001'),('p=0.029','p=0.027')]),
 ('Table~\\ref{tab:overall_effects} reports treatment effects',[('33.8','33.7'),('41.0','41.1'),('24.8','24.6'),('16.2','16.5'),('about 12 percentage','about 13 percentage')]),
 ('Bracelets increase take-up by 7.5',[('7.5','7.8'),('p<0.05','p<0.01'),('4.8','5.2'),('p=0.048','p=0.034'),('2.7','2.6')]),
 ('The distance interaction differs between',[('p=0.083','p=0.088')]),
 ('The reduced form evidence delivers',[('about 16 percentage','about 16.5 percentage')])])
edit_lines('online-appendix.tex',[
 ('Ink Far--Close interactions are 17.5',[('17.5 and 11.5','18.2 and 11.3')]),
 ('Calendar interaction is 2.4',[('2.4','2.8')])])
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Prepared',len(manifest),'files; manuscript not yet written.')
