import csv,collections
sums=collections.defaultdict(float);counts=collections.Counter()
with open('data/simulated-counterfactual-treatment-assignment.csv') as f:
 for r in csv.DictReader(f):
  if int(r['vill_seed'])<=100 and r['cluster.id']!='NA':
   k=r['cluster.id'];sums[k]+=float(r['dist']);counts[k]+=1
with open('build/introduction-pvalue-audit-20260929/expected-distance-100-seeds.csv','w') as f:
 w=csv.writer(f);w.writerow(['cluster.id','dist'])
 for k in sorted(sums,key=int):w.writerow([k,sums[k]/counts[k]])
print('Reconstructed',len(sums),'communities using village seeds 1:100.')
