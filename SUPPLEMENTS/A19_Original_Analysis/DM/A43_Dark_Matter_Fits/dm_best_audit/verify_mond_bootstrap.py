from pathlib import Path
import csv,json,hashlib
import numpy as np
B=Path(__file__).resolve().parent.parent;D=B/'dm_mond_repairs';O=Path(__file__).resolve().parent
rows=list(csv.DictReader((D/'input.csv').open()));group=np.array([r['group'] for r in rows]);pop=np.array([r['population'] for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]);J=json.loads((D/'results/results.json').read_text())
gs=np.unique(group);gi={g:i for i,g in enumerate(gs)};mult=np.random.default_rng(20261028).poisson(1,size=(5000,len(gs)));count=0
for split in ['cv_20261027','population_holdout']:
 b={};point={}
 for k in J['families']+['mixed_regular']:
  p=B/'dm_unified_regular/results/augmented'/f'{k}_{split}.csv' if k=='mixed_regular' else D/'results'/f'{k}_{split}.csv'
  rr=list(csv.DictReader(p.open()));pr=np.log10([float(r['gpred_m_s2'])/1e-10 for r in rr]);out=[];v=[]
  for popname in sorted(set(pop)):
   gg=sorted(set(group[pop==popname]));eg=np.array([np.mean((pr[(pop==popname)&(group==g)]-y[(pop==popname)&(group==g)])**2) for g in gg]);m=mult[:,[gi[g] for g in gg]];out.append(m@eg/m.sum(axis=1));v.append(eg.mean())
  b[k]=np.sqrt(np.mean(out,axis=0));point[k]=np.sqrt(np.mean(v))
 for key,d in J['bootstrap_'+split].items():
  k,base=key.split(' minus ');delta=b[k]-b[base]
  assert abs(point[k]-point[base]-d['point_difference_dex'])<1e-12
  assert max(abs(np.percentile(delta,[2.5,97.5])-d['percentile_95_dex']))<1e-12
  assert abs(np.mean(delta<0)-d['fraction_negative_descriptive_only'])<1e-12;count+=1
j=json.loads((O/'mond_verification.json').read_text());j['descriptive_bootstrap_comparisons_checked']=count;j['bootstrap_audit_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(O/'mond_verification.json').write_text(json.dumps(j,indent=2));print('PASS',count,'descriptive bootstrap comparisons')
