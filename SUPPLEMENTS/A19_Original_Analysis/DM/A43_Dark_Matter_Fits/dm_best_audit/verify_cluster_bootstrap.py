from pathlib import Path
import json,csv,hashlib
import numpy as np
B=Path(__file__).resolve().parent.parent;D=B/'dm_cluster_repair';O=Path(__file__).resolve().parent
rr=list(csv.DictReader((D/'predictions.csv').open()));DI=json.loads((D/'diagnostics.json').read_text());rng=np.random.default_rng(20261028);count=0
for split in ['host_holdout','source_CLASH','source_XCOP']:
 chosen=[r for r in rr if r['case']=='combined_forward_mid' and r['split']==split];hosts=sorted(set(r['system'] for r in chosen));sources=sorted(set(r['source'] for r in chosen));multiplicity=rng.poisson(1,(5000,len(hosts)));position={h:i for i,h in enumerate(hosts)};scores={}
 for k in ['simple_repair','simple_refit_a0','simple_frozen_gas']:
  pred=[r for r in chosen if r['model']==k];vals=[]
  for s in sources:
   hs=sorted(set(r['system'] for r in pred if r['source']==s));ids=[position[h] for h in hs];err=np.array([np.mean([np.log10(float(r['gpred_m_s2'])/float(r['gobs_m_s2']))**2 for r in pred if r['system']==h]) for h in hs]);m=multiplicity[:,ids];den=m.sum(axis=1);vals.append(np.divide(m@err,den,out=np.full(5000,np.nan),where=den>0))
  scores[k]=np.sqrt(np.mean(vals,axis=0))
 for base in ['simple_refit_a0','simple_frozen_gas']:
  dif=scores['simple_repair']-scores[base];valid=np.isfinite(dif);a=DI['bootstrap'][split+' repair minus '+base]
  assert sum(valid)==a['finite_draws'];assert max(abs(np.percentile(dif[valid],[2.5,97.5])-a['conditional_95_interval_dex']))<1e-12;assert abs(np.mean(dif[valid]<0)-a['fraction_better_descriptive'])<1e-12;count+=1
j=json.loads((O/'cluster_verification.json').read_text());j['descriptive_bootstrap_comparisons_checked']=count;j['bootstrap_audit_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(O/'cluster_verification.json').write_text(json.dumps(j,indent=2));print('PASS',count,'cluster bootstrap comparisons')
