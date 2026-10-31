"""Independent direct-formula and fold/score verification; no fitting."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
R=Path(__file__).resolve().parent;J=json.loads((R/'results/results.json').read_text());D=R.parent/'dm_broad_joint/combined_input.csv';rows=list(csv.DictReader(D.open()));group=np.array([r['group'] for r in rows]);pop=np.array([r['population'] for r in rows]);env=np.array([r.get('environment','') for r in rows]);POPS=sorted(set(pop))
assert hashlib.sha256(D.read_bytes()).hexdigest()==J['input_sha256'];assert hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest()==J['protocol_sha256']
def formula(k,c,g,r):
 x=np.log10(g/1e-10);u=np.log10(r);q=x+u;m=x+2*u
 if k=='RAR':a=1e-10*10**c[0]
 elif k=='potential_quad':a=1e-10*10**(c[0]+c[1]*q+c[2]*q*q)
 elif k=='mass_bounded':z=np.tanh((m-c[3])/c[4]);a=1e-10*10**(c[0]+c[1]*z+c[2]*z*z)
 else:
  alpha=.5 if k=='half_power' else c[3];beta=c[3] if k=='half_power' else c[4]
  s=q if k in ['power_sum','half_power'] else m
  t=q if k in ['power_sum','half_power','mixed_power'] else x+c[5]*u if k=='index_power' else m
  a=1e-10*(10**c[0]+10**(c[1]-alpha*s)+10**(c[2]+beta*t))
 if k=='mass_additive':return g+np.sqrt(g*a)
 if k=='mass_screened':return g+np.sqrt(g*a)*np.exp(-np.sqrt(g/(1e-10*10**c[0])))
 return g/(-np.expm1(-np.sqrt(g/a)))
def B(pred,obs,indices):
 pp=[]
 for p in sorted(set(pop[indices])):
  ii=np.where(pop[indices]==p)[0];gg=group[indices[ii]];pp.append(np.mean([np.mean(np.log10(pred[ii[gg==g]]/obs[ii[gg==g]])**2) for g in set(gg)]))
 return np.sqrt(np.mean(pp))
files=[];maxerr=0.
for k in J['models']:
 d=J['models'][k]
 for split in ['all','cv_20261027','cv_7331','cv_48105','population_holdout','host_MW','host_M31']:
  f=R/'results'/f'{k}_{split}.csv';rr=list(csv.DictReader(f.open()));ii=np.array([int(r['row']) for r in rr]);gb=np.array([float(r['gbar_m_s2']) for r in rr]);rad=np.array([float(r['radius_kpc']) for r in rr]);obs=np.array([float(r['gobs_m_s2']) for r in rr]);pred=np.array([float(r['gpred_m_s2']) for r in rr]);err=np.array([float(r['log10_residual']) for r in rr]);assert max(abs(np.log10(pred/obs)-err))<1e-12
  assert len(set(ii))==len(ii)
  for g in set(group[ii]):assert len({r['fold'] for r in rr if r['group']==g})==1
  if split=='all':cc=[d['coefficients']]*len(rr);met=d['training'];assert len(rr)==2910
  elif split.startswith('cv_'):
   seed=split.split('_')[1];cc=[d['CV'][seed]['coefficients'][int(r['fold'])] for r in rr];met=d['CV'][seed]['metrics'];assert len(rr)==2910
  elif split=='population_holdout':cc=[d[split]['coefficients'][r['population']] for r in rr];met=d[split]['metrics'];assert all(r['fold']==r['population'] for r in rr)
  else:
   h=split.split('_')[1];z=d['host_holdout'][h];cc=[z['coefficients']]*len(rr);met=z['metrics'];aliases={'MW':['Milky_Way','G/L'],'M31':['M31','A/L']};assert len(rr)=={'MW':26,'M31':18}[h];assert set(ii)==set(np.where((pop=='dSph')&np.isin(env,aliases[h]))[0]);assert z['training_rows']==2910-len(ii)
  direct=np.array([formula(k,c,g,r) for c,g,r in zip(cc,gb,rad)]);delta=max(abs(np.log10(direct/pred)));maxerr=max(maxerr,delta);assert delta<1e-10,(k,split,delta)
  assert abs(B(pred,obs,ii)-met['B_dex'])<1e-12,(k,split)
  files.append(f.name)
assert len(files)==70
out={'status':'passed','prediction_tables':len(files),'max_formula_difference_dex':maxerr,'input_and_protocol_hashes':'valid','full_host_exclusions':{'MW':26,'M31':18},'checks':['independent direct formulas','independent equal-population/equal-system score','residual consistency','whole groups kept in one fold','host alias exclusions','prediction row identities'],'files':files}
(R/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
