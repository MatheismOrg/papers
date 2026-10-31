from pathlib import Path
import csv,json,hashlib
import numpy as np
R=Path(__file__).resolve().parent;J=json.loads((R/'results/results.json').read_text())
assert hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest()==J['protocol_sha256']
for f,h in J['input_hashes'].items():assert hashlib.sha256((R/f).read_bytes()).hexdigest()==h
for f,h in J['source_hashes'].items():assert hashlib.sha256((R.parent/f).read_bytes()).hexdigest()==h

def direct(k,c,g,r):
 x=np.log10(g/1e-10);u=np.log10(r);q=x+u;m=x+2*u
 if k=='RAR':a=1e-10*10**c[0]
 elif k=='potential_quad':a=1e-10*10**(c[0]+c[1]*q+c[2]*q*q)
 else:
  ss=q if k=='power_sum' else m;neg=10**(c[1]-c[3]*ss);pos=10**(c[2]+c[4]*q)
  if k in ['mixed_gate','mixed_gate_range']:neg/=1+(g/1e-8)**2
  if k=='mixed_gate_range':pos/=1+(r/10**c[5])**2
  a=1e-10*(10**c[0]+neg+pos)
 return g/(-np.expm1(-np.sqrt(g/a)))

def B(pred,obs,indices,pop,group):
 pp=[]
 for p in sorted(set(pop[indices])):
  ii=np.where(pop[indices]==p)[0];gg=group[indices[ii]];pp.append(np.mean([np.mean(np.log10(pred[ii[gg==g]]/obs[ii[gg==g]])**2) for g in set(gg)]))
 return np.sqrt(np.mean(pp))

allfiles=[];maxerr=0.
for ds,j in J['models_by_dataset'].items():
 rows=list(csv.DictReader((R/f'{ds}_input.csv').open()));pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows]);source=np.array(['CLASH' if r['source']=='Tian2020_CDS_fig2' else r['source'] for r in rows]);N=len(rows)
 if ds=='augmented':
  original=list(csv.DictReader((R/'original_input.csv').open()));assert len(rows)==2931
  for a,b in zip(original,rows[:2910]):
   c=dict(a)
   if c['population']=='CLASH':c['population']='clusters'
   assert b==c
  assert len(set(group[pop=='clusters']))==27 and len(set(group))==285
  assert set(pop)=={'SPARC','ETG','dSph','MaNGA_BCG','clusters'}
 for k,d in j['models'].items():
  splits=['all','cv_20261027','cv_7331','cv_48105','population_holdout','host_MW','host_M31']+(['source_XCOP','source_CLASH'] if ds=='augmented' else [])
  for split in splits:
   f=R/'results'/ds/f'{k}_{split}.csv';rr=list(csv.DictReader(f.open()));ii=np.array([int(r['row']) for r in rr]);gb=np.array([float(r['gbar_m_s2']) for r in rr]);rad=np.array([float(r['radius_kpc']) for r in rr]);obs=np.array([float(r['gobs_m_s2']) for r in rr]);pred=np.array([float(r['gpred_m_s2']) for r in rr]);err=np.array([float(r['log10_residual']) for r in rr]);assert max(abs(np.log10(pred/obs)-err))<1e-12
   assert len(set(ii))==len(ii)
   for g in set(group[ii]):assert len({r['fold'] for r in rr if r['group']==g})==1
   if split=='all':cc=[d['coefficients']]*len(rr);met=d['training'];assert len(rr)==N
   elif split.startswith('cv_'):
    seed=split.split('_')[1];cc=[d['CV'][seed]['coefficients'][int(r['fold'])] for r in rr];met=d['CV'][seed]['metrics'];assert len(rr)==N
   elif split=='population_holdout':cc=[d[split]['coefficients'][r['population']] for r in rr];met=d[split]['metrics'];assert all(r['fold']==r['population'] for r in rr)
   elif split.startswith('host_'):
    h=split.split('_')[1];z=d['host_holdout'][h];cc=[z['coefficients']]*len(rr);met=z['metrics'];aliases={'MW':['Milky_Way','G/L'],'M31':['M31','A/L']};assert len(rr)=={'MW':26,'M31':18}[h];assert set(ii)==set(np.where((pop=='dSph')&np.isin(env,aliases[h]))[0]);assert z['training_rows']==N-len(ii)
   else:
    src=split.split('_')[1];z=d['source_holdout'][src];cc=[z['coefficients']]*len(rr);met=z['metrics'];assert len(rr)=={'XCOP':21,'CLASH':84}[src];assert set(ii)==set(np.where(source==src)[0]);assert z['training_rows']==N-len(ii)
   p2=np.array([direct(k,c,g,r) for c,g,r in zip(cc,gb,rad)]);delta=max(abs(np.log10(p2/pred)));maxerr=max(maxerr,delta);assert delta<1e-10
   assert abs(B(pred,obs,ii,pop,group)-met['B_dex'])<1e-12
   for src,ms in met['sources'].items():
    z=np.where(source[ii]==src)[0];assert len(z)==ms['rows'];assert abs(B(pred[z],obs[z],ii[z],pop,group)-ms['group_RMS_dex'])<1e-12
   allfiles.append(str(f.relative_to(R)))
assert len(allfiles)==96
out={'status':'passed','prediction_tables':len(allfiles),'max_formula_difference_dex':maxerr,'checks':['direct formulas','independent scores including source subsets','all original rows preserved except deliberate CLASH class rename','27 equal-weight cluster groups','96 tables and fold identities','26 MW and18 M31 host exclusions','21 X-COP and84 CLASH source exclusions','input/protocol/source hashes'],'files':allfiles}
(R/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='files'},indent=2))
