"""Audit stellar rules and apply the existing force blend to their new base."""
from pathlib import Path
import importlib.util,json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stellar_models',ROOT/'run.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
pop=s.pop;summary=json.loads((s.OUT/'summary.json').read_text());recs=[json.loads(p.read_text()) for p in sorted(s.OUT.glob('*_record.json'))]
assert len(recs)==20;maxerr=0.;bounds=[];failed=[]
for rec in recs:
 if rec['type']!='all':assert not set(rec['train_gi'])&set(rec['test_gi'])
 assert rec['fits']['stellar_curved']['loss']<=rec['fits']['stellar_power']['loss']+1e-12
 for m,fit in rec['fits'].items():
  if not fit.get('selected_success',True):failed.append([rec['name'],m])
  if m=='stellar_curved' and abs(fit['p'][4]-1)<1e-6:bounds.append(rec['name'])
for seed in [20261027,7331,48105]:
 seen=[j for r in recs if r['name'].startswith(f'cv_{seed}_') for j in r['test_gi']];assert sorted(seen)==list(range(pop.NG))
 for m in s.MODELS:
  d=pd.read_csv(s.OUT/f'{m}_cv_{seed}.csv');e=np.log10(d.gpred_m_s2/d.gobs_m_s2);rms=np.sqrt(pd.DataFrame({'g':d.system,'e2':e*e}).groupby('g').e2.mean().mean());maxerr=max(maxerr,abs(rms-summary['repeats'][str(seed)][m]['rms']))
assert not failed and maxerr<1e-10
assert np.max(abs(np.log10(s.moster(np.log10(s.AM))/s.D['ms'])))<1e-10
# New predictor is still everywhere finite and nonnegative, with a baryon floor.
for model,p in summary['full_coefficients'].items():
 mass=s.mass(model,p);assert np.all(mass>=s.D['mp']/s.FB) and np.all(np.isfinite(mass))
# Fixed published stellar relation stays monotone over the entire test domain.
xx=np.linspace(5,18,500);assert np.all(np.diff(s.moster(xx))>0)

ens=[]
for rec in recs:
 old=json.loads((pop.OLD/'results'/f"{rec['name']}_record.json").read_text());ra=pop.rar(old['models']['RAR']['fit']['p'][0])['log'];nf=s.predict('stellar_power',rec['fits']['stellar_power']['p']);tr=rec['train_gi'];te=rec['test_gi']
 def pred(w):return np.log10((1-w)*10**ra+w*10**nf)
 def obj(w):return float((((pred(w)-np.log10(s.D['go']))**2*s.MASK).sum(1)/s.COUNT)[tr].mean())
 z=minimize_scalar(obj,bounds=(0,1),method='bounded',options={'xatol':1e-12});w=min([0.,1.,z.x],key=obj)
 ens.append(dict(name=rec['name'],type=rec['type'],train_gi=tr,test_gi=te,weight=float(w),train_loss=obj(w),RAR_coefficients=old['models']['RAR']['fit']['p'],stellar_coefficients=rec['fits']['stellar_power']['p'],log_predictions=pred(w)[te].tolist(),optimizer_success=bool(z.success)))
(s.OUT/'stellar_blend_records.json').write_text(json.dumps(ens,indent=2))
cv={};mses=[]
for seed in [20261027,7331,48105]:
 pr=np.full(s.MASK.shape,np.nan)
 for rec in ens:
  if rec['name'].startswith(f'cv_{seed}_'):pr[rec['test_gi']]=rec['log_predictions']
 metric,mse=pop.evaluate(pr);cv[str(seed)]=metric;mses.append(mse)
 t=s.DF[['row','system','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();t['gpred_m_s2']=10**pop.flatten(pr);t['log10_residual']=pop.flatten(pr-np.log10(s.D['go']));t.to_csv(s.OUT/f'stellar_blend_cv_{seed}.csv',index=False)
means={k:float(np.mean([cv[str(seed)][k] for seed in [20261027,7331,48105]])) for k in ['rms','median_speed']};tails={}
for rec in ens:
 if rec['type']!='tail':continue
 j=rec['test_gi'];e=(np.array(rec['log_predictions'])-np.log10(s.D['go'][j]))*s.MASK[j];tails[rec['name']]=float(np.sqrt(np.mean(np.sum(e*e,axis=1)/s.COUNT[j])))
# All-data coefficients frozen before transfer. Retain subsets and conditional
# object bootstrap; host-correlated systematic errors are not modeled by it.
full=next(e for e in ens if e['type']=='all');dd=pd.read_csv(pop.OLD/'dwarfs_benchmark.csv');d=dict(r=dd.radius_kpc.to_numpy()[:,None],gb=dd.gbar_m_s2.to_numpy()[:,None],go=dd.gobs_m_s2.to_numpy()[:,None],star=dd.gbar_m_s2.to_numpy()[:,None],mp=dd.stellar_mass_Msun.to_numpy()[:,None],ms=dd.stellar_mass_Msun.to_numpy()[:,None],re=dd.projected_Re_circularized_kpc.to_numpy()[:,None]);pr=np.log10((1-full['weight'])*10**pop.rar(full['RAR_coefficients'][0],d)['log']+full['weight']*10**s.predict('stellar_power',full['stellar_coefficients'],d,s.inverse(d['ms'])))[:,0]
dt=pd.read_csv(s.OUT/'frozen_dwarf_transfer.csv');dt['stellar_blend']=10**pr
old=pd.read_csv(pop.OLD/'results/frozen_dwarf_transfer.csv');assert np.array_equal(old.system,dt.system)
for name in ['RAR','NFW','RAR_NFW_blend']:dt['old_'+name]=old[name]
dt.to_csv(s.OUT/'dwarf_comparison.csv',index=False)
models=s.MODELS+['stellar_blend','old_RAR','old_NFW','old_RAR_NFW_blend'];de={m:np.log10(dt[m]/dt.gobs_m_s2).to_numpy() for m in models}
subsets={}
for label,mask in [('all',np.ones(len(dt),bool)),('Milky_Way',dd.host.str.startswith('G')),('Andromeda',dd.host.str.startswith('A')),('classical',dd.population=='classical_dSph'),('ultrafaint',dd.population=='ultrafaint_dSph')]:
 mask=np.array(mask);subsets[label]=dict(n=int(sum(mask)),scores={m:float(np.sqrt(np.mean(de[m][mask]**2))) for m in models} if mask.any() else {})
ix=np.random.default_rng(2981043).integers(0,len(dt),(10000,len(dt)));dboot={m:np.sqrt((e[ix]**2).mean(1)) for m,e in de.items()};dcontr={}
for m in ['stellar_power','stellar_curved','stellar_blend']:
 dif=dboot[m]-dboot['old_NFW'];dcontr[m+' minus old_NFW']=dict(difference=float(np.sqrt(np.mean(de[m]**2))-np.sqrt(np.mean(de['old_NFW']**2))),CI95=np.percentile(dif,[2.5,97.5]).tolist())
# Correlated repeat bootstrap, not new data or model-search adjustment.
ix=np.random.default_rng(20261029).integers(0,pop.NG,(10000,pop.NG));bb=np.sqrt(np.array(mses)[:,ix].mean(2)).mean(0);contr={}
for name in ['RAR','NFW','RAR_NFW_blend']:
 vals=[]
 for seed in [20261027,7331,48105]:
  q=pd.read_csv(pop.OLD/f'results/{name}_cv_{seed}.csv');vals.append(q.groupby('system').log10_residual.apply(lambda z:np.mean(z*z)).reindex(s.NAMES).to_numpy())
 vv=np.array(vals);boot=np.sqrt(vv[:,ix].mean(2)).mean(0);diff=bb-boot;contr[name]=dict(difference=float(means['rms']-np.sqrt(vv.mean(1)).mean()),CI95=np.percentile(diff,[2.5,97.5]).tolist())
out=dict(audit=dict(max_metric_difference=maxerr,failed_selected=failed,curved_transition_upper_bound_cases=bounds,passed=True),blend=dict(means=means,repeats=cv,tails=tails,full_weight=full['weight'],contrasts=contr),dwarf_subsets=subsets,dwarf_contrasts=dcontr)
(s.OUT/'ensemble_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
