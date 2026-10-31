from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
import run as r
s=json.loads((r.OUT/'summary.json').read_text());records=[json.loads(p.read_text()) for p in sorted(r.OUT.glob('*_record.json'))]
assert len(records)==20
maxdiff=0.;selected=[];bounds=[]
for z in records:
 if z['type']!='all':assert not set(z['train_names'])&set(z['test_names'])
 for m,fit in z['fits'].items():
  if not fit['selected_success']:selected.append([z['name'],m])
  if m in r.MODELS:
   b,_=r.spec(m)
   for k,(lo,hi) in enumerate(b):
    if min(abs(fit['p'][k]-lo),abs(fit['p'][k]-hi))<1e-5:bounds.append([z['name'],m,k,fit['p'][k]])
for seed in [20261027,7331,48105]:
 seen=[]
 for z in records:
  if z['name'].startswith(f'cv_{seed}_'):seen+=z['test_names']
 assert sorted(seen)==r.NAMES
 for m in ['NFW_response_MSE']+r.MODELS:
  d=pd.read_csv(r.OUT/f'{m}_cv_{seed}.csv');e=np.log10(d.gpred_m_s2/d.gobs_m_s2)
  mse=pd.DataFrame({'g':d.system,'e2':e*e}).groupby('g').e2.mean();score=np.sqrt(mse.mean());maxdiff=max(maxdiff,abs(score-s['repeats'][str(seed)][m]['rms']))
  if m in r.MODELS:
   n=d.groupby('system').galaxy_mean_joint_nll.first().mean();maxdiff=max(maxdiff,abs(n-s['repeats'][str(seed)][m]['mean_nll']))
assert maxdiff<1e-11 and not selected
# Deterministic nonlinear Gaussian quadrature for halo state, not likelihood.
# This checks the first two moments of log acceleration for concentration and
# excess-mass variability; it does not validate every nuisance or the full density.
qs={}
for model in ['NFW_c','NFW_mc','NFW_response_mc']:
 p=s['full_coefficients'][model];pars=p[:4] if model=='NFW_response_mc' else p[:3];sm=0. if model=='NFW_c' else p[-1]
 center=r.nfw(pars);prev=None;ans=None
 for order in [20,32]:
  v,w=hermgauss(order);nodes=np.sqrt(2)*v;ww=np.outer(w,w).ravel()/np.pi
  mm,cc=np.meshgrid(nodes*sm,nodes*.11,indexing='ij');mm=mm.ravel()[:,None];cc=cc.ravel()[:,None]
  means=[];variances=[];linear=[];nominal=[]
  for j in range(r.NG):
   n=r.COUNT[j];d={k:arr[j:j+1,:n] if arr.ndim==2 and arr.shape[1]>1 else arr[j:j+1] for k,arr in r.D.items() if k in ['r','mp','ms','re','gb','star']}
   vals=r.nfw(pars,d,mm,cc)['log'];mu=ww@vals;var=ww@((vals-mu)**2)
   means.extend(mu);variances.extend(var);linear.extend((sm*center['JM'][j,:n])**2+(.11*center['JC'][j,:n])**2);nominal.extend(center['log'][j,:n])
  ans={k:np.array(v) for k,v in [('mean',means),('variance',variances),('linear',linear),('nominal',nominal)]}
  if prev is not None:convergence=float(max(abs(ans['mean']-prev['mean'])))
  prev=ans
 shift=ans['mean']-ans['nominal'];sdrel=np.sqrt(ans['variance']/ans['linear'])-1
 qs[model]=dict(mass_scatter_dex=sm,order_convergence_max_mean_dex=convergence,mean_shift_absolute_median_dex=float(np.median(abs(shift))),mean_shift_absolute_95th_dex=float(np.percentile(abs(shift),95)),mean_shift_absolute_max_dex=float(max(abs(shift))),relative_sd_error_median=float(np.median(sdrel)),relative_sd_error_95th_abs=float(np.percentile(abs(sdrel),95)),relative_sd_error_max_abs=float(max(abs(sdrel))))
 # Rows follow sorted NAMES, identical to D's row map; retain a full audit table.
 tt=pd.DataFrame({'system':np.repeat(r.NAMES,r.COUNT),'nonlinear_mean_logg':ans['mean'],'nominal_logg':ans['nominal'],'nonlinear_sd_dex':np.sqrt(ans['variance']),'linear_sd_dex':np.sqrt(ans['linear'])});tt.to_csv(r.OUT/f'{model}_nonlinear_moments.csv',index=False)
out=dict(max_recomputed_metric_difference=maxdiff,failed_selected=selected,bound_hits=bounds,nonlinear_halo_moment_checks=qs,passed_numerical_audit=True,interpretation='The moment check covers only halo scatter, not the complete nuisance-integrated likelihood. Large fitted scatter and linearization limit physical interpretation.')
(r.OUT/'audit_post.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
