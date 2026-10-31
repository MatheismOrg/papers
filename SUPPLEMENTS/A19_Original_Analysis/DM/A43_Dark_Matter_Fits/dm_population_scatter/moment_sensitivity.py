"""Post-audit nonlinear halo-moment sensitivity with coefficients frozen.

Gaussian moment matching is not the exact nonlinear mixture density. Distance,
inclination and stellar uncertainty remain first-order. No parameter is refit.
"""
import json
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
import run as r
v,w=hermgauss(32);nodes=np.sqrt(2)*v;weights=np.outer(w,w).ravel()/np.pi
out=[]
for path in sorted(r.OUT.glob('cv_*_record.json')):
 rec=json.loads(path.read_text());p=rec['fits']['NFW_mc']['p'];sm=p[4];m,c=np.meshgrid(nodes*sm,nodes*.11,indexing='ij');m=m.ravel()[:,None];c=c.ravel()[:,None]
 cent,var,U=r.distribution('NFW_mc',p);pred=[];sd=[];scores=[]
 for j in rec['test_gi']:
  n=r.COUNT[j];d={k:arr[j:j+1,:n] if arr.ndim==2 and arr.shape[1]>1 else arr[j:j+1] for k,arr in r.D.items() if k in ['r','mp','ms','re','gb','star']}
  vals=r.nfw(p[:3],d,m,c)['log'];mu=weights@vals;dev=vals-mu
  cov=dev.T@(weights[:,None]*dev);obs=U[j,:n,:3]
  total=np.diag(var[j,:n])+obs@obs.T+cov;sign,ld=np.linalg.slogdet(total);assert sign>0
  e=mu-np.log10(r.D['go'][j,:n]);score=.5*(n*np.log(2*np.pi)+ld+e@np.linalg.solve(total,e))/n
  pp=np.zeros(r.NMAX);ss=np.ones(r.NMAX);pp[:n]=mu;ss[:n]=np.sqrt(np.diag(total));pred.append(pp.tolist());sd.append(ss.tolist());scores.append(float(score))
 out.append(dict(name=rec['name'],test_gi=rec['test_gi'],log=pred,sd=sd,nll=scores))
(r.OUT/'moment_sensitivity_records.json').write_text(json.dumps(out,indent=2))
cv={};mse=[];nl=[];base=[]
for seed in [20261027,7331,48105]:
 logs=np.zeros(r.MASK.shape);sd=np.ones_like(logs);nll=np.zeros(r.NG)
 for rec in out:
  if rec['name'].startswith(f'cv_{seed}_'):logs[rec['test_gi']]=rec['log'];sd[rec['test_gi']]=rec['sd'];nll[rec['test_gi']]=rec['nll']
 cv[str(seed)],ee=r.evaluate(logs,sd,nll);mse.append(ee);nl.append(nll)
 tab=r.DF[['row','system','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();tab['gpred_m_s2']=10**r.flatten(logs);tab['predictive_sd_dex']=r.flatten(sd);tab['log10_residual']=r.flatten(logs-np.log10(r.D['go']));tab.to_csv(r.OUT/f'NFW_mc_moment_cv_{seed}.csv',index=False)
 t=pd.read_csv(r.OUT/f'RAR_offset_cv_{seed}.csv');base.append(t.groupby('system').galaxy_mean_joint_nll.first().reindex(r.NAMES).to_numpy())
means={k:float(np.mean([cv[str(seed)][k] for seed in [20261027,7331,48105]])) for k in cv['20261027']}
ix=np.random.default_rng(20261029).integers(0,r.NG,(10000,r.NG));diff=np.array(nl)-np.array(base);bs=diff[:,ix].mean((0,2));contrast=dict(difference=float(diff.mean()),CI95=np.percentile(bs,[2.5,97.5]).tolist())
s=dict(means=means,repeats=cv,minus_RAR_offset_logscore=contrast,scope='Frozen linear-fit coefficients, nonlinear first two halo moments, Gaussian approximation; remaining nuisance covariance is linearized. No exact nonlinear marginal likelihood or new sample.')
(r.OUT/'moment_sensitivity_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
