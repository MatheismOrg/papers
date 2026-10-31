"""Established feedback profiles under A43's shared-predictor protocol."""
from pathlib import Path
import os,json,hashlib,time,sys
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy.special import hyp2f1
from scipy.optimize import least_squares,minimize,minimize_scalar
from scipy.integrate import quad
from concurrent.futures import ProcessPoolExecutor,as_completed

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
G=6.67430e-11*1.988409870698051e30/3.085677581491367e19/1e6
UNIT=1e6/3.085677581491367e19;FB=.049/.3175;H=.671;H0=.0671
GYR=3.15576e16;KPC_KM=3.085677581491367e16
SEEDS=[20261027,7331,48105]
STARTS=[[11.5,.7,0],[12,1.2,.5],[10.5,.4,-.5]]
BOUNDS=([8,.1,-1.5],[14,2,1.5])
MODELS=['RAR','NFW','DC14','DC14_response','coreNFW_duration']
DF=pd.read_csv(ROOT/'input.csv');FEAT=pd.read_csv(ROOT/'galaxy_features.csv').set_index('system')
NAMES=sorted(DF.system.unique());GI={g:i for i,g in enumerate(NAMES)}
DATA={'id':DF.system.to_numpy(),'group':DF.group.to_numpy(),
 'gi':np.array([GI[g] for g in DF.system]),'r':DF.radius_kpc.to_numpy(),
 'gb':DF.gbar_m_s2.to_numpy(),'go':DF.gobs_m_s2.to_numpy(),
 'mp':np.array([FEAT.loc[g,'mass_proxy_msun'] for g in NAMES]),
 'ms':np.array([.5e9*FEAT.loc[g,'luminosity_3p6_1e9_lsun'] for g in NAMES]),
 're':np.array([FEAT.loc[g,'reff_kpc'] for g in NAMES])}
I=np.arange(len(DF))

def f(x):
 x=np.asarray(x,dtype=float);out=np.log1p(x)-x/(1+x);sm=x<.01
 out=np.where(sm,sum((-1)**j*(j-1)/j*x**j for j in range(2,13)),out)
 return out

def reference(p,d=DATA):
 mr=d['mp']/FB+10**p[0]*(d['mp']/1e10)**p[1]*(d['re']/3)**p[2]
 c=10**(.905-.101*np.log10(mr*H/1e12));R=np.cbrt(G*mr/(100*H0**2));rs=R/c
 u=np.full_like(c,np.log(1.4))
 for _ in range(8):
  z=c*np.exp(u);fz=f(z)
  val=np.log(fz/f(c))-np.log(92.8/200)-3*u
  deriv=z*z/((1+z)**2*fz)-3
  u-=val/deriv
 q=np.exp(u);cv=c*q;Rv=R*q;mh=(1-FB)*mr*f(cv)/f(c)
 return dict(mr=mr,c=c,R=R,rs=rs,cv=cv,Rv=Rv,mh=mh)

def shape(X):
 xx=np.clip(X,-4.1,-1.3)
 a=2.94-np.log10(10**(-1.08*(xx+2.33))+10**(2.29*(xx+2.33)))
 b=4.23+1.34*xx+.26*xx**2
 g=-.06+np.log10(10**(-.68*(xx+2.56))+10**(xx+2.56))
 cf=1+3e-5*np.exp(3.4*(xx+4.5))
 low=X<-4.1
 return np.where(low,1.,a),np.where(low,3.,b),np.where(low,1.,g),np.where(low,1.,cf)

def mass_integral(x,a,b,g):
 x,a,b,g=np.broadcast_arrays(x,a,b,g)
 p=(3-g)/a;q=(b-g)/a
 v=x**(3-g)/(3-g)*hyp2f1(p,q,p+1,-x**a)
 if not np.all(np.isfinite(v)&(v>0)):
  raise FloatingPointError('DC14 integral nonfinite/nonpositive')
 return v

def components(p,d=DATA):
 h=reference(p,d);gi=d['gi'];rr=np.minimum(d['r'],h['Rv'][gi])
 mn=h['mh'][gi]*f(rr/h['rs'][gi])/f(h['cv'][gi])
 X=np.log10(d['ms']/h['mh']);a,b,g,cf=shape(X)
 cdc=h['cv']*cf;rm=h['Rv']/cdc
 rs=rm*((b-2)/(2-g))**(1/a)
 denom=mass_integral(h['Rv']/rs,a,b,g)
 md=h['mh'][gi]*mass_integral(rr/rs[gi],a[gi],b[gi],g[gi])/denom[gi]
 low=X[gi]<-4.1;md=np.where(low,mn,md)
 h.update(X=X,a=a,b=b,g=g,cf=cf,dc_rs=rs)
 return mn,md,h

def predict(model,p,d=DATA):
 if model=='RAR':return d['gb']/(-np.expm1(-np.sqrt(d['gb']/(1e-10*10**p[0]))))
 if model in ('DC14','DC14_response'):
  mn,md,h=components(p,d)
  mass=md if model=='DC14' else (1-p[3])*mn+p[3]*md
 else:
  h=reference(p,d);gi=d['gi'];rr=np.minimum(d['r'],h['Rv'][gi]);mass=h['mh'][gi]*f(rr/h['rs'][gi])/f(h['cv'][gi])
  if model=='coreNFW_duration':
   mr_s=h['mh']*f(np.ones_like(h['cv']))/f(h['cv'])
   tdyn=2*np.pi*np.sqrt(h['rs']**3/(G*mr_s))*KPC_KM/GYR
   nn=np.tanh(.04*14/tdyn);rc=1.75*d['re']
   mass*=np.power(np.tanh(rr/rc[gi])/np.tanh(h['Rv'][gi]/rc[gi]),nn[gi])
 return d['gb']+G*mass/d['r']**2*UNIT

def weights(ii):
 gi=DATA['gi'][ii];counts=np.bincount(gi,minlength=len(NAMES));return 1/(np.count_nonzero(counts)*counts[gi])

def metrics(pred,ii):
 e=np.log10(pred/DATA['go'][ii]);w=weights(ii)
 return dict(rms=float(np.sqrt(w@e**2)),median_speed=float(np.median(np.abs(10**(.5*e)-1))*100),galaxies=int(len(np.unique(DATA['gi'][ii]))),rows=len(ii))

def fit(model,ii,base=None,dc=None):
 sw=np.sqrt(weights(ii));go=DATA['go'][ii]
 def resid(p):return sw*np.log10(predict(model,p)[ii]/go)
 def loss(p):return float(np.sum(resid(p)**2))
 if model=='RAR':
  z=minimize_scalar(lambda x:loss([x]),bounds=(-3,3),method='bounded',options={'xatol':1e-12});return dict(p=[float(z.x)],loss=float(z.fun),attempts=[dict(method='bounded',success=bool(z.success),nfev=z.nfev)])
 lo,hi=BOUNDS;starts=STARTS
 if model=='DC14_response':
  lo=lo+[0];hi=hi+[1];starts=[x+[v] for x in STARTS for v in [.1,.5,.9]]
  starts+=[base['p']+[.02],dc['p']+[.98]]
 attempts=[]
 for start in starts:
  z=least_squares(resid,start,bounds=(lo,hi),max_nfev=1500,ftol=1e-10,xtol=1e-10,gtol=1e-10)
  attempts.append(dict(method='least_squares',start=start,p=z.x.tolist(),loss=float(2*z.cost),success=bool(z.success),message=z.message,nfev=z.nfev,optimality=float(z.optimality),active_bounds=z.active_mask.tolist()))
 if model=='DC14_response':
  attempts.extend([dict(method='exact_NFW_endpoint',p=base['p']+[0.],loss=loss(base['p']+[0.]),success=True),dict(method='exact_DC14_endpoint',p=dc['p']+[1.],loss=loss(dc['p']+[1.]),success=True)])
 best=min(attempts,key=lambda a:a['loss'])
 if model in ('DC14','DC14_response'):
  z=minimize(loss,best['p'],method='Powell',bounds=list(zip(lo,hi)),options={'maxfev':4000,'xtol':1e-8,'ftol':1e-10})
  attempts.append(dict(method='Powell_polish',start=best['p'],p=z.x.tolist(),loss=float(z.fun),success=bool(z.success),message=z.message,nfev=z.nfev))
 best=min(attempts,key=lambda a:a['loss'])
 return dict(p=best['p'],loss=best['loss'],selected_method=best['method'],attempts=attempts)

def diagnostic(p):
 _,_,h=components(p)
 return dict(X_min=float(min(h['X'])),X_max=float(max(h['X'])),below_recipe=[NAMES[i] for i in np.flatnonzero(h['X']<-4.1)],upper_cap=[NAMES[i] for i in np.flatnonzero(h['X']>-1.3)],outside_simulated_halo_mass=[NAMES[i] for i in np.flatnonzero((h['mh']<9.4e9)|(h['mh']>7.2e11))],rvir_over_r200=[float(min(h['Rv']/h['R'])),float(max(h['Rv']/h['R']))],alpha_range=[float(min(h['a'])),float(max(h['a']))],gamma_range=[float(min(h['g'])),float(max(h['g']))])

def job_run(job):
 name,te,kind=job;te=np.array(te,dtype=int);tr=I if kind=='all' else I[~np.isin(DATA['id'],DATA['id'][te])]
 assert kind=='all' or not set(DATA['id'][tr])&set(DATA['id'][te])
 fits={};out={};start=time.time()
 for model in ['RAR','NFW','DC14','DC14_response','coreNFW_duration']:
  fitr=fit(model,tr,fits.get('NFW'),fits.get('DC14'));fits[model]=fitr;pred=predict(model,fitr['p'])[te]
  out[model]=dict(fit=fitr,metrics=metrics(pred,te),pred=pred.tolist())
  if model.startswith('DC14'):out[model]['diagnostic']=diagnostic(fitr['p'])
 ra=predict('RAR',fits['RAR']['p']);nf=predict('NFW',fits['NFW']['p']);w=weights(tr)
 def obj(x):return float(w@np.log10(((1-x)*ra[tr]+x*nf[tr])/DATA['go'][tr])**2)
 z=minimize_scalar(obj,bounds=(0,1),method='bounded',options={'xatol':1e-12});ww=min([0,1,z.x],key=obj);pred=((1-ww)*ra+ww*nf)[te]
 out['RAR_NFW_blend']=dict(fit=dict(weight=float(ww),loss=obj(ww),base_RAR=fits['RAR']['p'],base_NFW=fits['NFW']['p']),metrics=metrics(pred,te),pred=pred.tolist())
 return dict(name=name,type=kind,train=tr.tolist(),test=te.tolist(),models=out,seconds=time.time()-start)

def validation_jobs():
 jobs=[('all',I.tolist(),'all')]
 for seed in SEEDS:
  tab=pd.read_csv(ROOT.parent/f'dm_cosmological_halo/results/NFW_Planck_size_cv_{seed}.csv')
  assert np.array_equal(tab['row'],I) and np.array_equal(tab.system,DATA['id'])
  for fold in range(5):jobs.append((f'cv_{seed}_{fold}',I[tab.fold.to_numpy()==fold].tolist(),'cv'))
 for what,field in [('mass','mass_proxy_msun'),('surface','sigma_proxy_msun_kpc2')]:
  order=sorted(NAMES,key=lambda g:(float(FEAT.loc[g,field]),g))
  for side,ns in [('low',order[:37]),('high',order[-37:])]:jobs.append((f'{what}_{side}',I[np.isin(DATA['id'],ns)].tolist(),'tail'))
 return jobs

def math_checks():
 errs=[]
 for xx in np.linspace(-4.1,-1.3,13):
  aa,bb,gg,_=shape(np.array([xx]));a,b,g=aa[0],bb[0],gg[0]
  for x in np.logspace(-5,3,13):
   numeric=quad(lambda u:np.exp((3-g)*u-(b-g)/a*np.logaddexp(0,a*u)),np.log(x)-45,np.log(x),epsabs=1e-12,epsrel=2e-11)[0]
   exact=float(mass_integral(x,a,b,g));errs.append(abs(numeric-exact)/exact)
 assert max(errs)<1e-7,max(errs)
 p=[11.4672833579,.9937487126,-.6355007558];mn,md,h=components(p)
 assert np.all(md>0) and np.all(mn>0)
 nfw=predict('NFW',p);nested=predict('DC14_response',p+[0]);assert np.array_equal(nfw,nested)
 conversion=float(max(abs(f(h['cv'])/f(h['c'])/(h['Rv']/h['R'])**3-92.8/200)))
 return dict(integral_cases=len(errs),max_integral_relative_error=max(errs),reference_overdensity_error=conversion,exact_NFW_nested=True)

def summarize(res):
 res=sorted(res,key=lambda r:r['name']);models=list(res[0]['models']);repeats={};tables=[]
 for seed in SEEDS:
  rs=[r for r in res if r['name'].startswith(f'cv_{seed}_')];repeats[str(seed)]={}
  for model in models:
   pr=np.full(len(I),np.nan)
   for r in rs:pr[r['test']]=r['models'][model]['pred']
   assert np.isfinite(pr).all()
   repeats[str(seed)][model]=metrics(pr,I)
   tab=DF[['row','system','group','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();tab['gpred_m_s2']=pr;tab['log10_residual']=np.log10(pr/DATA['go']);tab.to_csv(OUT/f'{model}_cv_{seed}.csv',index=False)
 means={m:{k:float(np.mean([repeats[str(s)][m][k] for s in SEEDS])) for k in ['rms','median_speed']} for m in models}
 tails={r['name']:{m:r['models'][m]['metrics'] for m in models} for r in res if r['type']=='tail'}
 allr=next(r for r in res if r['type']=='all');dd=pd.read_csv(ROOT/'dwarfs_benchmark.csv');d=dict(id=dd.system.to_numpy(),gi=np.arange(len(dd)),r=dd.radius_kpc.to_numpy(),gb=dd.gbar_m_s2.to_numpy(),go=dd.gobs_m_s2.to_numpy(),mp=dd.stellar_mass_Msun.to_numpy(),ms=dd.stellar_mass_Msun.to_numpy(),re=dd.projected_Re_circularized_kpc.to_numpy())
 dwarf={};dtab=dd[['system','host','radius_kpc','gbar_m_s2','gobs_m_s2']].copy()
 for model in models:
  fr=allr['models'][model]['fit']
  if model=='RAR_NFW_blend':pr=(1-fr['weight'])*predict('RAR',fr['base_RAR'],d)+fr['weight']*predict('NFW',fr['base_NFW'],d)
  else:pr=predict(model,fr['p'],d)
  e=np.log10(pr/d['go']);dwarf[model]=dict(rms=float(np.sqrt(np.mean(e**2))),median_acceleration_error=float(np.median(abs(10**e-1))*100));dtab[model]=pr
 dtab.to_csv(OUT/'frozen_dwarf_transfer.csv',index=False)
 # Paired resampling retains a galaxy's predictions across all repeat partitions.
 gm={m:[] for m in models}
 for s in SEEDS:
  for m in models:
   t=pd.read_csv(OUT/f'{m}_cv_{s}.csv');gm[m].append(t.groupby('system').log10_residual.apply(lambda x:np.mean(x*x)).reindex(NAMES).to_numpy())
 rand=np.random.default_rng(20261029);ind=rand.integers(0,len(NAMES),size=(10000,len(NAMES)))
 boot={m:np.sqrt(np.stack(gm[m])[:,ind].mean(axis=2)).mean(axis=0) for m in models};contrasts={}
 for m in models:
  if m=='RAR':continue
  vals=boot[m]-boot['RAR'];contrasts[m]=dict(difference=means[m]['rms']-means['RAR']['rms'],CI95=np.percentile(vals,[2.5,97.5]).tolist())
 # Reproduce original stored NFW coefficients/predictions before comparisons.
 original={};maxpred=0
 for s in SEEDS:
  a=pd.read_csv(ROOT.parent/f'dm_cosmological_halo/results/NFW_Planck_size_cv_{s}.csv');b=pd.read_csv(OUT/f'NFW_cv_{s}.csv');maxpred=max(maxpred,float(max(abs(np.log10(a.gpred_m_s2/b.gpred_m_s2)))))
 original['maximum_absolute_log_prediction_difference']=maxpred
 failures=[dict(split=r['name'],model=m,attempt=a) for r in res for m in MODELS for a in r['models'][m]['fit'].get('attempts',[]) if not a.get('success',True)]
 return dict(means=means,repeats=repeats,tails=tails,dwarfs=dwarf,contrasts=contrasts,baseline_reproduction=original,failed_attempts=failures,all_coefficients={m:allr['models'][m]['fit'] for m in models},all_diagnostics={m:allr['models'][m].get('diagnostic') for m in ['DC14','DC14_response']})

if __name__=='__main__':
 checks=math_checks();(OUT/'math_checks.json').write_text(json.dumps(checks,indent=2));print('MATH',checks,flush=True)
 jobs=validation_jobs();res=[]
 with ProcessPoolExecutor(max_workers=4) as ex:
  pending={ex.submit(job_run,j):j[0] for j in jobs}
  for future in as_completed(pending):
   r=future.result();res.append(r);(OUT/(r['name']+'_record.json')).write_text(json.dumps(r,indent=2));print(r['name'],{m:round(z['metrics']['rms'],6) for m,z in r['models'].items()},'seconds',round(r['seconds'],1),flush=True)
 summary=summarize(res);(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print('SUMMARY',json.dumps(summary['means']),flush=True)
 hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='manifest.json' and '__pycache__' not in str(p)}
 (ROOT/'manifest.json').write_text(json.dumps(hashes,indent=2))
