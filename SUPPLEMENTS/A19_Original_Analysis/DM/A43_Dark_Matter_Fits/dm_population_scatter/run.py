"""A43: finite concentration response and linearized population prediction."""
from pathlib import Path
import os,io,json,zipfile,hashlib,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy.optimize import minimize,least_squares
from concurrent.futures import ProcessPoolExecutor,as_completed

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
OLD=ROOT.parent/'dm_feedback_profiles'
DF=pd.read_csv(OLD/'input.csv');FEAT=pd.read_csv(OLD/'galaxy_features.csv').set_index('system')
NAMES=sorted(DF.system.unique());NG=len(NAMES);GI={s:i for i,s in enumerate(NAMES)}
LN=np.log(10.);G=4.30091727003628e-6;UNIT=1e6/3.085677581491367e19;FB=.049/.3175;H0=.0671;H=.671
NMAX=int(DF.groupby('system').size().max());COUNT=np.array([sum(DF.system==g) for g in NAMES]);MASK=np.arange(NMAX)[None,:]<COUNT[:,None]
ROWS=[np.flatnonzero(DF.system.to_numpy()==g) for g in NAMES]

def data():
 d={key:np.ones((NG,NMAX)) for key in ['r','gb','go','star','ev','v']}
 cats={}
 for line in (ROOT.parent/'dm_fit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
  a=line.split()
  if len(a)==19 and a[1].isdigit():cats[a[0]]=[float(a[k]) for k in [2,3,5,6]]
 props=[];rr=[]
 with zipfile.ZipFile(ROOT.parent/'dm_fit/data/Rotmod_LTG.zip') as z:
  for j,name in enumerate(NAMES):
   raw=np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat')),ndmin=2)
   for k,row in enumerate(ROWS[j]):
    q=DF.iloc[row];matches=np.flatnonzero(np.isclose(raw[:,0],q.radius_kpc,rtol=0,atol=1e-9));assert len(matches)==1
    rad,v,ev,gas,disk,bulge,*_=raw[matches[0]];gb=(gas*abs(gas)+.5*disk*abs(disk)+.7*bulge*abs(bulge))/rad*UNIT;go=v*v/rad*UNIT
    assert np.isclose(gb,q.gbar_m_s2,rtol=1e-12) and np.isclose(go,q.gobs_m_s2,rtol=1e-12)
    d['r'][j,k]=rad;d['gb'][j,k]=gb;d['go'][j,k]=go;d['star'][j,k]=(.5*disk*abs(disk)+.7*bulge*abs(bulge))/rad*UNIT;d['v'][j,k]=v;d['ev'][j,k]=ev
    rr.append(dict(row=int(row),system=name,radius_kpc=rad,vobs_km_s=v,errv_km_s=ev,gstar_m_s2=d['star'][j,k],distance_Mpc=cats[name][0],distance_error_Mpc=cats[name][1],inclination_deg=cats[name][2],inclination_error_deg=cats[name][3]))
   props.append(cats[name])
 props=np.array(props);d['mp']=FEAT.loc[NAMES,'mass_proxy_msun'].to_numpy()[:,None];d['ms']=(.5e9*FEAT.loc[NAMES,'luminosity_3p6_1e9_lsun']).to_numpy()[:,None];d['re']=FEAT.loc[NAMES,'reff_kpc'].to_numpy()[:,None]
 d['sd']=props[:,1:2]/props[:,0:1]/LN;d['si']=2/np.tan(np.deg2rad(props[:,2:3]))*np.deg2rad(props[:,3:4])/LN
 d['sv']=2*d['ev']/d['v']/LN;d['mask']=MASK;d['count']=COUNT
 return d,pd.DataFrame(rr).sort_values('row')
D,RAW=data()

def f(x):
 out=np.log1p(x)-x/(1+x);small=x<.01
 if np.any(small):
  out[small]=sum((-1)**k*(k-1)/k*x[small]**k for k in range(2,13))
 return out

def nfw(p,d=D,delta_m=0.,delta_c=0.):
 A,B,C=p[:3];K=p[3] if len(p)>3 else 0.
 floor=d['mp']/FB;excess=10**A*(d['mp']/1e10)**B*(d['re']/3)**C*10**np.asarray(delta_m)
 mr=floor+excess;c=10**(.905-.101*np.log10(mr*H/1e12)+K+np.asarray(delta_c));R=np.cbrt(G*mr/(100*H0**2));x=np.minimum(d['r'],R)*c/R
 fc=f(c);fx=f(x);gh=G*(1-FB)*mr*fx/fc/d['r']**2*UNIT;total=d['gb']+gh
 hx=x*x/(1+x)**2/fx;hc=c*c/(1+c)**2/fc;inside=d['r']<R
 hm=np.where(inside,1-(.101+1/3)*hx+.101*hc,1.);hconc=np.where(inside,hx-hc,0.);hr=np.where(inside,hx,0.)-2
 frac=gh/total;dm=excess/mr
 JM=frac*hm*dm;JC=frac*hconc
 JD=1+frac*(hm*(2*floor+(2*B+C)*excess)/mr+hr)
 JY=d['star']/total+frac*hm*(floor+B*excess)/mr*d['ms']/d['mp']
 return dict(log=np.log10(total),JD=JD,JY=JY,JM=JM,JC=JC,mr=mr,c=c,R=R)

def rar(a,d=D):
 z=np.sqrt(d['gb']/(1e-10*10**a));den=-np.expm1(-z);t=d['gb']/den
 # z/expm1(z) is stable when z is large.
 rat=np.divide(z,np.expm1(np.minimum(z,700)))
 q=1-.5*rat
 return dict(log=np.log10(t),JD=np.ones_like(t),JY=q*d['star']/d['gb'])

MODELS=['RAR_nuisance','RAR_offset','NFW_c','NFW_mc','NFW_response_mc','Blend_offset']

def central(model,p,base=None,d=D):
 if model.startswith('RAR'):return rar(p[0],d)
 if model.startswith('NFW'):return nfw(p[:4] if model=='NFW_response_mc' else p[:3],d)
 r=rar(base['base_RAR'][0],d);n=nfw(base['base_NFW'],d);w=base['weight'];gr=10**r['log'];gn=10**n['log'];gg=(1-w)*gr+w*gn
 wn=w*gn/gg;wr=(1-w)*gr/gg
 return dict(log=np.log10(gg),JD=wr*r['JD']+wn*n['JD'],JY=wr*r['JY']+wn*n['JY'])

def noise(model,p):
 if model=='RAR_nuisance':return 10**p[1],0.,0.
 if model=='RAR_offset':return 10**p[1],p[2],0.
 if model=='NFW_c':return 10**p[3],0.,0.
 if model=='NFW_mc':return 10**p[3],0.,p[4]
 if model=='NFW_response_mc':return 10**p[4],0.,p[5]
 return 10**p[0],p[1],0.

def distribution(model,p,base=None,d=D):
 c=central(model,p,base,d);s,tau,sm=noise(model,p)
 cols=[c['JD']*d['sd'],np.broadcast_to(d['si'],c['log'].shape),.1*c['JY']]
 if model.startswith('NFW'):
  cols.append(.11*c['JC'])
  if model!='NFW_c':cols.append(sm*c['JM'])
 elif model in ['RAR_offset','Blend_offset']:cols.append(np.ones_like(c['log'])*tau)
 U=np.stack(cols,axis=-1)*d['mask'][...,None]
 v=np.where(d['mask'],d['sv']**2+s*s,1.)
 return c,v,U

def score_arrays(model,p,base=None,d=D):
 c,v,U=distribution(model,p,base,d);e=(c['log']-np.log10(d['go']))*d['mask'];iv=1/v
 mat=np.eye(U.shape[2])[None,:,:]+np.einsum('gnk,gn,gnl->gkl',U,iv,U)
 vec=np.einsum('gnk,gn->gk',U,e*iv)
 solution=np.linalg.solve(mat,vec[...,None])[...,0];sign,ld=np.linalg.slogdet(mat);assert np.all(sign>0)
 q=(e*e*iv).sum(1)-(vec*solution).sum(1)
 nll=.5*(d['count']*np.log(2*np.pi)+np.log(v).sum(1)+ld+q)/d['count']
 sd=np.sqrt(v+(U*U).sum(2));mask=d['mask'];cnt=d['count']
 return dict(log=c['log'],e=e,sd=sd,nll=nll,mse=(e*e).sum(1)/cnt,cover68=((abs(e)<=sd)*mask).sum(1)/cnt,cover95=((abs(e)<=1.95996398454*sd)*mask).sum(1)/cnt)

def spec(model):
 ab=[(8,14),(.1,2),(-1.5,1.5)];ss=(-2.30103,-.30103)
 ABC=[[11.5,.7,0],[12,1.2,.5],[10.5,.4,-.5]]
 if model=='RAR_nuisance':return [(-3,3),ss],[[0,-1.2],[.2,-1],[-.2,-.8]]
 if model=='RAR_offset':return [(-3,3),ss,(0,.5)],[[0,-1.2,.1],[.2,-1,.2],[-.2,-.8,.05]]
 if model=='NFW_c':return ab+[ss],[a+[-1.2] for a in ABC]
 if model=='NFW_mc':return ab+[ss,(0,.6)],[a+[-1.2,v] for a,v in zip(ABC,[.15,.35,.5])]
 if model=='NFW_response_mc':return ab+[(-.5,.5),ss,(0,.6)],[a+[k,-1.2,v] for a,k,v in zip(ABC,[0,-.15,.15],[.15,.35,.5])]
 return [ss,(0,.5)],[[-1.2,.1],[-1,.2],[-.8,.05]]

def fit(model,train,base):
 bounds,starts=spec(model)
 def objective(p):return float(score_arrays(model,p,base)['nll'][train].mean())
 attempts=[]
 for start in starts:
  z=minimize(objective,start,method='L-BFGS-B',bounds=bounds,options={'maxiter':1000,'ftol':2e-12,'gtol':2e-6,'maxls':40})
  attempts.append(dict(start=start,p=z.x.tolist(),loss=float(z.fun),success=bool(z.success),message=str(z.message),nit=int(z.nit),nfev=int(z.nfev),method='L-BFGS-B'))
 best=min(attempts,key=lambda x:x['loss'])
 if not best['success']:
  z=minimize(objective,best['p'],method='Powell',bounds=bounds,options={'maxiter':1000,'ftol':1e-10,'xtol':1e-7})
  attempts.append(dict(start=best['p'],p=z.x.tolist(),loss=float(z.fun),success=bool(z.success),message=str(z.message),nit=int(z.nit),nfev=int(z.nfev),method='Powell'))
  best=min(attempts,key=lambda x:x['loss'])
 return dict(p=best['p'],loss=best['loss'],selected_success=best['success'],attempts=attempts)

def point_fit(train,old):
 bounds=[(8,14),(.1,2),(-1.5,1.5),(-.5,.5)]
 def residual(p):
  e=(nfw(p)['log']-np.log10(D['go']))*MASK/np.sqrt(COUNT[:,None]);return e[train].ravel()/np.sqrt(len(train))
 attempts=[]
 for start in [[11.5,.7,0,0],[12,1.2,.5,-.15],[10.5,.4,-.5,.15]]:
  z=least_squares(residual,start,bounds=np.array(bounds).T,max_nfev=1500,ftol=1e-11,xtol=1e-11,gtol=1e-10)
  attempts.append(dict(start=start,p=z.x.tolist(),loss=float(2*z.cost),success=bool(z.success),nfev=z.nfev,message=z.message))
 p=old['models']['NFW']['fit']['p']+[0.];attempts.append(dict(p=p,loss=float(np.sum(residual(p)**2)),success=True,method='exact_old_endpoint'))
 best=min(attempts,key=lambda z:z['loss']);return dict(p=best['p'],loss=best['loss'],selected_success=best['success'],attempts=attempts)

def job(path):
 old=json.loads(Path(path).read_text());name=old['name'];train=sorted({GI[g] for g in DF.system.iloc[old['train']]});test=sorted({GI[g] for g in DF.system.iloc[old['test']]})
 assert old['type']=='all' or not set(train)&set(test)
 base=old['models']['RAR_NFW_blend']['fit'];st=time.time();fits={};preds={}
 pf=point_fit(train,old);fits['NFW_response_MSE']=pf;pp=nfw(pf['p'])['log'];preds['NFW_response_MSE']=dict(log=pp[test].tolist())
 for model in MODELS:
  fits[model]=fit(model,train,base);v=score_arrays(model,fits[model]['p'],base)
  preds[model]={k:np.asarray(val)[test].tolist() for k,val in v.items() if k!='e'}
 out=dict(name=name,type=old['type'],train_names=[NAMES[j] for j in train],test_names=[NAMES[j] for j in test],test_gi=test,base=base,fits=fits,predictions=preds,seconds=time.time()-st)
 return out

def flatten(a):
 b=np.zeros(len(DF))
 for j,ii in enumerate(ROWS):b[ii]=a[j,:len(ii)]
 return b

def evaluate(log,sd=None,nll=None):
 e=(log-np.log10(D['go']))*MASK;mse=(e*e).sum(1)/COUNT
 out=dict(rms=float(np.sqrt(mse.mean())),median_speed=float(np.median(abs(10**(.5*flatten(e))-1))*100))
 if sd is not None:out.update(mean_nll=float(nll.mean()),coverage68=float((((abs(e)<=sd)*MASK).sum(1)/COUNT).mean()),coverage95=float((((abs(e)<=1.95996398454*sd)*MASK).sum(1)/COUNT).mean()))
 return out,mse

def summarize(results):
 models=['NFW_response_MSE']+MODELS;cv={};mse={m:[] for m in models};nll={m:[] for m in MODELS}
 for seed in [20261027,7331,48105]:
  cv[str(seed)]={}
  rr=[r for r in results if r['name'].startswith(f'cv_{seed}_')]
  for model in models:
   logs=np.full(MASK.shape,np.nan);sds=logs.copy();scores=np.full(NG,np.nan)
   for r in rr:
    j=r['test_gi'];z=r['predictions'][model];logs[j]=z['log']
    if model in MODELS:sds[j]=z['sd'];scores[j]=z['nll']
   assert np.all(np.isfinite(logs))
   metric,ee=evaluate(logs,sds if model in MODELS else None,scores);cv[str(seed)][model]=metric;mse[model].append(ee)
   t=DF[['row','system','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();t['gpred_m_s2']=10**flatten(logs);t['log10_residual']=flatten(logs-np.log10(D['go']))
   if model in MODELS:
    nll[model].append(scores);t['predictive_sd_dex']=flatten(sds);t['galaxy_mean_joint_nll']=t.system.map({g:scores[j] for j,g in enumerate(NAMES)})
   t.to_csv(OUT/f'{model}_cv_{seed}.csv',index=False)
 means={m:{k:float(np.mean([cv[str(s)][m][k] for s in [20261027,7331,48105]])) for k in cv['20261027'][m]} for m in models}
 tails={}
 for r in results:
  if r['type']!='tail':continue
  tails[r['name']]={}
  for m in models:
   z=r['predictions'][m];j=r['test_gi'];e=(np.array(z['log'])-np.log10(D['go'][j]))*MASK[j];ee=(e*e).sum(1)/COUNT[j]
   tails[r['name']][m]=dict(rms=float(np.sqrt(ee.mean())))
   if m in MODELS:tails[r['name']][m]['nll']=float(np.mean(z['nll']))
 # Frozen central predictions only on dwarfs; disk-specific noise not transferred.
 full=next(r for r in results if r['type']=='all');dd=pd.read_csv(OLD/'dwarfs_benchmark.csv');
 d=dict(r=dd.radius_kpc.to_numpy()[:,None],gb=dd.gbar_m_s2.to_numpy()[:,None],go=dd.gobs_m_s2.to_numpy()[:,None],star=dd.gbar_m_s2.to_numpy()[:,None],mp=dd.stellar_mass_Msun.to_numpy()[:,None],ms=dd.stellar_mass_Msun.to_numpy()[:,None],re=dd.projected_Re_circularized_kpc.to_numpy()[:,None])
 dwarfs={};dt=dd[['system','host','radius_kpc','gbar_m_s2','gobs_m_s2']].copy()
 for m in models:
  p=full['fits'][m]['p'];a=nfw(p,d)['log'] if m=='NFW_response_MSE' else central(m,p,full['base'],d)['log'];err=a-np.log10(d['go'])
  dwarfs[m]=dict(rms=float(np.sqrt(np.mean(err*err))));dt[m]=10**a[:,0]
 dt.to_csv(OUT/'frozen_dwarf_transfer.csv',index=False)
 # Add archived point baselines for directly comparable central-prediction contrasts.
 for m,oldname in [('RAR_MSE','RAR'),('NFW_MSE','NFW'),('Blend_MSE','RAR_NFW_blend')]:
  mse[m]=[]
  for seed in [20261027,7331,48105]:
   t=pd.read_csv(OLD/f'results/{oldname}_cv_{seed}.csv');mse[m].append(t.groupby('system').log10_residual.apply(lambda x:np.mean(x*x)).reindex(NAMES).to_numpy())
 rand=np.random.default_rng(20261029);ix=rand.integers(0,NG,(10000,NG));bs={m:np.sqrt(np.array(v)[:,ix].mean(2)).mean(0) for m,v in mse.items()}
 boot={}
 for m in models:
  for baseline in ['RAR_MSE','Blend_MSE']:
   x=bs[m]-bs[baseline];point=np.mean(np.sqrt(np.mean(mse[m],axis=1)))-np.mean(np.sqrt(np.mean(mse[baseline],axis=1)))
   boot[m+' minus '+baseline]=dict(difference=float(point),CI95=np.percentile(x,[2.5,97.5]).tolist())
 likelihood={}
 for m in MODELS:
  if m=='RAR_offset':continue
  diff=np.array(nll[m])-np.array(nll['RAR_offset']);x=diff[:,ix].mean((0,2));likelihood[m+' minus RAR_offset']=dict(difference=float(diff.mean()),CI95=np.percentile(x,[2.5,97.5]).tolist())
 return dict(means=means,repeats=cv,tails=tails,dwarfs=dwarfs,RMS_contrasts=boot,logscore_contrasts=likelihood,full_coefficients={m:full['fits'][m]['p'] for m in models},failed_selected=[(r['name'],m) for r in results for m in models if not r['fits'][m]['selected_success']],optimizer_attempts=sum(len(r['fits'][m]['attempts']) for r in results for m in models))

if __name__=='__main__':
 RAW.to_csv(ROOT/'matched_observations.csv',index=False)
 paths=sorted((OLD/'results').glob('*_record.json'));assert len(paths)==20
 res=[]
 with ProcessPoolExecutor(max_workers=4) as ex:
  futures={ex.submit(job,str(p)):p.stem for p in paths}
  for fu in as_completed(futures):
   z=fu.result();res.append(z);(OUT/(z['name']+'_record.json')).write_text(json.dumps(z,indent=2));print(z['name'],'seconds',round(z['seconds'],1),'K_point',round(z['fits']['NFW_response_MSE']['p'][3],4),flush=True)
 s=summarize(res);(OUT/'summary.json').write_text(json.dumps(s,indent=2));print('SUMMARY',json.dumps(s['means']),flush=True)
