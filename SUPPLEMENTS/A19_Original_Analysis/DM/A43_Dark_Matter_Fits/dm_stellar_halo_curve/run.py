from pathlib import Path
import sys,os,json,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy.optimize import brentq,least_squares
from concurrent.futures import ProcessPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT.parent/'dm_population_scatter'))
import run as pop
# Loaded as script (__main__), import run resolves the sibling population module.
D=pop.D;NAMES=pop.NAMES;ROWS=pop.ROWS;COUNT=pop.COUNT;MASK=pop.MASK;DF=pop.DF;FB=pop.FB
MODELS=['Moster_fixed','Moster_scaled','stellar_power','stellar_curved']

def moster(logh):
 q=10**(logh-11.59);return 10**logh*2*.0351/(q**(-1.376)+q**.608)
def inverse(ms):return np.array([10**brentq(lambda x:np.log10(moster(x))-np.log10(z),5,18,xtol=1e-12) for z in np.asarray(ms).ravel()]).reshape(np.asarray(ms).shape)
AM=inverse(D['ms'])
def soft(x):return np.logaddexp(0,2*pop.LN*x)/(2*pop.LN)
def mass(model,p,d=D,am=AM):
 floor=d['mp']/FB;x=np.log10(d['ms']/1e10);z=np.log10(d['re']/3)
 if model=='Moster_fixed':return np.maximum(am,floor)
 if model=='Moster_scaled':return floor+10**(p[0]+p[1]*z)*np.maximum(am-floor,0)
 if model=='stellar_power':e=p[0]+p[1]*x+p[2]*z
 else:e=p[0]+p[1]*x+p[2]*(soft(x-p[4])-soft(-p[4]))+p[3]*z
 return floor+10**e

def predict(model,p,d=D,am=AM):
 mr=mass(model,p,d,am);c=10**(.905-.101*np.log10(mr*pop.H/1e12));R=np.cbrt(pop.G*mr/(100*pop.H0**2));x=np.minimum(d['r'],R)*c/R
 return np.log10(d['gb']+pop.G*(1-FB)*mr*pop.f(x)/pop.f(c)/d['r']**2*pop.UNIT)

def fit(model,tr,previous=None):
 if model=='Moster_fixed':return dict(p=[],attempts=[],loss=float(np.mean(np.sum(((predict(model,[])-np.log10(D['go']))*MASK)**2,axis=1)[tr]/COUNT[tr])))
 if model=='Moster_scaled':bounds=[(-2,2),(-1.5,1.5)];starts=[[0,0],[-.5,-.5],[.5,.5]]
 elif model=='stellar_power':bounds=[(8,14),(.1,2),(-1.5,1.5)];starts=[[11.5,.7,0],[12,1.2,.5],[10.5,.4,-.5]]
 else:
  bounds=[(8,14),(.1,2),(0,2),(-1.5,1.5),(-3,1)];starts=[[11.5,.5,.5,-.5,-.5],[11.5,.4,1,0,0],[12,.8,.2,.5,-1.5]]
  p=previous['p'];starts.append([p[0],p[1],.01,p[2],-1])
 def res(p):return (((predict(model,p)-np.log10(D['go']))*MASK/np.sqrt(COUNT[:,None]))[tr]/np.sqrt(len(tr))).ravel()
 aa=[]
 for st in starts:
  z=least_squares(res,st,bounds=np.array(bounds).T,max_nfev=2000,ftol=1e-10,xtol=1e-10,gtol=1e-10);aa.append(dict(start=st,p=z.x.tolist(),loss=float(2*z.cost),success=bool(z.success),nfev=int(z.nfev),message=z.message,active_bounds=z.active_mask.tolist()))
 if model=='stellar_curved':
  q=[p[0],p[1],0.,p[2],-1.];aa.append(dict(p=q,loss=float(np.sum(res(q)**2)),success=True,method='exact_power_endpoint'))
 if model=='Moster_scaled':aa.append(dict(p=[0.,0.],loss=float(np.sum(res([0.,0.])**2)),success=True,method='exact_published_endpoint'))
 best=min(aa,key=lambda z:z['loss']);return dict(p=best['p'],loss=best['loss'],selected_success=best['success'],attempts=aa)

def job(path):
 old=json.loads(Path(path).read_text());tr=sorted({pop.GI[g] for g in DF.system.iloc[old['train']]});te=sorted({pop.GI[g] for g in DF.system.iloc[old['test']]});assert old['type']=='all' or not set(tr)&set(te)
 fits={};preds={};st=time.time()
 for m in MODELS:
  fits[m]=fit(m,tr,fits.get('stellar_power'));preds[m]=predict(m,fits[m]['p'])[te].tolist()
 return dict(name=old['name'],type=old['type'],train_gi=tr,test_gi=te,fits=fits,predictions=preds,seconds=time.time()-st)

def summarize(rs):
 cv={};mse={m:[] for m in MODELS}
 for seed in [20261027,7331,48105]:
  cv[str(seed)]={}
  for m in MODELS:
   pr=np.full(MASK.shape,np.nan)
   for z in rs:
    if z['name'].startswith(f'cv_{seed}_'):pr[z['test_gi']]=z['predictions'][m]
   assert np.isfinite(pr).all();metric,ee=pop.evaluate(pr);cv[str(seed)][m]=metric;mse[m].append(ee)
   t=DF[['row','system','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();t['gpred_m_s2']=10**pop.flatten(pr);t['log10_residual']=pop.flatten(pr-np.log10(D['go']));t.to_csv(OUT/f'{m}_cv_{seed}.csv',index=False)
 means={m:{k:float(np.mean([cv[str(seed)][m][k] for seed in [20261027,7331,48105]])) for k in ['rms','median_speed']} for m in MODELS}
 tails={}
 for z in rs:
  if z['type']!='tail':continue
  tails[z['name']]={}
  for m in MODELS:
   j=z['test_gi'];e=(np.array(z['predictions'][m])-np.log10(D['go'][j]))*MASK[j];tails[z['name']][m]=float(np.sqrt(np.mean(np.sum(e*e,axis=1)/COUNT[j])))
 full=next(z for z in rs if z['type']=='all');dd=pd.read_csv(pop.OLD/'dwarfs_benchmark.csv');d=dict(r=dd.radius_kpc.to_numpy()[:,None],gb=dd.gbar_m_s2.to_numpy()[:,None],go=dd.gobs_m_s2.to_numpy()[:,None],mp=dd.stellar_mass_Msun.to_numpy()[:,None],ms=dd.stellar_mass_Msun.to_numpy()[:,None],re=dd.projected_Re_circularized_kpc.to_numpy()[:,None]);am=inverse(d['ms']);dt=dd[['system','host','radius_kpc','gbar_m_s2','gobs_m_s2']].copy();dwarfs={}
 for m in MODELS:
  pr=predict(m,full['fits'][m]['p'],d,am);e=pr-np.log10(d['go']);dwarfs[m]=float(np.sqrt(np.mean(e*e)));dt[m]=10**pr[:,0]
 dt.to_csv(OUT/'frozen_dwarf_transfer.csv',index=False)
 for alias,name in [('RAR','RAR'),('NFW','NFW'),('Blend','RAR_NFW_blend')]:
  mse[alias]=[]
  for seed in [20261027,7331,48105]:
   t=pd.read_csv(pop.OLD/f'results/{name}_cv_{seed}.csv');mse[alias].append(t.groupby('system').log10_residual.apply(lambda x:np.mean(x*x)).reindex(NAMES).to_numpy())
 ix=np.random.default_rng(20261029).integers(0,len(NAMES),(10000,len(NAMES)));boot={m:np.sqrt(np.array(v)[:,ix].mean(2)).mean(0) for m,v in mse.items()};contrasts={}
 for m in MODELS:
  for b in ['RAR','NFW','Blend']:
   diff=boot[m]-boot[b];contrasts[m+' minus '+b]=dict(difference=float(np.sqrt(np.array(mse[m]).mean(1)).mean()-np.sqrt(np.array(mse[b]).mean(1)).mean()),CI95=np.percentile(diff,[2.5,97.5]).tolist())
 return dict(means=means,repeats=cv,tails=tails,dwarfs=dwarfs,contrasts=contrasts,full_coefficients={m:full['fits'][m]['p'] for m in MODELS},failed_selected=[(z['name'],m) for z in rs for m in MODELS if not z['fits'][m].get('selected_success',True)],optimizer_attempts=sum(len(z['fits'][m]['attempts']) for z in rs for m in MODELS))

if __name__=='__main__':
 assert np.max(abs(np.log10(moster(np.log10(AM))/D['ms'])))<1e-10
 assert np.allclose(predict('Moster_fixed',[]),predict('Moster_scaled',[0,0]),rtol=0,atol=1e-13)
 rs=[]
 with ProcessPoolExecutor(max_workers=4) as ex:
  pending=[ex.submit(job,str(p)) for p in sorted((pop.OLD/'results').glob('*_record.json'))]
  for fu in as_completed(pending):
   z=fu.result();rs.append(z);(OUT/(z['name']+'_record.json')).write_text(json.dumps(z,indent=2));print(z['name'],round(z['seconds'],1),flush=True)
 s=summarize(rs);(OUT/'summary.json').write_text(json.dumps(s,indent=2));print('SUMMARY',json.dumps(s['means']),flush=True)
