from pathlib import Path
import csv,json,hashlib,time
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss
R=Path(__file__).resolve().parent;O=R/'results';O.mkdir(exist_ok=True)
HALF=1/np.sqrt(2**(2/3)-1);Q,W=leggauss(64);rr=HALF*(Q+1)/2;ww=W*HALF/2

def shifted(d):
 if d<1e-7:return .5
 A=(rr-d)**2+1
 return float(np.sum(ww*.5*rr/d*A**-1.5*(-np.expm1(-1.5*np.log1p(4*rr*d/A)))))

def Hdirect(z,kind):
 points=sorted(set([0.,80.]+[float(v) for v in [z*.01,z*.1,z,z*10,z*100,1.,5.,20.] if 0<v<80]))
 f=lambda t:(.5*(1+t) if kind=='NLG' else t)*np.exp(-t)*shifted(t/z)
 return sum(quad(f,lo,hi,epsabs=1e-16*min(1,z if kind=='NLG' else z*z),epsrel=2e-11,limit=150)[0] for lo,hi in zip(points[:-1],points[1:]))

def build_basis():
 grid=np.linspace(-8,8,401);out=[]
 for l in grid:out.append({'log10z':l,'NLG_H':Hdirect(10**l,'NLG'),'MOG_H':Hdirect(10**l,'MOG')})
 with (R/'basis.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 return out

def writecsv(path,rows):
 with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
FITLOG=[]
def main():
 basis=build_basis();H={k:PchipInterpolator([r['log10z'] for r in basis],np.log10([r[k+'_H'] for r in basis])) for k in ['NLG','MOG']}
 raw=list(csv.DictReader((R.parent/'dm_unified_regular/augmented_input.csv').open()));rows=[dict(row=i,**r) for i,r in enumerate(raw) if r['population']=='dSph'];assert len(rows)==44
 writecsv(R/'input.csv',rows);n=len(rows);I=np.arange(n);gb=np.array([float(r['gbar_m_s2']) for r in rows]);go=np.array([float(r['gobs_m_s2']) for r in rows]);rad=np.array([float(r['radius_kpc']) for r in rows]);a=rad/HALF;group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows]);y=np.log10(go);bound=[]
 for ml in [1,2,3]:
  D=go/(gb*ml/2)
  for i in I:bound.append({'system':rows[i]['system'],'stellar_ML':ml,'gobs_over_gbar':D[i],'necessary_alpha_any_spherical_profile':max(0.,(D[i]-1)/2),'necessary_alpha_Plummer':max(0.,D[i]-1),'exceeds_NLG_general_bound':bool(D[i]>1+2*10.94),'exceeds_MOG_general_bound':bool(D[i]>1+2*8.89)})
 writecsv(R/'mass_budget_bounds.csv',bound)
 models=['NLG_frozen','NLG_fit','MOG_frozen','MOG_fit','RAR_fixed','RAR_fit'];out={'protocol_sha256':sha(R/'PROTOCOL.txt'),'script_sha256':sha(Path(__file__)),'input_sha256':sha(R/'input.csv'),'basis_sha256':sha(R/'basis.csv'),'source_sha256':sha(R.parent/'dm_unified_regular/augmented_input.csv'),'fits':FITLOG,'models_by_ML':{},'basis_halfmass_radius_over_a':HALF}
 def pred(k,c,idx,ml):
  b=gb[idx]*ml/2
  if k.startswith('RAR'):
   aa=1.2e-10 if k=='RAR_fixed' else 1e-10*10**c[0];return np.log10(b/(-np.expm1(-np.sqrt(b/aa))))
  kind=k.split('_')[0]
  if k.endswith('frozen'):alpha,mu={'NLG':(10.94,.059),'MOG':(8.89,.042)}[kind]
  else:alpha,mu=10**c[0],10**c[1]
  z=np.log10(mu*a[idx]);assert np.all((z>=-8)&(z<=8)),z
  return np.log10(b*(1+2*alpha*10**H[kind](z)))
 def fit(k,idx,tag,ml):
  if not k.endswith('_fit'):return []
  starts=[[0],[1]] if k=='RAR_fit' else [[1,-1.2],[2,0],[3,-2],[1,2]];lo=[-3] if k=='RAR_fit' else [-3,-5];hi=[3] if k=='RAR_fit' else [5,5];res=[]
  for si,st in enumerate(starts):
   fit=least_squares(lambda c:(pred(k,c,idx,ml)-y[idx])/np.sqrt(len(idx)),st,bounds=(lo,hi),max_nfev=2000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
   FITLOG.append({'model':k,'stellar_ML':ml,'split':tag,'start':si,'start_coefficients':st,'lo':lo,'hi':hi,'training_rows':idx.tolist(),'coefficients':fit.x.tolist(),'cost':float(fit.cost),'success':bool(fit.success),'message':fit.message,'active_bounds':fit.active_mask.tolist(),'nfev':fit.nfev});res.append(fit)
  return min(res,key=lambda x:x.cost).x.tolist()
 def metrics(p,idx):
  e=p-y[idx];return {'RMS_dex':float(np.sqrt(np.mean(e*e))),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'rows':len(idx)}
 def save(k,p,idx,fold,ml,tag):
  records=[]
  for j,i in enumerate(idx):records.append({'input_row':int(i),'broad_row':rows[i]['row'],'system':rows[i]['system'],'group':group[i],'environment':env[i],'model':k,'stellar_ML':ml,'split':tag,'fold':fold[j],'radius_kpc':rad[i],'gbar_m_s2':gb[i]*ml/2,'gobs_m_s2':go[i],'gpred_m_s2':10**p[j],'log10_residual':p[j]-y[i]})
  savepath=O/f'ML{ml}_{k}_{tag}.csv';writecsv(savepath,records)
 for ml in [1,2,3]:
  ds={};out['models_by_ML'][str(ml)]=ds
  for k in models:
   cc=fit(k,I,'all',ml);pr=pred(k,cc,I,ml);d={'coefficients':cc,'all':metrics(pr,I),'CV':{},'host':{}};save(k,pr,I,['all']*n,ml,'all')
   for seed in [20261027,7331,48105]:
    rng=np.random.default_rng(seed);mp={g:j%5 for j,g in enumerate(rng.permutation(sorted(set(group))))};fd=np.array([mp[g] for g in group]);p=np.full(n,np.nan);cs=[]
    for f in range(5):
     tr=I[fd!=f];te=I[fd==f];c=fit(k,tr,f'cv_{seed}_{f}',ml);p[te]=pred(k,c,te,ml);cs.append(c)
    d['CV'][str(seed)]={'coefficients':cs,'metrics':metrics(p,I)};save(k,p,I,fd.tolist(),ml,f'cv_{seed}')
   for host,aliases in [('MW',['Milky_Way','G/L']),('M31',['M31','A/L'])]:
    te=I[np.isin(env,aliases)];tr=I[~np.isin(env,aliases)];assert len(te)=={'MW':26,'M31':18}[host];c=fit(k,tr,'host_'+host,ml);p=pred(k,c,te,ml);d['host'][host]={'coefficients':c,'training_rows':tr.tolist(),'metrics':metrics(p,te)};save(k,p,te,[host]*len(te),ml,'host_'+host)
   ds[k]=d;print(ml,k,'CV',d['CV']['20261027']['metrics']['RMS_dex'],'MW',d['host']['MW']['metrics']['RMS_dex'],'M31',d['host']['M31']['metrics']['RMS_dex'],'coeff',cc,flush=True)
 out['fit_sets']=sum(x['start']==0 for x in FITLOG);out['optimizer_attempts']=len(FITLOG);out['optimizer_failures']=[x for x in FITLOG if not x['success']];out['prediction_tables']=len(list(O.glob('*.csv')))
 out['bound_summary']={str(ml):{'NLG_general_violations':sum(r['stellar_ML']==ml and r['exceeds_NLG_general_bound'] for r in bound),'MOG_general_violations':sum(r['stellar_ML']==ml and r['exceeds_MOG_general_bound'] for r in bound),'max_alpha_required':max(r['necessary_alpha_any_spherical_profile'] for r in bound if r['stellar_ML']==ml)} for ml in [1,2,3]}
 (R/'results.json').write_text(json.dumps(out,indent=2));print('DONE',out['fit_sets'],out['optimizer_attempts'],out['prediction_tables'],out['bound_summary'],flush=True)
if __name__=='__main__':main()
