from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import csv,json,sys,time,hashlib
R=Path(__file__).resolve().parent;O=R/'results';O.mkdir(exist_ok=True);rng=np.random.default_rng(20261028)
CALLS=0;WARN=[]
def weights(g,idx,hosts=None):
 w=np.zeros(len(idx));cats=[None] if hosts is None else np.unique(hosts[idx])
 for h in cats:
  m=np.ones(len(idx),bool) if h is None else hosts[idx]==h
  for gg in np.unique(g[idx][m]):
   mm=m&(g[idx]==gg);w[mm]=1/(len(cats)*len(np.unique(g[idx][m]))*sum(mm))
 return w

def rar(g,a):return g/(-np.expm1(-np.sqrt(g/a)))
def simple(g,a,e):
 # Stable difference F(g+e)-F(e), F(z)=z/2+sqrt(z*z/4+a*z).
 A=np.sqrt(.25*(g+e)**2+a*(g+e));B=np.sqrt(.25*e**2+a*e)
 return g*(.5+(a+.25*g+.5*e)/(A+B))
def fit(fun,n,idx,w):
 global CALLS;CALLS+=1
 res=least_squares(lambda c:np.sqrt(w)*(fun(c,idx)-Y[idx]),np.zeros(n),bounds=([-3.]*n,[3.]*n),max_nfev=1500,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 if not res.success:WARN.append(res.message)
 return res.x

def metrics(pr,g,hosts=None):
 idx=np.arange(len(pr));d=pr-Y
 return {'RMS_dex':float(np.sqrt(np.sum(weights(g,idx,hosts)*d*d))),'median_acceleration_error_percent':float(np.median(abs(10**d-1))*100),'median_speed_error_percent':float(np.median(abs(10**(d/2)-1))*100)}
def evaluate(fun,n,g,fold,hosts=None):
 pr=np.zeros(len(g));coef={}
 for v in np.unique(fold):
  tr=np.where(fold!=v)[0];te=np.where(fold==v)[0];assert not set(g[tr])&set(g[te]);c=fit(fun,n,tr,weights(g,tr,hosts));pr[te]=fun(c,te);coef[str(v)]=c.tolist()
 return pr,{'coefficients':coef,'metrics':metrics(pr,g,hosts)}
def save(name,rows,pr,fold):
 data=[]
 for i,r in enumerate(rows):data.append(dict(r,fold=str(fold[i]),gpred_m_s2=10**pr[i],log10_residual=pr[i]-Y[i]))
 with (O/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

J={'protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'environment':{},'history':{}}
rows=list(csv.DictReader((R/'sparc_environment_input.csv').open()));g=np.array([r['system'] for r in rows]);G=np.array([float(r['gbar_m_s2']) for r in rows]);Y=np.log10([float(r['gobs_m_s2']) for r in rows]);D=np.array([float(r['distance_Mpc']) for r in rows]);ii=np.arange(len(rows));bands=np.digitize(D,[10,20,40,80]);gs=np.unique(g)
folds={}
for seed in [20261028,7331,48105]:
 perm=np.random.default_rng(seed).permutation(gs);mp={v:j%5 for j,v in enumerate(perm)};folds[str(seed)]=np.array([mp[v] for v in g])
folds['distance_bands']=bands
for scenario in ['MC','NC']:
 E=1.2e-10*10**np.array([float(r['log_eN_'+scenario]) for r in rows]);err=np.array([float(r['err_log_eN_'+scenario]) for r in rows])
 def pred(k,c,i):
  a=1e-10*10**c[0]
  if k=='RAR':p=rar(G[i],a)
  else:
   e=np.zeros(len(i)) if k=='simple' else E[i]*(10**c[1] if k=='EFE_scale' else 1)
   p=simple(G[i],a,e)
  return np.log10(p)
 jj={}
 for k in ['RAR','simple','EFE_fixed','EFE_scale']:
  n=2 if k=='EFE_scale' else 1;fn=lambda c,i:pred(k,c,i);co=fit(fn,n,ii,weights(g,ii));dd={'coefficients_all':co.tolist(),'validation':{}}
  for label,fd in folds.items():
   pr,v=evaluate(fn,n,g,fd);dd['validation'][label]=v;save('environment_'+scenario+'_'+k+'_'+label,rows,pr,fd)
  jj[k]=dd
  print('env',scenario,k,[(label,round(d['metrics']['RMS_dex'],5)) for label,d in dd['validation'].items()],co.tolist(),flush=True)
 J['environment'][scenario]=jj
 # Model-dependent field-uncertainty sensitivity, shared draw at all radii.
 baseE=E.copy();draws=[];k='EFE_fixed';fn=lambda c,i:pred(k,c,i)
 for it in range(100):
  zz={v:rng.normal() for v in gs};E=baseE*10**(err*np.array([zz[v] for v in g]));pp,v=evaluate(fn,1,g,bands);draws.append(v['metrics']['RMS_dex'])
 E=baseE;J['environment'][scenario]['field_error_sensitivity_100']={'RMS_percentile_5_50_95':np.percentile(draws,[5,50,95]).tolist(),'values':draws,'interpretation':'Independent per-galaxy log-field perturbations using quoted errors; not a joint field posterior; spatially correlated systematics omitted.'}
 (O/'results.json').write_text(json.dumps(J,indent=2))
rows=list(csv.DictReader((R/'history_input.csv').open()));g=np.array([r['system'] for r in rows]);H=np.array([r['host'] for r in rows]);G=np.array([float(r['gbar_m_s2']) for r in rows]);Y=np.log10([float(r['gobs_m_s2']) for r in rows]);T=np.array([float(r['log_tau90_yr']) for r in rows]);em=np.array([float(r['err_minus_total_dex']) for r in rows]);ep=np.array([float(r['err_plus_total_dex']) for r in rows]);d=np.log10([float(r['Dhost_kpc'])/100 for r in rows]);mass=np.log10([float(r['stellar_mass_Msun'])/1e6 for r in rows]);ii=np.arange(len(rows));fd=H
models={'RAR':(1,0),'history_5Gyr':(2,5),'history_2Gyr':(2,2),'history_10Gyr':(2,10),'distance':(2,0),'history_distance':(3,5),'mass':(2,0),'history_mass':(3,5)}
def histpred(k,c,i):
 aa=c[0]+np.zeros(len(i));tau=models[k][1]
 if tau:aa+=c[1]*np.exp(-np.minimum(10**(T[i]-9),13.8)/tau)
 if k=='distance':aa+=c[1]*d[i]
 if k=='mass':aa+=c[1]*mass[i]
 if k=='history_distance':aa+=c[2]*d[i]
 if k=='history_mass':aa+=c[2]*mass[i]
 return np.log10(rar(G[i],1e-10*10**aa))
for k,(n,tau) in models.items():
 fn=lambda c,i:histpred(k,c,i);co=fit(fn,n,ii,weights(g,ii,H));dd={'coefficients_all':co.tolist(),'validation':{}}
 for label,fd in [('host_holdout',H),('leave_one_dwarf_out',ii)]:
  pr,v=evaluate(fn,n,g,fd,H);v['host_metrics']={}
  for h in np.unique(H):
   m=H==h;e=pr[m]-Y[m];v['host_metrics'][h]={'RMS_dex':float(np.sqrt(np.mean(e*e))),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'n':int(sum(m))}
  dd['validation'][label]=v;save('history_'+k+'_'+label,rows,pr,fd)
 J['history'][k]=dd
 print('history',k,[(l,round(v['metrics']['RMS_dex'],5)) for l,v in dd['validation'].items()],dd['validation']['host_holdout']['coefficients'],flush=True)
(O/'results.json').write_text(json.dumps(J,indent=2))
# History-error sensitivity refits the baseline and the fixed 5-Gyr model each time.
baseT=T.copy();vals=[]
for it in range(500):
 z=rng.normal(size=len(T));T=np.minimum(baseT+z*np.where(z>=0,ep,em),np.log10(13.8e9));fn=lambda c,i:histpred('history_5Gyr',c,i);pr,v=evaluate(fn,2,g,H,H);vals.append(v['metrics']['RMS_dex'])
T=baseT;base=J['history']['RAR']['validation']['host_holdout']['metrics']['RMS_dex'];J['history_error_sensitivity']={'draws':500,'RMS_percentile_5_50_95':np.percentile(vals,[5,50,95]).tolist(),'fraction_improving_baseline':float(np.mean(np.array(vals)<base)),'values':vals,'interpretation':'Asymmetric independent perturbations of photometric log-age summaries; conditional sensitivity, not a posterior or significance from two hosts.'}
# Broad model reference refit on all five populations, excluding every dwarf of the target host.
sys.path.insert(0,str(R.parent/'dm_one_formula'));import search_one_formula as b
broad={};pr=np.zeros(len(rows))
for h in ['MW','M31']:
 keep=[not (r['population']=='dSph' and (r['environment'] in (['Milky_Way','G/L'] if h=='MW' else ['M31','A/L']))) for r in b.rows];tr=b.I[keep]
 withheld_names=set(g[H==h]);assert not withheld_names & set(b.rows[i]['system'] for i in tr)
 c=b.fit('power_sum',tr);inds=np.array([next(i for i,r in enumerate(b.rows) if r['system']==name) for name in g[H==h]]);pr[H==h]=b.predict('power_sum',c,inds)-10
 broad[h]={'coefficients':c.tolist(),'training_rows':len(tr),'excluded_rows':int(sum(~np.array(keep)))}
J['broad_power_sum_host_excluded']={'fits':broad,'metrics_on_15_matches':metrics(pr,g,H)};save('broad_power_sum_host_excluded',rows,pr,H)
J['fit_calls']=CALLS+b.CALLS;J['optimizer_warnings']=WARN+b.WARN
(O/'results.json').write_text(json.dumps(J,indent=2));print('DONE',J['fit_calls'],J['optimizer_warnings'],flush=True)
