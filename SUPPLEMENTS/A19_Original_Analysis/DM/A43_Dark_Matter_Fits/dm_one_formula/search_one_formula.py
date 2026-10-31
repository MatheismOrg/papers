from pathlib import Path
import csv,json,hashlib,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/'dm_broad_joint/combined_input.csv';OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader(DATA.open()));n=len(rows);I=np.arange(n)
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows]);u=np.log10([float(r['radius_kpc']) for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]);q=x+u
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r.get('environment','') for r in rows]);POPS=sorted(set(pop));LN=np.log(10.)
W=1/(1+10.**(2*(x-2)))
# Starts, lower/upper bounds. Each vector is in explicitly dimensionless logs.
SPEC={
 'RAR':([.1],[-3],[3]),
 'potential_quad':([.067,.35,.172],[-3,-2,-1],[3,2,1]),
 'poly2D':None,
 'piecewise_RAR':([.1]*5,[-3]*5,[3]*5),
 'Lelli_floor':([.1,-1.],[-3,-4],[3,2]),
 'Santos_log':None,
 'potential_floor':([.067,.35,.172,-1.],[-3,-2,-1,-4],[3,2,1,2]),
 'linear_potential_floor':([.1,.3,-1.],[-3,-2,-4],[3,2,2]),
 'radius_floor':([.1,-1.,.2],[-3,-4,-1],[3,2,1]),
 'potential_radius':([.067,.35,.172,0.],[-3,-2,-1,-2],[3,2,1,2]),
 'single_index':([.067,.35,.172,1.],[-3,-2,0.,-2],[3,2,1,3]),
 'power_sum':([-.1,-1.,-.5,.5,.5],[-3,-5,-5,0,0],[3,3,3,1,.95]),
 'potential_cubic':([.067,.35,.172,0.],[-3,-2,-1,-.5],[3,2,1,.5]),
 'effective_quad2D':([.067,.35,.35,.172,.344,.172],[-3,-2,-2,-1,-1,-1],[3,2,2,1,1,1])}
SHARED=[k for k in SPEC if k!='piecewise_RAR']
DEFS={
 'RAR':'g=F(gbar,a_ref*10^a), F(g,a)=g/[1-exp(-sqrt(g/a))]',
 'potential_quad':'loga=a+b*q+c*q^2; untapered previous model',
 'poly2D':'y=c0+c1*x+c2*u+c3*x^2+c4*x*u+c5*u^2; no positivity constraint',
 'piecewise_RAR':'one loga per population; prediction uses population label; NOT shared',
 'Lelli_floor':'g=F(gbar,a)+h*exp[-sqrt(gbar*a)/h]; loga=c0,logh=c1',
 'Santos_log':'log10(gobs in SI)=A/log10(gbar in SI)+B; coefficient order A,B',
 'potential_floor':'same floor; loga=W*(a+b*q+c*q^2), logh=c3',
 'linear_potential_floor':'same floor; loga=W*(a+b*q), logh=c2',
 'radius_floor':'same floor; loga=c0,logh=c1+c2*u',
 'potential_radius':'loga=W*(a+b*q+c*q^2+d*u)',
 'single_index':'loga=W*(a+b*s+c*s^2),s=x+k*u; coefficients a,b,c,k',
 'power_sum':'a_eff/a_ref=10^A+10^(B-alpha*q)+10^(C+beta*q); order A,B,C,alpha,beta',
 'potential_cubic':'loga=W*(a+b*q+c*q^2+d*q^3)',
 'effective_quad2D':'loga=W*(c0+c1*x+c2*u+c3*x^2+c4*x*u+c5*u^2)'}
def weights(idx):
 w=np.empty(len(idx));ps=np.unique(pop[idx])
 for p in ps:
  gs=np.unique(group[idx][pop[idx]==p])
  for g in gs:
   m=(pop[idx]==p)&(group[idx]==g);w[m]=1/(len(ps)*len(gs)*sum(m))
 return w

def rar(xx,aa):
 t=10.**np.clip(.5*(xx-aa),-290,290)
 return xx-np.log10(-np.expm1(-t))

def basefeature(idx):
 return np.column_stack([np.ones(len(idx)),x[idx],u[idx],x[idx]**2,x[idx]*u[idx],u[idx]**2])

def predict(k,c,idx):
 xx=x[idx];uu=u[idx];qq=q[idx];ww=W[idx]
 if k=='poly2D':return basefeature(idx)@c
 if k=='Santos_log':return c[0]/(xx-10)+c[1]+10
 if k=='piecewise_RAR':return rar(xx,np.array([c[POPS.index(p)] for p in pop[idx]]))
 if k in ['RAR','Lelli_floor','radius_floor']:aa=c[0]+xx*0
 elif k in ['potential_quad','potential_floor','potential_radius','potential_cubic']:
  aa=c[0]+c[1]*qq+c[2]*qq**2
  if k=='potential_radius':aa+=c[3]*uu
  if k=='potential_cubic':aa+=c[3]*qq**3
  if k!='potential_quad':aa*=ww
 elif k=='linear_potential_floor':aa=ww*(c[0]+c[1]*qq)
 elif k=='single_index':
  ss=xx+c[3]*uu;aa=ww*(c[0]+c[1]*ss+c[2]*ss**2)
 elif k=='power_sum':aa=logsumexp(LN*np.array([xx*0+c[0],c[1]-c[3]*qq,c[2]+c[4]*qq]),axis=0)/LN
 elif k=='effective_quad2D':aa=ww*(basefeature(idx)@c)
 else:raise ValueError(k)
 rr=rar(xx,aa)
 if k in ['Lelli_floor','potential_floor','linear_potential_floor','radius_floor']:
  hh=c[1]+c[2]*uu if k=='radius_floor' else c[-1]+xx*0
  ff=hh-10.**np.clip(.5*(xx+aa)-hh,-290,290)/LN
  rr=np.logaddexp(LN*rr,LN*ff)/LN
 return rr
WARN=[];CALLS=0

def fit(k,idx):
 global CALLS
 CALLS+=1;sw=np.sqrt(weights(idx))
 if k in ['poly2D','Santos_log']:
  A=basefeature(idx) if k=='poly2D' else np.column_stack([1/(x[idx]-10),np.ones(len(idx))]);yy=y[idx] if k=='poly2D' else y[idx]-10
  return np.linalg.lstsq(A*sw[:,None],yy*sw,rcond=None)[0]
 c0,lo,hi=[np.array(v,dtype=float) for v in SPEC[k]]
 starts=[c0]
 if k in ['power_sum','single_index']:
  c2=c0.copy()
  if k=='single_index':c2[-1]=.4
  else:c2[-2:]=[.95,.85]
  starts.append(c2)
 results=[least_squares(lambda c:sw*(predict(k,c,idx)-y[idx]),s,bounds=(lo,hi),max_nfev=1600,xtol=1e-9,ftol=1e-9,gtol=1e-9) for s in starts]
 r=min(results,key=lambda z:z.cost)
 if not r.success:WARN.append({'model':k,'n':len(idx),'message':r.message})
 return r.x

def B(pr,idx):return float(np.sqrt(np.sum(weights(idx)*(pr-y[idx])**2)))

def metrics(pr,idx):
 out={'B_dex':B(pr,idx),'populations':{}}
 for p in np.unique(pop[idx]):
  m=pop[idx]==p;e=pr[m]-y[idx[m]]
  out['populations'][p]={'median_acceleration_error_percent':float(np.median(abs(10**np.clip(e,-290,290)-1))*100),'RMS_dex':float(np.sqrt(np.mean(e*e))),'group_RMS_dex':B(pr[m],idx[m]),'median_speed_error_percent':float(np.median(abs(10**np.clip(e/2,-290,290)-1))*100),'rows':int(sum(m))}
 return out

def folds(idx,seed,nf=5):
 rng=np.random.default_rng(seed);mp={}
 for p in sorted(set(pop[idx])):
  gs=[g for g in np.unique(group[idx][pop[idx]==p]) if g not in mp];ct=np.array([sum(v==j for g,v in mp.items() if g in group[idx][pop[idx]==p]) for j in range(nf)])
  for g in rng.permutation(gs):
   j=int(np.argmin(ct));mp[g]=j;ct[j]+=1
 return np.array([mp[g] for g in group[idx]])

def savepred(k,pr,idx=I,fd=None):
 with (OUT/(k+'.csv')).open('w') as f:
  wr=csv.writer(f);wr.writerow(['row','population','system','group','fold','gbar_m_s2','radius_kpc','gobs_m_s2','gpred_m_s2','log10_residual'])
  for j,i in enumerate(idx):wr.writerow([i,pop[i],rows[i]['system'],group[i],'' if fd is None else int(fd[j]),10**x[i]*1e-10,10**u[i],10**y[i]*1e-10,10**np.clip(pr[j],-290,290)*1e-10,pr[j]-y[i]])

def main():
 result={'status':'Exploratory follow-up on previously inspected data. No claim of universal law or superiority to all published theories.','input_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'definitions':DEFS,'bounds_and_starts':SPEC,'population_order':POPS,'fixed_regulator':'W=1/[1+(gbar/1e-8)^2]','rows':n,'groups':len(set(group)),'models':{},'nested':{}}
 for k in SPEC:
  t=time.time();cc=fit(k,I);cv={};preds={}
  for seed in [20261027,7331,48105]:
   fd=folds(I,seed);pr=np.full(n,np.nan);co=[]
   for fold in range(5):
    tr=I[fd!=fold];te=I[fd==fold];assert not set(group[tr])&set(group[te]);c=fit(k,tr);pr[te]=predict(k,c,te);co.append(c.tolist())
   cv[str(seed)]={'metrics':metrics(pr,I),'coefficients':co};preds[seed]=pr;savepred(k+'_cv_'+str(seed),pr,fd=fd)
  d={'coefficients':cc.tolist(),'training':metrics(predict(k,cc,I),I),'CV':cv}
  if k!='piecewise_RAR':
   po=np.full(n,np.nan);pc={}
   for p in POPS:
    te=I[pop==p];tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];c=fit(k,tr);po[te]=predict(k,c,te);pc[p]=c.tolist()
   d['population_holdout']={'metrics':metrics(po,I),'coefficients':pc};savepred(k+'_population_holdout',po)
  result['models'][k]=d
  print(k,'CV',*[round(cv[str(s)]['metrics']['B_dex'],4) for s in [20261027,7331,48105]],'LOPO',round(d.get('population_holdout',{}).get('metrics',{}).get('B_dex',float('nan')),4),'coeff',cc.tolist(),'sec',round(time.time()-t,1),flush=True)
  (OUT/'results.json').write_text(json.dumps(result,indent=2))
 # Nested candidate selection uses all shared candidates; older baselines included.
 fd=folds(I,20261027);pr=np.full(n,np.nan);selections=[]
 for outer in range(5):
  tr=I[fd!=outer];te=I[fd==outer];inner=folds(tr,39173+outer,3);scores={}
  for k in SHARED:
   ip=np.full(len(tr),np.nan)
   for j in range(3):
    tri=tr[inner!=j];tei=tr[inner==j];c=fit(k,tri);ip[inner==j]=predict(k,c,tei)
   scores[k]=B(ip,tr)
  chosen=min(scores,key=scores.get);co=result['models'][chosen]['CV']['20261027']['coefficients'][outer];pr[te]=predict(chosen,co,te);selections.append({'outer_fold':outer,'chosen':chosen,'scores':scores})
  print('nested',outer,chosen,flush=True)
 result['nested']={'metrics':metrics(pr,I),'selection':selections};savepred('nested',pr,fd=fd)
 # Apply the same pre-specified follow-up host test to every shared candidate.
 result['dwarf_host_holdout']={}
 for host in ['MW','M31']:
  te=I[(pop=='dSph')&(env==host)]
  if not len(te):continue
  tr=np.setdiff1d(I,te);d={}
  for k in SHARED:
   c=fit(k,tr);pp=predict(k,c,te);d[k]={'coefficients':c.tolist(),'metrics':metrics(pp,te)}
  result['dwarf_host_holdout'][host]=d
 result['fit_calls']=CALLS;result['warnings']=WARN
 (OUT/'results.json').write_text(json.dumps(result,indent=2));print('DONE',CALLS,'fit calls; warnings',WARN,flush=True)
if __name__=='__main__':main()
