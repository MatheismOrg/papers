from pathlib import Path
import sys,json,hashlib,time,csv
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'dm_one_formula'))
import search_one_formula as b
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
SPEC={
 'softened_power':([-.14,-2.6,-.2,1.,.77,-3.],[-3,-5,-5,0,0,-5],[3,3,3,1,.95,2]),
 'bounded_U':([.1,.4,1.,0.,1.],[-3,-3,0,-4,.2],[3,3,6,4,5]),
 'half_power':([-.1,-2,-.2,.77],[-3,-5,-5,0],[3,3,3,.95])}
def predict(k,c,idx):
 x=b.x[idx];q=b.q[idx]
 if k=='softened_power':
  lsum=np.logaddexp(b.LN*q,b.LN*c[5])/b.LN
  aa=logsumexp(b.LN*np.array([q*0+c[0],c[1]-c[3]*lsum,c[2]+c[4]*q]),axis=0)/b.LN
 elif k=='half_power':
  aa=logsumexp(b.LN*np.array([q*0+c[0],c[1]-.5*q,c[2]+c[3]*q]),axis=0)/b.LN
 else:
  z=np.tanh((q-c[3])/c[4]);aa=c[0]+c[1]*z+c[2]*z*z
 return b.rar(x,aa)
WARN=[];CALLS=0

def fit(k,idx):
 global CALLS
 CALLS+=1
 start,lo,hi=map(lambda z:np.array(z,dtype=float),SPEC[k]);starts=[start]
 if k=='softened_power':
  for v in [-1.,-4.5]:
   z=start.copy();z[-1]=v;starts.append(z)
 elif k=='bounded_U':
  for v in [-1.,1.]:
   z=start.copy();z[3]=v;starts.append(z)
 sw=np.sqrt(b.weights(idx));rs=[least_squares(lambda c:sw*(predict(k,c,idx)-b.y[idx]),z,bounds=(lo,hi),max_nfev=2400,ftol=1e-9,xtol=1e-9,gtol=1e-9) for z in starts]
 rr=min(rs,key=lambda z:z.cost)
 if not rr.success:WARN.append({'model':k,'message':rr.message,'n':len(idx)})
 return rr.x

def save(name,pr,fd=None):
 old=b.OUT;b.OUT=OUT;b.savepred(name,pr,fd=fd);b.OUT=old

def bootstrap(preds):
 gs=np.unique(b.group);gi={g:i for i,g in enumerate(gs)};mat={}
 for k,pr in preds.items():
  per=[]
  for p in b.POPS:
   mask=b.pop==p;gg=np.unique(b.group[mask]);per.append((np.array([gi[g] for g in gg]),np.array([np.mean((pr[(b.pop==p)&(b.group==g)]-b.y[(b.pop==p)&(b.group==g)])**2) for g in gg])))
  mat[k]=per
 rng=np.random.default_rng(20261028);mult=rng.poisson(1,size=(5000,len(gs)));scores={}
 for k,per in mat.items():scores[k]=np.sqrt(sum((mult[:,inds]@err)/mult[:,inds].sum(axis=1) for inds,err in per)/len(per))
 out={}
 for k in preds:
  if k=='RAR':continue
  for base in ['RAR','potential_quad','power_sum']:
   if k==base:continue
   delta=scores[k]-scores[base];out[k+' minus '+base]={'point_difference_dex':b.B(preds[k],b.I)-b.B(preds[base],b.I),'percentile_95_dex':np.percentile(delta,[2.5,97.5]).tolist(),'fraction_negative_descriptive_only':float(np.mean(delta<0))}
 return out

def readpred(path):
 rr=list(csv.DictReader(path.open()));assert [int(r['row']) for r in rr]==list(b.I)
 return np.array([float(r['log10_residual']) for r in rr])+b.y

def main():
 j={'input_sha256':hashlib.sha256(b.DATA.read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((ROOT/'PROTOCOL.txt').read_bytes()).hexdigest(),'spec':SPEC,'models':{}}
 for k in SPEC:
  tt=time.time();co=fit(k,b.I);d={'coefficients':co.tolist(),'CV':{},'training':b.metrics(predict(k,co,b.I),b.I)}
  for seed in [20261027,7331,48105]:
   fd=b.folds(b.I,seed);pr=np.full(b.n,np.nan);cs=[]
   for fold in range(5):
    tr=b.I[fd!=fold];te=b.I[fd==fold];assert not set(b.group[tr])&set(b.group[te]);c=fit(k,tr);pr[te]=predict(k,c,te);cs.append(c.tolist())
   d['CV'][str(seed)]={'metrics':b.metrics(pr,b.I),'coefficients':cs};save(k+'_cv_'+str(seed),pr,fd)
  po=np.full(b.n,np.nan);cs={}
  for p in b.POPS:
   te=b.I[b.pop==p];groups=set(b.group[te]);tr=b.I[np.array([g not in groups for g in b.group])];c=fit(k,tr);po[te]=predict(k,c,te);cs[p]=c.tolist()
  d['population_holdout']={'metrics':b.metrics(po,b.I),'coefficients':cs};save(k+'_population_holdout',po)
  j['models'][k]=d
  print(k,'CV',*[round(d['CV'][str(s)]['metrics']['B_dex'],6) for s in [20261027,7331,48105]],'LOPO',round(b.B(po,b.I),6),'coef',co.tolist(),'seconds',round(time.time()-tt),flush=True)
  (OUT/'results.json').write_text(json.dumps(j,indent=2))
 for split in ['cv_20261027','population_holdout']:
  pp={k:readpred(b.OUT/(k+'_'+split+'.csv')) for k in ['RAR','potential_quad','power_sum']}
  pp.update({k:readpred(OUT/(k+'_'+split+'.csv')) for k in SPEC})
  j['bootstrap_'+split]=bootstrap(pp)
 j['fit_calls']=CALLS;j['warnings']=WARN
 (OUT/'results.json').write_text(json.dumps(j,indent=2));print('DONE',CALLS,'fits',WARN,flush=True)
if __name__=='__main__':main()
