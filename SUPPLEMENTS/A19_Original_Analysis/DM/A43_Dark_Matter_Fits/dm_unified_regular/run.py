from pathlib import Path
import csv,json,hashlib,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
R=Path(__file__).resolve().parent;D=R.parent/'dm_broad_joint/combined_input.csv';X=R.parent/'dm_fixed_radius/fixed_radius_inputs.csv';OUT=R/'results';OUT.mkdir(exist_ok=True);LN=np.log(10.)
NEW=['mixed_regular'];KINDS=NEW
SPEC={'mixed_regular':{'starts':[[-.1,-2.,-.2,.6,.7,-2.,2.5],[0.,-1.5,-.5,.3,.3,-1.,1.5],[-.4,-3.,-.7,.9,.9,-1.5,3.5]],'lo':[-3.,-5.,-5.,0.,0.,-2.,0.],'hi':[3.,3.,3.,.95,.95,0.,4.]}}
COMPARATORS=['RAR','potential_quad','power_sum','mixed_power','mixed_gate','mixed_gate_range']

def prepare():
 import shutil
 for name in ['original','augmented']:
  src=R.parent/'dm_unified_repair'/(name+'_input.csv');shutil.copy2(src,R/(name+'_input.csv'))
 return tuple(list(csv.DictReader((R/(name+'_input.csv')).open())) for name in ['original','augmented'])

def load(rows_,name):
 global rows,I,n,x,u,y,q,m,pop,group,env,source,POPS,O,DATASET
 DATASET=name;rows=rows_;n=len(rows);I=np.arange(n);O=OUT/name;O.mkdir(exist_ok=True)
 x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows]);u=np.log10([float(r['radius_kpc']) for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]);q=x+u;m=x+2*u
 pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r.get('environment','') for r in rows]);source=np.array(['CLASH' if r.get('source','')=='Tian2020_CDS_fig2' else r.get('source','') for r in rows]);POPS=sorted(set(pop))

def weights(idx):
 w=np.empty(len(idx));ps=np.unique(pop[idx])
 for p in ps:
  gs=np.unique(group[idx][pop[idx]==p])
  for g in gs:
   z=(pop[idx]==p)&(group[idx]==g);w[z]=1/(len(ps)*len(gs)*sum(z))
 return w

def rar(xx,aa):
 t=10**np.clip(.5*(xx-aa),-290,290);return xx-np.log10(-np.expm1(-t))

def logae(k,c,xx,uu):
 qq=xx+uu;mm=xx+2*uu
 neg=c[1]-c[3]*mm-np.logaddexp(0,LN*4*(c[5]-uu))/LN-np.logaddexp(0,LN*2*(xx-2))/LN
 pos=c[2]+c[4]*qq-np.logaddexp(0,LN*(uu-c[6]))/LN
 return logsumexp(LN*np.array([xx*0+c[0],neg,pos]),axis=0)/LN

def rawpredict(k,c,xx,uu):return rar(xx,logae(k,c,xx,uu))
def predict(k,c,idx):return rawpredict(k,c,x[idx],u[idx])
FIT_LOG=[]
def fit(k,idx,tag):
 sw=np.sqrt(weights(idx));s=SPEC[k];rs=[]
 for j,z in enumerate(s['starts']):
  rr=least_squares(lambda c:sw*(predict(k,c,idx)-y[idx]),np.array(z),bounds=(s['lo'],s['hi']),max_nfev=2000,ftol=1e-9,xtol=1e-9,gtol=1e-9)
  d={'dataset':DATASET,'model':k,'split':tag,'start_index':j,'train_rows':len(idx),'success':bool(rr.success),'message':rr.message,'nfev':rr.nfev,'cost':float(rr.cost),'optimality':float(rr.optimality),'coefficients':rr.x.tolist(),'active_bounds':rr.active_mask.tolist()};FIT_LOG.append(d);rs.append(rr)
 rr=min(rs,key=lambda z:z.cost)
 return rr.x

def folds(idx,seed,nf=5):
 rng=np.random.default_rng(seed);mp={}
 for p in sorted(set(pop[idx])):
  gs=[g for g in np.unique(group[idx][pop[idx]==p]) if g not in mp];ct=np.array([sum(v==j for g,v in mp.items() if g in group[idx][pop[idx]==p]) for j in range(nf)])
  for g in rng.permutation(gs):j=int(np.argmin(ct));mp[g]=j;ct[j]+=1
 return np.array([mp[g] for g in group[idx]])

def B(pr,idx):return float(np.sqrt(np.sum(weights(idx)*(pr-y[idx])**2)))
def submetrics(pr,idx):
 e=pr-y[idx];return {'RMS_dex':float(np.sqrt(np.mean(e*e))),'group_RMS_dex':B(pr,idx),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'median_speed_error_percent':float(np.median(abs(10**(e/2)-1))*100),'rows':len(idx),'groups':len(set(group[idx]))}
def metrics(pr,idx):
 d={'B_dex':B(pr,idx),'rows':len(idx),'groups':len(set(group[idx])),'populations':{},'sources':{}}
 for p in np.unique(pop[idx]):
  zz=pop[idx]==p;d['populations'][p]=submetrics(pr[zz],idx[zz])
 for s in np.unique(source[idx]):
  zz=source[idx]==s;d['sources'][s]=submetrics(pr[zz],idx[zz])
 return d

def save(name,pr,idx,fd=None):
 with (O/(name+'.csv')).open('w') as f:
  w=csv.writer(f);w.writerow(['row','population','system','group','environment','source','source_detail','fold','gbar_m_s2','radius_kpc','gobs_m_s2','gpred_m_s2','log10_residual'])
  for j,i in enumerate(idx):w.writerow([i,pop[i],rows[i]['system'],group[i],env[i],source[i],rows[i].get('source',''),'' if fd is None else str(fd[j]),10**x[i]*1e-10,10**u[i],10**y[i]*1e-10,10**pr[j]*1e-10,pr[j]-y[i]])

def bootstrap(preds):
 gs=np.unique(group);gi={g:i for i,g in enumerate(gs)};rng=np.random.default_rng(20261028);mult=rng.poisson(1,size=(5000,len(gs)));scores={}
 for k,pr in preds.items():
  pp=[]
  for p in POPS:
   gg=np.unique(group[pop==p]);inds=np.array([gi[g] for g in gg]);err=np.array([np.mean((pr[(pop==p)&(group==g)]-y[(pop==p)&(group==g)])**2) for g in gg]);pp.append((mult[:,inds]@err)/mult[:,inds].sum(axis=1))
  scores[k]=np.sqrt(np.mean(pp,axis=0))
 z={}
 for k in NEW:
  for base in COMPARATORS:
   e=scores[k]-scores[base];z[k+' minus '+base]={'point_difference_dex':B(preds[k],I)-B(preds[base],I),'percentile_95_dex':np.percentile(e,[2.5,97.5]).tolist(),'fraction_negative_descriptive_only':float(np.mean(e<0))}
 return z

def main():
 original,augmented=prepare();out={'status':'Adaptive stage3 single-family regularity repair; explicit imposed length prior; no independent confirmation.','protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'fit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_hashes':{str((R.parent/'dm_unified_repair'/fn).relative_to(R.parent)):hashlib.sha256((R.parent/'dm_unified_repair'/fn).read_bytes()).hexdigest() for fn in ['original_input.csv','augmented_input.csv']},'input_hashes':{fn:hashlib.sha256((R/fn).read_bytes()).hexdigest() for fn in ['original_input.csv','augmented_input.csv']},'spec':SPEC,'new_families':NEW,'comparators_reused':COMPARATORS,'models_by_dataset':{},'fit_log':FIT_LOG}
 for name,data in [('original',original),('augmented',augmented)]:
  load(data,name);preds={key:{} for key in ['cv_20261027','population_holdout']};models={};out['models_by_dataset'][name]={'rows':n,'systems':len(set(r['system'] for r in rows)),'groups':len(set(group)),'population_order':POPS,'models':models}
  for k in KINDS:
   t=time.time();co=fit(k,I,'all');d={'coefficients':co.tolist(),'training':metrics(predict(k,co,I),I),'CV':{},'host_holdout':{},'source_holdout':{}};save(k+'_all',predict(k,co,I),I)
   for seed in [20261027,7331,48105]:
    fd=folds(I,seed);pr=np.full(n,np.nan);cs=[]
    for fold in range(5):
     tr=I[fd!=fold];te=I[fd==fold];assert not set(group[tr])&set(group[te]);cc=fit(k,tr,f'cv_{seed}_{fold}');pr[te]=predict(k,cc,te);cs.append(cc.tolist())
    d['CV'][str(seed)]={'metrics':metrics(pr,I),'coefficients':cs};save(k+'_cv_'+str(seed),pr,I,fd)
    if seed==20261027:preds['cv_20261027'][k]=pr
   pr=np.full(n,np.nan);cs={}
   for p in POPS:
    te=I[pop==p];tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];assert not set(group[tr])&tg;cc=fit(k,tr,'population_'+p);pr[te]=predict(k,cc,te);cs[p]=cc.tolist()
   d['population_holdout']={'metrics':metrics(pr,I),'coefficients':cs};save(k+'_population_holdout',pr,I,pop);preds['population_holdout'][k]=pr
   for h in ['MW','M31']:
    host_aliases={'MW':['Milky_Way','G/L'],'M31':['M31','A/L']};te=I[(pop=='dSph')&np.isin(env,host_aliases[h])];assert len(te)=={'MW':26,'M31':18}[h];tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];cc=fit(k,tr,'host_'+h);pr=predict(k,cc,te);d['host_holdout'][h]={'excluded_rows':len(te),'training_rows':len(tr),'coefficients':cc.tolist(),'metrics':metrics(pr,te)};save(k+'_host_'+h,pr,te,np.array([h]*len(te)))
   if name=='augmented':
    for src in ['XCOP','CLASH']:
     te=I[source==src];assert len(te)=={'XCOP':21,'CLASH':84}[src],(src,len(te));tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];cc=fit(k,tr,'source_'+src);pr=predict(k,cc,te);d['source_holdout'][src]={'excluded_rows':len(te),'training_rows':len(tr),'coefficients':cc.tolist(),'metrics':metrics(pr,te)};save(k+'_source_'+src,pr,te,np.array([src]*len(te)))
   models[k]=d;(OUT/'results.json').write_text(json.dumps(out,indent=2));print(name,k,'CV',*[round(d['CV'][str(s)]['metrics']['B_dex'],6) for s in [20261027,7331,48105]],'LOPO',round(d['population_holdout']['metrics']['B_dex'],6),'coeff',co.tolist(),'sec',round(time.time()-t,1),flush=True)
  prior=json.loads((R.parent/'dm_unified_repair/results/results.json').read_text())
  assert out['input_hashes'][name+'_input.csv']==prior['input_hashes'][name+'_input.csv']
  for split,pp in preds.items():
   current=list(csv.DictReader((O/('mixed_regular_'+split+'.csv')).open()))
   for base in COMPARATORS:
    f=R.parent/'dm_unified_repair/results'/name/(base+'_'+split+'.csv');rr=list(csv.DictReader(f.open()))
    assert len(rr)==n
    assert all(a['row']==b['row'] and a['population']==b['population'] and a['group']==b['group'] and a['fold']==b['fold'] for a,b in zip(current,rr))
    pp[base]=np.log10([float(r['gpred_m_s2'])/1e-10 for r in rr])
   out['models_by_dataset'][name]['bootstrap_'+split]=bootstrap(pp)
  out['models_by_dataset'][name]['comparator_metrics_reused']={k:prior['models_by_dataset'][name]['models'][k] for k in COMPARATORS}
  out['models_by_dataset'][name]['comparison_source_sha256']=hashlib.sha256((R.parent/'dm_unified_repair/results/results.json').read_bytes()).hexdigest()
 out['optimizer_calls']=len(FIT_LOG);out['fit_calls']=sum(1 for v in FIT_LOG if v['start_index']==0);out['optimizer_failures']=[v for v in FIT_LOG if not v['success']]
 (OUT/'results.json').write_text(json.dumps(out,indent=2));print('DONE',out['fit_calls'],out['optimizer_calls'],'failures',len(out['optimizer_failures']),flush=True)
if __name__=='__main__':main()
