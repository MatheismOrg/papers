from pathlib import Path
import csv,json,hashlib,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
R=Path(__file__).resolve().parent;D=R.parent/'dm_broad_joint/combined_input.csv';O=R/'results';O.mkdir(exist_ok=True)
rows=list(csv.DictReader(D.open()));I=np.arange(len(rows));n=len(rows);LN=np.log(10.)
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows]);u=np.log10([float(r['radius_kpc']) for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]);q=x+u;m=x+2*u
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r.get('environment','') for r in rows]);POPS=sorted(set(pop))
NEW=['mass_power','mass_bounded','mass_additive','mixed_power','index_power','mass_screened'];BASE=['RAR','potential_quad','power_sum','half_power'];KINDS=BASE+NEW
MP=[[-.1,-2.,-1.,.6,.5],[0.,-1.5,-1.,.3,.3],[-.4,-3.,-2.,.9,.8]]
SPEC={
'RAR':{'starts':[[.1]],'lo':[-3.],'hi':[3.]},
'potential_quad':{'starts':[[.067,.35,.172]],'lo':[-3.,-2.,-1.],'hi':[3.,2.,1.]},
'power_sum':{'starts':[[-.1,-1.,-.5,.5,.5],[-.1,-1.,-.5,.95,.85]],'lo':[-3.,-5.,-5.,0.,0.],'hi':[3.,3.,3.,1.,.95]},
'half_power':{'starts':[[-.1,-2.,-.2,.77],[0.,-1.,-.5,.5]],'lo':[-3.,-5.,-5.,0.],'hi':[3.,3.,3.,.95]},
'mass_bounded':{'starts':[[.1,.3,1.,0.,2.],[0.,.1,2.,-1.,3.],[-.4,.5,2.,1.,1.]],'lo':[-3.,-3.,0.,-6.,.2],'hi':[3.,3.,6.,6.,8.]}}
for k in ['mass_power','mass_additive','mass_screened']:SPEC[k]={'starts':MP,'lo':[-3.,-5.,-5.,0.,0.],'hi':[3.,3.,3.,.95,.95]}
SPEC['mixed_power']={'starts':[[-.1,-2.,-.2,.6,.7],[0.,-1.5,-.5,.3,.3],[-.4,-3.,-.7,.9,.9]],'lo':[-3.,-5.,-5.,0.,0.],'hi':[3.,3.,3.,.95,.95]}
SPEC['index_power']={'starts':[[-.1,-2.,-.2,.6,.7,1.],[-.1,-2.,-.7,.5,.5,2.],[-.4,-3.,-.7,.9,.9,.5]],'lo':[-3.,-5.,-5.,0.,0.,0.],'hi':[3.,3.,3.,.95,.95,2.]}

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
 if k=='RAR':return np.zeros_like(xx)+c[0]
 if k=='potential_quad':return c[0]+c[1]*qq+c[2]*qq*qq
 if k=='mass_bounded':
  z=np.tanh((mm-c[3])/c[4]);return c[0]+c[1]*z+c[2]*z*z
 if k=='power_sum':ss,tt,al,be=qq,qq,c[3],c[4]
 elif k=='half_power':ss,tt,al,be=qq,qq,.5,c[3]
 else:
  ss=mm;tt=qq if k=='mixed_power' else xx+c[5]*uu if k=='index_power' else mm;al,be=c[3],c[4]
 return logsumexp(LN*np.array([xx*0+c[0],c[1]-al*ss,c[2]+be*tt]),axis=0)/LN

def rawpredict(k,c,xx,uu):
 aa=logae(k,c,xx,uu)
 if k=='mass_additive':return np.logaddexp(LN*xx,LN*.5*(xx+aa))/LN
 if k=='mass_screened':
  extra=.5*(xx+aa)-10**np.clip(.5*(xx-c[0]),-290,290)/LN;return np.logaddexp(LN*xx,LN*extra)/LN
 return rar(xx,aa)

def predict(k,c,idx):return rawpredict(k,c,x[idx],u[idx])
FIT_LOG=[]
def fit(k,idx,tag):
 sw=np.sqrt(weights(idx));s=SPEC[k];rs=[]
 for j,z in enumerate(s['starts']):
  rr=least_squares(lambda c:sw*(predict(k,c,idx)-y[idx]),np.array(z),bounds=(s['lo'],s['hi']),max_nfev=1600,ftol=1e-9,xtol=1e-9,gtol=1e-9)
  d={'model':k,'split':tag,'start_index':j,'train_rows':len(idx),'success':bool(rr.success),'message':rr.message,'nfev':rr.nfev,'cost':float(rr.cost),'optimality':float(rr.optimality),'coefficients':rr.x.tolist(),'active_bounds':rr.active_mask.tolist()};FIT_LOG.append(d);rs.append(rr)
 rr=min(rs,key=lambda z:z.cost)
 return rr.x

def folds(idx,seed,nf=5):
 rng=np.random.default_rng(seed);mp={}
 for p in sorted(set(pop[idx])):
  gs=[g for g in np.unique(group[idx][pop[idx]==p]) if g not in mp];ct=np.array([sum(v==j for g,v in mp.items() if g in group[idx][pop[idx]==p]) for j in range(nf)])
  for g in rng.permutation(gs):j=int(np.argmin(ct));mp[g]=j;ct[j]+=1
 return np.array([mp[g] for g in group[idx]])

def B(pr,idx):return float(np.sqrt(np.sum(weights(idx)*(pr-y[idx])**2)))
def metrics(pr,idx):
 d={'B_dex':B(pr,idx),'rows':len(idx),'groups':len(set(group[idx])),'populations':{}}
 for p in np.unique(pop[idx]):
  ii=idx[pop[idx]==p];e=pr[pop[idx]==p]-y[ii];d['populations'][p]={'RMS_dex':float(np.sqrt(np.mean(e*e))),'group_RMS_dex':B(pr[pop[idx]==p],ii),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'median_speed_error_percent':float(np.median(abs(10**(e/2)-1))*100),'rows':len(ii)}
 return d

def save(name,pr,idx,fd=None):
 with (O/(name+'.csv')).open('w') as f:
  w=csv.writer(f);w.writerow(['row','population','system','group','environment','fold','gbar_m_s2','radius_kpc','gobs_m_s2','gpred_m_s2','log10_residual'])
  for j,i in enumerate(idx):w.writerow([i,pop[i],rows[i]['system'],group[i],env[i],'' if fd is None else str(fd[j]),10**x[i]*1e-10,10**u[i],10**y[i]*1e-10,10**pr[j]*1e-10,pr[j]-y[i]])

def bootstrap(preds):
 gs=np.unique(group);gi={g:i for i,g in enumerate(gs)};rng=np.random.default_rng(20261028);mult=rng.poisson(1,size=(5000,len(gs)));scores={}
 for k,pr in preds.items():
  pp=[]
  for p in POPS:
   gg=np.unique(group[pop==p]);inds=np.array([gi[g] for g in gg]);err=np.array([np.mean((pr[(pop==p)&(group==g)]-y[(pop==p)&(group==g)])**2) for g in gg]);pp.append((mult[:,inds]@err)/mult[:,inds].sum(axis=1))
  scores[k]=np.sqrt(np.mean(pp,axis=0))
 z={}
 for k in NEW:
  for base in ['RAR','potential_quad','power_sum','half_power']:
   e=scores[k]-scores[base];z[k+' minus '+base]={'point_difference_dex':B(preds[k],I)-B(preds[base],I),'percentile_95_dex':np.percentile(e,[2.5,97.5]).tolist(),'fraction_negative_descriptive_only':float(np.mean(e<0))}
 return z

def limits(k,c):
 # Stable numerical spot checks supplement analytic proofs, not observations.
 uu=np.array([0.,0.,0.,0.]);xx=np.array([-20.,-40.,20.,40.]);yy=rawpredict(k,c,xx,uu)
 # An isolated source with fixed m=0, gbar/a_ref=10^(m-2u).
 ur=np.array([8.,9.,12.,13.]);xr=-2*ur;yr=rawpredict(k,c,xr,ur)
 return {'fixed_radius_log_gbar':xx.tolist(),'fixed_radius_log_gpred':yy.tolist(),'relative_high_g_excess':(10**np.clip(yy[2:]-xx[2:],-290,290)-1).tolist(),'point_mass_far_u':ur.tolist(),'point_mass_far_logg':yr.tolist(),'far_slopes':[(yr[1]-yr[0]),(yr[3]-yr[2])]}

def main():
 out={'status':'Exploratory constrained follow-up on repeatedly inspected data; no universal theory or novelty claim.','input_sha256':hashlib.sha256(D.read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'spec':SPEC,'new_families':NEW,'comparator_families':BASE,'population_order':POPS,'models':{},'fit_log':FIT_LOG}
 preds={key:{} for key in ['cv_20261027','population_holdout']}
 for k in KINDS:
  t=time.time();co=fit(k,I,'all');d={'coefficients':co.tolist(),'training':metrics(predict(k,co,I),I),'limits':limits(k,co),'CV':{},'host_holdout':{}}
  save(k+'_all',predict(k,co,I),I)
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
   host_aliases={'MW':['Milky_Way','G/L'],'M31':['M31','A/L']}
   te=I[(pop=='dSph')&np.isin(env,host_aliases[h])];assert len(te)=={'MW':26,'M31':18}[h];tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];cc=fit(k,tr,'host_'+h);pr=predict(k,cc,te);d['host_holdout'][h]={'excluded_rows':len(te),'training_rows':len(tr),'coefficients':cc.tolist(),'metrics':metrics(pr,te)};save(k+'_host_'+h,pr,te,np.array([h]*len(te)))
  out['models'][k]=d
  (O/'results.json').write_text(json.dumps(out,indent=2));print(k,'CV',*[round(d['CV'][str(s)]['metrics']['B_dex'],6) for s in [20261027,7331,48105]],'LOPO',round(d['population_holdout']['metrics']['B_dex'],6),'SPARCspeed',round(d['CV']['20261027']['metrics']['populations']['SPARC']['median_speed_error_percent'],3),'coeff',co.tolist(),'seconds',round(time.time()-t,1),flush=True)
 for split,pp in preds.items():out['bootstrap_'+split]=bootstrap(pp)
 out['optimizer_calls']=len(FIT_LOG);out['fit_calls']=sum(1 for v in FIT_LOG if v['start_index']==0);out['optimizer_failures']=[v for v in FIT_LOG if not v['success']]
 (O/'results.json').write_text(json.dumps(out,indent=2));print('DONE',out['fit_calls'],out['optimizer_calls'],'optimizer failures',len(out['optimizer_failures']),flush=True)
if __name__=='__main__':main()
