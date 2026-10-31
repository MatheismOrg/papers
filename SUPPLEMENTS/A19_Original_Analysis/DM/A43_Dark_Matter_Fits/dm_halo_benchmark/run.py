"""Frozen protocol: contiguous radial prediction, not unseen-galaxy prediction."""
from pathlib import Path
import csv, io, json, hashlib, zipfile, time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
R=Path(__file__).resolve().parent;O=R/'results';O.mkdir(exist_ok=True)
IN=R.parent/'dm_unified_regular/augmented_input.csv';RAW=R.parent/'dm_fit/data/Rotmod_LTG.zip'
LN=np.log(10.);KPC=3.085677581491367e19;A0=np.log10(1.2)
rows=list(csv.DictReader(IN.open()));N=len(rows);I=np.arange(N)
pop=np.array([z['population'] for z in rows]);gal=np.array([z['system'] for z in rows]);group=np.array([z['group'] for z in rows]);rad=np.array([float(z['radius_kpc']) for z in rows]);gb=np.array([float(z['gbar_m_s2'])/1e-10 for z in rows]);go=np.array([float(z['gobs_m_s2'])/1e-10 for z in rows]);x=np.log10(gb);u=np.log10(rad);y=np.log10(go)
SG=sorted(set(gal[pop=='SPARC']));ELIG=[g for g in SG if sum((pop=='SPARC')&(gal==g))>=9]
gas=np.full(N,np.nan);star=np.full(N,np.nan);vobs=np.full(N,np.nan);verr=np.full(N,np.nan)
with zipfile.ZipFile(RAW) as z:
 for g in SG:
  arr=np.loadtxt(io.BytesIO(z.read(g+'_rotmod.dat')),ndmin=2)
  for i in I[(pop=='SPARC')&(gal==g)]:
   mm=np.isclose(arr[:,0],rad[i],rtol=0,atol=1e-9);assert sum(mm)==1,(i,g,rad[i])
   rr,v,e,vg,vd,vb,*_=arr[mm][0];gas[i]=vg*abs(vg)/rr*1e6/KPC/1e-10;star[i]=(.5*vd*abs(vd)+.7*vb*abs(vb))/rr*1e6/KPC/1e-10;vobs[i]=v;verr[i]=e
   assert np.isclose(gb[i],gas[i]+star[i],rtol=1e-12);assert np.isclose(go[i],v*v/rr*1e6/KPC/1e-10,rtol=1e-12)
lo_ml={}
for g in SG:
 ii=I[(pop=='SPARC')&(gal==g)];assert np.all(star[ii]>=0)
 lo_ml[g]=max(.2,float(np.max(np.divide(-gas[ii],star[ii],out=np.full(len(ii),-np.inf),where=star[ii]>0)))+1e-6)
 assert lo_ml[g]<1
with (R/'matched_components.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['row','galaxy','radius_kpc','gas_accel_over_aref','stellar_accel_over_aref','gbar_m_s2','gobs_m_s2','vobs_km_s','errv_km_s','stellar_multiplier_lower_bound'])
 for i in I[pop=='SPARC']:w.writerow([i,gal[i],rad[i],gas[i],star[i],gb[i]*1e-10,go[i]*1e-10,vobs[i],verr[i],lo_ml[gal[i]]])
(R/'augmented_input.csv').write_bytes(IN.read_bytes())

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def rar(xx,aa):
 t=10**np.clip(.5*(xx-aa),-290,290);return xx-np.log10(-np.expm1(-t))

def regular(c,xx,uu):
 inv=c[1]-c[3]*(xx+2*uu)-np.logaddexp(0,LN*4*(c[5]-uu))/LN-np.logaddexp(0,LN*2*(xx-2))/LN
 pos=c[2]+c[4]*(xx+uu)-np.logaddexp(0,LN*(uu-c[6]))/LN
 aa=logsumexp(LN*np.array([np.zeros_like(xx)+c[0],inv,pos]),axis=0)/LN
 return rar(xx,aa)

REGSPEC={'lo':[-3.,-5.,-5.,0.,0.,-2.,0.],'hi':[3.,3.,3.,.95,.95,0.,4.], 'starts':[[-.1,-2.,-.2,.6,.7,-2.,2.5],[0.,-1.5,-.5,.3,.3,-1.,1.5],[-.4,-3.,-.7,.9,.9,-1.5,3.5]]}

def massf(kind,z):
 z=np.array(z,float)
 if kind=='NFW':
  ans=np.log1p(z)-z/(1+z);s=sum(((-1)**n)*(n-1)/n*z**n for n in range(2,13));return np.where(z<.01,s,ans)
 ans=.5*np.log1p(z)+.25*np.log1p(z*z)-.5*np.arctan(z)
 s=sum(z**(4*k+3)/(4*k+3)-z**(4*k+4)/(4*k+4) for k in range(4));return np.where(z<.03,s,ans)

def predict(kind,c,ii):
 if kind=='RAR_fixed':return rar(x[ii],A0)
 if kind in ['RAR_shared','RAR_individual']:return rar(x[ii],c[0])
 if kind=='RAR_stellar':return rar(np.log10(gas[ii]+10**c[0]*star[ii]),A0)
 if kind.startswith('mixed_regular'):return regular(c,x[ii],u[ii])
 h=kind.split('_')[0];xx=rad[ii]/10**c[1];dm=10**c[0]*massf(h,xx)/(xx*xx*massf(h,1.))
 bar=gas[ii]+10**c[2]*star[ii] if kind.endswith('_stellar') else gb[ii]
 assert np.all(bar>0)
 return np.log10(bar+dm)

def weights(ii,broad=False):
 w=np.empty(len(ii));ps=np.unique(pop[ii]) if broad else ['SPARC']
 for p in ps:
  gg=np.unique(group[ii][pop[ii]==p])
  for g in gg:
   m=(pop[ii]==p)&(group[ii]==g);w[m]=1/(len(ps)*len(gg)*sum(m))
 return w

def metric(pr,ii):
 ee=pr-y[ii];gg=np.unique(gal[ii]);rms=[np.mean(ee[gal[ii]==g]**2) for g in gg]
 return {'galaxy_RMS_dex':float(np.sqrt(np.mean(rms))),'row_RMS_dex':float(np.sqrt(np.mean(ee*ee))),'median_acceleration_error_percent':float(np.median(abs(10**ee-1))*100),'median_speed_error_percent':float(np.median(abs(10**(ee/2)-1))*100),'rows':len(ii),'galaxies':len(gg)}

LOG=[];SELECTED=[]
def fit(kind,ii,split,g='shared',broad=False):
 if kind.startswith('mixed_regular'):spec=REGSPEC
 elif kind.startswith('RAR'):
  spec={'lo':[-3.],'hi':[3.],'starts':[[A0],[-.5],[.5]]}
  if kind=='RAR_stellar':
   low=np.log10(lo_ml[g]);hi=np.log10(2.);spec={'lo':[low],'hi':[hi],'starts':[[np.clip(np.log10(lam),low+1e-9,hi-1e-9)] for lam in [.5,1,1.5]]}
 else:
  mid=np.log10(np.median(rad[ii]));delta=max(1e-6,float(np.median(np.maximum(go[ii]-gb[ii],1e-6))));spec={'lo':[-6.,-2.],'hi':[5.,3.],'starts':[]}
  for off,lam in zip([-.5,0,.5],[.5,1,1.5]):
   rs=np.clip(mid+off,-1.99,2.99);xx=10**(mid-rs);aa=np.log10(delta/(massf(kind.split('_')[0],xx)/(xx*xx*massf(kind.split('_')[0],1.))));spec['starts'].append([np.clip(aa,-5.99,4.99),rs])
  if kind.endswith('_stellar'):
   low=np.log10(lo_ml[g]);hi=np.log10(2.);spec['lo'].append(low);spec['hi'].append(hi)
   for st,lam in zip(spec['starts'],[.5,1,1.5]):st.append(float(np.clip(np.log10(lam),low+1e-9,hi-1e-9)))
 sw=np.sqrt(weights(ii,broad) if g=='shared' else np.ones(len(ii))/len(ii));candidates=[]
 tag=f'{split}:{kind}:{g}'
 for j,start in enumerate(spec['starts']):
  rr=least_squares(lambda c:sw*(predict(kind,c,ii)-y[ii]),start,bounds=(spec['lo'],spec['hi']),max_nfev=2500,ftol=1e-10,xtol=1e-10,gtol=1e-10)
  record={'id':len(LOG),'tag':tag,'split':split,'model':kind,'galaxy':g,'start':j,'initial':start,'bounds':[spec['lo'],spec['hi']],'coefficients':rr.x.tolist(),'active_bounds':rr.active_mask.tolist(),'success':bool(rr.success),'message':str(rr.message),'nfev':int(rr.nfev),'cost':float(rr.cost),'optimality':float(rr.optimality),'train_rows':ii.tolist(),'weighting':'equal_population_galaxy_point' if broad else 'equal_galaxy_point' if g=='shared' else 'equal_point'};LOG.append(record);candidates.append(record)
 best=min(candidates,key=lambda d:d['cost']);SELECTED.append(best['id']);return np.array(best['coefficients']),best['id']

PERGAL=['RAR_individual','NFW_fixed','Burkert_fixed','RAR_stellar','NFW_stellar','Burkert_stellar']
SHARED=['RAR_shared','mixed_regular_disk','mixed_regular_broad']
KINDS=['RAR_fixed']+SHARED+PERGAL
NPARAM={'RAR_fixed':'0','RAR_shared':'1 shared','mixed_regular_disk':'7 shared','mixed_regular_broad':'7 shared','RAR_individual':'1 per galaxy (104 total)','NFW_fixed':'2 per galaxy (208 total)','Burkert_fixed':'2 per galaxy (208 total)','RAR_stellar':'1 per galaxy (104 total)','NFW_stellar':'3 per galaxy (312 total)','Burkert_stellar':'3 per galaxy (312 total)'}

def save_prediction(split,kind,pr,te,fitid):
 with (O/f'{split}_{kind}.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['row','system','group','radius_kpc','gbar_m_s2','gobs_m_s2','gpred_m_s2','log10_residual','split','model','selected_fit_id'])
  for i,p,cid in zip(te,pr,fitid):w.writerow([i,gal[i],group[i],rad[i],gb[i]*1e-10,go[i]*1e-10,10**p*1e-10,p-y[i],split,kind,cid])

def bootstrap(preds,te):
 gg=np.unique(gal[te]);rng=np.random.default_rng(20261028);idx=rng.integers(0,len(gg),size=(5000,len(gg)));scores={}
 for k,pr in preds.items():
  e=pr-y[te];eg=np.array([np.mean(e[gal[te]==g]**2) for g in gg]);scores[k]=np.sqrt(np.mean(eg[idx],axis=1))
 out={}
 for k in KINDS:
  for b in ['RAR_shared','NFW_fixed','Burkert_fixed']:
   if k==b:continue
   dif=scores[k]-scores[b];out[k+' minus '+b]={'difference_dex':metric(preds[k],te)['galaxy_RMS_dex']-metric(preds[b],te)['galaxy_RMS_dex'],'percentile_95_dex':np.percentile(dif,[2.5,97.5]).tolist()}
 return out

def main():
 out={'status':'Exploratory contiguous radial prediction inside observed galaxies; no unseen-system or cosmological comparison.','protocol_sha256':sha(R/'PROTOCOL.txt'),'script_sha256':sha(Path(__file__)),'input_sha256':sha(IN),'raw_rotmod_sha256':sha(RAW),'matched_components_sha256':sha(R/'matched_components.csv'),'eligibility':{'galaxies':ELIG,'count':len(ELIG),'rows':int(sum(np.isin(gal,ELIG)&(pop=='SPARC'))),'excluded_fewer_than_9_radii':[g for g in SG if g not in ELIG]},'parameter_counts':NPARAM,'stellar_bounds_by_galaxy':{g:[lo_ml[g],2.] for g in ELIG},'splits':{},'fit_log':LOG,'selected_ids':SELECTED}
 for split in ['outer_third','middle_third']:
  t=time.time();tests=[];trg={};teg={}
  for g in ELIG:
   ii=I[(pop=='SPARC')&(gal==g)];ii=ii[np.argsort(rad[ii],kind='stable')];assert np.all(np.diff(rad[ii])>0);k=len(ii)//3;te=ii[-k:] if split=='outer_third' else ii[k:2*k];tr=np.array([i for i in ii if i not in te]);assert len(tr)>=6 and len(te)>=3;trg[g]=tr;teg[g]=te;tests.extend(te.tolist())
  te=np.array(sorted(tests));tr=np.array([i for i in I if i not in set(te)]);trdisk=tr[pop[tr]=='SPARC'];preds={};dd={'test_rows':te.tolist(),'disk_train_rows':trdisk.tolist(),'broad_train_rows':tr.tolist(),'galaxy_train_rows':{g:z.tolist() for g,z in trg.items()},'galaxy_test_rows':{g:z.tolist() for g,z in teg.items()},'models':{}}
  out['splits'][split]=dd
  for kind in KINDS:
   pr=np.full(len(te),np.nan);fitids=np.full(len(te),-1);co={};trainpred=[];trainidx=[]
   if kind=='RAR_fixed':
    pr=predict(kind,[],te);tr_eval=np.concatenate(list(trg.values()));trainpred=predict(kind,[],tr_eval);trainidx=tr_eval
   elif kind in SHARED:
    train=tr if kind=='mixed_regular_broad' else trdisk;c,cid=fit(kind,train,split,broad=kind=='mixed_regular_broad');pr=predict(kind,c,te);fitids[:]=cid;co={'shared':c.tolist()};tr_eval=np.concatenate(list(trg.values()));trainpred=predict(kind,c,tr_eval);trainidx=tr_eval
   else:
    for g in ELIG:
     c,cid=fit(kind,trg[g],split,g);zz=gal[te]==g;pr[zz]=predict(kind,c,te[zz]);fitids[zz]=cid;co[g]=c.tolist();trainpred.extend(predict(kind,c,trg[g]).tolist());trainidx.extend(trg[g].tolist())
   assert np.all(np.isfinite(pr));preds[kind]=pr;trainidx=np.array(trainidx,int);trainpred=np.array(trainpred,float)
   dd['models'][kind]={'test':metric(pr,te),'train_eligible_radii':metric(trainpred,trainidx),'coefficients':co,'selected_fit_ids':sorted(set(fitids.tolist())),'by_galaxy':{g:metric(pr[gal[te]==g],te[gal[te]==g]) for g in ELIG}}
   save_prediction(split,kind,pr,te,fitids);(O/'results.json').write_text(json.dumps(out,indent=2));print(split,kind,dd['models'][kind]['test'],flush=True)
  dd['paired_galaxy_bootstrap']=bootstrap(preds,te);print('SPLITDONE',split,'seconds',time.time()-t,flush=True)
 out['fit_sets']=len(SELECTED);out['optimizer_calls']=len(LOG);out['optimizer_failures']=[d['id'] for d in LOG if not d['success']];out['selected_optimizer_failures']=[i for i in SELECTED if not LOG[i]['success']]
 (O/'results.json').write_text(json.dumps(out,indent=2));print('DONE',out['fit_sets'],out['optimizer_calls'],'failed',out['optimizer_failures'],'selected_failed',out['selected_optimizer_failures'],flush=True)
if __name__=='__main__':main()
