from pathlib import Path
import csv,json,hashlib,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
R=Path(__file__).resolve().parent;OUT=R/'results';OUT.mkdir(exist_ok=True)
INPUT=R.parent/'dm_unified_regular/augmented_input.csv';CAT=R.parent/'dm_fit/data/SPARC_Lelli2016c.mrt'
LN=np.log(10.);I=None;LOG=[]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writecsv(p,rr):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)

def prepare():
 global rows,I,N,group,gal,rad,gb,go,x,u,y,m,h,galfeatures,GS
 cats={}
 for line in CAT.read_text().splitlines():
  try:
   a=line.split()
   if len(a)!=19 or not a[1].isdigit():continue
   name=a[0];lum=float(a[7]);reff=float(a[9]);mhi=float(a[13])
  except (ValueError,IndexError):continue
  cats[name]={'system':name,'luminosity_3p6_1e9_lsun':lum,'mhi_1e9_msun':mhi,'reff_kpc':reff,'mass_proxy_msun':1e9*(.5*lum+1.33*mhi),'sigma_proxy_msun_kpc2':1e9*(.5*lum+1.33*mhi)/(2*np.pi*reff**2)}
 rows=[]
 for j,z in enumerate(csv.DictReader(INPUT.open())):
  if z['population']!='SPARC':continue
  f=cats[z['system']];assert f['mass_proxy_msun']>0 and f['reff_kpc']>0
  rows.append(dict(row=len(rows),broad_row=j,**z,**{k:v for k,v in f.items() if k!='system'}))
 N=len(rows);assert N==2700;I=np.arange(N);gal=np.array([r['system'] for r in rows]);group=np.array([r['group'] for r in rows]);GS=sorted(set(gal));assert len(GS)==149
 assert len(set(group))==149
 galfeatures={g:cats[g] for g in GS};writecsv(R/'galaxy_features.csv',list(galfeatures.values()));writecsv(R/'input.csv',rows)
 rad=np.array([float(r['radius_kpc']) for r in rows]);gb=np.array([float(r['gbar_m_s2'])/1e-10 for r in rows]);go=np.array([float(r['gobs_m_s2'])/1e-10 for r in rows]);x=np.log10(gb);u=np.log10(rad);y=np.log10(go)
 m=np.log10([r['mass_proxy_msun']/1e10 for r in rows]);h=np.log10([r['reff_kpc']/3 for r in rows])

KINDS=['RAR_fixed','RAR','mixed_regular','NFW_mass','Burkert_mass','NFW_size','Burkert_size']
SPEC={'RAR':{'lo':[-3.],'hi':[3.],'starts':[[.1]]},'mixed_regular':{'lo':[-3.,-5.,-5.,0.,0.,-2.,0.],'hi':[3.,3.,3.,.95,.95,0.,4.], 'starts':[[-.1,-2.,-.2,.6,.7,-2.,2.5],[0.,-1.5,-.5,.3,.3,-1.,1.5],[-.4,-3.,-.7,.9,.9,-1.5,3.5]]}}
for k in ['NFW_mass','Burkert_mass']:SPEC[k]={'lo':[-3,-1.5,-1,-1.5],'hi':[2,1.5,2,1.5],'starts':[[-.5,0,1,.4],[0,.4,.5,.2],[-1,-.2,1.5,.6]]}
for k in ['NFW_size','Burkert_size']:SPEC[k]={'lo':[-3,-1.5,-1.5,-1,-1.5,-1.5],'hi':[2,1.5,1.5,2,1.5,1.5],'starts':[[-.5,0,0,1,.4,0],[0,.4,-.5,.5,.2,.5],[-1,-.2,.5,1.5,.6,-.3]]}

def massf(kind,z):
 z=np.asarray(z,dtype=float);ans=np.empty_like(z);small=z<.01
 if kind=='NFW':
  zs=z[small];ans[small]=sum(((-1)**n)*(n-1)/n*zs**n for n in range(2,13));zl=z[~small];ans[~small]=np.log1p(zl)-zl/(1+zl)
 else:
  zs=z[small];ans[small]=sum(zs**(4*k+3)/(4*k+3)-zs**(4*k+4)/(4*k+4) for k in range(4));zl=z[~small];ans[~small]=.5*np.log1p(zl)+.25*np.log1p(zl*zl)-.5*np.arctan(zl)
 return ans

def rar(xx,aa):
 t=10**np.clip(.5*(xx-aa),-290,290);return xx-np.log10(-np.expm1(-t))
def halopars(k,c,mm,hh):
 if k.endswith('mass'):return c[0]+c[1]*mm,c[2]+c[3]*mm
 return c[0]+c[1]*mm+c[2]*hh,c[3]+c[4]*mm+c[5]*hh
def rawpredict(k,c,xx,uu,mm,hh):
 if k=='RAR_fixed':return rar(xx,np.log10(1.2))
 if k=='RAR':return rar(xx,c[0])
 if k=='mixed_regular':
  inv=c[1]-c[3]*(xx+2*uu)-np.logaddexp(0,LN*4*(c[5]-uu))/LN-np.logaddexp(0,LN*2*(xx-2))/LN
  pos=c[2]+c[4]*(xx+uu)-np.logaddexp(0,LN*(uu-c[6]))/LN
  aa=logsumexp(LN*np.array([np.zeros_like(xx)+c[0],inv,pos]),axis=0)/LN;return rar(xx,aa)
 la,lr=halopars(k,c,mm,hh);z=10**(uu-lr);kind=k.split('_')[0];dm=10**la*massf(kind,z)/(z*z*massf(kind,1.))
 return np.log10(10**xx+dm)
def predict(k,c,ii):return rawpredict(k,c,x[ii],u[ii],m[ii],h[ii])
def weights(ii):
 w=np.empty(len(ii));gg=np.unique(group[ii])
 for g in gg:
  z=group[ii]==g;w[z]=1/(len(gg)*sum(z))
 return w
def metrics(pr,ii):
 e=pr-y[ii];return {'galaxy_RMS_dex':float(np.sqrt(np.sum(weights(ii)*e*e))),'row_RMS_dex':float(np.sqrt(np.mean(e*e))),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'median_speed_error_percent':float(np.median(abs(10**(.5*e)-1))*100),'rows':len(ii),'galaxies':len(set(gal[ii]))}
def fit(k,ii,tag):
 if k=='RAR_fixed':return [],None
 sw=np.sqrt(weights(ii));spec=SPEC[k];rrs=[];ids=[]
 for j,start in enumerate(spec['starts']):
  rr=least_squares(lambda c:sw*(predict(k,c,ii)-y[ii]),start,bounds=(spec['lo'],spec['hi']),max_nfev=3000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
  ids.append(len(LOG));LOG.append({'model':k,'split':tag,'start_index':j,'train_rows':ii.tolist(),'train_galaxies':len(set(gal[ii])),'coefficients':rr.x.tolist(),'success':bool(rr.success),'message':rr.message,'cost':float(rr.cost),'optimality':float(rr.optimality),'nfev':rr.nfev,'active_bounds':rr.active_mask.tolist()});rrs.append(rr)
 j=int(np.argmin([v.cost for v in rrs]));return rrs[j].x,ids[j]
def folds(seed):
 rng=np.random.default_rng(seed);mp={g:j%5 for j,g in enumerate(rng.permutation(sorted(set(group))))};return np.array([mp[g] for g in group])
def save(k,tag,pr,ii,fd=None):
 rr=[]
 for j,i in enumerate(ii):rr.append({'row':int(i),'broad_row':rows[i]['broad_row'],'system':gal[i],'group':group[i],'fold':'' if fd is None else str(fd[j]),'radius_kpc':rad[i],'gbar_m_s2':gb[i]*1e-10,'gobs_m_s2':go[i]*1e-10,'mass_proxy_msun':rows[i]['mass_proxy_msun'],'reff_kpc':rows[i]['reff_kpc'],'gpred_m_s2':10**pr[j]*1e-10,'log10_residual':pr[j]-y[i]})
 writecsv(OUT/(k+'_'+tag+'.csv'),rr)
def score(pr,ii):return {'metrics':metrics(pr,ii),'by_galaxy':{g:metrics(pr[gal[ii]==g],ii[gal[ii]==g]) for g in sorted(set(gal[ii]))}}
def boot(pred,ii):
 gg=sorted(set(gal[ii]));rng=np.random.default_rng(20261029);mult=rng.poisson(1,size=(5000,len(gg)));den=mult.sum(1);assert np.all(den>0);scores={}
 for k,p in pred.items():
  ee=np.array([np.mean((p[gal[ii]==g]-y[ii][gal[ii]==g])**2) for g in gg]);scores[k]=np.sqrt(mult@ee/den)
 out={}
 for k in pred:
  if k=='RAR':continue
  z=scores[k]-scores['RAR'];out[k+' minus RAR']={'difference_dex':metrics(pred[k],ii)['galaxy_RMS_dex']-metrics(pred['RAR'],ii)['galaxy_RMS_dex'],'percentile_95_dex':np.percentile(z,[2.5,97.5]).tolist()}
 return out

def main():
 prepare();tails={}
 for name,field in [('mass','mass_proxy_msun'),('surface','sigma_proxy_msun_kpc2')]:
  order=sorted(GS,key=lambda g:(galfeatures[g][field],g))
  for side,gg in [('low',order[:37]),('high',order[-37:])]:tails[name+'_'+side]=I[np.isin(gal,gg)]
 result={'status':'Exploratory shared photometry-to-halo relations; no per-galaxy fitted parameters; no all-scale or novelty claim.','input_source_sha256':sha(INPUT),'catalog_sha256':sha(CAT),'input_sha256':sha(R/'input.csv'),'features_sha256':sha(R/'galaxy_features.csv'),'protocol_sha256':sha(R/'PROTOCOL.txt'),'script_sha256':sha(Path(__file__)),'rows':N,'galaxies':len(GS),'models':{},'spec':SPEC,'fit_log':LOG,'tails':{k:{'rows':ii.tolist(),'galaxies':sorted(set(gal[ii])),'feature_min':{f:min(galfeatures[g][f] for g in set(gal[ii])) for f in ['mass_proxy_msun','sigma_proxy_msun_kpc2']},'feature_max':{f:max(galfeatures[g][f] for g in set(gal[ii])) for f in ['mass_proxy_msun','sigma_proxy_msun_kpc2']}} for k,ii in tails.items()}}
 preds={f'cv_{s}':{} for s in [20261027,7331,48105]};preds.update({k:{} for k in tails})
 for k in KINDS:
  st=time.time();co,fitid=fit(k,I,'all');pr=predict(k,co,I);save(k,'all',pr,I);d={'coefficients':list(co),'selected_all_fit':fitid,'all':score(pr,I),'CV':{},'tail_holdout':{}}
  for seed in [20261027,7331,48105]:
   fd=folds(seed);pr=np.full(N,np.nan);cs=[];selected=[]
   for f in range(5):
    tr=I[fd!=f];te=I[fd==f];assert not set(gal[tr])&set(gal[te]);cc,fid=fit(k,tr,f'cv_{seed}_{f}');pr[te]=predict(k,cc,te);cs.append(list(cc));selected.append(fid)
   d['CV'][str(seed)]={**score(pr,I),'coefficients':cs,'selected_fit_ids':selected};save(k,f'cv_{seed}',pr,I,fd);preds[f'cv_{seed}'][k]=pr
  for name,te in tails.items():
   tr=I[~np.isin(gal,gal[te])];assert not set(gal[tr])&set(gal[te]);cc,fid=fit(k,tr,name);pr=predict(k,cc,te);d['tail_holdout'][name]={**score(pr,te),'coefficients':list(cc),'selected_fit_id':fid};save(k,name,pr,te,np.array([name]*len(te)));preds[name][k]=pr
  if k.startswith(('NFW','Burkert')):
   la,lr=halopars(k,co,m,h);a=co[1];b=co[3] if k.endswith('mass') else co[4];exponent=a+2*b if b>0 else a if k.startswith('NFW') else a-b
   d['physical_diagnostics']={'implied_A_m_s2_minmax':(10**np.array([la.min(),la.max()])*1e-10).tolist(),'implied_rs_kpc_minmax':(10**np.array([lr.min(),lr.max()])).tolist(),'mass_zero_limit_fixed_Reff_and_r_halo_power':float(exponent),'mass_zero_limit_has_log_factor':bool(b>0),'halo_force_vanishes_in_formal_mass_zero_limit':bool(exponent>0),'outer_untruncated_mass_diverges':True,'cosmological_prior_imposed':False}
  result['models'][k]=d;(OUT/'results.json').write_text(json.dumps(result,indent=2));print(k,'CV',*[round(d['CV'][str(s)]['metrics']['galaxy_RMS_dex'],6) for s in [20261027,7331,48105]],'speed',round(d['CV']['20261027']['metrics']['median_speed_error_percent'],4),'tails',*[round(d['tail_holdout'][z]['metrics']['galaxy_RMS_dex'],6) for z in tails],'coef',list(co),'sec',round(time.time()-st,1),flush=True)
 result['bootstrap']={k:boot(v,I if k.startswith('cv_') else tails[k]) for k,v in preds.items()}
 result['fit_sets']=sum(z['start_index']==0 for z in LOG);result['optimizer_attempts']=len(LOG);result['optimizer_failures']=[z for z in LOG if not z['success']];result['prediction_tables']=len(list(OUT.glob('*.csv')))
 (OUT/'results.json').write_text(json.dumps(result,indent=2));print('DONE',result['fit_sets'],result['optimizer_attempts'],'failures',len(result['optimizer_failures']),flush=True)
if __name__=='__main__':main()
