"""Frozen two-family Dutton-Maccio c(M) test; see PROTOCOL.txt."""
from pathlib import Path
import csv, json, hashlib, time, shutil
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
COMPARE=ROOT.parent/'dm_shared_halo/results'
G=6.67430e-11; MSUN=1.988409870698051e30; KPC=3.085677581491367e19
H=0.671; H0=67.1e3/(1000*KPC); FB=0.0490/0.3175
SEEDS=[20261027,7331,48105]
SPEC={
 'NFW_Planck_mass':{'lo':[8.,.1],'hi':[14.,2.],'starts':[[11.5,.7],[12.,1.2],[10.5,.4]],'parameter_names':['A','B']},
 'NFW_Planck_size':{'lo':[8.,.1,-1.5],'hi':[14.,2.,1.5],'starts':[[11.5,.7,0.],[12.,1.2,.5],[10.5,.4,-.5]],'parameter_names':['A','B','C']}}
LOG=[]
rows=list(csv.DictReader((ROOT/'input.csv').open()))
N=len(rows); I=np.arange(N)
gal=np.array([r['system'] for r in rows]); group=np.array([r['group'] for r in rows])
GS=sorted(set(gal)); assert N==2700 and len(GS)==149 and len(set(group))==149
rad=np.array([float(r['radius_kpc']) for r in rows])
gb=np.array([float(r['gbar_m_s2']) for r in rows])
go=np.array([float(r['gobs_m_s2']) for r in rows]); y=np.log10(go/1e-10)
mp=np.array([float(r['mass_proxy_msun']) for r in rows]); re=np.array([float(r['reff_kpc']) for r in rows])
m=np.log10(mp/1e10); hh=np.log10(re/3.)
gf={r['system']:r for r in csv.DictReader((ROOT/'galaxy_features.csv').open())}
assert set(gf)==set(GS) and np.all(gb>0) and np.all(go>0) and np.all(rad>0)
assert np.array_equal(np.arange(N),np.array([int(r['row']) for r in rows]))

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writecsv(p,rr):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def massf(z):
 z=np.asarray(z,dtype=float); out=np.empty_like(z); small=z<.01
 zz=z[small];out[small]=sum(((-1)**n)*(n-1)/n*zz**n for n in range(2,13))
 zz=z[~small];out[~small]=np.log1p(zz)-zz/(1+zz)
 return out
def halo(k,co,ii):
 exponent=co[0]+co[1]*m[ii]
 if k.endswith('size'): exponent=exponent+co[2]*hh[ii]
 mr=mp[ii]/FB+10**exponent
 c=10**(.905-.101*np.log10(mr*H/1e12))
 r200=(G*MSUN*mr/(100*H0**2))**(1/3)/KPC
 rs=r200/c; z=np.minimum(rad[ii]/rs,c)
 gd=G*MSUN*(1-FB)*mr*massf(z)/(massf(c)*(rad[ii]*KPC)**2)
 return mr,c,r200,rs,gd
def predict(k,co,ii):return np.log10((gb[ii]+halo(k,co,ii)[-1])/1e-10)
def weights(ii):
 w=np.empty(len(ii));gg=np.unique(group[ii])
 for g in gg:
  z=group[ii]==g;w[z]=1/(len(gg)*sum(z))
 return w
def metrics(pr,ii):
 e=pr-y[ii]
 return {'galaxy_RMS_dex':float(np.sqrt(np.sum(weights(ii)*e**2))),
         'row_RMS_dex':float(np.sqrt(np.mean(e**2))),
         'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),
         'median_speed_error_percent':float(np.median(abs(10**(.5*e)-1))*100),
         'rows':len(ii),'galaxies':len(set(gal[ii]))}
def score(pr,ii):return {'metrics':metrics(pr,ii),'by_galaxy':{g:metrics(pr[gal[ii]==g],ii[gal[ii]==g]) for g in sorted(set(gal[ii]))}}
def fit(k,ii,tag):
 sw=np.sqrt(weights(ii));sp=SPEC[k]; rr=[];ids=[]
 for j,start in enumerate(sp['starts']):
  z=least_squares(lambda co:sw*(predict(k,co,ii)-y[ii]),start,bounds=(sp['lo'],sp['hi']),max_nfev=3000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
  ids.append(len(LOG));LOG.append({'model':k,'split':tag,'start_index':j,'train_rows':ii.tolist(),'train_galaxies':len(set(gal[ii])),'coefficients':z.x.tolist(),'success':bool(z.success),'message':z.message,'cost':float(z.cost),'optimality':float(z.optimality),'nfev':z.nfev,'active_bounds':z.active_mask.tolist()});rr.append(z)
 q=int(np.argmin([z.cost for z in rr]));return rr[q].x,ids[q]
def folds(seed):
 rng=np.random.default_rng(seed);fm={g:j%5 for j,g in enumerate(rng.permutation(sorted(set(group))))}
 return np.array([fm[g] for g in group])
def diagnostics(k,co,ii):
 mr,c,r200,rs,gd=halo(k,co,ii); outside=(mr*H<1e10)|(mr*H>1e15)
 def mm(a):return [float(np.min(a)),float(np.max(a))]
 return {'Mref_msun_minmax':mm(mr),'c200_minmax':mm(c),'r200_kpc_minmax':mm(r200),'rs_kpc_minmax':mm(rs),'proxy_retention_fraction_minmax':mm(mp[ii]/(FB*mr)),
  'reference_cM_calibration_extrapolated_galaxies':sorted(set(gal[ii][outside])),
  'reference_cM_calibration_extrapolated_rows':int(sum(outside)),
  'rows_outside_reference_r200':int(sum(rad[ii]>r200)),
  'galaxies_outside_reference_r200':sorted(set(gal[ii][rad[ii]>r200])),
  'maximum_observed_radius_over_reference_r200':float(np.max(rad[ii]/r200)),
  'min_halo_acceleration_m_s2':float(np.min(gd)),
  'formal_mass_zero_halo_acceleration_vanishes':bool(co[1]>0),
  'finite_total_dark_mass':True,'deterministic_cosmological_structural_prior':True,'full_LCDM_prediction':False}
def save(k,tag,pr,ii,fd,co_by_fold):
 tab=[]
 for f,co in co_by_fold.items():
  jj=np.flatnonzero(fd==f); mr,c,r200,rs,gd=halo(k,co,ii[jj])
  for z,j in enumerate(jj):
   i=ii[j]
   tab.append({'row':int(i),'broad_row':rows[i]['broad_row'],'system':gal[i],'group':group[i],'fold':str(fd[j]),'radius_kpc':rad[i],'gbar_m_s2':gb[i],'gobs_m_s2':go[i],'mass_proxy_msun':mp[i],'reff_kpc':re[i],'gpred_m_s2':10**pr[j]*1e-10,'log10_residual':pr[j]-y[i],
    'Mref_msun':mr[z],'c200':c[z],'reference_r200_kpc':r200[z],'rs_kpc':rs[z],'halo_g_m_s2':gd[z],'reference_halo_truncated_here':bool(rad[i]>r200[z]),'cM_mass_extrapolation':bool(mr[z]*H<1e10 or mr[z]*H>1e15)})
 tab.sort(key=lambda r:r['row']);writecsv(OUT/(k+'_'+tag+'.csv'),tab)
def getcomparison(model,tag,ii,fd):
 p=COMPARE/(model+'_'+tag+'.csv'); zz=list(csv.DictReader(p.open()))
 assert len(zz)==len(ii)
 for z,i,f in zip(zz,ii,fd):
  assert int(z['row'])==i and z['system']==gal[i] and z['group']==group[i]
  assert z['fold']==str(f)
  for key,arr in [('radius_kpc',rad),('gbar_m_s2',gb),('gobs_m_s2',go),('mass_proxy_msun',mp),('reff_kpc',re)]:
   assert np.isclose(float(z[key]),arr[i],rtol=1e-12,atol=0), (model,tag,key,i)
 # Persist exact comparator tables rather than only the path/hash.
 d=ROOT/'comparators';d.mkdir(exist_ok=True);shutil.copyfile(p,d/p.name)
 return np.log10(np.array([float(z['gpred_m_s2']) for z in zz])/1e-10),{'file':str(p),'copied_file':str((d/p.name).relative_to(ROOT)),'sha256':sha(p),'verified_rows_folds_features':True}
def paired_bootstrap(pred,ii,reference):
 gg=sorted(set(gal[ii]));rng=np.random.default_rng(20261029);multi=rng.poisson(1,size=(5000,len(gg)));den=multi.sum(axis=1);assert np.all(den>0)
 scores={}
 for k,p in pred.items():
  e=np.array([np.mean((p[gal[ii]==g]-y[ii][gal[ii]==g])**2) for g in gg]);scores[k]=np.sqrt(multi@e/den)
 out={}
 for k in SPEC:
  d=scores[k]-scores[reference]
  out[k+' minus '+reference]={'difference_dex':metrics(pred[k],ii)['galaxy_RMS_dex']-metrics(pred[reference],ii)['galaxy_RMS_dex'],'percentile_95_dex':np.percentile(d,[2.5,97.5]).tolist()}
 return out

def main():
 tails={}
 for nm,field in [('mass','mass_proxy_msun'),('surface','sigma_proxy_msun_kpc2')]:
  order=sorted(GS,key=lambda g:(float(gf[g][field]),g))
  for side,gg in [('low',order[:37]),('high',order[-37:])]:tails[nm+'_'+side]=I[np.isin(gal,gg)]
 allsplits={f'cv_{s}':(I,folds(s)) for s in SEEDS}
 allsplits.update({name:(ii,np.array([name]*len(ii))) for name,ii in tails.items()})
 result={'status':'Exploratory deterministic published-cM-prior baseline; no full LCDM, all-scale, or novelty claim.',
  'input_sha256':sha(ROOT/'input.csv'),'features_sha256':sha(ROOT/'galaxy_features.csv'),'protocol_sha256':sha(ROOT/'PROTOCOL.txt'),'script_sha256':sha(__file__),
  'rows':N,'galaxies':len(GS),'constants':{'G_SI':G,'Msun_kg':MSUN,'kpc_m':KPC,'h':H,'H0_SI':H0,'fb':FB,'cM_intercept':.905,'cM_slope':-.101},'spec':SPEC,'models':{},'fit_log':LOG,
  'tails':{name:{'rows':ii.tolist(),'galaxies':sorted(set(gal[ii]))} for name,ii in tails.items()},'comparator_sources':{},'comparison_metrics':{},'bootstrap':{}}
 preds={tag:{} for tag in allsplits}
 for k in SPEC:
  st=time.time();co,fid=fit(k,I,'all');pr=predict(k,co,I)
  save(k,'all',pr,I,np.array(['']*N),{'':co})
  d={'coefficients':co.tolist(),'selected_all_fit':fid,'all':score(pr,I),'physical_diagnostics':diagnostics(k,co,I),'CV':{},'tail_holdout':{}}
  for seed in SEEDS:
   fd=folds(seed);pr=np.full(N,np.nan);cs=[];sel=[];diags={}
   for f in range(5):
    tr=I[fd!=f];te=I[fd==f];assert not set(gal[tr])&set(gal[te]);cc,fitid=fit(k,tr,f'cv_{seed}_{f}');pr[te]=predict(k,cc,te);cs.append(cc.tolist());sel.append(fitid);diags[str(f)]=diagnostics(k,cc,te)
   d['CV'][str(seed)]={**score(pr,I),'coefficients':cs,'selected_fit_ids':sel,'test_physical_diagnostics':diags}
   save(k,f'cv_{seed}',pr,I,fd,{i:cs[i] for i in range(5)});preds[f'cv_{seed}'][k]=pr
  for nm,te in tails.items():
   tr=I[~np.isin(gal,gal[te])];assert not set(gal[tr])&set(gal[te]);cc,fitid=fit(k,tr,nm);pr=predict(k,cc,te)
   d['tail_holdout'][nm]={**score(pr,te),'coefficients':cc.tolist(),'selected_fit_id':fitid,'test_physical_diagnostics':diagnostics(k,cc,te)}
   save(k,nm,pr,te,np.array([nm]*len(te)),{nm:cc});preds[nm][k]=pr
  result['models'][k]=d
  (OUT/'results.json').write_text(json.dumps(result,indent=2))
  print(k,'CV',[round(d['CV'][str(s)]['metrics']['galaxy_RMS_dex'],6) for s in SEEDS],'speed',round(d['CV'][str(SEEDS[0])]['metrics']['median_speed_error_percent'],5),'tails',[round(d['tail_holdout'][s]['metrics']['galaxy_RMS_dex'],6) for s in tails],'coef',co.tolist(),'seconds',round(time.time()-st,2),flush=True)
 for tag,(ii,fd) in allsplits.items():
  result['comparator_sources'][tag]={};result['comparison_metrics'][tag]={}
  for k in ['RAR','RAR_fixed','NFW_mass','NFW_size','Burkert_mass','Burkert_size','mixed_regular']:
   pp,source=getcomparison(k,tag,ii,fd);preds[tag][k]=pp;result['comparator_sources'][tag][k]=source;result['comparison_metrics'][tag][k]=metrics(pp,ii)
  for k in SPEC:result['comparison_metrics'][tag][k]=metrics(preds[tag][k],ii)
  result['bootstrap'][tag]={}
  for ref in ['RAR','NFW_mass','NFW_size']:result['bootstrap'][tag].update(paired_bootstrap(preds[tag],ii,ref))
 result['fit_sets']=sum(x['start_index']==0 for x in LOG);result['optimizer_attempts']=len(LOG);result['optimizer_failures']=[x for x in LOG if not x['success']]
 result['new_prediction_tables']=len(list(OUT.glob('*.csv')))
 (OUT/'results.json').write_text(json.dumps(result,indent=2))
 summary=[]
 for split,mm in result['comparison_metrics'].items():
  for model,met in mm.items():summary.append({'split':split,'model':model,**met})
 writecsv(OUT/'metrics_summary.csv',summary)
 (ROOT/'manifest.json').write_text(json.dumps({'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='manifest.json'}},indent=2))
 print('DONE',result['fit_sets'],'sets',result['optimizer_attempts'],'attempts',len(result['optimizer_failures']),'failures',flush=True)

if __name__=='__main__':main()
