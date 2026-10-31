"""Independent shared-halo audit; does not import analysis prediction code."""
from pathlib import Path
import csv,json,hashlib,re,collections
import numpy as np
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
B=Path(__file__).resolve().parent.parent;D=B/'dm_shared_halo';O=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csvs(p):return list(csv.DictReader(p.open()))
J=json.loads((D/'results/results.json').read_text());rows=csvs(D/'input.csv');features=csvs(D/'galaxy_features.csv');broad=csvs(B/'dm_unified_regular/augmented_input.csv')
assert J['input_sha256']==sha(D/'input.csv') and J['features_sha256']==sha(D/'galaxy_features.csv')
assert J['catalog_sha256']==sha(B/'dm_fit/data/SPARC_Lelli2016c.mrt') and J['input_source_sha256']==sha(B/'dm_unified_regular/augmented_input.csv')
assert J['protocol_sha256']==sha(D/'PROTOCOL.txt') and J['script_sha256']==sha(D/'run.py')
# Independent lexical parse: accept only data rows following final separator;
# enumerate the documented named fields before selecting predictor columns.
text=(B/'dm_fit/data/SPARC_Lelli2016c.mrt').read_text();body=text.rsplit('--------------------------------------------------------------------------------',1)[1]
names=['name','T','distance','distance_error','method','inclination','inclination_error','L','L_error','Reff','SB_eff','Rdisk','SB_disk','MHI','RHI','Vflat','Vflat_error','Q','ref']
cat={}
for line in body.strip().splitlines():
 a=re.findall(r'\S+',line);assert len(a)==len(names);v=dict(zip(names,a));cat[v['name']]=v
assert len(cat)==175
ff={v['system']:v for v in features}
for g,v in ff.items():
 q=cat[g];lum,hi,r=map(float,[q['L'],q['MHI'],q['Reff']]);mass=(.5*lum+1.33*hi)*1e9
 for field,val in [('luminosity_3p6_1e9_lsun',lum),('mhi_1e9_msun',hi),('reff_kpc',r),('mass_proxy_msun',mass),('sigma_proxy_msun_kpc2',mass/(2*np.pi*r*r))]:assert np.isclose(float(v[field]),val,rtol=2e-15,atol=0),(g,field)
for z in rows:
 q=broad[int(z['broad_row'])];assert q['system']==z['system'] and q['population']=='SPARC'
 for field in ['radius_kpc','gbar_m_s2','gobs_m_s2']:assert float(z[field])==float(q[field])
 for field in ff[z['system']]:
  if field!='system':assert z[field]==ff[z['system']][field]
N=len(rows);I=np.arange(N);gal=np.array([z['system'] for z in rows]);group=np.array([z['group'] for z in rows]);gb=np.array([float(z['gbar_m_s2'])/1e-10 for z in rows]);r=np.array([float(z['radius_kpc']) for z in rows]);target=np.log10([float(z['gobs_m_s2'])/1e-10 for z in rows]);m=np.log10([float(z['mass_proxy_msun'])/1e10 for z in rows]);h=np.log10([float(z['reff_kpc'])/3 for z in rows]);assert N==2700 and len(set(gal))==149
nodes,w=leggauss(32);t=(nodes+1)/2;w=w/2
def massfn(k,z):
 z=np.atleast_1d(np.asarray(z,float));out=np.empty_like(z);tiny=z<.03
 # Independent integration of the density, instead of analysis small-z series.
 zz=z[tiny,None];tt=t[None,:]
 if k=='NFW':out[tiny]=z[tiny]**2*np.sum(w*tt/(1+zz*tt)**2,axis=1);v=z[~tiny];out[~tiny]=np.log1p(v)-v/(1+v)
 else:out[tiny]=z[tiny]**3*np.sum(w*tt**2/((1+zz*tt)*(1+(zz*tt)**2)),axis=1);v=z[~tiny];out[~tiny]=np.log((1+v)**2*(1+v*v))/4-np.arctan(v)/2
 return out
def pred(k,c,ii):
 b=gb[ii];rr=r[ii]
 if k in ['RAR_fixed','RAR']:a=1.2 if k=='RAR_fixed' else 10**c[0]
 elif k=='mixed_regular':
  S=rr**4/(rr**4+10**(4*c[5]));W=1/(1+(b/100)**2)
  a=10**c[0]+10**c[1]*(b*rr**2)**(-c[3])*S*W+10**c[2]*(b*rr)**c[4]/(1+rr/10**c[6])
 else:
  if k.endswith('_mass'):A=10**(c[0]+c[1]*m[ii]);rs=10**(c[2]+c[3]*m[ii])
  else:A=10**(c[0]+c[1]*m[ii]+c[2]*h[ii]);rs=10**(c[3]+c[4]*m[ii]+c[5]*h[ii])
  typ=k.split('_')[0];z=rr/rs;return np.log10(b+A*massfn(typ,z)/(z*z*massfn(typ,[1])[0]))
 return np.log10(b/(-np.expm1(-np.sqrt(b/a))))
def weights(ii):
 gs=np.unique(group[ii]);w=np.zeros(len(ii))
 for g in gs:z=group[ii]==g;w[z]=1/(len(gs)*sum(z))
 return w
def metrics(p,ii):
 e=p-target[ii];return {'galaxy_RMS_dex':np.sqrt(sum(weights(ii)*e*e)),'row_RMS_dex':np.sqrt(np.mean(e*e)),'median_acceleration_error_percent':100*np.median(abs(10**e-1)),'median_speed_error_percent':100*np.median(abs(10**(.5*e)-1)),'rows':len(ii),'galaxies':len(set(gal[ii]))}
def checkmetrics(p,ii,d):
 for k,v in metrics(p,ii).items():assert np.isclose(v,d[k],rtol=1e-10,atol=1e-12),(k,v,d[k])
def folds(seed):
 gs=sorted(set(group));rng=np.random.default_rng(seed);permuted=rng.permutation(gs);parts=[set(permuted[j::5]) for j in range(5)]
 return np.array([next(j for j in range(5) if g in parts[j]) for g in group])
tails={}
for name,f in [('mass','mass_proxy_msun'),('surface','sigma_proxy_msun_kpc2')]:
 gs=sorted(ff,key=lambda g:(float(ff[g][f]),g))
 for suffix,names_ in [('low',gs[:37]),('high',gs[-37:])]:
  name_=name+'_'+suffix;tails[name_]=I[np.isin(gal,names_)];assert set(J['tails'][name_]['galaxies'])==set(names_) and J['tails'][name_]['rows']==tails[name_].tolist()
maxpred=maxcost=0.;tablecount=predictioncount=0;allpred={};trainrecords={}
for k,d in J['models'].items():
 variants=[('all',I,d['coefficients'],d['all'],d['selected_all_fit'])]
 for seed,z in d['CV'].items():variants.append(('cv_'+seed,I,z['coefficients'],z,z['selected_fit_ids']))
 for tag,z in d['tail_holdout'].items():variants.append((tag,tails[tag],z['coefficients'],z,z['selected_fit_id']))
 for tag,ii,co,score,selected in variants:
  cc=csvs(D/'results'/f'{k}_{tag}.csv');assert [int(z['row']) for z in cc]==ii.tolist();out=np.log10([float(z['gpred_m_s2'])/1e-10 for z in cc]);expected=np.empty(len(ii))
  if tag.startswith('cv_'):
   fd=folds(int(tag[3:]));assert np.array_equal(fd,np.array([int(z['fold']) for z in cc]))
   for fold in range(5):
    te=I[fd==fold];tr=I[fd!=fold];assert not set(gal[te])&set(gal[tr]);expected[fd==fold]=pred(k,co[fold],te)
    if selected[fold] is not None:trainrecords[selected[fold]]=(k,f'{tag}_{fold}',tr,co[fold])
  else:
   expected=pred(k,co,ii);tr=I if tag=='all' else I[~np.isin(gal,gal[ii])]
   if tag!='all':assert not set(gal[tr])&set(gal[ii])
   if selected is not None:trainrecords[selected]=(k,tag,tr,co)
  assert np.max(abs(out-expected))<1e-10,(k,tag,np.max(abs(out-expected)))
  maxpred=max(maxpred,float(np.max(abs(out-expected))));checkmetrics(out,ii,score['metrics'])
  for g in set(gal[ii]):zz=gal[ii]==g;checkmetrics(out[zz],ii[zz],score['by_galaxy'][g])
  assert np.allclose([float(z['log10_residual']) for z in cc],out-target[ii],rtol=1e-11,atol=1e-12)
  tablecount+=1;predictioncount+=len(ii);allpred[k,tag]=(ii,out)
attempts=collections.defaultdict(list)
for no,f in enumerate(J['fit_log']):
 k=f['model'];ii=np.array(f['train_rows']);tag=f['split'];spec=J['spec'][k];coef=np.array(f['coefficients']);assert f['success'];assert np.all(coef>=np.array(spec['lo'])-1e-10) and np.all(coef<=np.array(spec['hi'])+1e-10)
 if tag=='all':expected=I
 elif tag.startswith('cv_'):seed,fold=map(int,tag[3:].split('_'));expected=I[folds(seed)!=fold]
 else:expected=I[~np.isin(gal,gal[tails[tag]])]
 assert np.array_equal(ii,expected);e=pred(k,coef,ii)-target[ii];cost=.5*np.sum(weights(ii)*e*e);assert abs(cost-f['cost'])<1e-11;maxcost=max(maxcost,abs(cost-f['cost']));attempts[k,tag].append((no,f['cost']))
for idx,(k,tag,ii,coef) in trainrecords.items():
 f=J['fit_log'][idx];assert f['model']==k and f['split']==tag and np.array_equal(ii,f['train_rows']) and np.array_equal(coef,f['coefficients']);assert f['cost']==min(z[1] for z in attempts[k,tag])
quaderror=0.
for typ in ['NFW','Burkert']:
 for z in np.geomspace(1e-6,1e4,41):
  integ=(lambda s:s/(1+s)**2) if typ=='NFW' else (lambda s:s*s/((1+s)*(1+s*s)))
  q=quad(integ,0,z,epsabs=1e-24,epsrel=2e-11,limit=300)[0];f=massfn(typ,[z])[0];quaderror=max(quaderror,abs(f/q-1));assert abs(f/q-1)<2e-10
bootcount=0
for tag,comp in J['bootstrap'].items():
 ii,base=allpred['RAR',tag];gs=sorted(set(gal[ii]));rng=np.random.default_rng(20261029);mult=rng.poisson(1,size=(5000,len(gs)));den=mult.sum(1);err0=np.array([np.mean((base[gal[ii]==g]-target[ii][gal[ii]==g])**2) for g in gs]);s0=np.sqrt(mult@err0/den)
 for label,v in comp.items():
  k=label.split(' minus ')[0];_,pr=allpred[k,tag];err=np.array([np.mean((pr[gal[ii]==g]-target[ii][gal[ii]==g])**2) for g in gs]);s=np.sqrt(mult@err/den);ci=np.percentile(s-s0,[2.5,97.5]);assert np.allclose(ci,v['percentile_95_dex'],rtol=1e-9,atol=1e-12);assert np.isclose(metrics(pr,ii)['galaxy_RMS_dex']-metrics(base,ii)['galaxy_RMS_dex'],v['difference_dex'],atol=1e-12);bootcount+=1
report={'status':'passed','tables':tablecount,'prediction_rows':predictioncount,'fit_sets':len(attempts),'optimizer_attempts':len(J['fit_log']),'bootstrap_comparisons':bootcount,'raw_photometric_galaxies_matched':len(ff),'max_prediction_discrepancy_dex':maxpred,'max_objective_cost_discrepancy':maxcost,'density_integral_checks':82,'max_density_integral_relative_discrepancy':quaderror,'all_group_and_feature_tail_exclusions_checked':True,'no_Vflat_prediction_feature':True,'result_sha256':sha(D/'results/results.json'),'audit_script_sha256':sha(Path(__file__)),'limits':['Whole-galaxy transfer is exploratory reuse, not an independent survey.','Mass proxy uses common stellarM/L and stellar light radius; no full baryonic geometry or nuisance likelihood.','Affine halo-state laws do not reproduce full abundance-matching/CDM priors or unify every scale.']}
assert tablecount==56 and len(attempts)==120 and len(J['fit_log'])==320
(O/'shared_halo_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
