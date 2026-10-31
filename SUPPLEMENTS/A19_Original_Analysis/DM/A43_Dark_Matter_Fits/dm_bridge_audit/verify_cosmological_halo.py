"""Independent deterministic concentration-mass halo audit."""
from pathlib import Path
import csv,json,hashlib,collections
import numpy as np
from numpy.polynomial.legendre import leggauss
B=Path(__file__).resolve().parent.parent;D=B/'dm_cosmological_halo';O=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open()))
J=json.loads((D/'results/results.json').read_text());rows=read(D/'input.csv');features=read(D/'galaxy_features.csv')
for key,name in [('protocol','PROTOCOL.txt'),('script','run.py'),('input','input.csv'),('features','galaxy_features.csv')]:assert J[key+'_sha256']==sha(D/name)
for name in ['input.csv','galaxy_features.csv']:assert sha(D/name)==sha(B/'dm_shared_halo'/name)
# Raw photometry provenance is independently checked by verify_shared_halo.py;
# these are exact immutable copies, not a second join using kinematic targets.
N=len(rows);I=np.arange(N);gal=np.array([r['system'] for r in rows]);group=np.array([r['group'] for r in rows]);rad=np.array([float(r['radius_kpc']) for r in rows]);gb=np.array([float(r['gbar_m_s2']) for r in rows]);go=np.array([float(r['gobs_m_s2']) for r in rows]);mp=np.array([float(r['mass_proxy_msun']) for r in rows]);re=np.array([float(r['reff_kpc']) for r in rows]);target=np.log10(go/1e-10);ff={r['system']:r for r in features}
G=6.67430e-11;MSUN=1.988409870698051e30;KPC=3.085677581491367e19;h=.671;H0=67.1*1000/(1e3*KPC);fb=.049/.3175
for key,val in {'G_SI':G,'Msun_kg':MSUN,'kpc_m':KPC,'h':h,'H0_SI':H0,'fb':fb,'cM_intercept':.905,'cM_slope':-.101}.items():assert J['constants'][key]==val
qn,qw=leggauss(40);tt=(qn+1)/2;qw=qw/2
def f(z):
 z=np.asarray(z);out=np.empty_like(z);small=z<.1
 # Integrate dimensionless NFW density rather than using analysis series.
 out[small]=z[small]**2*np.sum(qw*tt/(1+z[small,None]*tt)**2,axis=1)
 out[~small]=np.log1p(z[~small])-1+1/(1+z[~small]);return out
def physical(k,c,ii):
 M=mp[ii]/fb+10**c[0]*(mp[ii]/1e10)**c[1]*(re[ii]/3)**(c[2] if k.endswith('size') else 0)
 conc=10**.905*(M*h/1e12)**(-.101)
 # Derive r200 from mean-density definition, independent algebraic route.
 rho_critical=3*H0*H0/(8*np.pi*G);r200=(3*M*MSUN/(4*np.pi*200*rho_critical))**(1/3)/KPC;rs=r200/conc
 enclosed=(1-fb)*M*f(np.minimum(rad[ii]/rs,conc))/f(conc)
 acc=G*MSUN*enclosed/(rad[ii]*KPC)**2
 return M,conc,r200,rs,acc
def pred(k,c,ii):return np.log10((gb[ii]+physical(k,c,ii)[-1])/1e-10)
def weights(ii):
 groups,counts=np.unique(group[ii],return_counts=True);cc=dict(zip(groups,counts));return np.array([1/(len(groups)*cc[g]) for g in group[ii]])
def metrics(p,ii):
 e=p-target[ii];return {'galaxy_RMS_dex':np.sqrt(np.dot(weights(ii),e*e)),'row_RMS_dex':np.sqrt(np.mean(e*e)),'median_acceleration_error_percent':100*np.median(abs(10**e-1)),'median_speed_error_percent':100*np.median(abs(10**(.5*e)-1)),'rows':len(ii),'galaxies':len(set(gal[ii]))}
def checkmetrics(p,ii,d):
 for key,value in metrics(p,ii).items():assert np.isclose(value,d[key],rtol=1e-10,atol=1e-12),(key,value,d[key])
def checkdiag(k,c,ii,d):
 M,conc,r200,rs,acc=physical(k,c,ii);outside=(M*h<1e10)|(M*h>1e15);over=rad[ii]>r200
 for key,value in [('Mref_msun_minmax',M),('c200_minmax',conc),('r200_kpc_minmax',r200),('rs_kpc_minmax',rs),('proxy_retention_fraction_minmax',mp[ii]/(fb*M))]:assert np.allclose(d[key],[min(value),max(value)],rtol=1e-12,atol=0)
 assert d['reference_cM_calibration_extrapolated_rows']==sum(outside) and d['reference_cM_calibration_extrapolated_galaxies']==sorted(set(gal[ii][outside]))
 assert d['rows_outside_reference_r200']==sum(over) and d['galaxies_outside_reference_r200']==sorted(set(gal[ii][over]))
 assert abs(d['maximum_observed_radius_over_reference_r200']-max(rad[ii]/r200))<1e-12 and np.isclose(d['min_halo_acceleration_m_s2'],min(acc),rtol=1e-12,atol=0)
 assert d['formal_mass_zero_halo_acceleration_vanishes']==bool(c[1]>0) and d['finite_total_dark_mass'] and d['deterministic_cosmological_structural_prior'] and not d['full_LCDM_prediction']
 assert np.all(mp[ii]/(fb*M)<=1) and np.all(acc>0)
def folds(seed):
 order=np.random.default_rng(seed).permutation(sorted(set(group)));return np.array([next(j for j in range(5) if g in order[j::5]) for g in group])
tails={}
for name,field in [('mass','mass_proxy_msun'),('surface','sigma_proxy_msun_kpc2')]:
 gs=sorted(ff,key=lambda g:(float(ff[g][field]),g))
 for tag,gg in [('low',gs[:37]),('high',gs[-37:])]:tails[name+'_'+tag]=I[np.isin(gal,gg)];assert J['tails'][name+'_'+tag]['rows']==tails[name+'_'+tag].tolist()
fitgroups=collections.defaultdict(list);maxcost=0.
for i,d in enumerate(J['fit_log']):
 ii=np.array(d['train_rows']);c=np.array(d['coefficients']);k=d['model'];spec=J['spec'][k];assert d['success'] and np.all(c>=np.array(spec['lo'])-1e-10) and np.all(c<=np.array(spec['hi'])+1e-10)
 value=np.dot(weights(ii),(pred(k,c,ii)-target[ii])**2)/2;maxcost=max(maxcost,abs(value-d['cost']));assert abs(value-d['cost'])<1e-11;fitgroups[k,d['split']].append(i)
def checkfit(fid,k,tag,ii,c):
 d=J['fit_log'][fid];assert d['model']==k and d['split']==tag and d['train_rows']==ii.tolist() and d['coefficients']==c
 assert d['cost']==min(J['fit_log'][i]['cost'] for i in fitgroups[k,tag]) and len(fitgroups[k,tag])==3
allpred={};count=0;maxpred=0.;maxphysical=0.;tables=0
for k,d in J['models'].items():
 variants=[('all',I,d['coefficients'],d,d['selected_all_fit'])]
 variants += [('cv_'+s,I,v['coefficients'],v,v['selected_fit_ids']) for s,v in d['CV'].items()]
 variants += [(tag,tails[tag],v['coefficients'],v,v['selected_fit_id']) for tag,v in d['tail_holdout'].items()]
 for tag,ii,c,score,fid in variants:
  tab=read(D/'results'/f'{k}_{tag}.csv');assert [int(r['row']) for r in tab]==ii.tolist();p=np.log10([float(r['gpred_m_s2'])/1e-10 for r in tab]);expected=np.empty(len(ii));phys=np.empty((5,len(ii)))
  if tag.startswith('cv_'):
   fd=folds(int(tag[3:]));assert np.array_equal([int(r['fold']) for r in tab],fd)
   for fold in range(5):
    tr=I[fd!=fold];te=I[fd==fold];assert not set(gal[tr])&set(gal[te]);checkfit(fid[fold],k,tag+'_'+str(fold),tr,c[fold]);expected[fd==fold]=pred(k,c[fold],te);phys[:,fd==fold]=physical(k,c[fold],te);checkdiag(k,c[fold],te,score['test_physical_diagnostics'][str(fold)])
  else:
   tr=I if tag=='all' else I[~np.isin(gal,gal[ii])];checkfit(fid,k,tag,tr,c);expected=pred(k,c,ii);phys[:]=physical(k,c,ii)
   if tag=='all':checkdiag(k,c,ii,score['physical_diagnostics']);score=score['all']
   else:assert not set(gal[tr])&set(gal[ii]);checkdiag(k,c,ii,score['test_physical_diagnostics'])
  err=max(abs(p-expected));maxpred=max(maxpred,float(err));assert err<1e-10;checkmetrics(p,ii,score['metrics'])
  for g in set(gal[ii]):z=gal[ii]==g;checkmetrics(p[z],ii[z],score['by_galaxy'][g])
  for j,r in enumerate(tab):
   i=ii[j];assert r['system']==gal[i] and r['group']==group[i] and int(r['broad_row'])==int(rows[i]['broad_row']);assert abs(float(r['log10_residual'])-(p[j]-target[i]))<1e-12
   for pidx,name in enumerate(['Mref_msun','c200','reference_r200_kpc','rs_kpc','halo_g_m_s2']):
    err=abs(float(r[name])/phys[pidx,j]-1);maxphysical=max(maxphysical,err);assert err<1e-10
   assert (r['reference_halo_truncated_here']=='True')==(rad[i]>phys[2,j]);assert (r['cM_mass_extrapolation']=='True')==bool(phys[0,j]*h<1e10 or phys[0,j]*h>1e15)
  allpred[k,tag]=(ii,p);count+=len(ii);tables+=1
copies=0;bootcount=0
for tag,d in J['comparator_sources'].items():
 ii=I if tag.startswith('cv_') else tails[tag]
 for k,s in d.items():
  pth=D/s['copied_file'];assert sha(pth)==s['sha256']==sha(B/'dm_shared_halo/results'/pth.name);tab=read(pth);assert [int(r['row']) for r in tab]==ii.tolist();p=np.log10([float(r['gpred_m_s2'])/1e-10 for r in tab]);allpred[k,tag]=(ii,p);checkmetrics(p,ii,J['comparison_metrics'][tag][k]);copies+=1
 for k in J['models']:checkmetrics(allpred[k,tag][1],ii,J['comparison_metrics'][tag][k])
for tag,d in J['bootstrap'].items():
 ii=I if tag.startswith('cv_') else tails[tag];gs=sorted(set(gal[ii]));mult=np.random.default_rng(20261029).poisson(1,size=(5000,len(gs)));den=mult.sum(1);scores={}
 for label,s in d.items():
  left,right=label.split(' minus ')
  for k in [left,right]:
   if k not in scores:
    p=allpred[k,tag][1];err=np.array([np.mean((p[gal[ii]==g]-target[ii][gal[ii]==g])**2) for g in gs]);scores[k]=np.sqrt(mult@err/den)
  ci=np.percentile(scores[left]-scores[right],[2.5,97.5]);assert np.allclose(ci,s['percentile_95_dex'],rtol=1e-10,atol=1e-12);dif=metrics(allpred[left,tag][1],ii)['galaxy_RMS_dex']-metrics(allpred[right,tag][1],ii)['galaxy_RMS_dex'];assert abs(dif-s['difference_dex'])<1e-12;bootcount+=1
assert tables==J['new_prediction_tables']==16 and len(fitgroups)==J['fit_sets']==40 and len(J['fit_log'])==J['optimizer_attempts']==120 and J['optimizer_failures']==[] and copies==49
# Explicit tests beyond observed range validate the claimed finite-mass branch.
far=np.geomspace(1e-3,1e4,101);c=12.;massfraction=f(np.minimum(far,c))/f(np.array([c]))[0];assert np.all(massfraction<=1+1e-13) and np.allclose(massfraction[far>=c],1)
report={'status':'PASS','independent_method':'Immutable photometric-input identity; independent density integral and critical-density radius construction; all forecasts, structural states, diagnostics, fits, minima, exclusions, comparator copies and bootstraps. No analysis-code imports.','prediction_tables':tables,'prediction_rows':count,'fit_sets':len(fitgroups),'optimizer_attempts':len(J['fit_log']),'comparator_tables_verified':copies,'bootstrap_comparisons':bootcount,'max_prediction_discrepancy_dex':maxpred,'max_structural_state_relative_discrepancy':maxphysical,'max_objective_discrepancy':maxcost,'truncation_check_radii':len(far),'no_Vflat_prediction_feature':True,'limits':['Mref is dark-matter-only-equivalent reference M200c; dark component is (1-fb)Mref. This is an explicit deterministic convention, not an empirical halo mass measurement.','Mass mapping and concentration are conditional on common photometric stellar M/L and an approximate HI-plus-stellar mass proxy.','Published concentration-mass calibration is extrapolated for some low-mass objects; scatter and baryonic halo response are omitted.','Fixed Planck cosmology and concentration-mass relation supply structural priors, not a full LambdaCDM prediction.','No observed radius tests the imposed outer truncation; finite total mass is a declared extrapolation.','The three-coefficient size-model RMS advantage is small, changes sign across CV seeds and has descriptive intervals spanning zero; median speed error is worse than RAR.','No per-galaxy target-informed halo parameter is fitted in held-out forecasts; affine mapping coefficients are shared.'],'result_sha256':sha(D/'results/results.json'),'audit_script_sha256':sha(Path(__file__))}
(O/'cosmological_halo_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
