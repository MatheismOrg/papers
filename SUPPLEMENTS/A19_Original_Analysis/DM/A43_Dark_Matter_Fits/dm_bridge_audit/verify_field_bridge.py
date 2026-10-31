"""Independent positive-kernel/Plummer audit; never imports the fit code."""
from pathlib import Path
import csv,json,hashlib,collections
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
B=Path(__file__).resolve().parent.parent;D=B/'dm_field_bridge';O=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open()))
J=json.loads((D/'results.json').read_text());rows=read(D/'input.csv');basis=read(D/'basis.csv')
for key,name in [('protocol','PROTOCOL.txt'),('script','run.py'),('input','input.csv'),('basis','basis.csv')]:assert J[key+'_sha256']==sha(D/name)
assert J['source_sha256']==sha(B/'dm_unified_regular/augmented_input.csv')
broad=read(B/'dm_unified_regular/augmented_input.csv');assert len(rows)==44
for row in rows:
 orig=broad[int(row['row'])];assert orig['population']=='dSph'
 for key in orig:assert row[key]==orig[key],(row['system'],key)
R=1/np.sqrt(2**(2/3)-1)
assert abs(R-J['basis_halfmass_radius_over_a'])<1e-14
def F(d):
 # Exact translated-ball Plummer mass, rationalized to avoid cancellation.
 A=np.sqrt(1+(R-d)**2);BB=np.sqrt(1+(R+d)**2);C=d*d+1-R*R
 return 4*R**3/(A*BB*(A+BB)*(A*BB+C))
def F_shell(d):
 # Independent source-centered spherical shells and solid-angle intersection.
 if d==0:return R**3/(1+R*R)**1.5
 s0=max(0.,R-d);full=s0**3/(1+s0*s0)**1.5
 lo=abs(d-R);hi=d+R
 part=quad(lambda s:3*s*s/(1+s*s)**2.5*.5*(1-(s*s+d*d-R*R)/(2*s*d)),lo,hi,epsabs=1e-20,epsrel=2e-12)[0]
 return full+part
shift_errors=[]
for d in np.r_[0.,np.geomspace(1e-6,1e3,51)]:
 exact=F(d);direct=F_shell(d);shift_errors.append(abs(exact/direct-1));assert np.isclose(exact,direct,rtol=5e-9,atol=1e-16)
def H(z,kind):
 # Integrate full infinite support using the analytic F, unlike the analysis'
 # finite t interval and 64-point radial quadrature. Scale tiny-z integrals.
 if z<=1:
  if kind=='NLG':fun=lambda d:.5*(1+z*d)*np.exp(-z*d)*F(d);scale=z
  else:fun=lambda d:d*np.exp(-z*d)*F(d);scale=z*z
  pieces=[0.,1.,10.,np.inf]
 else:
  fun=lambda t:(.5*(1+t) if kind=='NLG' else t)*np.exp(-t)*F(t/z);scale=1.;pieces=[0.,1.,10.,np.inf]
 return scale*sum(quad(fun,a,b,epsabs=1e-13,epsrel=3e-11,limit=200)[0] for a,b in zip(pieces[:-1],pieces[1:]))
norms={k:quad(lambda t:(.5*(1+t) if k=='NLG' else t)*np.exp(-t),0,np.inf,epsabs=1e-12)[0] for k in ['NLG','MOG']}
assert all(abs(v-1)<1e-12 for v in norms.values())
logs=np.array([float(z['log10z']) for z in basis]);inter={k:PchipInterpolator(logs,np.log10([float(z[k+'_H']) for z in basis])) for k in ['NLG','MOG']}
basis_errors={};mid_errors={};proofcheck=0
for kind in inter:
 direct=np.array([H(10**x,kind) for x in logs]);saved=np.array([float(z[kind+'_H']) for z in basis]);err=np.max(abs(saved/direct-1));basis_errors[kind]=float(err)
 # Frozen grid's adaptive tiny-z tail is accurate to about 1.4 ppm for MOG.
 assert err<2e-6,(kind,err)
 assert np.all(direct>0) and np.all(direct<=.5*(1+1e-10)) and np.all(np.diff(direct)>=-1e-12)
 mid=(logs[1:]+logs[:-1])/2;truth=np.array([H(10**x,kind) for x in mid]);err=np.max(abs(10**inter[kind](mid)/truth-1));mid_errors[kind]=float(err)
 assert err<5e-5,(kind,err)
 proofcheck+=len(direct)+len(truth)
N=len(rows);I=np.arange(N);gb=np.array([float(r['gbar_m_s2']) for r in rows]);go=np.array([float(r['gobs_m_s2']) for r in rows]);rad=np.array([float(r['radius_kpc']) for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows]);target=np.log10(go)
frozen={'NLG':(10.94,.059),'MOG':(8.89,.042)}
def pred(k,c,ii,ml,direct=False):
 b=gb[ii]*ml/2
 if k.startswith('RAR'):
  a0=1.2e-10 if k=='RAR_fixed' else 1e-10*10**c[0]
  return np.log10(b/(1-np.exp(-np.sqrt(b/a0))))
 kind=k.split('_')[0];alpha,mu=frozen[kind] if k.endswith('_frozen') else 10**np.array(c)
 z=mu*rad[ii]/R;h=np.array([H(v,kind) for v in z]) if direct else 10**inter[kind](np.log10(z))
 return np.log10(b*(1+2*alpha*h))
def metric(p,ii):
 e=p-target[ii];return dict(RMS_dex=float(np.sqrt(np.mean(e*e))),median_acceleration_error_percent=float(100*np.median(abs(10**e-1))),rows=len(ii))
def folds(seed):
 perm=np.random.default_rng(seed).permutation(sorted(set(group)));return np.array([next(f for f in range(5) if g in perm[f::5]) for g in group])
bound=read(D/'mass_budget_bounds.csv');assert len(bound)==132
summ={}
for ml in [1,2,3]:
 DD=go/(gb*ml/2);vv=[v for v in bound if int(v['stellar_ML'])==ml];assert [v['system'] for v in vv]==[v['system'] for v in rows]
 for i,v in enumerate(vv):
  for k,value in [('gobs_over_gbar',DD[i]),('necessary_alpha_any_spherical_profile',max(0.,(DD[i]-1)/2)),('necessary_alpha_Plummer',max(0.,DD[i]-1))]:assert np.isclose(float(v[k]),value,rtol=1e-14)
  for kind,(alpha,mu) in frozen.items():assert (v['exceeds_'+kind+'_general_bound']=='True')==(DD[i]>1+2*alpha)
 ss={k+'_general_violations':int(sum(DD>1+2*v[0])) for k,v in frozen.items()};ss['max_alpha_required']=float(max(np.maximum(0,(DD-1)/2)));assert ss==J['bound_summary'][str(ml)];summ[str(ml)]=dict(ss,**{k+'_Plummer_violations':int(sum(DD>1+v[0])) for k,v in frozen.items()})
attempts=collections.defaultdict(list);maxcost=0.
for fit in J['fits']:
 k=fit['model'];ml=fit['stellar_ML'];ii=np.array(fit['training_rows']);c=np.array(fit['coefficients']);assert fit['success']
 assert np.all(c>=np.array(fit['lo'])-1e-10) and np.all(c<=np.array(fit['hi'])+1e-10)
 value=np.mean((pred(k,c,ii,ml)-target[ii])**2)/2;maxcost=max(maxcost,abs(value-fit['cost']));assert abs(value-fit['cost'])<1e-10
 attempts[ml,k,fit['split']].append(fit)
def checkfit(k,ml,tag,train,coef):
 if not k.endswith('_fit'):assert coef==[];return
 fits=attempts[ml,k,tag];assert len(fits)==(2 if k=='RAR_fit' else 4)
 assert all(f['training_rows']==train.tolist() for f in fits)
 best=min(fits,key=lambda f:f['cost']);assert np.array_equal(coef,best['coefficients'])
seen=0;count=0;maxpred=0.;maxdirect=0.;all_mu_bounds=[]
for mlstr,models in J['models_by_ML'].items():
 ml=int(mlstr)
 for k,data in models.items():
  variants=[('all',I,data['coefficients'],data['all'])]+[('cv_'+s,I,d['coefficients'],d['metrics']) for s,d in data['CV'].items()]+[('host_'+s,I[np.isin(env,['Milky_Way','G/L'] if s=='MW' else ['M31','A/L'])],d['coefficients'],d['metrics']) for s,d in data['host'].items()]
  for tag,ii,coef,score in variants:
   tab=read(D/'results'/f'ML{ml}_{k}_{tag}.csv');assert [int(r['input_row']) for r in tab]==ii.tolist();actual=np.log10([float(r['gpred_m_s2']) for r in tab]);expected=np.empty(len(ii))
   if tag.startswith('cv_'):
    fd=folds(int(tag[3:]));assert np.array_equal(fd,[int(r['fold']) for r in tab])
    for f in range(5):
     train=I[fd!=f];test=I[fd==f];assert not set(group[train])&set(group[test]);checkfit(k,ml,tag+'_'+str(f),train,coef[f]);expected[fd==f]=pred(k,coef[f],test,ml)
   else:
    train=I if tag=='all' else I[~np.isin(group,group[ii])];checkfit(k,ml,tag,train,coef);expected=pred(k,coef,ii,ml)
    if tag.startswith('host_'):
     assert not set(group[train])&set(group[ii]);assert data['host'][tag[5:]]['training_rows']==train.tolist();assert len(ii)==(26 if tag=='host_MW' else 18)
    if tag=='all' and not k.startswith('RAR'):
     exact=pred(k,coef,ii,ml,True);maxdirect=max(maxdirect,float(max(abs(exact-actual))))
     if k.endswith('_fit'):all_mu_bounds.append({'ML':ml,'model':k,'mu_kpc_inverse':10**coef[1],'at_upper_bound':abs(coef[1]-5)<1e-5})
   error=float(max(abs(expected-actual)));maxpred=max(maxpred,error);assert error<1e-10
   for key,val in metric(actual,ii).items():assert np.isclose(val,score[key],rtol=1e-11,atol=1e-12)
   for j,r in enumerate(tab):
    i=ii[j];assert r['system']==rows[i]['system'] and int(r['broad_row'])==int(rows[i]['row']);assert np.isclose(float(r['gbar_m_s2']),gb[i]*ml/2,rtol=1e-15,atol=0)
    assert abs(float(r['log10_residual'])-(actual[j]-target[i]))<1e-12
   seen+=1;count+=len(ii)
assert seen==J['prediction_tables']==108 and len(attempts)==J['fit_sets']==162 and len(J['fits'])==J['optimizer_attempts']==540
assert J['optimizer_failures']==[] and maxdirect<2e-5
selected_kernels=[min(fs,key=lambda f:f['cost']) for (ml,k,tag),fs in attempts.items() if not k.startswith('RAR')]
assert len(selected_kernels)==108 and all(abs(f['coefficients'][1]-5)<1e-4 for f in selected_kernels)
alpha_limit_error=0.
for f in selected_kernels:
 ii=np.array(f['training_rows']);ml=f['stellar_ML'];lim=10**np.mean(np.log10(go[ii]/(gb[ii]*ml/2)))-1
 alpha_limit_error=max(alpha_limit_error,abs(10**f['coefficients'][0]/lim-1))
assert alpha_limit_error<2e-6
report={'status':'PASS','independent_method':'Exact rationalized translated Plummer sphere; independent infinite-support scalar integrals; source joins, all predictions, objectives, minima, metrics and complete galaxy/host exclusions. No analysis-code imports.','prediction_tables':seen,'prediction_rows':count,'fit_sets':len(attempts),'optimizer_attempts':len(J['fits']),'kernel_normalizations':norms,'translated_sphere_checks':len(shift_errors),'max_translated_sphere_relative_discrepancy':max(shift_errors),'independent_basis_integrals':proofcheck,'basis_max_relative_error':basis_errors,'midpoint_interpolation_max_relative_error':mid_errors,'max_prediction_discrepancy_dex':maxpred,'max_all_fit_direct_integral_discrepancy_dex':maxdirect,'max_objective_discrepancy':maxcost,'bound_summary':summ,'all_fitted_mu_at_upper_bound':all_mu_bounds,'selected_kernel_fit_sets_at_mu_upper_bound':len(selected_kernels),'max_selected_alpha_relative_difference_from_analytic_infinite_mu_limit':alpha_limit_error,'limits':['General 1+2alpha bound is necessary only for isolated spherical finite nonnegative baryon source at its baryon half-mass radius.','Sharper 1+alpha bound needs a radially nonincreasing source density such as Plummer.','Dwarf accelerations use a half-light estimator and assumed stellar mass-to-light ratios; binaries, anisotropy, tides and host field are not modeled.','Finite satellite host density and external-field contributions are absent: this does not exclude every formulation of NLG or MOG.','Stationary kernels provide no fitted time history or evidence for discrete spacetime.','mu saturates the fitted bound and the short-range limit is effectively a common multiplicative mass enhancement, not a resolved finite length measurement.'],'result_sha256':sha(D/'results.json'),'audit_script_sha256':sha(Path(__file__))}
(O/'field_bridge_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
