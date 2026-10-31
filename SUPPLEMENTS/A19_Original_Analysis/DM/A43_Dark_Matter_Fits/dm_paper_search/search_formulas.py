"""Small exploratory search of shared DM phenomenology; no new observations.

Run with numpy, scipy, matplotlib installed. Input defaults to ../dm_fit/results.
The seven candidate forms and bounds below are fixed before this script is run.
All samples were previously inspected, so nested CV is internal validation only.
"""
from pathlib import Path
import csv, json, hashlib, argparse
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
P=argparse.ArgumentParser();P.add_argument('--input',type=Path,default=ROOT.parent/'dm_fit/results')
args=P.parse_args();OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader((args.input/'held_out_predictions.csv').open()))
cl=json.loads((args.input/'cluster_test.json').read_text())['clusters']
groups=np.array([r['galaxy'] for r in rows]+[r['cluster'] for r in cl])
x=np.log10(np.array([float(r['gbar_m_s2']) for r in rows]+[r['gbar'] for r in cl])/1e-10)
y=np.log10(np.array([float(r['gobs_m_s2']) for r in rows]+[r['gobs'] for r in cl])/1e-10)
z=np.log10(np.array([float(r['radius_kpc']) for r in rows]+[r['r500_kpc'] for r in cl]))
cluster=np.arange(len(x))>=len(rows)
fold=np.array([int(r['fold']) for r in rows]+[i%5 for i in np.random.default_rng(20261027).permutation(len(cl))])
allidx=np.arange(len(y))
KINDS=['RAR','quadratic','quadratic_radius','RAR_radius','RAR_exponent','RAR_radius_exponent','RAR_two_populations']
DESCRIPTIONS={
 'RAR':'y=x-log10(1-exp(-10**((x-a)/2))); a=log10(a0/1e-10)',
 'quadratic':'y=c0+c1*x+c2*x*x',
 'quadratic_radius':'y=c0+c1*x+c2*z+c3*x*x+c4*x*z+c5*z*z',
 'RAR_radius':'y=x-log10(1-exp(-10**((x-a-b*z)/2))); a0(r)=1e-10*10**a*(r/kpc)**b',
 'RAR_exponent':'y=x-log10(1-exp(-10**(p*(x-a)))); p is a free low-acceleration exponent',
 'RAR_radius_exponent':'y=x-log10(1-exp(-10**(p*(x-a-b*z))))',
 'RAR_two_populations':'Standard RAR with a=agal for galaxies and a=acluster for clusters; explicitly piecemeal, not a shared function of continuous observables'}

def design(kind,idx):
    return np.column_stack([np.ones(len(idx)),x[idx],x[idx]**2] if kind=='quadratic' else [np.ones(len(idx)),x[idx],z[idx],x[idx]**2,x[idx]*z[idx],z[idx]**2])

def predict(kind,c,idx):
    if kind.startswith('quadratic'):return design(kind,idx)@c
    a=np.full(len(idx),c[0]); p=.5
    if kind in ['RAR_radius','RAR_radius_exponent']:a+=c[1]*z[idx]
    if kind=='RAR_exponent':p=c[1]
    if kind=='RAR_radius_exponent':p=c[2]
    if kind=='RAR_two_populations':a=np.where(cluster[idx],c[1],c[0])
    return x[idx]-np.log10(-np.expm1(-10.**(p*(x[idx]-a))))

def weights(idx,mode):
    if mode=='equal_point':return np.ones(len(idx))
    w=np.array([1/np.sum(groups[idx]==g) for g in groups[idx]])
    if mode=='equal_scale':
        for iscl in [False,True]:
            m=cluster[idx]==iscl
            if m.any():w[m]/=np.sum(w[m])
    return w

FITCOUNT=0
def fit(kind,idx,mode):
    global FITCOUNT
    FITCOUNT+=1
    w=weights(idx,mode);sw=np.sqrt(w)
    if kind.startswith('quadratic'):
        return np.linalg.lstsq(design(kind,idx)*sw[:,None],y[idx]*sw,rcond=None)[0]
    spec={
      'RAR':([0.1],[-3],[2]),
      'RAR_radius':([0,.1],[-3,-2],[2,2]),
      'RAR_exponent':([0,.5],[-3,.1],[2,.9]),
      'RAR_radius_exponent':([0,.1,.5],[-3,-2,.1],[2,2,.9]),
      'RAR_two_populations':([0,.5],[-3,-3],[2,2])}
    init,lo,hi=spec[kind]
    fun=lambda c:sw*(predict(kind,c,idx)-y[idx])
    r=least_squares(fun,init,bounds=(lo,hi),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=1500)
    assert r.success,(kind,r.message)
    return r.x

def mse_by_system(e,idx):
    return {g:float(np.mean(e[groups[idx]==g]**2)) for g in np.unique(groups[idx])}

def score(pred,idx):
    # Equal populations, and equal systems within each population.
    scores=[]
    for iscl in [False,True]:
        m=cluster[idx]==iscl
        if m.any():scores.append(np.mean(list(mse_by_system((pred-y[idx])[m],idx[m]).values())))
    return float(np.mean(scores))

def metrics(pred,idx=allidx):
    e=pred-y[idx];g=~cluster[idx];c=cluster[idx]
    d={'balanced_system_rmse_dex':float(np.sqrt(score(pred,idx)))}
    if g.any():
        d.update(galaxy_point_rmse_dex=float(np.sqrt(np.mean(e[g]**2))),
          galaxy_equal_system_rmse_dex=float(np.sqrt(np.mean(list(mse_by_system(e[g],idx[g]).values())))),
          galaxy_median_abs_speed_error=float(np.median(abs(10**(e[g]/2)-1))),
          galaxy_median_of_system_median_abs_speed_error=float(np.median([np.median(abs(10**(e[groups[idx]==s]/2)-1)) for s in np.unique(groups[idx][g])])),
          galaxy_fraction_below_baryons=float(np.mean(pred[g]<x[idx][g])))
    if c.any():
        d.update(cluster_rmse_dex=float(np.sqrt(np.mean(e[c]**2))),
          cluster_median_abs_acceleration_error=float(np.median(abs(10**e[c]-1))),
          cluster_max_abs_acceleration_error=float(np.max(abs(10**e[c]-1))),
          cluster_signed_log_residuals=e[c].tolist())
    return d

def save_predictions(name,pred):
    with (OUT/f'{name}_predictions.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(['system','cluster','fold','x','z','y','prediction_y','log10_residual'])
        wr.writerows(zip(groups,cluster,fold,x,z,y,pred,pred-y))

summary={'status':'Exploratory internal validation. Samples and the radius gap were already examined before this candidate search. No claim of all-DM validation.',
  'critical_radius_caveat':'Cluster R500 is derived jointly with mass: gobs(R500)=(4*pi*G/3)*500*rho_crit(z)*R500. Radius thus contains target information. Cross-validation by system does not remove this definitional dependency. A radius-only null is mandatory; low errors for radius-augmented fits do not establish physical unification.',
  'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.input/'held_out_predictions.csv',args.input/'cluster_test.json']},
  'definitions':DESCRIPTIONS,'parameter_bounds':{'a':[-3,2],'b':[-2,2],'p':[.1,.9]},
  'outer_validation':'Existing five whole-system folds; inner selection uses three stratified system folds per outer training set, seed 49173+outer_fold.',
  'selection_score':'Mean of galaxy equal-system log-MSE and cluster log-MSE. No target-system parameters are refit.',
  'weighting':'equal_system gives each system weight 1; equal_scale normalizes these weights separately within galaxies and clusters, giving each population equal total weight.',
  'counts':{'galaxies':149,'galaxy_points':2700,'clusters':7,'cluster_points':7},
  'ranges':{'galaxy_radius_kpc':[float(10**z[~cluster].min()),float(10**z[~cluster].max())],
   'cluster_radius_kpc':[float(10**z[cluster].min()),float(10**z[cluster].max())],
   'galaxy_gbar_m_s2':[float(1e-10*10**x[~cluster].min()),float(1e-10*10**x[~cluster].max())],
   'cluster_gbar_m_s2':[float(1e-10*10**x[cluster].min()),float(1e-10*10**x[cluster].max())]},'weightings':{}}

for mode in ['equal_system','equal_scale']:
    models={};predictions={};foldcoeff={}
    for kind in KINDS:
        pred=np.full(len(y),np.nan);coefs=[]
        for k in range(5):
            tr=allidx[fold!=k];te=allidx[fold==k]
            assert not set(groups[tr])&set(groups[te])
            c=fit(kind,tr,mode);pred[te]=predict(kind,c,te);coefs.append(c.tolist())
        assert np.isfinite(pred).all()
        cf=fit(kind,allidx,mode)
        models[kind]={'coefficients':cf.tolist(),'cv':metrics(pred),'outer_fold_coefficients':coefs}
        predictions[kind]=pred;foldcoeff[kind]=coefs
        save_predictions(f'{mode}_{kind}',pred)
    nested=np.full(len(y),np.nan);selection=[]
    for k in range(5):
        train=allidx[fold!=k];test=allidx[fold==k]
        rng=np.random.default_rng(49173+k);mapping={}
        for iscl in [False,True]:
            ug=np.unique(groups[train][cluster[train]==iscl])
            for j,g in enumerate(rng.permutation(ug)):mapping[g]=j%3
        infold=np.array([mapping[g] for g in groups[train]])
        scores={}
        for kind in KINDS:
            ip=np.full(len(train),np.nan)
            for j in range(3):
                tr=train[infold!=j];te=train[infold==j]
                assert not set(groups[tr])&set(groups[te])
                c=fit(kind,tr,mode);ip[infold==j]=predict(kind,c,te)
            scores[kind]=score(ip,train)
        chosen=min(KINDS,key=lambda q:scores[q])
        nested[test]=predictions[chosen][test]
        selection.append({'outer_fold':k,'selected':chosen,'inner_balanced_mse_dex2':scores})
    save_predictions(f'{mode}_nested_selection',nested)
    shared_nested=np.full(len(y),np.nan);shared_selections=[]
    for s in selection:
        k=s['outer_fold'];scores=s['inner_balanced_mse_dex2']
        shared=min([q for q in KINDS if q!='RAR_two_populations'],key=lambda q:scores[q])
        shared_nested[fold==k]=predictions[shared][fold==k]
        shared_selections.append({'outer_fold':k,'selected':shared})
    save_predictions(f'{mode}_shared_only_nested_selection',shared_nested)
    summary['weightings'][mode]={'models':models,'nested_selection':{'cv':metrics(nested),'selections':selection},'shared_only_nested_selection':{'cv':metrics(shared_nested),'selections':shared_selections}}
    print(mode,[(k,round(v['cv']['galaxy_median_abs_speed_error']*100,3),round(v['cv']['cluster_median_abs_acceleration_error']*100,3),round(v['cv']['balanced_system_rmse_dex'],4)) for k,v in models.items()],flush=True)
    print('nested',metrics(nested),[s['selected'] for s in selection],flush=True)

# Galaxy-only experiment: same folds and point weights as the earlier published result.
# Equal-system weighting is supplied separately to expose the training-weight change.
gi=allidx[~cluster];galonly={}
for kind in ['RAR','quadratic','quadratic_radius','RAR_radius','RAR_exponent','RAR_radius_exponent']:
    pred=np.zeros(len(gi));coeff=[]
    for k in range(5):
        tr=gi[fold[gi]!=k];te=gi[fold[gi]==k]
        cf=fit(kind,tr,'equal_system');pred[fold[gi]==k]=predict(kind,cf,te);coeff.append(cf.tolist())
    cf=fit(kind,gi,'equal_system')
    cp=predict(kind,cf,allidx[cluster])
    galonly[kind]={'coefficients':cf.tolist(),'galaxy_cv':metrics(pred,gi),'zero_cluster_training_transfer':metrics(cp,allidx[cluster])}
summary['galaxy_only_equal_system_training']=galonly

# Repeat galaxy-only candidates using exactly the original point-weighted objective.
galpoint={}
for kind in ['RAR','quadratic','quadratic_radius','RAR_radius','RAR_exponent','RAR_radius_exponent']:
    pred=np.zeros(len(gi))
    for k in range(5):
        tr=gi[fold[gi]!=k];te=gi[fold[gi]==k]
        cf=fit(kind,tr,'equal_point');pred[fold[gi]==k]=predict(kind,cf,te)
    cf=fit(kind,gi,'equal_point')
    galpoint[kind]={'coefficients':cf.tolist(),'galaxy_cv':metrics(pred,gi),'zero_cluster_training_transfer':metrics(predict(kind,cf,allidx[cluster]),allidx[cluster])}
    with (OUT/f'galaxy_only_equal_point_{kind}_predictions.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(['system','fold','x','z','y','prediction_y','log10_residual'])
        wr.writerows(zip(groups[gi],fold[gi],x[gi],z[gi],y[gi],pred,pred-y[gi]))
summary['galaxy_only_equal_point_training']=galpoint

# Diagnostic nulls on the seven R500 estimates. The first ignores baryons entirely;
# the second exploits the overdensity-radius definition and is not a new gravity law.
ci=allidx[cluster];nulls={}
for name,slope in [('constant_acceleration',0),('radius_definition_g_equals_C_r',1)]:
    pred=np.full(len(ci),np.nan);coefs=[]
    for k in range(5):
        tr=ci[fold[ci]!=k];te=ci[fold[ci]==k]
        a=float(np.mean(y[tr]-slope*z[tr]));pred[fold[ci]==k]=a+slope*z[te];coefs.append(a)
    a=float(np.mean(y[ci]-slope*z[ci]));nulls[name]={'intercept_y_minus_slope_z':a,'slope':slope,'outer_fold_intercepts':coefs,'cv':metrics(pred,ci)}
    with (OUT/f'cluster_null_{name}_predictions.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(['system','fold','x','z','y','prediction_y','log10_residual'])
        wr.writerows(zip(groups[ci],fold[ci],x[ci],z[ci],y[ci],pred,pred-y[ci]))
summary['cluster_nulls']=nulls

# Paired system bootstrap of the fixed candidate comparison, not a post-selection p-value.
# It resamples already cross-validated system residuals; training uncertainty is omitted.
base=list(csv.DictReader((OUT/'equal_system_RAR_predictions.csv').open()))
alt=list(csv.DictReader((OUT/'equal_system_RAR_radius_predictions.csv').open()))
eb=np.array([float(r['log10_residual']) for r in base]);ea=np.array([float(r['log10_residual']) for r in alt])
ug=np.unique(groups[~cluster]);ub=np.array([np.mean(eb[groups==g]**2) for g in ug]);ua=np.array([np.mean(ea[groups==g]**2) for g in ug])
rng=np.random.default_rng(123456);draw=rng.integers(0,len(ug),(10000,len(ug)))
delta=np.sqrt(np.mean(ua[draw],axis=1))-np.sqrt(np.mean(ub[draw],axis=1))
summary['paired_system_bootstrap_RAR_radius_minus_RAR']={'quantity':'Galaxy equal-system CV RMS difference in dex, joint equal-system training; conditional on existing CV predictions, excludes retraining uncertainty and model selection',
 'point':float(np.sqrt(ua.mean())-np.sqrt(ub.mean())),'percentile_95_interval':np.quantile(delta,[.025,.975]).tolist()}
summary['fit_calls']=FITCOUNT
(OUT/'search_results.json').write_text(json.dumps(summary,indent=2))
print('saved',OUT/'search_results.json','fits',FITCOUNT,flush=True)
