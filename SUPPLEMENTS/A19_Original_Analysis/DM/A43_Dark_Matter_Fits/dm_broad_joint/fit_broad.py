"""Exploratory equal-population shared-formula comparison. See PROTOCOL.txt.

The input combined_input.csv must have population, system, group, radius_kpc,
gbar_m_s2, gobs_m_s2; optional subpopulation, environment and source fields.
Only physically specified positive radii are admitted by the input builder.
"""
from pathlib import Path
import json,csv,argparse
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
P=argparse.ArgumentParser();P.add_argument('--quick',action='store_true')
P.add_argument('--input',type=Path,default=ROOT/'combined_input.csv')
P.add_argument('--metadata',type=Path,default=ROOT/'input_metadata.json')
P.add_argument('--output',type=Path,default=ROOT/'results');args=P.parse_args()
OUT=args.output;OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader(args.input.open()))
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows])
y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows])
u=np.log10([float(r['radius_kpc']) for r in rows])
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows])
subpop=np.array([r.get('subpopulation',r['population']) for r in rows])
env=np.array([r.get('environment','') for r in rows])
I=np.arange(len(rows));POPS=sorted(set(pop));a_fixed=np.log10(1.1603274962269242e-10/1e-10)
SHARED=['RAR','RAR_exponent','RAR_radius','RAR_radius_quad','RAR_potential_quad','RAR_density_quad','quadratic_2D']
KINDS=SHARED+['RAR_piecewise']

def rawweights(idx):
    w=np.empty(len(idx))
    for p in np.unique(pop[idx]):
        mp=pop[idx]==p;gs=np.unique(group[idx][mp])
        for g in gs:
            m=mp&(group[idx]==g);w[m]=1/(len(gs)*np.sum(m))
    return w

def features(kind,idx):
    if kind=='quadratic_2D':return np.column_stack([x[idx]*0+1,x[idx],u[idx],x[idx]**2,x[idx]*u[idx],u[idx]**2])
    q={'RAR_radius':u,'RAR_radius_quad':u,'RAR_potential_quad':x+u,'RAR_density_quad':x-u}[kind][idx]
    return np.column_stack([q*0+1,q] if kind=='RAR_radius' else [q*0+1,q,q*q])

def rar(xv,a,p=.5):
    logt=p*(xv-a)
    # Stable log(1-exp(-t)); clipping merely avoids exponential overflow.
    t=10**np.clip(logt,-290,290)
    return xv-np.log10(-np.expm1(-t))

def predict(kind,c,idx):
    if kind=='quadratic_2D':return features(kind,idx)@c
    if kind=='RAR':return rar(x[idx],c[0])
    if kind=='RAR_exponent':return rar(x[idx],c[0],c[1])
    if kind=='RAR_piecewise':
        return rar(x[idx],np.array([c[POPS.index(p)] for p in pop[idx]]))
    return rar(x[idx],features(kind,idx)@c)

COUNTS=0;FIT_WARNINGS=[]
def fit(kind,idx):
    global COUNTS
    COUNTS+=1;sw=np.sqrt(rawweights(idx))
    if kind=='quadratic_2D':
        A=features(kind,idx);return np.linalg.lstsq(A*sw[:,None],y[idx]*sw,rcond=None)[0]
    if kind=='RAR':c0=[.1];lo=[-3];hi=[3]
    elif kind=='RAR_exponent':c0=[.1,.5];lo=[-3,.1];hi=[3,.9]
    elif kind=='RAR_radius':c0=[.1,0];lo=[-3,-2];hi=[3,2]
    elif kind=='RAR_piecewise':c0=[.1]*len(POPS);lo=[-3]*len(POPS);hi=[3]*len(POPS)
    else:c0=[.1,0,0];lo=[-3,-2,-1];hi=[3,2,1]
    r=least_squares(lambda c:sw*(predict(kind,c,idx)-y[idx]),c0,bounds=(lo,hi),max_nfev=800,xtol=1e-9,ftol=1e-9,gtol=1e-9)
    if not r.success:FIT_WARNINGS.append(dict(kind=kind,n=len(idx),message=r.message))
    return r.x

def score(pred,idx):
    e=pred-y[idx]
    return float(np.sum(rawweights(idx)*e*e)/len(np.unique(pop[idx])))

def met(pred,idx):
    e=pred-y[idx]
    return {'rows':len(idx),'systems':len({rows[i]['system'] for i in idx}),'validation_groups':len(np.unique(group[idx])),
      'rms_log10_acceleration_dex':float(np.sqrt(np.mean(e*e))),
      'equal_system_rms_dex':float(np.sqrt(np.mean([np.mean(e[group[idx]==g]**2) for g in np.unique(group[idx])]))),
      'median_absolute_acceleration_fraction':float(np.median(abs(10**np.clip(e,-300,300)-1))),
      'median_observed_over_predicted':float(np.median(10**np.clip(-e,-300,300))),
      'median_absolute_speed_fraction_under_same_radius':float(np.median(abs(10**np.clip(e/2,-300,300)-1))),
      'median_signed_log10_residual':float(np.median(e)),
      'fraction_predicted_below_baryons':float(np.mean(pred<x[idx]))}

def metrics(pred,idx=I):
    out={'balanced_population_system_rms_dex':float(np.sqrt(score(pred,idx))),'populations':{}}
    for p in np.unique(pop[idx]):
        m=pop[idx]==p;out['populations'][p]=met(pred[m],idx[m])
    out['dwarf_subpopulations']={}
    for s in ['classical_dSph','ultrafaint_dSph']:
        m=subpop[idx]==s
        if m.any():out['dwarf_subpopulations'][s]=met(pred[m],idx[m])
    return out

def makefold(idx,n,seed):
    rng=np.random.default_rng(seed);mapping={}
    # Shared physical systems receive the same fold even across different populations.
    for p in sorted(np.unique(pop[idx])):
        gs=[g for g in np.unique(group[idx][pop[idx]==p]) if g not in mapping]
        counts=np.array([sum(v==j for g,v in mapping.items() if g in group[idx][pop[idx]==p]) for j in range(n)])
        for g in rng.permutation(gs):
            j=int(np.argmin(counts));mapping[g]=j;counts[j]+=1
    return np.array([mapping[g] for g in group[idx]])

fold=makefold(I,5,20261027)
def savepred(name,pred,idx=I,folds=None):
    with (OUT/(name+'_predictions.csv')).open('w') as f:
        w=csv.writer(f);w.writerow(['population','subpopulation','system','group','fold','radius_kpc','gbar_m_s2','gobs_m_s2','pred_gobs_m_s2','log10_residual'])
        for j,i in enumerate(idx):w.writerow([pop[i],subpop[i],rows[i]['system'],group[i],int(fold[i]) if folds is None else folds[j],10**u[i],1e-10*10**x[i],1e-10*10**y[i],float(1e-10*10**np.clip(pred[j],-300,300)),float(pred[j]-y[i])])

meta=json.loads(args.metadata.read_text())
result={'status':'Exploratory heterogeneous comparison. Raw observations and derived mass estimators retain different systematic assumptions. No universal-DM or memory detection claim.',
'inputs':meta,'variables':{'x':'log10(gbar/[1e-10 m/s^2])','u':'log10(radius/kpc)','y':'log10(gobs/[1e-10 m/s^2])'},
'training':'Equal population total weights; equal validation-group weights within each population; equal radial point weights within group. dSph subtypes share one training population. Two nearby MaNGA BCGs conservatively share a group.',
'validation':'Five physical-system folds,3innerfolds for shared-familyselection; seed20261027,innerseed39173+outerfold. Known same systems across sources are grouped. Not all dwarfs in samehost are grouped; hostenvironment sensitivity separately.',
'parameter_bounds':{'a':[-3,3],'b':[-2,2],'c':[-1,1],'p':[.1,.9]},
'model_definitions':{'RAR':'g=gbar/[1-exp(-sqrt(gbar/a0))]','RAR_exponent':'sqrt replaced with (gbar/a0)^p','RAR_radius':'log10(a_eff/1e-10)=a+b*u','RAR_radius_quad':'log10(a_eff/1e-10)=a+b*u+c*u^2','RAR_potential_quad':'same quadratic in q=x+u; potential proxy only','RAR_density_quad':'same quadratic in q=x-u; density proxy only','quadratic_2D':'y=c0+c1*x+c2*u+c3*x*x+c4*x*u+c5*u*u','RAR_piecewise':'standard RAR with one loga0 per population, population order given'},
'piecewise_population_order':POPS,
'frozen_SPARC_RAR':metrics(rar(x,a_fixed)),
'baryons_only':metrics(x),
'models':{},'leave_population_out':{}}
savepred('frozen_SPARC_RAR',rar(x,a_fixed));savepred('baryons_only',x)
CV={}
for kind in KINDS:
    pred=np.full(len(y),np.nan);foldcoef=[]
    for k in range(5):
        tr=I[fold!=k];te=I[fold==k];assert not set(group[tr])&set(group[te])
        c=fit(kind,tr);pred[te]=predict(kind,c,te);foldcoef.append(c.tolist())
    assert np.isfinite(pred).all();c=fit(kind,I);CV[kind]=pred
    result['models'][kind]={'all_data_coefficients':c.tolist(),'all_data_training_metrics':metrics(predict(kind,c,I)),
      'cross_validation':metrics(pred),'fold_coefficients':foldcoef}
    savepred(kind,pred)
    print(kind,'balanced RMS',round(result['models'][kind]['cross_validation']['balanced_population_system_rms_dex'],4),'median%',{p:round(v['median_absolute_acceleration_fraction']*100,1) for p,v in result['models'][kind]['cross_validation']['populations'].items()},flush=True)

if not args.quick:
    nested=np.full(len(y),np.nan);selections=[]
    for k in range(5):
        train=I[fold!=k];test=I[fold==k];innerfold=makefold(train,3,39173+k);scores={}
        for kind in SHARED:
            ip=np.full(len(train),np.nan)
            for j in range(3):
                tr=train[innerfold!=j];te=train[innerfold==j]
                assert not set(group[tr])&set(group[te])
                c=fit(kind,tr);ip[innerfold==j]=predict(kind,c,te)
            scores[kind]=score(ip,train)
        chosen=min(scores,key=scores.get);nested[test]=CV[chosen][test]
        selections.append({'outer_fold':k,'chosen':chosen,'inner_balanced_MSE_dex2':scores})
        print('nested fold',k,chosen,flush=True)
    result['nested_shared_selection']={'metrics':metrics(nested),'selection':selections};savepred('nested_shared',nested)

for p in POPS:
    te=I[pop==p];testgroups=set(group[te]);tr=I[np.array([g not in testgroups for g in group])]
    assert not set(group[tr])&testgroups
    d={}
    for kind in SHARED:
        c=fit(kind,tr);pr=predict(kind,c,te);d[kind]={'coefficients':c.tolist(),'metrics':met(pr,te)}
        savepred('leave_'+p+'_out_'+kind,pr,te)
    result['leave_population_out'][p]={'training_rows':len(tr),'test_rows':len(te),'models':d}
    print('population holdout',p,flush=True)

# Sensitivity to withholding all dSphs in one host environment simultaneously.
result['dwarf_host_holdout']={}
for h in sorted(set(env[pop=='dSph'])):
    if not h:continue
    te=I[(pop=='dSph')&(env==h)]
    if len(te)<3:continue
    tr=I[~((pop=='dSph')&(env==h))]
    d={}
    for kind in SHARED:
        c=fit(kind,tr);pr=predict(kind,c,te);d[kind]={'coefficients':c.tolist(),'metrics':met(pr,te)}
    result['dwarf_host_holdout'][h]={'test_rows':len(te),'models':d}

result['fit_calls']=COUNTS;result['fit_warnings']=FIT_WARNINGS
(OUT/'broad_results.json').write_text(json.dumps(result,indent=2))
print('complete',len(rows),'rows',len(set(group)),'systems',COUNTS,'fits',flush=True)
