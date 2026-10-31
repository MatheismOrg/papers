from pathlib import Path
import json,csv,hashlib
import numpy as np
R=Path(__file__).resolve().parent
rr=list(csv.DictReader((R/'predictions.csv').open()));inp=list(csv.DictReader((R/'inputs.csv').open()));q=json.loads((R/'quadrature.json').read_text());result=json.loads((R/'results.json').read_text())
G=6.67430e-11;MSUN=1.988409870698051e30;KPC=3.085677581491367e19;A0=1.2e-10
def pred(g,nu):
    if nu=='simple':return .5*g+np.sqrt(.25*g*g+A0*g)
    return g/(-np.expm1(-np.sqrt(g/A0)))
maxquad=0.;maxcore=0.;records=[]
for row in rr:
    nu,kind=row['model'].split('_',1)
    if kind not in ['repair','frozen_gas']:continue
    i=int(row['row']);dat=inp[i];z=q[i];s=np.array(z['radius_nodes_kpc']);w=np.array(z['mass_weights_msun']);s2=np.array(z['check_nodes_kpc']);w2=np.array(z['check_weights_msun']);eta,rc=(10.,400.) if kind=='frozen_gas' else 10**np.array(json.loads(row['coefficients_json']))
    b=float(row['gbar_m_s2']);fac=G*MSUN/(float(row['radius_kpc'])*KPC)**2
    weighted=np.sum(w*np.exp(-s/rc));wp=np.sum(w2*np.exp(-s2/rc));p1=pred(b+fac*eta*weighted,nu);p2=pred(b+fac*eta*wp,nu);change=abs(np.log10(p2/p1));maxquad=max(maxquad,change)
    r0=float(dat['inner_radius_kpc']);m0=float(dat['inner_mass_msun']);span=0.
    if r0>0:
        inner=np.sum(w[s<r0]*np.exp(-s[s<r0]/rc));lo=weighted-inner+m0*np.exp(-r0/rc);hi=weighted-inner+m0
        plo=pred(b+fac*eta*lo,nu);phi=pred(b+fac*eta*hi,nu);span=float(np.log10(phi/plo));maxcore=max(maxcore,span)
    records.append({'case':row['case'],'model':row['model'],'split':row['split'],'row':i,'refined_quadrature_change_dex':float(change),'unresolved_inner_mass_full_bound_span_dex':span})
with (R/'numerical_diagnostics.csv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(records[0]));wr.writeheader();wr.writerows(records)

boot={};rng=np.random.default_rng(20261028)
for split in ['host_holdout','source_CLASH','source_XCOP']:
    chosen=[r for r in rr if r['case']=='combined_forward_mid' and r['split']==split]
    systems=sorted({r['system'] for r in chosen});sources=sorted({r['source'] for r in chosen});mult=rng.poisson(1,(5000,len(systems)));smap={s:i for i,s in enumerate(systems)}
    scores={}
    for model in ['simple_repair','simple_refit_a0','simple_frozen_gas']:
        z=[r for r in chosen if r['model']==model];bys=[]
        for src in sources:
            gs=sorted({r['system'] for r in z if r['source']==src});idx=[smap[g] for g in gs];e=np.array([np.mean([float(r['log10_residual'])**2 for r in z if r['system']==g]) for g in gs]);den=mult[:,idx].sum(axis=1);bys.append(np.divide(mult[:,idx]@e,den,out=np.full(5000,np.nan),where=den>0))
        scores[model]=np.sqrt(np.mean(bys,axis=0))
    for base in ['simple_refit_a0','simple_frozen_gas']:
        diff=scores['simple_repair']-scores[base];good=np.isfinite(diff)
        boot[split+' repair minus '+base]={'conditional_95_interval_dex':np.percentile(diff[good],[2.5,97.5]).tolist(),'finite_draws':int(sum(good)),'fraction_better_descriptive':float(np.mean(diff[good]<0))}
out={'status':'passed','maximum_refined_quadrature_change_dex':float(maxquad),'maximum_unresolved_core_bound_span_dex':float(maxcore),'checks':len(records),'bootstrap':boot,'bootstrap_caveat':'Paired source-stratified whole-host Poisson resampling of frozen predictions only. Does not include refitting, source modeling, observational covariance or candidate selection uncertainty.','parameter_bounds_hit':sum(any(x!=0 for x in f['active_bounds']) for f in result['fits']),'optimizer_failures':result['optimizer_failures']}
(R/'diagnostics.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
