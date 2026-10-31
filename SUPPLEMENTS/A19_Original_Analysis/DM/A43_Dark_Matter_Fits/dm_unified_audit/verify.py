from pathlib import Path
import csv, json, hashlib, math
import numpy as np

ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'dm_unified_limits'
INPUT=ROOT.parent/'dm_broad_joint/combined_input.csv'
J=json.loads((SRC/'results/results.json').read_text())
ROWS=list(csv.DictReader(INPUT.open()))
assert hashlib.sha256(INPUT.read_bytes()).hexdigest()==J['input_sha256']
assert hashlib.sha256((SRC/'PROTOCOL.txt').read_bytes()).hexdigest()==J['protocol_sha256']

def direct(model,co,gb,rk):
    # Independent dimensionful evaluation of the published family definitions.
    c=np.asarray(co,float); gb=np.asarray(gb,float); rk=np.asarray(rk,float)
    s=gb/1e-10; mass=s*rk*rk; potential=s*rk
    if model=='RAR': ae=gb*0+1e-10*10**c[0]
    elif model=='potential_quad':
        q=np.log10(potential);ae=1e-10*10**(c[0]+c[1]*q+c[2]*q*q)
    elif model=='mass_bounded':
        z=np.tanh((np.log10(mass)-c[3])/c[4]);ae=1e-10*10**(c[0]+c[1]*z+c[2]*z*z)
    elif model=='power_sum':
        ae=1e-10*(10**c[0]+10**c[1]*potential**(-c[3])+10**c[2]*potential**c[4])
    elif model=='half_power':
        ae=1e-10*(10**c[0]+10**c[1]*potential**(-.5)+10**c[2]*potential**c[3])
    else:
        pos=potential if model=='mixed_power' else s*rk**c[5] if model=='index_power' else mass
        ae=1e-10*(10**c[0]+10**c[1]*mass**(-c[3])+10**c[2]*pos**c[4])
    if model=='mass_additive': gp=gb+np.sqrt(gb*ae)
    elif model=='mass_screened':gp=gb+np.sqrt(gb*ae)*np.exp(-np.sqrt(gb/(1e-10*10**c[0])))
    else:gp=gb/(-np.expm1(-np.sqrt(gb/ae)))
    return gp,ae

def metric(rs):
    # Equal populations, then equal groups, then equal radial observations.
    bypop={}
    for p in sorted({r['population'] for r in rs}):
        subset=[r for r in rs if r['population']==p]
        gg=sorted({r['group'] for r in subset})
        vals=[np.mean([r['resid']**2 for r in subset if r['group']==g]) for g in gg]
        bypop[p]=float(np.mean(vals))
    return math.sqrt(float(np.mean(list(bypop.values())))),bypop

def close(a,b,tol=2e-12):
    assert abs(float(a)-float(b))<tol,(a,b)

checks=[]; maxerr=0.; predicted={}; sample_rows={}
for f in sorted((SRC/'results').glob('*.csv')):
    if not any(f.stem.startswith(k+'_') for k in J['models']):continue
    model=next(k for k in sorted(J['models'],key=len,reverse=True) if f.stem.startswith(k+'_'))
    suffix=f.stem[len(model)+1:];m=J['models'][model]
    rows=list(csv.DictReader(f.open())); ids=[int(r['row']) for r in rows]
    assert len(ids)==len(set(ids))
    if suffix=='all':rec=m['training'];coeff={None:m['coefficients']}
    elif suffix.startswith('cv_'):rec=m['CV'][suffix[3:]]['metrics'];coeff={str(i):c for i,c in enumerate(m['CV'][suffix[3:]]['coefficients'])}
    elif suffix=='population_holdout':rec=m[suffix]['metrics'];coeff=m[suffix]['coefficients']
    else:
        h=suffix[5:];rec=m['host_holdout'][h]['metrics'];coeff={h:m['host_holdout'][h]['coefficients']}
        aliases={'MW':{'Milky_Way','G/L'},'M31':{'M31','A/L'}}[h]
        expected={i for i,r in enumerate(ROWS) if r['population']=='dSph' and r['environment'] in aliases}
        assert set(ids)==expected
        assert len(ids)=={'MW':26,'M31':18}[h]
    if not suffix.startswith('host_'):assert set(ids)==set(range(len(ROWS)))
    for r in rows:
        org=ROWS[int(r['row'])]
        for name in ['population','system','group','environment']:assert r[name]==org[name]
        for name in ['gbar_m_s2','radius_kpc','gobs_m_s2']:close(float(r[name])/float(org[name]),1.)
        key=None if suffix=='all' else r['fold']
        gp,_=direct(model,coeff[key],float(r['gbar_m_s2']),float(r['radius_kpc']))
        err=abs(math.log10(float(gp)/float(r['gpred_m_s2'])));maxerr=max(maxerr,err)
        assert err<2e-11,(f.name,r['row'],err)
        r['resid']=math.log10(float(r['gpred_m_s2'])/float(r['gobs_m_s2']))
        close(r['resid'],r['log10_residual'])
    for g in {r['group'] for r in rows}:assert len({r['fold'] for r in rows if r['group']==g})==1
    if suffix=='population_holdout':assert all(r['fold']==r['population'] for r in rows)
    score,popmse=metric(rows);close(score,rec['B_dex'])
    assert len(rows)==rec['rows'] and len({r['group'] for r in rows})==rec['groups']
    for p,v in rec['populations'].items():
        e=np.array([r['resid'] for r in rows if r['population']==p]);assert len(e)==v['rows']
        close(np.sqrt(np.mean(e*e)),v['RMS_dex']);close(np.sqrt(popmse[p]),v['group_RMS_dex'])
        close(np.median(abs(10**e-1))*100,v['median_acceleration_error_percent'],1e-10)
        close(np.median(abs(10**(.5*e)-1))*100,v['median_speed_error_percent'],1e-10)
    checks.append({'file':f.name,'rows':len(rows),'B_dex':score})
    predicted[model,suffix]=np.array([r['resid'] for r in rows])
    sample_rows[model,suffix]=rows

# Verify optimizer records against coefficient selections and training exclusions.
fits={}
for v in J['fit_log']:fits.setdefault((v['model'],v['split']),[]).append(v)
assert len(fits)==J['fit_calls'] and len(J['fit_log'])==J['optimizer_calls']
selected_failures=[]
for (model,split),runs in fits.items():
    best=min(runs,key=lambda v:v['cost']);c=best['coefficients']
    if not best['success']:selected_failures.append((model,split))
    if split=='all':saved=J['models'][model]['coefficients'];tr=list(range(len(ROWS)))
    elif split.startswith('cv_'):
        _,seed,fold=split.split('_');saved=J['models'][model]['CV'][seed]['coefficients'][int(fold)]
        rr=sample_rows[model,'cv_'+seed];tr=[int(r['row']) for r in rr if r['fold']!=fold]
        test_groups={r['group'] for r in rr if r['fold']==fold}
        assert not test_groups.intersection(ROWS[i]['group'] for i in tr)
    elif split.startswith('population_'):
        p=split[len('population_'):];saved=J['models'][model]['population_holdout']['coefficients'][p]
        test_groups={r['group'] for r in ROWS if r['population']==p};tr=[i for i,r in enumerate(ROWS) if r['group'] not in test_groups]
        assert p not in {ROWS[i]['population'] for i in tr}
    else:
        h=split[5:];saved=J['models'][model]['host_holdout'][h]['coefficients']
        te=sample_rows[model,'host_'+h];test_groups={r['group'] for r in te};tr=[i for i,r in enumerate(ROWS) if r['group'] not in test_groups]
        assert len(tr)=={'MW':2884,'M31':2892}[h]
    assert np.max(abs(np.array(c)-np.array(saved)))<1e-13
    assert all(len(tr)==v['train_rows'] for v in runs)
    # Cost evaluated directly with the retained coefficients and reconstructed training set.
    rr=[]
    for i in tr:
        r=ROWS[i];gp,_=direct(model,c,float(r['gbar_m_s2']),float(r['radius_kpc']))
        rr.append(dict(r,resid=math.log10(float(gp)/float(r['gobs_m_s2']))))
    score,_=metric(rr);close(.5*score*score,best['cost'])
assert not selected_failures

# Recompute paired group-Poisson percentile intervals without importing fit code.
gg=sorted({r['group'] for r in ROWS});gi={g:i for i,g in enumerate(gg)}
mult=np.random.default_rng(20261028).poisson(1,size=(5000,len(gg)))
bootstrap_checks=0
for suffix in ['cv_20261027','population_holdout']:
    scores={}
    for model in J['models']:
        rs=sample_rows[model,suffix];bypop=[]
        for p in sorted({r['population'] for r in rs}):
            group_mse={g:np.mean([r['resid']**2 for r in rs if r['population']==p and r['group']==g]) for g in sorted({r['group'] for r in rs if r['population']==p})}
            ww=mult[:,[gi[g] for g in group_mse]];assert np.all(ww.sum(axis=1)>0)
            bypop.append(np.sum(ww*np.array(list(group_mse.values())),axis=1)/ww.sum(axis=1))
        scores[model]=np.sqrt(np.mean(bypop,axis=0))
    for key,v in J['bootstrap_'+suffix].items():
        model,base=key.split(' minus ');d=scores[model]-scores[base]
        assert np.max(abs(np.percentile(d,[2.5,97.5])-v['percentile_95_dex']))<2e-12
        close(np.mean(d<0),v['fraction_negative_descriptive_only']);bootstrap_checks+=1

# Finite isolated-Solar check: illustrative model predictions, no ephemeris fit.
G=6.67430e-11;MSUN=1.98847e30;KPC=3.0856775814913673e19;AU=149597870700.
solar=[]
for model in J['new_families']:
    for rau in [1.,9.58,30.]:
        r=rau*AU;gb=G*MSUN/r**2;gp,ae=direct(model,J['models'][model]['coefficients'],gb,r/KPC)
        solar.append({'model':model,'radius_AU':rau,'gbar_m_s2':gb,'a_eff_m_s2':float(ae),'gpred_m_s2':float(gp),'relative_excess':float(gp/gb-1),'extra_acceleration_m_s2':float(gp-gb)})
with (ROOT/'solar_stress.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(solar[0]));w.writeheader();w.writerows(solar)

result={'status':'passed','prediction_tables_checked':len(checks),'optimizer_fit_sets_checked':len(fits),'optimizer_attempts_checked':len(J['fit_log']),'optimizer_failures_recorded':len(J['optimizer_failures']),'selected_failed_fits':selected_failures,'bootstrap_intervals_checked':bootstrap_checks,'max_prediction_discrepancy_dex':maxerr,'input_sha256':J['input_sha256'],'protocol_sha256':J['protocol_sha256'],'checks':checks,'solar_constants':{'G_m3_kg_s2':G,'M_sun_kg':MSUN,'kpc_m':KPC,'AU_m':AU},'important_scientific_limitation':'Five of six new all-data fits have large finite isolated-Solar deviations despite asymptotic Newtonian recovery; see solar_stress.csv. No ephemeris, external-field quadrupole, or covariant-gravity likelihood is evaluated.'}
(ROOT/'verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
