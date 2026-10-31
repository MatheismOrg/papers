from pathlib import Path
import csv,json,hashlib
import numpy as np

R=Path(__file__).resolve().parent;C=R.parent/'dm_unified_limits/results'
data=R.parent/'dm_broad_joint/combined_input.csv';rows=list(csv.DictReader(data.open()));J=json.loads((C/'results.json').read_text())
assert hashlib.sha256(data.read_bytes()).hexdigest()==J['input_sha256']
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows]);I=np.arange(len(rows))
xx=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows]);uu=np.log10([float(r['radius_kpc']) for r in rows]);yy=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows])

def direct(k,c,x,u):
    q=x+u;m=x+2*u
    if k=='RAR':ae=10**c[0]+np.zeros_like(x)
    elif k=='potential_quad':ae=10**(c[0]+c[1]*q+c[2]*q*q)
    elif k=='mass_bounded':
        z=np.tanh((m-c[3])/c[4]);ae=10**(c[0]+c[1]*z+c[2]*z*z)
    elif k=='power_sum':ae=10**c[0]+10**(c[1]-c[3]*q)+10**(c[2]+c[4]*q)
    elif k=='half_power':ae=10**c[0]+10**(c[1]-.5*q)+10**(c[2]+c[3]*q)
    else:
        z=q if k=='mixed_power' else x+c[5]*u if k=='index_power' else m
        ae=10**c[0]+10**(c[1]-c[3]*m)+10**(c[2]+c[4]*z)
    gb=10**x
    if k=='mass_additive':pred=gb+np.sqrt(gb*ae)
    elif k=='mass_screened':pred=gb+np.sqrt(gb*ae)*np.exp(-np.sqrt(gb/(10**c[0])))
    else:pred=gb/(-np.expm1(-np.sqrt(gb/ae)))
    return np.log10(pred)

def B(pr,idx):
    return np.sqrt(np.mean([np.mean([np.mean((pr[(pop[idx]==p)&(group[idx]==g)]-yy[idx][(pop[idx]==p)&(group[idx]==g)])**2) for g in set(group[idx][pop[idx]==p])]) for p in set(pop[idx])]))

tabs=[];maxerr=0
for file in sorted(C.glob('*.csv')):
    model=next((k for k in J['models'] if file.name.startswith(k+'_')),None)
    if model is None:continue
    records=list(csv.DictReader(file.open()));ind=np.array([int(r['row']) for r in records]);pr=np.log10([float(r['gpred_m_s2'])/1e-10 for r in records]);fold=np.array([r['fold'] for r in records]);d=J['models'][model]
    tail=file.stem[len(model)+1:];expect=np.empty(len(ind))
    assert len(ind)>0 and len(set(ind))==len(ind)
    if tail=='all':
        expect=direct(model,d['coefficients'],xx[ind],uu[ind]);met=d['training']
    else:
        for f in set(fold):
            mask=fold==f;te=ind[mask]
            if tail.startswith('cv_'):
                seed=tail[3:];co=d['CV'][seed]['coefficients'][int(f)];tr=np.setdiff1d(I,te);met=d['CV'][seed]['metrics']
            elif tail=='population_holdout':
                assert set(pop[te])=={f};co=d['population_holdout']['coefficients'][f];tr=I[pop!=f];met=d['population_holdout']['metrics']
            else:
                assert tail=='host_'+f
                names=['Milky_Way','G/L'] if f=='MW' else ['M31','A/L'];assert set(te)==set(I[(pop=='dSph')&np.isin(env,names)]);assert len(te)==(26 if f=='MW' else 18)
                co=d['host_holdout'][f]['coefficients'];tr=np.setdiff1d(I,te);met=d['host_holdout'][f]['metrics']
            assert not set(group[tr])&set(group[te]);expect[mask]=direct(model,co,xx[te],uu[te])
    err=float(np.max(abs(expect-pr)));maxerr=max(maxerr,err);assert err<1e-11
    assert abs(B(pr,ind)-met['B_dex'])<1e-11
    assert max(abs(pr-yy[ind]-np.array([float(r['log10_residual']) for r in records])))<1e-12
    tabs.append(file.name)

# Illustrative point-source probes, not an ephemeris likelihood or covariant test.
GM_sun=1.32712440018e20;AU=149597870700.;KPC=3.0856775814913673e19
au=np.array([.387,1.,5.2,10.,30.,100.]);r=au*AU;gb=GM_sun/r**2;x=np.log10(gb/1e-10);u=np.log10(r/KPC)
probes=[]
for model,d in J['models'].items():
    pred=direct(model,d['coefficients'],x,u);ratio=10**(pred-x)
    for j,v in enumerate(au):probes.append({'model':model,'radius_AU':v,'radius_kpc':r[j]/KPC,'gbar_m_s2':gb[j],'gpred_m_s2':10**pred[j]*1e-10,'predicted_over_Newtonian':ratio[j],'absolute_anomaly_m_s2':gb[j]*np.expm1((pred[j]-x[j])*np.log(10.))})
with (R/'solar_point_source_probes.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(probes[0]));w.writeheader();w.writerows(probes)

# Check fixed-radius monotonicity on a stated grid; no global proof follows.
gridx=np.linspace(-15,15,601);gridu=np.linspace(-12,6,73);monotonic={}
for model,d in J['models'].items():
    vals=direct(model,d['coefficients'],gridx[:,None],gridu[None,:]);der=np.diff(vals,axis=0)/(gridx[1]-gridx[0]);z=np.unravel_index(np.argmin(der),der.shape)
    monotonic[model]={'minimum_dlogg_dloggbar':float(der[z]),'grid_x':float(gridx[z[0]]),'grid_u':float(gridu[z[1]]),'negative_steps':int(np.sum(der<0))}

out={'status':'passed','input_sha256':J['input_sha256'],'candidate_results_sha256':hashlib.sha256((C/'results.json').read_bytes()).hexdigest(),'prediction_tables_checked':len(tabs),'maximum_direct_equation_error_dex':maxerr,'checks':['all-data and held-out predictions','equal population/group RMS','residuals','group disjointness','26/18 host membership'],'solar_probe_constants':{'GM_sun_m3_s2':GM_sun,'AU_m':AU,'kpc_m':KPC},'solar_interpretation':'Illustrative isolated-point-source predictions against Newtonian baryonic acceleration. No actual ephemeris, external field, relativistic or Cassini likelihood fitted.','monotonicity_grid':{'x_range':[-15,15],'x_points':601,'u_range':[-12,6],'u_points':73,'results':monotonic},'files':tabs}
(R/'candidate_verification.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['files','monotonicity_grid']},indent=2))
for z in probes:
    if z['model']=='mixed_power':print(z)
