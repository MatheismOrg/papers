"""Exploratory fixed-physical-radius cluster check; not independent new observations.

Default inputs are the sibling dm_fit directory. Install numpy scipy astropy.
The 0.75 projected-stellar mass conversion is an approximation away from R500.
No-star and full-projected-star calculations bracket that contribution, not all
stellar-systematic uncertainty. Profile values at three radii are correlated.
"""
from pathlib import Path
import argparse, csv, json, sys, hashlib
import numpy as np
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser(); p.add_argument('--data-root', type=Path, default=HERE.parent/'dm_fit')
args=p.parse_args(); ROOT=args.data_root.resolve()
sys.path.insert(0,str(ROOT/'deps'))
from astropy.io import fits

G=6.67430e-11; MSUN=1.988409870698051e30; KPC=3.085677581491367e19
rows=[]; skipped=[]
def interp(r,rs,m):
    assert np.min(rs)<=r<=np.max(rs)
    assert np.all(np.diff(rs)>0) and np.all(np.array(m)>0)
    return float(np.exp(np.interp(np.log(r),np.log(rs),np.log(m))))

for path in sorted((ROOT/'data/xcop').glob('*_mstar.fits')):
    name=path.stem.replace('_mstar','')
    with fits.open(path) as star, fits.open(path.with_name(name+'_fgas_profile.fits')) as gas, fits.open(path.with_name(name+'_hydro_mass.fits')) as hydro:
        hs=star['MSTAR_SMOOTHED'].data; hg=gas['FGAS']; hm=hydro['HYDRO_MASS'].data
        rs=hs['RADIUS']; rg=hg.data['RADIUS']*hg.header['R500']; rh=hm['RADIUS']
        for r in [500.,750.,1000.]:
            if not all(np.min(rr)<=r<=np.max(rr) for rr in [rs,rg,rh]):
                skipped.append([name,r]); continue
            mg=interp(r,rg,hg.data['MGAS']); ms=interp(r,rs,hs['MSTAR'])
            fac=G*MSUN/(r*KPC)**2
            rows.append(dict(cluster=name,radius_kpc=r,mgas_msun=mg,mstar_projected_msun=ms,
                gbar=fac*(mg+.75*ms),gbar_no_stars=fac*mg,gbar_full_projected_stars=fac*(mg+ms),
                gobs_forward=fac*interp(r,rh,hm['M_FORW']),
                gobs_nfw=fac*interp(r,rh,hm['M_NFW']),
                gobs_einasto=fac*interp(r,rh,hm['M_EIN'])))

gr=list(csv.DictReader((ROOT/'results/held_out_predictions.csv').open()))
fit=json.loads((ROOT/'results/fit_results.json').read_text())
prior=json.loads((ROOT/'results/two_scale_results.json').read_text())
prior_coef=np.array(prior['models']['quadratic_acceleration_radius']['coefficients'])
gal_a0=fit['models']['RAR']['coefficients'][0]
gg=np.array([r['galaxy'] for r in gr]); gx=np.log10([float(r['gbar_m_s2'])/1e-10 for r in gr])
gy=np.log10([float(r['gobs_m_s2'])/1e-10 for r in gr]); gz=np.log10([float(r['radius_kpc']) for r in gr])
gw=np.array([1/np.sum(gg==g) for g in gg])
groups=np.array([r['cluster'] for r in rows]); z=np.log10([r['radius_kpc'] for r in rows])
cw=np.array([1/np.sum(groups==g) for g in groups]); names=np.unique(groups)
def design(x,z): return np.column_stack([x*0+1,x,z,x*x,x*z,z*z])
gX=design(gx,gz)
def rar(x,la0):
    u=10**((x-la0)/2)
    return x-np.log10(-np.expm1(-u))
def rar_fit(x,y,w):
    return float(minimize_scalar(lambda a: np.sum(w*(rar(x,a)-y)**2),bounds=(-2,2),method='bounded').x)
def linear_fit(X,y,w):return np.linalg.lstsq(X*np.sqrt(w[:,None]), y*np.sqrt(w),rcond=None)[0]
def metric(pred,y):
    e=pred-y
    return dict(rmse_dex=float(np.sqrt(np.mean(e*e))),median_absolute_acceleration_fraction=float(np.median(abs(10**e-1))),
        maximum_absolute_acceleration_fraction=float(np.max(abs(10**e-1))),
        median_observed_over_predicted=float(np.median(10**(-e))),
        per_cluster={str(g):dict(median_absolute_acceleration_fraction=float(np.median(abs(10**e[groups==g]-1))),
             rmse_dex=float(np.sqrt(np.mean(e[groups==g]**2)))) for g in names})

all_results={}; predictions=[]
for target in ['forward','nfw','einasto']:
    for baryon in ['gbar','gbar_no_stars','gbar_full_projected_stars']:
        x=np.log10([r[baryon]/1e-10 for r in rows]); y=np.log10([r['gobs_'+target]/1e-10 for r in rows]); X=design(x,z)
        pp={'galaxy_RAR_transfer':rar(x,np.log10(gal_a0/1e-10)),
            'previous_R500_quadratic_transfer':X@prior_coef}
        for name in ['cluster_RAR_LOCO','joint_quadratic_LOCO','radius_power_law_LOCO']:
            pp[name]=np.zeros(len(rows))
        fold_params=[]
        for hold in names:
            test=groups==hold; tr=~test
            assert not set(groups[tr])&set(groups[test])
            la0=rar_fit(x[tr],y[tr],cw[tr]); pp['cluster_RAR_LOCO'][test]=rar(x[test],la0)
            coef=linear_fit(np.vstack([gX,X[tr]]),np.r_[gy,y[tr]],np.r_[gw,cw[tr]])
            pp['joint_quadratic_LOCO'][test]=X[test]@coef
            R=np.column_stack([np.ones(len(z)),z]); c=linear_fit(R[tr],y[tr],cw[tr])
            pp['radius_power_law_LOCO'][test]=R[test]@c
            fold_params.append(dict(held_out_cluster=str(hold),cluster_RAR_a0_m_s2=float(1e-10*10**la0),joint_quadratic_coefficients=coef.tolist(),radius_power_coefficients=c.tolist()))
        key=target+'__'+baryon
        all_results[key]=dict(models={k:metric(v,y) for k,v in pp.items()},fold_parameters=fold_params,
            full_cluster_RAR_a0_m_s2=float(1e-10*10**rar_fit(x,y,cw)))
        for i,r in enumerate(rows):
            predictions.append(dict(case=key,cluster=r['cluster'],radius_kpc=r['radius_kpc'],gbar_m_s2=1e-10*10**x[i],gobs_m_s2=1e-10*10**y[i],**{k:float(1e-10*10**v[i]) for k,v in pp.items()}))

result=dict(status='Exploratory stress check on previously inspected clusters, not external validation.',
    sampling='Fixed physical radii 500, 750, 1000 kpc, chosen independently of each cluster overdensity radius. No radius extrapolation.',
    target='Primary: archived forward hydrostatic mass profile M_FORW; NFW and Einasto profiles are sensitivity cases. All are model-derived.',
    baryons='Gas mass profile plus 0.75 projected stellar mass; the 0.75 factor is not derived at these radii. Also test zero and full projected stellar mass. These are sensitivity cases, not confidence bounds.',
    validation='Leave one whole cluster out for fitted models. Joint polynomial trains on all 149 galaxies and six clusters; only cluster prediction is evaluated. Every system has total weight one. Three radii per cluster are correlated and are not independent observations.',
    prior_transfer='Previous R500 polynomial is frozen, but its training already included the same seven clusters at R500. Its physical-radius transfer is a within-cluster radial stress test, not held-out-cluster validation.',
    counts=dict(clusters=len(names),radii_per_cluster=3,cluster_rows=len(rows),galaxies=len(np.unique(gg)),galaxy_rows=len(gr)),
    skipped=skipped,galaxy_RAR_a0_m_s2=gal_a0,previous_R500_coefficients=prior_coef.tolist(),
    input_rows=rows,cases=all_results,
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data/xcop').glob('*.fits')})
(HERE/'stress_test_results.json').write_text(json.dumps(result,indent=2))
for filename,records in [('fixed_radius_inputs.csv',rows),('fixed_radius_predictions.csv',predictions)]:
    with (HERE/filename).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader();w.writerows(records)
print(json.dumps({k:v for k,v in result.items() if k in ['counts','skipped']},indent=2))
for k,v in all_results['forward__gbar']['models'].items():
    print(k, {a:b for a,b in v.items() if a!='per_cluster'})
print('Cluster RAR a0',all_results['forward__gbar']['full_cluster_RAR_a0_m_s2'])
