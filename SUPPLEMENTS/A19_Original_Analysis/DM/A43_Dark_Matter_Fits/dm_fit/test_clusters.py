from pathlib import Path
import sys,json,csv,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'deps'))
from astropy.io import fits
ROOT=Path(__file__).resolve().parent
fit=json.loads((ROOT/'results/fit_results.json').read_text())
G=6.67430e-11; MSUN=1.988409870698051e30; KPC=3.085677581491367e19
def interp(r,rs,m):return float(np.exp(np.interp(np.log(r),np.log(rs),np.log(m))))
rows=[]
for path in sorted((ROOT/'data/xcop').glob('*_mstar.fits')):
    name=path.stem.replace('_mstar','')
    with fits.open(path) as star,fits.open(path.with_name(name+'_fgas_profile.fits')) as gas,fits.open(path.with_name(name+'_hydro_mass.fits')) as hydro:
        h=gas['FGAS'];r=h.header['R500']; mt=h.header['M500']*1e14
        rs=star['MSTAR_SMOOTHED'].data['RADIUS'];ms=star['MSTAR_SMOOTHED'].data['MSTAR']
        rg=h.data['RADIUS']*r;mg=h.data['MGAS']
        assert rs.min()<r<rs.max() and rg.min()<r<rg.max()
        mstar_proj=interp(r,rs,ms)
        # Ghizzardi et al. (2021), section 4.1: deprojection factor at R500.
        mstar=.75*mstar_proj;mgas=interp(r,rg,mg)
        gb=G*(mstar+mgas)*MSUN/(r*KPC)**2;go=G*mt*MSUN/(r*KPC)**2
        x=np.log10(gb/1e-10)
        poly=1e-10*10**np.polynomial.polynomial.polyval(x,fit['models']['poly2']['coefficients'])
        a0=fit['models']['RAR']['coefficients'][0];rar=gb/(-np.expm1(-np.sqrt(gb/a0)))
        d=hydro['HYDRO_MASS'].data
        goforw=G*interp(r,d['RADIUS'],d['M_FORW'])*MSUN/(r*KPC)**2
        goein=G*interp(r,d['RADIUS'],d['M_EIN'])*MSUN/(r*KPC)**2
        gb_no_deproj=G*(mstar_proj+mgas)*MSUN/(r*KPC)**2
        rar_no_deproj=gb_no_deproj/(-np.expm1(-np.sqrt(gb_no_deproj/a0)))
        rows.append(dict(cluster=name,r500_kpc=r,mtot_msun=mt,mgas_msun=mgas,mstar_msun=mstar,gbar=gb,gobs=go,
            polynomial_prediction=poly,rar_prediction=rar,ratio_poly=go/poly,ratio_rar=go/rar,
            ratio_rar_forward=goforw/rar,ratio_rar_einasto=goein/rar,ratio_rar_without_stellar_deprojection=go/rar_no_deproj))
res={'scope':'Seven clusters with public stellar profiles; one comparison at R500 each; no cluster refitting',
     'method':'Published hydrostatic M500 plus gas mass and 0.75 times projected stellar mass at R500. Log-log interpolation. HSE and stellar-mass assumptions retained; not seven new raw acceleration measurements.',
     'clusters':rows,'median_observed_over_polynomial':float(np.median([r['ratio_poly'] for r in rows])),
     'median_observed_over_rar':float(np.median([r['ratio_rar'] for r in rows])),
     'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data/xcop').glob('*.fits')}}
(ROOT/'results/cluster_test.json').write_text(json.dumps(res,indent=2))
with (ROOT/'results/cluster_test.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps({k:v for k,v in res.items() if k!='source_hashes'},indent=2))
