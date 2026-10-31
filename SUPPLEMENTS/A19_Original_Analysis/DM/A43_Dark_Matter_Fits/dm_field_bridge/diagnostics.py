"""Post-fit numerical checks and descriptive diagnostics; does not tune any fit."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
import run

R=Path(__file__).resolve().parent
def readcsv(p):return list(csv.DictReader(p.open()))
def exact_shifted(d):
    r=run.HALF
    A=np.sqrt(1+(r-d)**2); B=np.sqrt(1+(r+d)**2)
    return 4*r**3/(A*B*(A+B)*(A*B+d*d+1-r*r))
def exact_H(z,kind):
    # Log-distance integration resolves both source and kernel scales.
    def f(l):
        d=np.exp(l);t=z*d
        p=.5*(1+t) if kind=='NLG' else t
        return z*d*p*np.exp(-t)*exact_shifted(d)
    lo=-40-max(0,np.log(z)); hi=np.log(90/z)
    cuts=sorted(set([lo,hi]+[x for x in [-20.,-10.,0.,10.,20.,-np.log(z)] if lo<x<hi]))
    return sum(quad(f,a,b,epsabs=1e-18*min(1,z if kind=='NLG' else z*z),epsrel=1e-10,limit=200)[0] for a,b in zip(cuts[:-1],cuts[1:]))

def main():
    data=readcsv(R/'input.csv');basis=readcsv(R/'basis.csv')
    H={k:PchipInterpolator([float(r['log10z']) for r in basis],np.log10([float(r[k+'_H']) for r in basis])) for k in ['NLG','MOG']}
    dvals=np.logspace(-9,9,361)
    shifted_errors=[abs(run.shifted(d)/exact_shifted(d)-1) for d in dvals]
    rng=np.random.default_rng(61029)
    zs=np.r_[np.logspace(-8,8,33),10**rng.uniform(-8,8,100)]
    checks=[]
    for kind in ['NLG','MOG']:
        for z in zs:
            direct=exact_H(z,kind);interp=10**H[kind](np.log10(z))
            checks.append({'kernel':kind,'z':float(z),'independent_H':direct,'interpolated_H':float(interp),'relative_error':float(interp/direct-1)})
    run.writecsv(R/'quadrature_checks.csv',checks)
    dwarf={r['system']:r for r in readcsv(R.parent/'dm_broad_dwarfs/results/dwarfs_benchmark.csv')}
    subset={'all':[r['system'] for r in data],
            'classical':[r['system'] for r in data if r['subpopulation']=='classical_dSph'],
            'ultrafaint':[r['system'] for r in data if r['subpopulation']=='ultrafaint_dSph'],
            'previously_defined_higher_precision':[r['system'] for r in data if max(float(dwarf[r['system']]['sigma_err_plus_km_s']),float(dwarf[r['system']]['sigma_err_minus_km_s']))/float(dwarf[r['system']]['sigma_los_km_s'])<=.25 and dwarf[r['system']]['possible_tidal_or_classification_caution']=='False']}
    bounds=readcsv(R/'mass_budget_bounds.csv');bdiag={}
    for ml in [1,2,3]:
        bdiag[str(ml)]={}
        for name,names in subset.items():
            rows=[r for r in bounds if int(r['stellar_ML'])==ml and r['system'] in names]
            D=np.array([float(r['gobs_over_gbar']) for r in rows]);mx=max(rows,key=lambda r:float(r['necessary_alpha_any_spherical_profile']))
            bdiag[str(ml)][name]={'count':len(rows),'NLG_general_violations':int(sum(D>1+2*10.94)),'MOG_general_violations':int(sum(D>1+2*8.89)),
                 'NLG_Plummer_violations':int(sum(D>1+10.94)),'MOG_Plummer_violations':int(sum(D>1+8.89)),
                 'median_D':float(np.median(D)),'largest_required_alpha_general':float(mx['necessary_alpha_any_spherical_profile']),
                 'largest_required_alpha_Plummer':float(mx['necessary_alpha_Plummer']),'largest_alpha_system':mx['system']}
    source_errors=[]
    for r in data:
        src=dwarf[r['system']];rh=float(r['radius_kpc']);gb=float(r['gbar_m_s2'])
        M=2*gb*(rh*3.085677581491367e19)**2/(6.67430e-11*1.988409870698051e30)
        source_errors.append(abs(M/float(src['stellar_mass_Msun'])-1))
    result=json.loads((R/'results.json').read_text())
    c=result['models_by_ML']['2']['NLG_fit']['coefficients']
    allalpha=10**c[0]
    # At mu -> infinity both kernels are alpha delta, giving g=(1+alpha)g_b.
    D=np.array([float(r['gobs_m_s2'])/float(r['gbar_m_s2']) for r in data])
    analytic_alpha=10**np.mean(np.log10(D))-1
    fitsets={}
    for f in result['fits']:
        key=(f['model'],f['stellar_ML'],f['split'])
        if key not in fitsets or f['cost']<fitsets[key]['cost']:fitsets[key]=f
    hit_mu=[f for f in fitsets.values() if f['model']!='RAR_fit' and f['coefficients'][1]>4.999]
    summary={'purpose':'Post-fit checks; no parameter or model tuning',
        'shifted_Plummer_max_relative_error_361_distances':max(shifted_errors),
        'kernel_normalizations':{kind:quad(lambda t:(.5*(1+t) if kind=='NLG' else t)*np.exp(-t),0,np.inf)[0] for kind in ['NLG','MOG']},
        'basis_max_relative_error_133_independent_points_per_kernel':{kind:max(abs(x['relative_error']) for x in checks if x['kernel']==kind) for kind in ['NLG','MOG']},
        'basis_bounds':{kind:[min(float(r[kind+'_H']) for r in basis),max(float(r[kind+'_H']) for r in basis)] for kind in ['NLG','MOG']},
        'reconstructed_total_stellar_mass_max_relative_error':max(source_errors),
        'refitted_NLG_alpha_ML2':allalpha,'constant_boost_analytic_optimum_alpha_ML2':analytic_alpha,
        'refitted_NLG_alpha_relative_difference_from_constant_boost':allalpha/analytic_alpha-1,
        'kernel_fit_sets_at_mu_upper_bound':len(hit_mu),'kernel_fit_sets_total':sum(k[0]!='RAR_fit' for k in fitsets),
        'bound_subsets':bdiag,
        'limitations':'All bound counts use central estimates and selected stellar M/L; they are not significance levels. Higher-precision criterion was defined in the earlier dwarf curation. Isolation, spherical equilibrium and Wolf acceleration reconstruction are conditional assumptions. No host source or Jeans profiles are fitted.'}
    (R/'diagnostics.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
