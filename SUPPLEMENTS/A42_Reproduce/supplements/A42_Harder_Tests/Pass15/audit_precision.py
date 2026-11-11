"""Separate precision, regulator, and massless-limit checks for pass 15."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from lattice_stress import simulate,subtract,dispersion
from continuum_reference import simulate as continuum,geometry

ROOT=Path(__file__).resolve().parent

def main():
    rows=json.loads((ROOT/'lattice_results.json').read_text())['rows']
    refinements=[];reference_checks=[]
    for kind in ['tanh','bump']:
        old=np.load(ROOT/f'continuum_{kind}.npz')
        stats,d=continuum(kind=kind,cutoff=256.,order=96,samples=2401,tol=3e-13)
        delta={v:float(np.max(abs(old[v]-d[v]))) for v in ['rho_ren','p_ren']}
        reference_checks.append(dict(kind=kind,changes='K=128 to256, order48 to96, tolerance2e-12 to3e-13',differences=delta))
        assert max(delta.values())<3e-11,delta
        for stencil,spacing in [('nearest',.05),('improved',.2),('improved',.05),('centered',.05)]:
            r=next(r for r in rows if r['kind']==kind and r['stencil']==stencil and r['spacing']==spacing)
            old=np.load(ROOT/r['history_file'])
            stats,d=simulate(spacing,stencil,kind,order=64,tol=2e-13)
            delta={v:float(np.max(abs(old[v]-d[v]))) for v in ['rho','p','rho2','p2']}
            refinements.append(dict(kind=kind,stencil=stencil,spacing=spacing,differences=delta))
            assert max(delta.values())<5e-10,delta
    # A denser time grid tests the work integral, separately from spatial grids.
    work=[]
    for stencil in ['nearest','improved','centered']:
        for samples in [1201,2401,4801]:
            stats,d=simulate(.2,stencil,samples=samples)
            work.append(dict(stencil=stencil,samples=samples,max_abs_work=stats['max_abs_work']))
        assert work[-1]['max_abs_work']<5e-11,work[-1]
    # Independent adaptive BZ quadrature at eta=0 in the pulse: H=0 and R!=0.
    # Split intervals around each IR zero; no target coefficient in the integral.
    coefficient=[];a,ap,H,Hp=geometry(0.,kind='bump');target=Hp/(12*np.pi*a*a)
    for stencil in ['nearest','improved','centered']:
        spacing=.4;top=np.pi/spacing
        for mass in [.1,.01,.001]:
            points=sorted(set([0.,mass*a,10*mass*a,top/2,top-10*mass*a,top-mass*a,top]))
            points=[p for p in points if 0<=p<=top]
            def integrand(k):
                M2=(mass*a)**2;w2=dispersion(k,spacing,stencil)+M2
                q2=subtract(w2,M2,H,Hp)[2]
                return mass*mass*q2/np.pi
            value=sum(quad(integrand,l,r,epsabs=2e-13,epsrel=2e-11,limit=200)[0]
                      for l,r in zip(points[:-1],points[1:]))
            coefficient.append(dict(stencil=stencil,spacing=spacing,mass=mass,
                                    anomaly_coefficient=value/target))
        expected=2 if stencil=='centered' else 1
        assert abs(coefficient[-1]['anomaly_coefficient']-expected)<1e-6,coefficient[-1]
    out=dict(reference_checks=reference_checks,lattice_refinements=refinements,
             work_refinement=work,massless_coefficients=coefficient,all_checks_passed=True)
    (ROOT/'audit_precision.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
