"""Independent formula audits for the principal quadratic KDK generator.

This script does not import the pass16 evolution implementation. It checks
static maps/vacua, exact false-particle formulas, metric pressure, and the
continuous-shadow second-order adiabatic subtraction. The latter is a
specified extension to the step sequence, not an exact local discrete metric
stress tensor. Units hbar=c=1; variables ordered (q,p); signature (-,+).
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.linalg import logm

ROOT=Path(__file__).resolve().parent
MU=1/np.sqrt(np.pi)

def geometry(eta,kind='tanh',tau=.6):
    z=np.tanh(eta/tau);s=1-z*z
    if kind=='tanh':return 1+.5*z,.5*s/tau,-s*z/tau**2
    if kind=='bump':return 1+.5*s,-s*z/tau,s*(3*z*z-1)/tau**2
    raise ValueError(kind)

def params(kappa2,a,h,mu=MU):
    w=np.sqrt(kappa2+(mu*a)**2);y=h*w/2;d=1-y*y
    if d<=0:raise ValueError('No stable principal Gaussian vacuum: h*omega must be <2')
    theta=2*np.arcsin(y);Omega=theta/h;r=w*np.sqrt(d)
    wa=mu*mu*a/w;ya=h*wa/2
    Oa=wa/np.sqrt(d)
    b=wa/w-y*ya/d
    # Independent a derivative through w, y and d.
    waa=mu*mu/w-mu**4*a*a/w**3;yaa=h*waa/2
    ba=waa/w-wa*wa/w**2-(ya*ya+y*yaa)/d-2*y*y*ya*ya/d**2
    return dict(w=w,y=y,d=d,theta=theta,Omega=Omega,r=r,Oa=Oa,b=b,ba=ba)

def terms(kappa2,eta,h,kind='tanh',mu=MU):
    a,ap,app=geometry(eta,kind);v=params(kappa2,a,h,mu)
    f=v['b']*ap/2
    fp=(v['ba']*ap*ap+v['b']*app)/2
    Op=v['Oa']*ap;O=v['Omega']
    n2=f*f/(4*O*O)
    re2=fp/(4*O*O)-f*Op/(4*O**3)
    E2=O*n2
    P2=(O-a*v['Oa'])*n2-a*O*v['b']*re2
    return E2,P2,n2,re2

def static_audits():
    rows=[];J=np.array([[0.,1.],[-1.,0.]])
    for w,h in [(0.3,.2),(1.,.4),(2.,.5),(3.,.6)]:
        y=h*w/2;d=1-y*y;theta=2*np.arcsin(y);r=w*np.sqrt(d);O=theta/h
        kick=np.array([[1.,0.],[-h*w*w/2,1.]])
        drift=np.array([[1.,h],[0.,1.]])
        direct=kick@drift@kick
        rotation=np.array([[np.cos(theta),np.sin(theta)/r],[-r*np.sin(theta),np.cos(theta)]])
        C=np.diag([1/(2*r),r/2]);G=np.diag([O*r,O/r])
        reconstructed=-J@logm(direct)/h
        initial=np.diag([1/(2*w),w/2]);occupations=[]
        for count in [1,2,7,31]:
            S=np.linalg.matrix_power(direct,count);evolved=S@initial@S.T
            number=(w*evolved[0,0]+evolved[1,1]/w-1)/2
            target=np.sin(count*theta)**2*y**4/(4*d)
            occupations.append(dict(steps=count,direct=float(number),formula=float(target),error=float(abs(number-target))))
        result=dict(omega=w,h=h,y=y,
            map_error=float(np.max(abs(direct-rotation))),
            vacuum_invariance_error=float(np.max(abs(direct@C@direct.T-C))),
            determinant_error=float(abs(np.linalg.det(C)-.25)),
            logarithm_generator_error=float(np.max(abs(reconstructed-G))),
            one_step_particle_number=occupations[0]['direct'],one_step_y6=float(y**6),
            occupations=occupations)
        assert result['map_error']<1e-12
        assert result['vacuum_invariance_error']<1e-12
        assert result['logarithm_generator_error']<2e-12
        assert max(x['error'] for x in occupations)<2e-12
        assert abs(result['one_step_particle_number']-y**6)<1e-12
        rows.append(result)
    return rows

def derivative_audits():
    rows=[]
    for kappa2,a,h in [(.1,.7,.2),(3.,1.,.25),(25.,1.4,.3)]:
        v=params(kappa2,a,h);da=1e-4*a
        def fd(fn):
            return (fn(a-2*da)-8*fn(a-da)+8*fn(a+da)-fn(a+2*da))/(12*da)
        error_O=abs(fd(lambda x:params(kappa2,x,h)['Omega'])-v['Oa'])
        error_b=abs(fd(lambda x:np.log(params(kappa2,x,h)['r']))-v['b'])
        error_ba=abs(fd(lambda x:params(kappa2,x,h)['b'])-v['ba'])
        q=.73;p=-.91
        def Ephys(scale):
            z=params(kappa2,scale,h)
            return z['Omega']*(z['r']*q*q+p*p/z['r'])/(2*scale)
        Heta=v['Omega']*(v['r']*q*q+p*p/v['r'])/2
        derivative=(v['Oa']*(v['r']*q*q+p*p/v['r'])
                   +v['Omega']*v['b']*(v['r']*q*q-p*p/v['r']))/2
        pressure=(Heta-a*derivative)/a**2
        difference=abs(-fd(Ephys)-pressure)
        row=dict(kappa2=kappa2,a=a,h=h,Omega_derivative_error=float(error_O),
                 log_r_derivative_error=float(error_b),b_derivative_error=float(error_ba),
                 metric_pressure_error=float(difference))
        assert max(error_O,error_b,error_ba,difference)<2e-10,row
        rows.append(row)
    return rows

def ward_audits():
    rows=[]
    for kind in ['tanh','bump']:
        for kappa2 in [0.,.1,2.,30.]:
            for eta in [-.7,0.,.3]:
                h=.2;dt=3e-5;a,ap,app=geometry(eta,kind)
                derivative=(terms(kappa2,eta-2*dt,h,kind)[0]-8*terms(kappa2,eta-dt,h,kind)[0]
                            +8*terms(kappa2,eta+dt,h,kind)[0]-terms(kappa2,eta+2*dt,h,kind)[0])/(12*dt)
                E,P,_,_=terms(kappa2,eta,h,kind);expected=ap/a*(E-P)
                err=abs(derivative-expected)
                assert err<2e-9,(kind,kappa2,eta,err)
                rows.append(dict(kind=kind,kappa2=kappa2,eta=eta,h=h,absolute_error=float(err)))
    return rows

def continuum_time_limit():
    rows=[]
    # kappa is kept fixed; this is h->0, independent of spatial refinement.
    for kind in ['tanh','bump']:
        for kappa2,eta in [(0.2,-.4),(1.,0.),(5.,.3)]:
            a,ap,app=geometry(eta,kind);w=np.sqrt(kappa2+MU*MU*a*a)
            wp=MU*MU*a*ap/w
            wpp=MU*MU*(ap*ap+a*app)/w-MU**4*a*a*ap*ap/w**3
            q2=wpp/(8*w**4)-3*wp*wp/(16*w**5)
            E=wp*wp/(16*w**3);P=E-MU*MU*a*a*q2
            errors=[]
            for h in [.08,.04,.02,.01]:
                e,p,_,_=terms(kappa2,eta,h,kind)
                errors.append(dict(h=h,E2_error=float(abs(e-E)),P2_error=float(abs(p-P))))
            # Both E and P can vanish specially; use their vector norm.
            norms=[np.hypot(z['E2_error'],z['P2_error']) for z in errors]
            ratios=[norms[i]/norms[i+1] for i in range(3)]
            assert all(3.9<r<4.1 for r in ratios),(kind,kappa2,eta,ratios)
            rows.append(dict(kind=kind,kappa2=kappa2,eta=eta,errors=errors,refinement_ratios=ratios))
    return rows

def first_step_density():
    rows=[]
    for nu in [.25,.5,.75]:
        for epsilon in [.4,.2,.1]:
            # Massless nearest-neighbor lattice, a=1: g=omega*epsilon.
            # The initial state is the continuous-time Hamiltonian vacuum.
            energy=quad(lambda x:2*np.sin(x/2)*(nu*np.sin(x/2))**6,
                        0,np.pi,epsabs=1e-14,epsrel=1e-12)[0]/(np.pi*epsilon**2)
            density=quad(lambda x:(nu*np.sin(x/2))**6,
                         0,np.pi,epsabs=1e-14,epsrel=1e-12)[0]/(np.pi*epsilon)
            target=64*nu**6/(35*np.pi*epsilon**2)
            target_n=5*nu**6/(16*epsilon)
            assert abs(energy-target)<2e-13
            assert abs(density-target_n)<2e-13
            rows.append(dict(nu=nu,epsilon=epsilon,h=nu*epsilon,energy_density=energy,
                analytic_energy_density=target,particle_density=density,
                analytic_particle_density=target_n))
    return rows

def main():
    out=dict(static=static_audits(),parameter_and_metric_derivatives=derivative_audits(),
        second_order_ward=ward_audits(),continuous_time_limit=continuum_time_limit(),
        massless_nearest_first_step=first_step_density(),all_checks_passed=True,
        limitations=[
          'The invariant Gaussian vacuum selects the positive principal quadratic generator; it is not a global minimum of many-body quasienergy modulo 2pi/h.',
          'The principal log generator is generally spatially nonlocal despite the local factorization of the step.',
          'The second-order subtraction is derived for continuous shadow-generator evolution and is an explicitly chosen extension when applied to frozen midpoint steps.',
          'The exact static false-particle formulas are occupation relative to the continuous-time vacuum, not actual cosmological production.'])
    (ROOT/'audit_shadow.json').write_text(json.dumps(out,indent=2)+'\n')
    summary=dict(all_checks_passed=True,
        max_map_error=max(r['map_error'] for r in out['static']),
        max_invariant_vacuum_error=max(r['vacuum_invariance_error'] for r in out['static']),
        max_logarithm_error=max(r['logarithm_generator_error'] for r in out['static']),
        max_metric_pressure_error=max(r['metric_pressure_error'] for r in out['parameter_and_metric_derivatives']),
        max_second_order_ward_error=max(r['absolute_error'] for r in out['second_order_ward']),
        output=str((ROOT/'audit_shadow.json').resolve()))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
