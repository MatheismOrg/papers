"""Independent direct-coordinate and metric-variation audits of pass 16.

The reference implementation evolves Bogoliubov amplitudes in a moving basis.
This audit instead propagates the normalized complex q,p mode coordinates by
the three local shears, and differentiates the prescribed shadow energy at
fixed q,p to obtain pressure. It also measures the frozen invariant-vacuum
residual numerically, and distinguishes symplecticity from stability.
"""
from pathlib import Path
import json
import numpy as np
import step_scalar as s

ROOT=Path(__file__).resolve().parent

def independent_parameters(K,a,h,mu):
    omega=np.sqrt(K+mu*mu*a*a)
    margin=1-h*h*omega*omega/4
    if np.min(margin)<=0:
        raise ValueError('The frozen oscillator has no stable Gaussian vacuum.')
    r=omega*np.sqrt(margin)
    Omega=2*np.arcsin(h*omega/2)/h
    return omega,r,Omega

def direct_coordinate_checks():
    rows=[];eps=.4;nu=.5;mu=s.MASS;tau=.6;order=24
    for stencil in ['nearest','improved']:
        for kind in ['tanh','bump']:
            stats,history=s.simulate(eps,nu,stencil,kind,order=order)
            h=stats['h'];t=history['eta']
            k,weights=s.nodes(eps,mu,order);K=s.dispersion(k,eps,stencil)
            a=s.geometry(t[0],tau,kind)[0]
            _,r,_=independent_parameters(K,a,h,mu)
            q=1/np.sqrt(2*r)+0j;p=-1j*np.sqrt(r/2)
            for tn in t[:-1]:
                amid=s.geometry(tn+h/2,tau,kind)[0]
                potential=K+mu*mu*amid*amid
                p-=h*potential*q/2
                q+=h*p
                p-=h*potential*q/2
            a=s.geometry(t[-1],tau,kind)[0]
            _,r,Omega=independent_parameters(K,a,h,mu)
            B=np.sqrt(r/2)*q-1j*p/np.sqrt(2*r)
            rho=float(weights@(Omega*abs(B)**2)/a**2)

            # Hold the evolved q,p fixed while varying the external scale.
            # Subtracting the frozen vacuum energy is compatible with this
            # variation by the Hellmann-Feynman identity for that Gaussian.
            def physical_excess_energy(scale):
                _,rr,OO=independent_parameters(K,scale,h,mu)
                return float(weights@(OO*(abs(p)**2/rr+rr*abs(q)**2-1)/(2*scale)))
            da=a*1e-4
            pressure=-(physical_excess_energy(a-2*da)-8*physical_excess_energy(a-da)
                       +8*physical_excess_energy(a+da)-physical_excess_energy(a+2*da))/(12*da)
            rho_error=abs(rho-history['rho0'][-1])
            pressure_error=abs(pressure-history['p0'][-1])
            wronskian=float(np.max(abs(q*np.conj(p)-np.conj(q)*p-1j)))
            assert rho_error<2e-12,(stencil,kind,rho_error)
            assert pressure_error<2e-10,(stencil,kind,pressure_error)
            assert wronskian<2e-12,(stencil,kind,wronskian)
            rows.append(dict(stencil=stencil,kind=kind,spacing=eps,nu=nu,h=h,
                nodes=len(k),steps=stats['steps'],mass=float(mu),
                initial_eta=float(t[0]),final_eta=float(t[-1]),
                rho_direct=rho,rho_reference=float(history['rho0'][-1]),
                rho_absolute_error=float(rho_error),
                pressure_metric_difference=pressure,pressure_reference=float(history['p0'][-1]),
                pressure_absolute_error=float(pressure_error),
                wronskian_error=wronskian,metric_difference_increment=da))
    return rows

def frozen_vacuum_checks():
    rows=[];eps=.4;nu=.5;h=eps*nu;mu=s.MASS
    for stencil in ['nearest','improved']:
        k,weights=s.nodes(eps,mu,24);K=s.dispersion(k,eps,stencil)
        omega,r,Omega=independent_parameters(K,1.,h,mu)
        q=1/np.sqrt(2*r)+0j;p=-1j*np.sqrt(r/2)
        for _ in range(100):
            p-=h*omega*omega*q/2
            q+=h*p
            p-=h*omega*omega*q/2
        B=np.sqrt(r/2)*q-1j*p/np.sqrt(2*r)
        occupation=float(np.max(abs(B)**2))
        covariance_error=float(max(np.max(abs(abs(q)**2-1/(2*r))),
                                   np.max(abs(abs(p)**2-r/2)),
                                   np.max(abs((q*np.conj(p)).real))))
        assert occupation<2e-25,occupation
        assert covariance_error<2e-12,covariance_error
        rows.append(dict(stencil=stencil,spacing=eps,nu=nu,h=h,steps=100,
            max_invariant_occupation=occupation,
            excitation_density=float(weights@(Omega*abs(B)**2)),
            max_covariance_error=covariance_error))
    return rows

def unstable_symplectic_control():
    # A canonical map need not have a stable invariant Gaussian vacuum.
    omega=1.;h=2.1
    d=1-h*h*omega*omega/4
    S=np.array([[1-h*h*omega*omega/2,h],
                [-h*omega*omega*d,1-h*h*omega*omega/2]])
    J=np.array([[0.,1.],[-1.,0.]])
    symplectic_error=float(np.max(abs(S@J@S.T-J)))
    spectral_radius=float(np.max(abs(np.linalg.eigvals(S))))
    rejected=False
    try:s.frozen(np.array([0.]),1.,h,mass=omega)
    except ValueError:rejected=True
    assert symplectic_error<1e-12 and spectral_radius>1. and rejected
    return dict(omega=omega,h=h,h_omega=h*omega,stability_margin=d,
        symplectic_error=symplectic_error,determinant=float(np.linalg.det(S)),
        spectral_radius=spectral_radius,stable_vacuum_rejected=rejected,
        interpretation='Canonical/quantum-unitary shears can form an unstable squeezing map; stability is a distinct requirement.')

def main():
    out=dict(direct_coordinates=direct_coordinate_checks(),
        frozen_vacuum=frozen_vacuum_checks(),
        unstable_symplectic_control=unstable_symplectic_control(),
        all_checks_passed=True,
        limitations=[
          'Pressure here is the explicitly selected principal-shadow-generator metric derivative, not a proof of finite-step covariant gravitational stress.',
          'This checks the step implementation and unrenormalized excitation observables; it does not independently derive the source subtraction.',
          'The invariant vacuum is a selected pure Gaussian state, not an absolute minimum of quasienergy modulo 2*pi/h.'])
    target=ROOT/'audit_steps.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
