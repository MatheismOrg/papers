"""Independent real-space operator and finite Gaussian-step audits for pass 15.

The finite periodic ring here checks the operators and their step implementation.
It is not the infinite-chain continuum stress benchmark, and no particle creation
is inferred from finite-step error. Only covariance/map convergence is tested.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from lattice_stress import dispersion,subtract,geometry

ROOT=Path(__file__).resolve().parent
STENCILS=('nearest','improved','centered')
MASS=1/np.sqrt(np.pi)

def matrix(N,epsilon,stencil):
    coefficients={'nearest':[(1,1.)],
                  'improved':[(1,4/3),(2,-1/12)],
                  'centered':[(2,.25)]}[stencil]
    K=np.zeros((N,N))
    for j in range(N):
        for shift,weight in coefficients:
            v=np.zeros(N);v[j]=-1;v[(j+shift)%N]=1
            K+=weight*np.outer(v,v)/epsilon**2
    return K

def matrix_checks():
    rows=[]
    for N in (16,32):
        epsilon=.2;k=2*np.pi*np.fft.fftfreq(N,d=epsilon)
        for stencil in STENCILS:
            K=matrix(N,epsilon,stencil);values=np.linalg.eigvalsh(K)
            expected=np.sort(dispersion(k,epsilon,stencil))
            zero_count=int(np.sum(abs(values)<1e-9))
            out=dict(N=N,spacing=epsilon,stencil=stencil,
                eigenspectrum_error=float(np.max(abs(values-expected))),
                min_eigenvalue=float(values.min()),zero_count=zero_count,
                symmetry_error=float(np.max(abs(K-K.T))))
            assert out['eigenspectrum_error']<2e-11,out
            assert values.min()>-2e-11,out
            assert zero_count==(2 if stencil=='centered' else 1),out
            if stencil=='centered':
                evens=np.arange(0,N,2);odds=np.arange(1,N,2)
                small=matrix(N//2,2*epsilon,'nearest')
                out['parity_offblock_error']=float(np.max(abs(K[np.ix_(evens,odds)])))
                out['parity_block_error']=float(max(np.max(abs(K[np.ix_(evens,evens)]-small)),
                                                   np.max(abs(K[np.ix_(odds,odds)]-small))))
                assert out['parity_offblock_error']==0 and out['parity_block_error']==0
            rows.append(out)
    return rows

def pressure_checks():
    rng=np.random.default_rng(4215);N=16;epsilon=.2;L=N*epsilon
    q=rng.normal(size=N);p=rng.normal(size=N);rows=[]
    for stencil in STENCILS:
        K=matrix(N,epsilon,stencil)
        for a in (.7,1.,1.4):
            def physical_energy(scale):
                return (p@p+q@K@q+(MASS*scale)**2*(q@q))/(2*scale)
            h=a*1e-4
            derivative=(physical_energy(a-2*h)-8*physical_energy(a-h)
                       +8*physical_energy(a+h)-physical_energy(a+2*h))/(12*h)
            from_variation=-derivative/L
            from_operator=(p@p+q@K@q-(MASS*a)**2*(q@q))/(2*a*a*L)
            rel=abs(from_variation-from_operator)/max(1.,abs(from_operator))
            assert rel<1e-10
            rows.append(dict(stencil=stencil,a=a,operator_pressure=float(from_operator),
               metric_variation_pressure=float(from_variation),relative_error=float(rel)))
    return rows

def wkb_ward_checks():
    rows=[];epsilon=.2
    for stencil in STENCILS:
        for eta in (-.7,0.,.3):
            for k in (0.,.1,1.,5.):
                K=dispersion(k,epsilon,stencil)
                def terms(t):
                    a,ap,H,Hp=geometry(t,.6,'tanh')
                    e,p,q=subtract(K+(MASS*a)**2,(MASS*a)**2,H,Hp)
                    return e,p,H
                h=3e-5
                derivative=(terms(eta-2*h)[0]-8*terms(eta-h)[0]+8*terms(eta+h)[0]-terms(eta+2*h)[0])/(12*h)
                e,p,H=terms(eta);expected=H*(e-p)
                error=abs(derivative-expected)
                assert error<2e-9
                rows.append(dict(stencil=stencil,eta=eta,k=k,
                    derivative=float(derivative),ward_rhs=float(expected),absolute_error=float(error)))
    return rows

def doubled_history_checks():
    rows=[]
    for kind in ('tanh','bump'):
        for epsilon in (.4,.2,.1,.05,.025):
            path_c=ROOT/f'{kind}_centered_{epsilon:g}.npz'.replace('.','p').replace('pnpz','.npz')
            path_n=ROOT/f'{kind}_nearest_{2*epsilon:g}.npz'.replace('.','p').replace('pnpz','.npz')
            if not(path_c.exists() and path_n.exists()):
                raise FileNotFoundError(f'Required completed history: {path_c} and {path_n}')
            c=np.load(path_c);n=np.load(path_n)
            assert np.max(abs(c['eta']-n['eta']))<1e-12
            errors={key:float(np.max(abs(c[key]-2*n[key])))
                    for key in ('rho','p','rho0','p0','rho2','p2','compensated_trace')}
            assert max(errors.values())<2e-9,errors
            rows.append(dict(kind=kind,centered_spacing=epsilon,nearest_spacing=2*epsilon,errors=errors))
    return rows

def gaussian_steps():
    N=16;epsilon=.2;t0=-1.2;tf=1.2;tau=.6
    eye=np.eye(N);zero=np.zeros((N,N));J=np.block([[zero,eye],[-eye,zero]])
    rows=[]
    for stencil in STENCILS:
        K=matrix(N,epsilon,stencil)
        def V(t):
            a=1+.5*np.tanh(t/tau)
            return K+MASS**2*a*a*eye
        evals,U=np.linalg.eigh(V(t0));freq=np.sqrt(evals)
        Q=(U*(.5/freq))@U.T;P=(U*(.5*freq))@U.T
        C0=np.block([[Q,zero],[zero,P]])
        def rhs(t,y):
            S=y.reshape(2*N,2*N)
            return np.vstack((S[N:],-V(t)@S[:N])).ravel()
        sol=solve_ivp(rhs,(t0,tf),np.eye(2*N).ravel(),method='DOP853',
                      rtol=2e-12,atol=2e-14,t_eval=[tf])
        assert sol.success
        reference=sol.y[:,-1].reshape(2*N,2*N)
        Cref=reference@C0@reference.T
        steps=[]
        for count in (128,256,512,1024):
            h=(tf-t0)/count;S=np.eye(2*N)
            for j in range(count):
                potential=V(t0+(j+.5)*h)
                # All q-only mass and bond potentials commute, so this kick
                # is exactly the product of finite-range Gaussian bond gates.
                momentum=S[N:]-.5*h*potential@S[:N]
                position=S[:N]+h*momentum
                momentum-=.5*h*potential@position
                S=np.vstack((position,momentum))
            C=S@C0@S.T
            err=float(np.linalg.norm(C-Cref)/np.linalg.norm(Cref))
            symplectic=float(np.max(abs(S@J@S.T-J)))
            pure=float(np.max(abs(C@J@C-.25*J)))
            steps.append(dict(steps=count,dt=h,covariance_relative_frobenius_error=err,
                map_max_absolute_error=float(np.max(abs(S-reference))),
                symplectic_error=symplectic,pure_covariance_error=pure))
            assert symplectic<2e-11
            assert pure<2e-10
        ratios=[steps[i]['covariance_relative_frobenius_error']/steps[i+1]['covariance_relative_frobenius_error'] for i in range(3)]
        assert all(3.7<r<4.3 for r in ratios),ratios
        rows.append(dict(stencil=stencil,N=N,spacing=epsilon,comoving_length=N*epsilon,
            mass=float(MASS),t_initial=t0,t_final=tf,tau=tau,
            initial_state='finite-ring instantaneous Gaussian ground state at t_initial',
            reference='independent adaptive real-space symplectic map',
            reference_symplectic_error=float(np.max(abs(reference@J@reference.T-J))),
            runs=steps,covariance_error_refinement_ratios=ratios))
    return rows

def main():
    out=dict(real_space_spectrum=matrix_checks(),metric_pressure=pressure_checks(),
        wkb_mode_ward=wkb_ward_checks(),doubled_histories=doubled_history_checks(),
        finite_gaussian_steps=gaussian_steps(),all_checks_passed=True,
        limitations=[
          'The periodic ring is an independent operator/gate audit, not the infinite-chain stress calculation.',
          'Finite-step covariance errors are numerical errors; they are not physical particle production.',
          'The full spatial sweeps use continuous-time adaptive evolution; these local steps were audited only on the stated finite rings.',
          'Gaussian oscillator Hilbert spaces are untruncated; this does not construct a finite-alphabet fundamental TM.'])
    path=ROOT/'audit_operators.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),all_checks_passed=True,
       max_spectrum_error=max(r['eigenspectrum_error'] for r in out['real_space_spectrum']),
       max_metric_pressure_relative_error=max(r['relative_error'] for r in out['metric_pressure']),
       max_wkb_ward_error=max(r['absolute_error'] for r in out['wkb_mode_ward']),
       max_doubling_error=max(max(r['errors'].values()) for r in out['doubled_histories']),
       gaussian_steps=out['finite_gaussian_steps']),indent=2))

if __name__=='__main__':main()
