"""Independent continuum target: massless Schwinger fermions vs free boson.

Zero shifted electric flux is imposed at both walls. In the original staggered
charge convention this means B=-1/4, since E_tilde_i=E_i-(-1)^i/4.
The improved massless diagonal equals J sum(E_tilde_i^2) plus a scalar.
This removes a boundary displacement source without fitting any outcome.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
import lattice_scaling as s

def continuum_modes(modes,length=8.,coupling=1.,tau=.15,tol=1e-11):
    lo=-10*tau;hi=10*tau;k=np.arange(1,modes+1)*np.pi/length
    mu=coupling/np.sqrt(np.pi)
    wi=np.sqrt(k*k+(mu*s.profile(lo,tau))**2)
    wf=np.sqrt(k*k+(mu*s.profile(hi,tau))**2)
    u=1/np.sqrt(2*wi);v=-1j*wi*u
    def rhs(eta,y):return np.r_[y[modes:],-(k*k+(mu*s.profile(eta,tau))**2)*y[:modes]]
    sol=solve_ivp(rhs,(lo,hi),np.r_[u,v],method='DOP853',rtol=tol,atol=tol/50,t_eval=[hi])
    assert sol.success
    u=sol.y[:modes,-1];v=sol.y[modes:,-1]
    beta=np.sqrt(wf/2)*u-1j*v/np.sqrt(2*wf)
    occupation=abs(beta)**2
    wronskian=float(np.max(abs(u*np.conj(v)-np.conj(u)*v-1j)))
    survival=float(np.exp(-.5*np.log1p(occupation).sum()))
    energy=float(wf@occupation)
    return dict(modes=modes,ground_survival=survival,final_excess_energy=energy,
        proper_excess_energy_density=energy/(length*s.profile(hi,tau)**2),
        wronskian_error=wronskian,occupations=occupation.tolist(),frequencies=wf.tolist())

def main():
    root=Path(__file__).parent;refs=[continuum_modes(n) for n in [20,40,80,160]]
    assert abs(refs[-1]['ground_survival']-refs[-2]['ground_survival'])<1e-9
    assert abs(refs[-1]['final_excess_energy']-refs[-2]['final_excess_energy'])<1e-9
    print(json.dumps(dict(continuum=refs[-1])),flush=True)
    rows=[];identity=[]
    for N in [6,8,10,12,14,16,18,20]:
        b=s.make_basis(N,boundary=-.25);a=8/(N+1);J=a/2
        tilde=b['e']-((-1.)**np.arange(N-1))/4
        diag=J*b['V']-a/8*b['Q'];balanced=J*(tilde**2).sum(axis=1)
        identity_error=float(np.ptp(diag-balanced));identity.append(dict(N=N,constant_shift=float((diag-balanced)[0]),error=identity_error))
        assert identity_error<1e-12
        z,state=s.evolve(b,length=8.,mass=0.,coupling=1.,tau=.15,improved=True)
        z['survival_error_vs_boson']=z['ground_survival']-refs[-1]['ground_survival']
        z['excess_energy_relative_error']=z['final_excess_energy']/refs[-1]['final_excess_energy']-1
        rows.append(z);print(json.dumps(z),flush=True)
        result=dict(parameters=dict(length=8.,mass=0.,coupling=1.,tau=.15,boundary=-.25,
                    physical_shifted_boundary=0.,wall_offset=1),continuum=refs,rows=rows,identity=identity,
                    complete=N==20)
        (root/'bosonization_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
