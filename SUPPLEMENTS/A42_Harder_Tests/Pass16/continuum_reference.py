"""Pass 14: stress-energy tests for a real scalar on an infinite spatial line.

Units hbar=c=1, metric a(eta)^2(-deta^2+dx^2), R=2 H'/a^2.
This is a supplied-background continuum benchmark, not a gravity solver.
The covariance representation avoids subtracting UV zero-point energies.
"""
from pathlib import Path
import json, argparse, time
import numpy as np
from scipy.integrate import solve_ivp, cumulative_simpson
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent

def geometry(t, tau=.6, kind='tanh'):
    z=np.tanh(np.asarray(t)/tau); s=1-z*z
    if kind=='tanh':
        a=1+.5*z; ap=.5*s/tau; app=-s*z/tau**2
    elif kind=='bump':
        a=1+.5*s; ap=-s*z/tau; app=s*(3*z*z-1)/tau**2
    else: raise ValueError(kind)
    H=ap/a; Hp=app/a-H*H
    return a,ap,H,Hp

def quadrature(cutoff=64., order=16, mass=1/np.sqrt(np.pi)):
    # Resolve IR near the mass scale as well as oscillatory UV pressure.
    edges=np.unique(np.r_[0.,mass*np.array([.0625,.125,.25,.5,1,2]),
                            np.array([.25,.5,1,2,4,8,16,32,64,128])])
    edges=edges[edges<cutoff];edges=np.r_[edges,cutoff]
    z,w=leggauss(order)
    k=np.concatenate([l+(z+1)*(r-l)/2 for l,r in zip(edges[:-1],edges[1:])])
    weights=np.concatenate([w*(r-l)/(2*np.pi) for l,r in zip(edges[:-1],edges[1:])])
    return k,weights

def simulate(mass=1/np.sqrt(np.pi), tau=.6, kind='tanh', cutoff=64., order=16,
             samples=2401, tol=2e-10, start=12.):
    begun=time.time(); k,w=quadrature(cutoff,order,mass); N=len(k)
    times=np.linspace(-start*tau,start*tau,samples)
    def rhs(t,y):
        a,ap,_,_=geometry(t,tau,kind); omega=np.sqrt(k*k+(mass*a)**2)
        f=mass*mass*a*ap/(2*omega**2)
        n,re,im=y.reshape(3,N)
        return np.r_[2*f*re,2*omega*im+f*(2*n+1),-2*omega*re]
    sol=solve_ivp(rhs,(times[0],times[-1]),np.zeros(3*N),method='DOP853',
                  t_eval=times,rtol=tol,atol=tol/100)
    assert sol.success
    n,re,im=sol.y.reshape(3,N,len(times));a,ap,H,Hp=geometry(times,tau,kind)
    M2=(mass*a)**2;omega=np.sqrt(k[:,None]**2+M2)
    rho0=np.einsum('k,kt->t',w,omega*n)/a**2
    pgas=np.einsum('k,kt->t',w,k[:,None]**2*n/omega)/a**2
    coherence=-np.einsum('k,kt->t',w,M2*re/omega)/a**2
    p0=pgas+coherence
    phi0=np.einsum('k,kt->t',w,(n+re)/omega)
    rho2=H*H/(24*np.pi*a*a);p2=(H*H-2*Hp)/(24*np.pi*a*a)
    rho=rho0-rho2;p=p0-p2;anomaly=Hp/(12*np.pi*a*a)
    def work(r,p):
        return a*r-a[0]*r[0]+cumulative_simpson(p*ap,x=times,initial=0)
    ward0=work(rho0,p0);ward=work(rho,p);wardgas=work(rho0,pgas)
    ward_repairgas=work(rho,pgas-p2)
    trace0=p0-rho0+mass*mass*phi0
    trace=p-rho+mass*mass*phi0
    determinant=(n+.5)**2-re**2-im**2-.25
    j=len(times)//2
    fields=dict(eta=times,a=a,rho_instant=rho0,p_instant=p0,p_particle_gas=pgas,
                pressure_coherence=coherence,rho_ren=rho,p_ren=p,
                phi2_instant=phi0,anomaly=anomaly,trace_ren=p-rho,
                trace_identity_lhs=trace,
                ward_instant=ward0,ward_ren=ward,ward_gas=wardgas)
    stats=dict(parameters=dict(mass=float(mass),tau=tau,kind=kind,cutoff=cutoff,
               order=order,samples=samples,tol=tol,start=start,nodes=N),
        center={key:float(value[j]) for key,value in fields.items()},
        max_abs_ward_instant=float(np.max(abs(ward0))),
        max_abs_ward_ren=float(np.max(abs(ward))),
        max_abs_ward_particle_gas=float(np.max(abs(wardgas))),
        max_abs_ward_particle_gas_with_anomaly=float(np.max(abs(ward_repairgas))),
        max_trace_identity_error=float(np.max(abs(trace-anomaly))),
        max_abs_instant_mass_compensated_trace=float(np.max(abs(trace0))),
        max_anomaly=float(np.max(abs(anomaly))),
        max_abs_coherence_pressure=float(np.max(abs(coherence))),
        max_covariance_determinant_error=float(np.max(abs(determinant))),
        min_occupation=float(np.min(n)),
        final_rho_instant=float(rho0[-1]),final_p_instant=float(p0[-1]),
        nfev=sol.nfev,seconds=time.time()-begun)
    return stats,fields

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');args=parser.parse_args()
    configs=[dict(cutoff=32.,order=12,samples=1201),
             dict(cutoff=64.,order=24,samples=2401),
             dict(cutoff=128.,order=48,samples=4801),
             dict(cutoff=128.,order=48,samples=4801,tol=2e-12),
             dict(kind='bump',cutoff=64.,order=24,samples=2401),
             dict(kind='bump',cutoff=128.,order=48,samples=4801),
             dict(mass=.1,cutoff=64.,order=24,samples=2401),
             dict(mass=.03,cutoff=64.,order=24,samples=2401),
             dict(mass=.01,cutoff=64.,order=24,samples=2401),
             dict(mass=.003,cutoff=64.,order=24,samples=2401)]
    if args.quick:configs=configs[:1]
    rows=[]
    for i,config in enumerate(configs):
        stats,fields=simulate(**config);rows.append(stats)
        np.savez_compressed(ROOT/f'history_{i:02d}.npz',**fields)
        print(json.dumps(dict(run=i,**stats)),flush=True)
        (ROOT/'stress_results.json').write_text(json.dumps(dict(rows=rows,complete=i==len(configs)-1),indent=2)+'\n')

if __name__=='__main__':main()
