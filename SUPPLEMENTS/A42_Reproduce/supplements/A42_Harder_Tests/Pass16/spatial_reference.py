"""Three spatially local harmonic lattices on a supplied 1+1 FLRW metric.

Same vacuum and second-order modewise adiabatic subtraction for all stencils.
No continuum geometric subtraction is inserted into the lattice calculation.
The infinite lattice is evaluated in its Brillouin zone; not a finite ring.
Time evolution is adaptive continuous-time Hamiltonian evolution, with a
separate finite-ring local-gate audit in audit_operators.py.
"""
from pathlib import Path
import json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp,cumulative_simpson
from continuum_reference import geometry,simulate as continuum

ROOT=Path(__file__).resolve().parent
MASS=1/np.sqrt(np.pi)

def dispersion(k,spacing,stencil):
    x=np.asarray(k)*spacing
    if stencil=='nearest':return 4*np.sin(x/2)**2/spacing**2
    if stencil=='improved':
        s=np.sin(x/2)**2
        return 4*s*(1+s/3)/spacing**2
    if stencil=='centered':return np.sin(x)**2/spacing**2
    raise ValueError(stencil)

def nodes(spacing,mass=MASS,order=32):
    top=np.pi/spacing
    # Resolve both k=0 and the BZ edge, where the centered rule has a doubler.
    raw=np.unique(np.r_[0.,mass*np.array([.03125,.0625,.125,.25,.5,1,2]),
                           [.25,.5,1,2,4,8,16,32,64,128,256]])
    half=raw[raw<top/2]
    edges=np.unique(np.r_[half,top/2,top-half])
    z,w=leggauss(order)
    k=np.concatenate([l+(z+1)*(r-l)/2 for l,r in zip(edges[:-1],edges[1:])])
    weights=np.concatenate([w*(r-l)/(2*np.pi) for l,r in zip(edges[:-1],edges[1:])])
    return k,weights

def subtract(w2,M2,H,Hp):
    w=np.sqrt(w2)
    wp=M2*H/w
    wpp=M2*(Hp+2*H*H)/w-M2*M2*H*H/w**3
    q2=wpp/(8*w**4)-3*wp*wp/(16*w**5)
    e2=wp*wp/(16*w**3)
    p2=e2-M2*q2
    return e2,p2,q2

def simulate(spacing=.2,stencil='nearest',kind='tanh',mass=MASS,tau=.6,
             samples=2401,order=32,tol=2e-11,start=12.):
    begun=time.time();k,weights=nodes(spacing,mass,order);K=dispersion(k,spacing,stencil);N=len(k)
    t=np.linspace(-start*tau,start*tau,samples)
    def rhs(eta,y):
        a,ap,H,Hp=geometry(eta,tau,kind);w=np.sqrt(K+(mass*a)**2)
        f=mass**2*a*ap/(2*w*w);n,re,im=y.reshape(3,N)
        return np.r_[2*f*re,2*w*im+f*(2*n+1),-2*w*re]
    sol=solve_ivp(rhs,(t[0],t[-1]),np.zeros(3*N),method='DOP853',
                  t_eval=t,rtol=tol,atol=tol/100)
    assert sol.success
    n,re,im=sol.y.reshape(3,N,len(t));a,ap,H,Hp=geometry(t,tau,kind)
    M2=(mass*a)**2;w2=K[:,None]+M2;w=np.sqrt(w2)
    rho0=(weights@(w*n))/a**2
    p0=(weights@((K[:,None]*n-M2*re)/w))/a**2
    q0=weights@((n+re)/w)
    e2,pressure2,q2=subtract(w2,M2,H,Hp)
    rho2=(weights@e2)/a**2;p2=(weights@pressure2)/a**2
    rho=rho0-rho2;p=p0-p2
    work=a*rho-a[0]*rho[0]+cumulative_simpson(ap*p,x=t,initial=0)
    compensated=p-rho+mass**2*q0
    target=Hp/(12*np.pi*a*a)
    d=dict(eta=t,a=a,rho=rho,p=p,rho0=rho0,p0=p0,q0=q0,
           rho2=rho2,p2=p2,compensated_trace=compensated,anomaly_target=target,work=work)
    j=len(t)//2
    stats=dict(spacing=spacing,stencil=stencil,kind=kind,mass=float(mass),tau=tau,
      samples=samples,order=order,nodes=N,tol=tol,start=start,
      max_abs_work=float(np.max(abs(work))),
      max_covariance_error=float(np.max(abs((n+.5)**2-re**2-im**2-.25))),
      center={key:float(v[j]) for key,v in d.items()},
      center_anomaly_ratio=float(compensated[j]/target[j]),
      trace_sup_error=float(np.max(abs(compensated-target))),
      trace_relative_sup_error=float(np.max(abs(compensated-target))/np.max(abs(target))),
      final_excitation_density=float(rho0[-1]),nfev=sol.nfev,seconds=time.time()-begun)
    return stats,d

def main():
    references={}
    for kind in ['tanh','bump']:
        stats,d=continuum(kind=kind,cutoff=128.,order=48,samples=2401,tol=2e-12)
        np.savez_compressed(ROOT/f'continuum_{kind}.npz',**d)
        references[kind]=stats
    rows=[]
    for kind in ['tanh','bump']:
        ref=np.load(ROOT/f'continuum_{kind}.npz')
        for stencil in ['nearest','improved','centered']:
            for spacing in [.8,.4,.2,.1,.05,.025]:
                stats,d=simulate(spacing,stencil,kind)
                for name,key in [('rho','rho_ren'),('p','p_ren')]:
                    error=float(np.max(abs(d[name]-ref[key])))
                    stats[name+'_sup_error']=error
                    stats[name+'_relative_sup_error']=error/float(np.max(abs(ref[key])))
                name=f'{kind}_{stencil}_{spacing:g}'.replace('.','p')
                np.savez_compressed(ROOT/f'{name}.npz',**d)
                stats['history_file']=name+'.npz';rows.append(stats)
                print(json.dumps(stats),flush=True)
                (ROOT/'lattice_results.json').write_text(json.dumps(dict(references=references,rows=rows,
                    complete=len(rows)==36),indent=2)+'\n')

if __name__=='__main__':main()
