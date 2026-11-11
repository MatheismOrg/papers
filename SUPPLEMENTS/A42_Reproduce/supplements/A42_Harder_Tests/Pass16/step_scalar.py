"""Pass 16: local KDK step, invariant Gaussian vacuum and shadow-source convention.

The frozen step's positive quadratic generator uses theta in (0,pi). Its
metric derivative and continuous-shadow adiabatic subtraction are explicitly
specified observables; they are not asserted to be a local covariant finite-h
stress tensor. All joint runs use h/epsilon fixed, no separate tiny-time ODE.
"""
from pathlib import Path
import json,time,argparse
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline
from continuum_reference import geometry,simulate as continuum
from spatial_reference import nodes,dispersion

ROOT=Path(__file__).resolve().parent;MASS=1/np.sqrt(np.pi)

def frozen(K,a,h,mass=MASS):
    w2=K+(mass*a)**2;w=np.sqrt(w2);d=1-h*h*w2/4
    if np.min(d)<=0:raise ValueError('No stable invariant Gaussian vacuum: h*omega >= 2')
    r=w*np.sqrt(d);theta=2*np.arcsin(h*w/2);Omega=theta/h
    b=mass**2*a*(1/w2-h*h/(4*d))
    ba=mass**2*(1/w2-h*h/(4*d))-2*mass**4*a*a/w2**2-h**4*mass**4*a*a/(8*d*d)
    Omega_a=mass**2*a/r
    return r,theta,Omega,b,ba,Omega_a,d

def source(A,B,K,weights,eta,h,kind,mass=MASS,tau=.6):
    a,ap,H,Hp=geometry(eta,tau,kind);app=a*(Hp+H*H)
    r,theta,O,b,ba,Oa,d=frozen(K,a,h,mass)
    n=abs(B)**2;cre=(A*np.conj(B)).real
    E=O*n;P=(O-a*Oa)*n-a*O*b*cre
    f=b*ap/2;fp=(ba*ap*ap+b*app)/2;Op=Oa*ap
    n2=f*f/(4*O*O);c2=fp/(4*O*O)-f*Op/(4*O**3)
    E2=O*n2;P2=(O-a*Oa)*n2-a*O*b*c2
    rho0=float(weights@E/a**2);p0=float(weights@P/a**2)
    rho2=float(weights@E2/a**2);p2=float(weights@P2/a**2)
    compensator=float(weights@(a*Oa*n+a*O*b*cre)/a**2)
    return dict(rho=rho0-rho2,p=p0-p2,rho0=rho0,p0=p0,rho2=rho2,p2=p2,
                mass_response=compensator,compensated_trace=rho2-p2,
                covariance_error=float(np.max(abs(abs(A)**2-abs(B)**2-1))),
                min_stability_margin=float(np.min(d)))

def change_basis(A,B,new_r,old_r):
    chi=.5*np.log(new_r/old_r);c=np.cosh(chi);s=np.sinh(chi)
    return c*A+s*B,s*A+c*B

def simulate(spacing=.2,nu=.5,stencil='nearest',kind='tanh',mass=MASS,tau=.6,order=32):
    begun=time.time();lo=-12*tau;hi=12*tau;h=nu*spacing;count=round((hi-lo)/h)
    assert abs(count*h-(hi-lo))<1e-11
    t=np.linspace(lo,hi,count+1);h=(hi-lo)/count
    k,weights=nodes(spacing,mass,order);K=dispersion(k,spacing,stencil)
    a=geometry(lo,tau,kind)[0];r=frozen(K,a,h,mass)[0]
    A=np.ones(len(k),complex);B=np.zeros(len(k),complex)
    histories=[source(A,B,K,weights,lo,h,kind,mass,tau)]
    for j in range(count):
        amid=geometry(t[j]+h/2,tau,kind)[0]
        rmid,theta,*_=frozen(K,amid,h,mass)
        A,B=change_basis(A,B,rmid,r)
        A*=np.exp(-1j*theta);B*=np.exp(1j*theta)
        anew=geometry(t[j+1],tau,kind)[0];rnew=frozen(K,anew,h,mass)[0]
        A,B=change_basis(A,B,rnew,rmid);r=rnew
        histories.append(source(A,B,K,weights,t[j+1],h,kind,mass,tau))
    d={key:np.array([x[key] for x in histories]) for key in histories[0]}
    a,ap,H,Hp=geometry(t,tau,kind)
    work=a*d['rho']-a[0]*d['rho'][0]+cumulative_simpson(ap*d['p'],x=t,initial=0)
    d.update(eta=t,a=a,work=work)
    stats=dict(spacing=spacing,nu=nu,h=h,steps=count,stencil=stencil,kind=kind,mass=float(mass),tau=tau,
      order=order,nodes=len(k),max_abs_work=float(np.max(abs(work))),
      max_covariance_error=float(np.max(d['covariance_error'])),
      min_stability_margin=float(np.min(d['min_stability_margin'])),
      final_excitation_density=float(d['rho0'][-1]),
      center={key:float(value[len(t)//2]) for key,value in d.items()},seconds=time.time()-begun)
    return stats,d

def static_controls():
    rows=[]
    for stencil in ['nearest','improved']:
        for spacing in [.4,.2,.1,.05,.025]:
            for nu in [.25,.5,.75]:
                h=nu*spacing;k,w=nodes(spacing,order=64);K=dispersion(k,spacing,stencil)
                omega=np.sqrt(K+MASS*MASS);r,theta,Omega,*_=frozen(K,1.,h)
                y=h*omega/2;count=round(1.2/h)
                n1=y**6
                nT=y**4*np.sin(count*theta)**2/(4*(1-y*y))
                d=r/omega
                invariant_mismatch=(d+1/d-2)/4
                u=1/np.sqrt(2*r);v=-1j*r*u
                c=1-h*h*omega*omega/2
                unew=c*u+h*v;vnew=-h*r*r*u+c*v
                matched_beta=np.sqrt(r/2)*unew-1j*vnew/np.sqrt(2*r)
                rows.append(dict(stencil=stencil,spacing=spacing,nu=nu,
                    one_step_apparent_energy=float(w@(omega*n1)),
                    fixed_time=1.2,steps=count,fixed_time_apparent_energy=float(w@(omega*nT)),
                    invariant_vacuum_energy_in_continuum_basis=float(w@(omega*invariant_mismatch)),
                    matched_invariant_max_occupation=float(np.max(abs(matched_beta)**2)),
                    matched_invariant_excitation_density=float(w@(Omega*abs(matched_beta)**2)),
                    massless_nearest_one_step_scaled_energy_limit=64*nu**6/(35*np.pi)))
    return rows

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    kinds=['tanh'] if args.pilot else ['tanh','bump']
    refs={};interp={}
    for kind in kinds:
        stats,d=continuum(kind=kind,cutoff=128.,order=48,samples=5761,tol=2e-12)
        np.savez_compressed(ROOT/f'continuum_{kind}.npz',**d);refs[kind]=stats
        interp[kind]={key:CubicSpline(d['eta'],d[key]) for key in ['rho_ren','p_ren','anomaly']}
    rows=[]
    configs=[(kind,stencil,nu,eps) for kind in kinds for stencil in ['nearest','improved']
             for nu in [.25,.5,.75] for eps in [.4,.2,.1,.05,.025]]
    if args.pilot:configs=[('tanh','nearest',.5,eps) for eps in [.4,.2,.1,.05,.025]]
    for kind,stencil,nu,eps in configs:
        stats,d=simulate(eps,nu,stencil,kind)
        for key,refkey in [('rho','rho_ren'),('p','p_ren'),('compensated_trace','anomaly')]:
            target=interp[kind][refkey](d['eta']);err=float(np.max(abs(d[key]-target)))
            stats[key+'_sup_error']=err;stats[key+'_relative_sup_error']=err/float(np.max(abs(target)))
        name=f'{kind}_{stencil}_nu{nu:g}_eps{eps:g}'.replace('.','p')+'.npz'
        np.savez_compressed(ROOT/name,**d);stats['history_file']=name;rows.append(stats)
        print(json.dumps(stats),flush=True)
        (ROOT/'step_results.json').write_text(json.dumps(dict(references=refs,rows=rows,
                complete=len(rows)==len(configs),pilot=args.pilot),indent=2)+'\n')
    (ROOT/'static_results.json').write_text(json.dumps(dict(rows=static_controls()),indent=2)+'\n')

if __name__=='__main__':main()
