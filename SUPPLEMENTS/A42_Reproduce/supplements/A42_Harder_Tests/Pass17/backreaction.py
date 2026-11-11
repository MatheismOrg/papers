"""Pass 17: fixed-cutoff 3+1 Gaussian scalar / homogeneous GR benchmark.

Units hbar=c=M_P=1 unless configured. Infinite oscillator Hilbert spaces are
represented by Gaussian covariances, not a finite-bit fundamental automaton.
Spatial nearest-neighbor cubic torus; fixed regulator, no UV-renormalization claim.
A symmetric composition of exact gravity, kinetic, potential and clock flows.
No Friedmann projection after initialization. All canonical lapse derivatives
are included through the gauge-specific Hamiltonian a**s * C.
"""
from pathlib import Path
from dataclasses import dataclass,asdict
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent

@dataclass
class Model:
    n: int=4
    spacing: float=1.
    mass: float=.4
    planck2: float=1.
    vacuum_density: float=0.
    squeeze: float=0.
    angle: float=.6
    @property
    def volume(self):return (self.n*self.spacing)**3
    def spectrum(self):
        k=2*np.pi*np.arange(self.n)/self.n
        one=4*np.sin(k/2)**2/self.spacing**2
        full=(one[:,None,None]+one[None,:,None]+one[None,None,:]).ravel()
        vals,counts=np.unique(np.round(full,12),return_counts=True)
        return vals,counts

def stress(a,x,z,model):
    K,d=model.spectrum();V=model.volume
    kinetic=d@z/(2*V*a**6)
    gradient=d@(K*x)/(2*V*a**2)
    mass=model.mass**2*(d@x)/(2*V)
    rho=kinetic+gradient+mass+model.vacuum_density
    pressure=kinetic-gradient/3-mass-model.vacuum_density
    return rho,pressure

def initial(model,s):
    K,d=model.spectrum();w=np.sqrt(K+model.mass**2)
    c=np.cos(model.angle);ss=np.sin(model.angle);ep=np.exp(2*model.squeeze);em=1/ep
    x=(c*c*ep+ss*ss*em)/(2*w)
    z=w*(ss*ss*ep+c*c*em)/2
    y=np.full_like(w,c*ss*(ep-em)/2)
    rho,_=stress(1.,x,z,model);H=np.sqrt(rho/(3*model.planck2))
    b=(3-s)/2;P=-6*model.planck2*model.volume*H/b
    return 1.,P,x,y,z

def observables(X,P,x,y,z,model,s,t):
    b=(3-s)/2;a=X**(1/b)
    H=-b*P*a**(b-3)/(6*model.planck2*model.volume)
    rho,p=stress(a,x,z,model)
    return np.array([t,a,H,rho,p,(3*model.planck2*H*H-rho)/rho,
                     np.max(np.abs(x*z-y*y-.25))])

def step(state,h,model,s=0,feedback=True):
    X,P,x,y,z,t=state;x=x.copy();y=y.copy();z=z.copy()
    b=(3-s)/2;gamma=(s+1)/b;delta=(s+3)/b
    K,d=model.spectrum();V=model.volume
    t+=h/2*X**(s/b)
    X-=h/2*b*b*P/(6*model.planck2*V)
    if X<=0:raise ValueError('Reached singular/negative scale variable; no bounce prescription')
    # Half kinetic flow, both field shear and canonical geometry impulse.
    alpha=h/(2*X*X)
    if feedback:P+=h/2*(d@z)/X**3
    x,y=x+2*alpha*y+alpha*alpha*z,y+alpha*z
    W=K*X**gamma+model.mass**2*X**delta
    if feedback:
        P-=h/2*(d@((gamma*K*X**(gamma-1)+delta*model.mass**2*X**(delta-1))*x))
        P-=h*V*model.vacuum_density*delta*X**(delta-1)
    beta=h*W
    z,y=z-2*beta*y+beta*beta*x,y-beta*x
    alpha=h/(2*X*X)
    if feedback:P+=h/2*(d@z)/X**3
    x,y=x+2*alpha*y+alpha*alpha*z,y+alpha*z
    X-=h/2*b*b*P/(6*model.planck2*V)
    if X<=0:raise ValueError('Reached singular/negative scale variable; no bounce prescription')
    t+=h/2*X**(s/b)
    return X,P,x,y,z,t

def simulate(model,s,h,end=2.,feedback=True):
    X,P,x,y,z=initial(model,s);state=(X,P,x,y,z,0.)
    hist=[observables(*state[:5],model,s,state[5])]
    cov=[np.array([x,y,z])]
    while state[-1]<end:
        state=step(state,h,model,s,feedback)
        hist.append(observables(*state[:5],model,s,state[-1]))
        cov.append(np.array(state[2:5]))
    arr=np.array(hist)
    return dict(zip(['time','a','H','rho','pressure','constraint','purity_error'],arr.T)),np.array(cov)

def reference(model,end=2.5,tol=2e-12):
    X,P,x,y,z=initial(model,0);K,d=model.spectrum();n=len(K)
    H=observables(X,P,x,y,z,model,0,0)[2]
    init=np.r_[1.,H,x,y,z]
    def rhs(t,v):
        a,H=v[:2];x,y,z=v[2:2+n],v[2+n:2+2*n],v[2+2*n:]
        rho,p=stress(a,x,z,model);W=a*K+a**3*model.mass**2
        return np.r_[a*H,-(rho+p)/(2*model.planck2),2*y/a**3,z/a**3-W*x,-2*W*y]
    sol=solve_ivp(rhs,(0,end),init,rtol=tol,atol=tol/20,method='DOP853',dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def evaluate_ref(sol,t,model):
    v=sol.sol(t);n=len(model.spectrum()[0]);a,H=v[:2]
    x=v[2:2+n].T;y=v[2+n:2+2*n].T;z=v[2+2*n:].T
    rho,p=np.array([stress(aa,xx,zz,model) for aa,xx,zz in zip(a,x,z)]).T
    return dict(time=t,a=a,H=H,rho=rho,pressure=p),np.stack([x,y,z],axis=1)

def main():
    cases={'vacuum':Model(),'squeezed':Model(squeeze=.4),'vacuum_Lambda':Model(vacuum_density=.05)}
    rows=[];clock_rows=[];refs={};history_cache={}
    for name,model in cases.items():
        sol=reference(model);dense,cc=evaluate_ref(sol,np.linspace(0,2.,2001),model)
        np.savez_compressed(ROOT/f'reference_{name}.npz',**dense,cov=cc)
        refs[name]=dict(model=asdict(model),final={k:float(v[-1]) for k,v in dense.items()},
                       initial={k:float(v[0]) for k,v in dense.items()})
        for feedback in [True,False]:
            # Wrong-feedback control only vacuum; all steps and clocks retained.
            if not feedback and name!='vacuum':continue
            for s in [0,1]:
                for h in [.08,.04,.02,.01,.005]:
                    hist,cov=simulate(model,s,h,feedback=feedback)
                    mask=hist['time']<=2+1e-12
                    tt=hist['time'][mask];ref,cr=evaluate_ref(sol,tt,model)
                    row=dict(case=name,lapse_power=s,h=h,feedback=feedback,steps=len(hist['time'])-1)
                    for k in ['a','H','rho','pressure']:
                        row[k+'_relative_sup_error']=float(np.max(abs(hist[k][mask]-ref[k]))/np.max(abs(ref[k])))
                    row['cov_relative_sup_error']=float(np.max(abs(cov[mask]-cr))/np.max(abs(cr)))
                    row['max_constraint']=float(np.max(abs(hist['constraint'][mask])))
                    row['max_purity_error']=float(np.max(hist['purity_error'][mask]))
                    row['final_a']=float(CubicSpline(hist['time'],hist['a'])(2.))
                    row['final_H']=float(CubicSpline(hist['time'],hist['H'])(2.))
                    filename=f'{name}_s{s}_h{h:g}_feedback{int(feedback)}'.replace('.','p')+'.npz'
                    row['history_file']=filename
                    np.savez_compressed(ROOT/filename,**hist,cov=cov)
                    rows.append(row);history_cache[name,s,h,feedback]=(hist,cov)
                    print(json.dumps(row),flush=True)
        grid=np.linspace(0,2.,1001)
        for h in [.08,.04,.02,.01,.005]:
            h0,c0=history_cache[name,0,h,True];h1,c1=history_cache[name,1,h,True]
            rr=dict(case=name,h=h)
            for key in ['a','H','rho','pressure']:
                v0=CubicSpline(h0['time'],h0[key])(grid);v1=CubicSpline(h1['time'],h1[key])(grid)
                rr[key+'_clock_relative_sup_difference']=float(np.max(abs(v0-v1))/np.max(abs(v0)))
            clock_rows.append(rr)
    result=dict(description=__doc__,references=refs,rows=rows,clock_rows=clock_rows,complete=True)
    (ROOT/'backreaction_results.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
