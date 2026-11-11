"""Independent canonical-ADM / lattice / reversibility audits for pass 17.

Only import public main step/reference to compare with separately implemented
canonical (a,p_a) equations and a real-space lattice matrix. No code edits.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from backreaction import Model,initial,step,reference,observables
from gravity_formula_audit import check as formula_check
ROOT=Path(__file__).resolve().parent

def full_spectrum(model):
    n=model.n;L=np.zeros((n,n))
    for i in range(n):
        L[i,i]+=2/model.spacing**2
        L[i,(i+1)%n]-=1/model.spacing**2
        L[i,(i-1)%n]-=1/model.spacing**2
    I=np.eye(n)
    K=np.kron(np.kron(L,I),I)+np.kron(np.kron(I,L),I)+np.kron(np.kron(I,I),L)
    return np.linalg.eigvalsh(K)

def canonical_reference(model,end=2.,tol=1e-13):
    kval=full_spectrum(model);n=len(kval);w=np.sqrt(kval+model.mass**2)
    # Rotate a pure squeezed covariance in dimensionless oscillator coordinates.
    R=np.array([[np.cos(model.angle),np.sin(model.angle)],[-np.sin(model.angle),np.cos(model.angle)]])
    # Main convention has positive cross covariance for positive angle, so use R^T.
    D=np.diag([np.exp(2*model.squeeze),np.exp(-2*model.squeeze)])/2
    G=R.T@D@R
    qq=G[0,0]/w;qp=np.full(n,G[0,1]);pp=G[1,1]*w
    V=model.volume;M=model.planck2;m=model.mass;Lam=model.vacuum_density
    E0=np.sum(pp+(kval+m*m)*qq)/2+V*Lam
    pa=-np.sqrt(12*M*V*E0)
    z0=np.r_[1.,pa,qq,qp,pp]
    def rhs(t,z):
        a,pa=z[:2];qq=z[2:2+n];qp=z[2+n:2+2*n];pp=z[2+2*n:]
        W=a*kval+a**3*m*m
        adot=-pa/(6*M*V*a)
        padot=-pa*pa/(12*M*V*a*a)+1.5*np.sum(pp)/a**4-.5*np.sum((kval+3*a*a*m*m)*qq)-3*V*Lam*a*a
        return np.r_[adot,padot,2*qp/a**3,pp/a**3-W*qq,-2*W*qp]
    sol=solve_ivp(rhs,(0,end),z0,method='DOP853',rtol=tol,atol=tol/50,dense_output=True)
    assert sol.success
    return sol,kval

def main():
    result={'formula_audit':formula_check(),'canonical_references':[],'reversibility':[],'pure_lambda':[],'wrong_feedback_exact':[]}
    for model in [Model(),Model(squeeze=.4),Model(vacuum_density=.05)]:
        sol,K=canonical_reference(model);target=reference(model,end=2.,tol=2e-13)
        tt=np.linspace(0,2.,401);z=sol.sol(tt);u=target.sol(tt)
        H=-z[1]/(6*model.planck2*model.volume*z[0]**2)
        kval,d=model.spectrum();expected=np.repeat(kval,d)
        x,p,qq,qp,pp=initial(model,0)
        initial_constraint=observables(x,p,qq,qp,pp,model,0,0)[5]
        result['canonical_references'].append(dict(squeeze=model.squeeze,Lambda=model.vacuum_density,
            spectrum_max_error=float(np.max(abs(K-expected))),real_mode_count=len(K),sum_degeneracy=int(d.sum()),
            a_sup_error=float(np.max(abs(z[0]-u[0]))),H_sup_error=float(np.max(abs(H-u[1]))),
            initial_relative_constraint=float(abs(initial_constraint)),initial_purity_error=float(np.max(abs(qq*pp-qp*qp-.25)))))
        for s in [0,1]:
            st=(*initial(model,s),0.)
            final=step(step(st,.037,model,s),-.037,model,s)
            err=max(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))) for a,b in zip(st,final))
            result['reversibility'].append(dict(squeeze=model.squeeze,Lambda=model.vacuum_density,lapse_power=s,max_absolute_error=err))
    model=Model();st0=(*initial(model,0),0.)
    H0=observables(*st0[:5],model,0,0.)[2]
    for s in [0,1]:
        for h in [.04,.01]:
            st=(*initial(model,s),0.)
            while st[-1]<2.:
                st=step(st,h,model,s,feedback=False)
            a=st[0]**(2/(3-s));t=st[-1]
            exact=(1+1.5*H0*t)**(2/3) if s==0 else np.sqrt(1+2*H0*t)
            result['wrong_feedback_exact'].append(dict(lapse_power=s,h=h,proper_time=t,a_error=float(abs(a-exact))))
    result['wrong_feedback_clock_gap_at_t2']=float(1-np.sqrt(1+4*H0)/(1+3*H0)**(2/3))
    class EmptyModel(Model):
        def spectrum(self):return np.array([],dtype=float),np.array([],dtype=int)
    model=EmptyModel(vacuum_density=.12)
    exact_H=np.sqrt(model.vacuum_density/(3*model.planck2))
    for s in [0,1]:
        previous=None
        for h in [.04,.02,.01,.005]:
            st=(*initial(model,s),0.)
            while st[-1]<1.:
                st=step(st,h,model,s)
            b=(3-s)/2;a=st[0]**(1/b);H=-b*st[1]*a**(b-3)/(6*model.planck2*model.volume)
            err=abs(a-np.exp(exact_H*st[-1]));Herr=abs(H-exact_H)
            result['pure_lambda'].append(dict(lapse_power=s,h=h,proper_time=st[-1],a_error=err,H_error=Herr,error_ratio=None if previous is None else previous/err))
            previous=err
    result['all_checks_passed']=(result['formula_audit']['all_checks_passed']
       and all(r['a_sup_error']<1e-10 and r['H_sup_error']<1e-10 and r['spectrum_max_error']<1e-10 and r['real_mode_count']==r['sum_degeneracy']==64 and r['initial_relative_constraint']<1e-14 and r['initial_purity_error']<1e-14 for r in result['canonical_references'])
       and all(r['max_absolute_error']<1e-11 for r in result['reversibility'])
       and all(r['a_error']<1e-11 for r in result['wrong_feedback_exact'])
       and all(r['error_ratio'] is None or 3.8<r['error_ratio']<4.2 for r in result['pure_lambda']))
    (ROOT/'audit_backreaction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));assert result['all_checks_passed']
if __name__=='__main__':main()
