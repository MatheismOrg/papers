"""Independent Pass 17 audit of metric responses of local KDK gates.

No imports from the main backreaction implementation. Canonical vector z=(q,p),
quadratic observables H=z^T G z/2, hbar=1. This audits an exactly specified
finite-step response, not a GR constraint-preservation theorem.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import logm, expm
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent
J=np.array([[0.,1.],[-1.,0.]])

def diff5(fun,x,delta=1e-4):
    return (fun(x-2*delta)-8*fun(x-delta)+8*fun(x+delta)-fun(x+2*delta))/(12*delta)

def step(t,v,h):
    K=np.array([[1.,0.],[-h*v/2,1.]])
    D=np.array([[1.,h*t],[0.,1.]])
    return K@D@K

def response(t,v,ta,va,h):
    e=np.array([1.,0.]); midp=np.array([-h*v/2,1.]); endq=np.array([1-h*h*t*v/2,h*t])
    return va/2*(np.outer(e,e)+np.outer(endq,endq))+ta*np.outer(midp,midp)

def coeff(a,k2=1.3,mass=.7):
    return a**-3, a*k2+a**3*mass**2

def coeff_derivative(a,k2=1.3,mass=.7):
    return -3*a**-4,k2+3*a*a*mass**2

def shadow_G(t,v,h):
    return (-J@logm(step(t,v,h))/h).real

def action(q0,q1,a,h):
    t,v=coeff(a)
    return (q1-q0)**2/(2*h*t)-h*v*(q0*q0+q1*q1)/4

def gate_audit():
    rows=[]
    for a in [.7,1.,1.4]:
        for h in [.025,.1,.3]:
            t,v=coeff(a);ta,va=coeff_derivative(a)
            S=step(t,v,h);G=response(t,v,ta,va,h)
            Sn=diff5(lambda aa:step(*coeff(aa),h),a)
            Gn=-J@np.linalg.solve(S,Sn)/h
            GF=shadow_G(t,v,h)
            GFa=diff5(lambda aa:shadow_G(*coeff(aa),h),a)
            nodes,weights=leggauss(32);Ga=np.zeros((2,2))
            for x,w in zip(nodes,weights):
                B=expm(J@GF*((x+1)*h/2))
                Ga+=w/2*(B.T@GFa@B)
            action_error=0.
            for z in [np.array([.4,.7]),np.array([-.8,.2]),np.array([.2,-.9])]:
                q1=float((S@z)[0]);r=.5*float(z@G@z)
                rn=-float(diff5(lambda aa:action(z[0],q1,aa,h),a))/h
                action_error=max(action_error,abs(r-rn))
            rows.append(dict(a=a,h=h,matrix_response_error=float(np.max(abs(G-Gn))),
                action_response_error=action_error,
                shadow_average_error=float(np.max(abs(G-Ga))),
                unaveraged_shadow_difference=float(np.max(abs(G-GFa)))))
    return rows

def lattice_audit():
    n=10;h=.5;K=2.49*np.eye(n)
    for j in range(n):
        K[j,(j+1)%n]=-1.;K[j,(j-1)%n]=-1.
    I=np.eye(n);Z=np.zeros((n,n));JJ=np.block([[Z,I],[-I,Z]])
    Q=np.zeros((n,n));Q[0,0]=1.
    def S_for(lam):
        M=K+lam*Q;C=I-h*h*M/2
        return np.block([[C,h*I],[-h*M@(I-h*h*M/4),C]])
    S=S_for(0.);e=np.eye(2*n)[0];u=S[0]
    G=.5*(np.outer(e,e)+np.outer(u,u))
    Gn=-JJ@np.linalg.solve(S,diff5(S_for,0.))/h
    GF=-JJ@logm(S).real/h
    GFa=diff5(lambda x:-JJ@logm(S_for(x)).real/h,0.)
    nodes,weights=leggauss(40);Ga=np.zeros_like(G)
    for x,w in zip(nodes,weights):
        B=expm(JJ@GF*((x+1)*h/2));Ga+=w/2*(B.T@GFa@B)
    allowed={0,1,n-1,n};mask=np.ones_like(G,dtype=bool)
    for i in allowed:
        for j in allowed:mask[i,j]=False
    return dict(sites=n,h=h,local_support_indices=sorted(allowed),
        local_response_outside_support=float(np.max(abs(G[mask]))),
        shadow_derivative_outside_local_support=float(np.max(abs(GFa[mask]))),
        matrix_response_error=float(np.max(abs(G-Gn))),
        shadow_average_error=float(np.max(abs(G-Ga))))

def lapse_control():
    rows=[]
    for h in [.4,.2,.1,.05,.025]:
        # Both t and v scale with the lapse N; response is at N=1.
        t,v=1.,1.3;S=step(t,v,h);G=response(t,v,t,v,h);GF=shadow_G(t,v,h)
        # A squeezed pure Gaussian, det C=1/4. The lapse response need not be
        # conserved by the same step even though the shadow quadratic is.
        C=np.diag([.7,1/(4*.7)]);start=.5*np.trace(G@C);fs=.5*np.trace(GF@C)
        err=0.;ferr=0.
        for _ in range(200):
            C=S@C@S.T
            err=max(err,abs(.5*np.trace(G@C)-start));ferr=max(ferr,abs(.5*np.trace(GF@C)-fs))
        rows.append(dict(h=h,lapse_response_change=float(err),shadow_energy_change=float(ferr)))
    return rows

def tkt_impulse_audit():
    """Exact frozen-X T/2,V,T/2 response used inside the main G/2,T/2,V,T/2,G/2 map."""
    rows=[]
    for lapse_power in [0,1]:
        b=(3-lapse_power)/2
        gamma=(lapse_power+1)/b;delta=(lapse_power+3)/b
        def coeff_x(X):
            return X**-2,1.3*X**gamma+.7**2*X**delta
        def derivative_x(X):
            return -2*X**-3,1.3*gamma*X**(gamma-1)+.7**2*delta*X**(delta-1)
        def tkt(X,h):
            t,v=coeff_x(X)
            D=np.array([[1.,h*t/2],[0.,1.]])
            K=np.array([[1.,0.],[-h*v,1.]])
            return D@K@D
        for X in [.7,1.,1.4]:
            for h in [.025,.1,.3]:
                t,v=coeff_x(X);tx,vx=derivative_x(X)
                halfq=np.array([1.,h*t/2]);endp=np.array([-h*v,1-h*h*t*v/2]);e=np.array([0.,1.])
                G=tx/2*(np.outer(e,e)+np.outer(endp,endp))+vx*np.outer(halfq,halfq)
                S=tkt(X,h);Gn=-J@np.linalg.solve(S,diff5(lambda xx:tkt(xx,h),X))/h
                C=np.array([[.7,.2],[.2,(.25+.2**2)/.7]])
                D=np.array([[1.,h*t/2],[0.,1.]]);K=np.array([[1.,0.],[-h*v,1.]])
                impulse=-h*tx*C[1,1]/4
                C1=D@C@D.T;impulse-=h*vx*C1[0,0]/2
                C2=K@C1@K.T;impulse-=h*tx*C2[1,1]/4
                target=-h/2*np.trace(G@C)
                rows.append(dict(lapse_power=lapse_power,X=X,h=h,
                    matrix_response_error=float(np.max(abs(G-Gn))),
                    geometry_impulse_error=float(abs(impulse-target))))
    return rows

def main():
    gates=gate_audit();lattice=lattice_audit();lapse=lapse_control()
    tkt=tkt_impulse_audit()
    result=dict(description=__doc__,gate_rows=gates,lattice=lattice,lapse_negative_control=lapse,tkt_impulse_rows=tkt)
    result['all_checks_passed']=bool(
        max(r['matrix_response_error'] for r in gates)<1e-8 and
        max(r['action_response_error'] for r in gates)<1e-8 and
        max(r['shadow_average_error'] for r in gates)<1e-8 and
        lattice['matrix_response_error']<1e-8 and lattice['shadow_average_error']<1e-8 and
        lattice['local_response_outside_support']==0. and
        lapse[0]['lapse_response_change']>1e-4 and
        max(r['shadow_energy_change'] for r in lapse)<1e-12 and
        max(r['matrix_response_error'] for r in tkt)<1e-8 and
        max(r['geometry_impulse_error'] for r in tkt)<1e-12)
    (ROOT/'audit_gate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));assert result['all_checks_passed']

if __name__=='__main__':main()
