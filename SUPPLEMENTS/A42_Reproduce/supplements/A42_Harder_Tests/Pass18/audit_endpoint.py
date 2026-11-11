"""Independent exact RST endpoint-continuation diagnostics.

All formulas follow the RST shock solution in hep-th/9206070 plus explicit
chain-rule calculations. No claim of new physical dynamics or singularity cure.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import lambertw
from decimal import Decimal as D, getcontext
ROOT=Path(__file__).resolve().parent

def shock(u,v,kappa=4.,lam=1.,m=2.,u0=1.):
    a=kappa/4;D=m/(lam**3*u0)
    return -lam**2*u*(v+D)+m/lam-a*np.log(-lam**2*u*v)

def vacuum(u,v,kappa=4.,lam=1.,m=2.,u0=1.):
    a=kappa/4;D=m/(lam**3*u0);X=-lam**2*u*(v+D)
    return X-a*np.log(X)

def splice(kappa,lam,m,u0):
    a=kappa/4;D=m/(lam**3*u0);alpha=m/(a*lam)
    us=a*lam*u0/m*np.expm1(alpha)
    vs=-D/(-np.expm1(-alpha))
    jump=a*(1/vs-1/(vs+D))
    # T_vv = -jump * delta(v-vs), E= integral du T_uu.
    # du/dv = -1/[lambda(v+D)] for affine outgoing retarded coordinate.
    # Thus E = [-lambda(vs+D)] * (-jump).
    energy=lam*(vs+D)*jump
    target=-a*lam*(-np.expm1(-alpha))
    return dict(kappa=kappa,lam=lam,m=m,u0=u0,us=us,vs=vs,
       continuity_error=max(abs(shock(u,vs,kappa,lam,m,u0)-vacuum(u,vs,kappa,lam,m,u0)) for u in [1.1*us,2*us,5*us]),
       derivative_jump=jump,coordinate_delta_source=-jump,
       legacy_matter_source_energy=energy,RST1992_source_energy=target,energy_error=abs(energy-target))

def observer_curvature(u,X0,kappa=4.,lam=1.,m=2.,u0=1.):
    a=kappa/4;D=m/(lam**3*u0);W0=X0-a*np.log(X0)
    def f(Y): return Y+m/lam-a*np.log(lam*lam*D*u+Y)-W0
    hi=max(10.,10*a*np.log1p(lam*lam*D*u)+10*abs(W0)+10*m/lam)
    Y=brentq(f,0.,hi,xtol=1e-13)
    product=-lam*lam*(Y-a)*(1-a/(lam*lam*D*u+Y))
    gp=-1/(2*(X0-a));gpp=X0/(2*(X0-a)**3)
    R=8*X0*(-lam*lam*gp+gpp*product)
    return R,Y,product

rows=[splice(k,lam,m,u0) for k,lam,m,u0 in [(4.,1.,.1,1.),(4.,1.,1.,1.),(4.,1.,2.,1.),(8.,.7,3.,1.2),(12.,2.,5.,.6)]]
obs=[]
for X0 in [2.,4.,10.,100.]:
    r1=observer_curvature(1e8,X0)[0];r2=observer_curvature(1e16,X0)[0]
    slope=(r2-r1)/(np.log(1e16)-np.log(1e8))
    target=-4*X0*X0/(X0-1)**3
    obs.append(dict(X0=X0,R_u1e8=r1,R_u1e16=r2,slope=slope,
                    target=target,relative_error=abs(slope/target-1)))
# In a nonempty fully modified patch use a smooth bilinear local solution Z:
# Z(u,v)=c-u*v, Z_uv=-1. Discrete mixed derivative is exact for every h.
# Choosing c and patch center makes Z<delta throughout every tested cell.
floors=[]
for h in [.1,.05,.025,.0125]:
    u0=v0=1.;c=.5;delta=.2
    Z=np.array([[c-(u0+i*h)*(v0+j*h) for j in range(2)] for i in range(2)])
    def mixed(F):return (F[1,1]-F[1,0]-F[0,1]+F[0,0])/(h*h)
    floor=np.maximum(Z,delta);reflect=delta+abs(Z-delta)
    floors.append(dict(h=h,original_residual=mixed(Z)+1,
        floor_residual=mixed(floor)+1,reflection_residual=mixed(reflect)+1))
branches=[]
for Z in [-.1,-.001,1e-8,.01,1.,10.]:
    argument=-np.exp(-1-Z)
    weak=-lambertw(argument,-1);strong=-lambertw(argument,0)
    branches.append(dict(Z=Z,weak_real=float(weak.real),weak_imag=float(weak.imag),
        strong_real=float(strong.real),strong_imag=float(strong.imag),
        both_real=bool(abs(weak.imag)<1e-12 and abs(strong.imag)<1e-12)))
# A concrete smooth fold reconstruction, not an asserted consistent new action:
# replace Z by softplus_delta(Z), then use the weak inverse. It is regular at
# every finite W and tends back to the original branch for Z >> delta.
# Direct high-precision mixed differentiation audits the broad chain-rule claim.
getcontext().prec=75
delta=D('.2')
def smooth_X(W):
    e=((W-1)/delta).exp();S=delta*(1+e).ln()
    X=1+(2*S).sqrt()+S
    for _ in range(40):
        step=(X-X.ln()-1-S)/(1-1/X)
        X-=step
        if abs(step)<D('1e-68'):break
    return X,e/(1+e)
def smooth_g(W):
    X,_=smooth_X(W)
    return -X.ln()/2
def curved(W0,u):
    Y=(2*u).ln()+W0-2
    for _ in range(20):
        step=(Y+2-(2*u+Y).ln()-W0)/(1-1/(2*u+Y));Y-=step
        if abs(step)<D('1e-68'):break
    v=-2-Y/u
    def phi(uu,vv):return smooth_g(-uu*(vv+2)+2-(-uu*vv).ln())
    hu=u*D('1e-7');hv=D('1e-7')/u
    duv=(phi(u+hu,v+hv)-phi(u+hu,v-hv)-phi(u-hu,v+hv)+phi(u-hu,v-hv))/(4*hu*hv)
    return 8*(-2*phi(u,v)).exp()*duv
smooth=[]
for W0 in [D('.5'),D('1.'),D('2.'),D('5.')]:
    r1=curved(W0,D('1e8'));r2=curved(W0,D('1e16'))
    slope=(r2-r1)/D('1e8').ln()
    X,L=smooth_X(W0)
    gsecond=X*L*L/(2*(X-1)**3)-L*(1-L)/(2*(X-1)*delta)
    target=-8*X*gsecond
    smooth.append(dict(W0=float(W0),delta=float(delta),
        g=float(smooth_g(W0)),g_second=float(gsecond),
        direct_curvature_u1e8=float(r1),direct_curvature_u1e16=float(r2),
        direct_slope=float(slope),predicted_slope=float(target),
        relative_error=float(abs(slope/target-1))))
result=dict(splice_rows=rows,observer_rows=obs,bulk_residual_rows=floors,
     inverse_branch_rows=branches,smooth_reconstruction_rows=smooth,
     interpretations=[
      'The historical RST vacuum splice requires the derived negative null-constraint matter/source impulse. This is not identified here with a complete corrected Bondi endpoint budget; that requires matching Polyakov boundary data and charges too.',
      'Floor/reflection preserve neither the unchanged mixed field equation nor its source; residuals do not vanish with grid refinement.',
      'For rho=phi=g(W), the unspliced late observer curvature has logarithmic coefficient -8*a*lambda^2*exp(-2*g)*g_second. A smooth local reconstruction alone cannot cure all observers while retaining g~-.5logW at weak coupling.',
      'This is conditional screening, not a no-go theorem for modified gravity or a demonstrated singularity resolution.'
     ])
result['all_checks_passed']=(max(r['continuity_error'] for r in rows)<1e-11 and max(r['energy_error'] for r in rows)<1e-12 and max(r['relative_error'] for r in obs)<1e-6 and all(abs(r['floor_residual']-1)<1e-12 and abs(r['reflection_residual']-2)<1e-10 for r in floors) and all(r['both_real']==(r['Z']>=0) for r in branches) and max(r['relative_error'] for r in smooth)<1e-6)
(ROOT/'audit_endpoint.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
