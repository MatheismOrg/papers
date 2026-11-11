import json,math,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq,minimize_scalar
from model import F,Fp,critical_dust,turning_points,y_from_a,rhs,quantum_stress,discrete_step

ROOT=Path(__file__).resolve().parent
NS,NW,NV=4,45,12
neff=NS+11*NW/2+62*NV
sm_lam=neff/(360*math.pi)
lambda_star=8/(3*math.sqrt(3))
config=dict(lambdas=[0.,.5,sm_lam,1.25,1.5],dust_ratios=[1.1,2.,10.],
            state_radiation=[0.,.1,1.],failure_lambda=1.6,
            sm_proxy=dict(real_scalars=NS,weyl_fermions=NW,vectors=NV,alpha_over_planck_area=1),
            renormalization='total Box R coefficient zero',main_state='radiation integration constant zero')
(ROOT/'inputs.json').write_text(json.dumps(config,indent=2)+'\n')
rows=[];allcurves={};stresschecks=[]

def evolve(c,lam,rad=0.,label=''):
    amin,amax=turning_points(c,lam,rad)
    a0=(amin+amax)/2
    y0=y_from_a(a0,c,lam,rad);v0=-math.sqrt(a0*a0*y0-1)
    def bounce(t,z):return z[1]
    bounce.direction=1;bounce.terminal=False
    def back(t,z):return z[0]-a0
    back.direction=1;back.terminal=True
    sol=solve_ivp(lambda t,z:rhs(t,z,c,lam,rad),(0,100*(amax+1)),[a0,v0],
                  method='DOP853',rtol=2e-11,atol=2e-12,events=[bounce,back],dense_output=True,max_step=.05)
    assert sol.success and len(sol.t_events[0])==1 and len(sol.t_events[1])==1
    t=np.linspace(0,sol.t[-1],1201);a,v=sol.sol(t);y=(1+v*v)/a**2
    fp=Fp(y,lam);density=c/a**3+rad/a**4
    D=y-(3*c/a**3+4*rad/a**4)/(2*fp)
    kappa=3*lam/(8*math.pi)
    state_density=3*rad/(8*math.pi*a**4)
    rho,p=quantum_stress(y,D,kappa,state_density)
    H=v/a;ydot=2*H*(D-y)
    rhodot=2*kappa*y*ydot-4*H*state_density
    ward=rhodot+3*H*(rho+p)
    trace=-rho+3*p
    # Curvature invariant expression independent of the compact 24*y*D.
    riem=12*(y*y+D*D);ric=12*(D*D+D*y+y*y);R=6*(D+y)
    euler=riem-4*ric+R*R
    anomaly=-kappa*euler/6
    quantum_fraction=np.divide(lam*y*y,y/np.sqrt(1-y*y),out=np.zeros_like(y),where=y!=0)
    btime=float(sol.t_events[0][0]);ab,vb=sol.y_events[0][0]
    out=dict(label=label,lam=lam,c=c,radiation=rad,amin_exact=amin,amax_exact=amax,
             amin_evolved=float(ab),bounce_time=btime,return_time=float(sol.t[-1]),
             bounce_radius_error=abs(float(ab)-amin),min_response=float(np.min(fp)),
             max_constraint_error=float(np.max(abs(F(y,lam)-density)/(1+abs(density)))),
             max_ward_residual=float(np.max(abs(ward)/(1+abs(rho)+abs(p)))),
             max_anomaly_residual=float(np.max(abs(trace-anomaly)/(1+abs(trace)))),
             max_kretschmann=float(np.max(riem)),max_abs_pressure=float(np.max(abs(p))),
             max_quantum_fraction=float(np.max(quantum_fraction)),
             velocity_return_error=abs(float(sol.y[1,-1])+v0))
    allcurves[label]=dict(t=t,a=a,v=v,y=y,D=D,rho_q=rho,p_q=p,kretschmann=riem,response=fp)
    return out

for lam in config['lambdas']:
    cc,xc=critical_dust(lam)
    for ratio in config['dust_ratios']:
        rows.append(evolve(ratio*cc,lam,label=f'lambda_{lam:.8f}_dust_{ratio}'))
state_rows=[]
for rad in config['state_radiation']:
    state_rows.append(evolve(2*critical_dust(sm_lam)[0],sm_lam,rad,label=f'state_{rad}'))

# Compare the same classical dust constant with and without quantum feedback.
fixed_c=2*critical_dust(0.)[0]
fixed_mass=[evolve(fixed_c,lam,label=f'fixed_c_lambda_{lam:.8f}') for lam in [0.,sm_lam]]

# Exact response bound: min_y 1/[2 y(1-y^2)^(3/2)].
opt=minimize_scalar(lambda y:1/(2*y*(1-y*y)**1.5),bounds=(.001,.999),method='bounded',options={'xatol':1e-14})
# Sufficient bound |D|<=1 for dust and rad=0; not claimed for general matter.
curv=minimize_scalar(lambda y:(2-y+3*y**3)/(y*(y+4)*(1-y*y)**1.5),bounds=(.001,.999),method='bounded',options={'xatol':1e-14})
# Low-curvature Einstein + anomaly comparison, fixed lam; no refitting.
low=[]
for density in [1e-3,1e-4,1e-5]:
    y=brentq(lambda y:F(y,sm_lam)-density,0,.01,xtol=1e-18)
    grq=2*density/(1+math.sqrt(1-4*sm_lam*density))
    low.append(dict(density=density,full_y=y,einstein_anomaly_y=grq,relative_difference=abs(y/grq-1),
                    normalized_leading_difference=(y-grq)/density**3))

# Failure control: finite-time encounter with F'=0 on the contracting branch.
lam=config['failure_lambda'];c=2.
yfold=brentq(lambda y:Fp(y,lam),1e-6,.5,xtol=1e-14)
afold=(c/F(yfold,lam))**(1/3);vfold=-math.sqrt(afold*afold*yfold-1)
ystart=yfold/2
def dtdy(y):
    a=(c/F(y,lam))**(1/3)
    H=-math.sqrt(y-1/a**2)
    return -Fp(y,lam)/(3*H*F(y,lam))
t_to_fold,terr=quad(dtdy,ystart,yfold,epsabs=1e-11,epsrel=1e-11)
approach=[]
for eps in [1e-3,1e-4,1e-5,1e-6]:
    y=yfold-eps;D=y-3*F(y,lam)/(2*Fp(y,lam))
    approach.append(dict(delta_y=eps,response=Fp(y,lam),D=D,K=12*(y*y+D*D)))
failure=dict(lam=lam,c=c,yfold=yfold,afold=afold,vfold=vfold,time_from_y_half=t_to_fold,
             quadrature_error=terr,approach=approach)

# One shared discrete time step, checked on a full contraction/bounce segment.
reference=allcurves['state_0.0'];T=float(reference['t'][-1]);a0=float(reference['a'][0]);v0=float(reference['v'][0])
c=2*critical_dust(sm_lam)[0]
discrete=[]
for steps in [300,600,1200]:
    dt=T/steps;a,v=a0,v0;maxres=0.;has_neg=False;has_pos=False
    for i in range(steps):
        a,v=discrete_step(a,v,dt,c,sm_lam)
        y=(1+v*v)/(a*a)
        maxres=max(maxres,abs(F(y,sm_lam)-c/a**3))
        has_neg|=v<0;has_pos|=v>0
    error=float(np.linalg.norm([a-reference['a'][-1],v-reference['v'][-1]]))
    discrete.append(dict(steps=steps,dt=dt,final_error=error,max_constraint_drift=maxres,bounced=bool(has_neg and has_pos)))

# External-geometry state checks. de Sitter is not a vacuum solution of
# the total NPQT equations at these weak lambdas without other sources.
kappa=neff/(960*math.pi**2)
limits=[]
for y in [0.,.01,.2,.9]:
    rho,p=quantum_stress(y,y,kappa)
    limits.append(dict(Phi=y,D=y,rho=rho,p=p,trace=-rho+3*p,expected_trace=-4*kappa*y*y))

# A proven uniform curvature bound uses F'>=q(1-y^2)^(-3/2),
# F<=y/sqrt(1-y^2), q=1-lambda/lambda_star. For zero extra radiation,
# D>=y(1-b+b*y^2), b=3/(2q); its exact minimum is -m below.
q=1-sm_lam/lambda_star;b=3/(2*q)
m=2*(b-1)**1.5/(3*math.sqrt(3*b))
uniform_bound=dict(q=q,b=b,min_D_lower_bound=-m,max_D_upper_bound=1.,
                   dimensionless_K_upper_bound=12*(1+max(1.,m*m)))

# Independent finite differences of rho along a high-accuracy orbit.
reference=allcurves['state_0.0'];T=float(reference['t'][-1]);c=2*critical_dust(sm_lam)[0]
sol=solve_ivp(lambda t,z:rhs(t,z,c,sm_lam),(0,T),[reference['a'][0],reference['v'][0]],
              method='DOP853',rtol=5e-13,atol=5e-14,dense_output=True,max_step=.01)
def values(t):
    a,v=sol.sol(t);y=(1+v*v)/a**2;acc=rhs(t,[a,v],c,sm_lam)[1];D=acc/a
    rho,p=quantum_stress(y,D,kappa)
    # Full angular/acceleration equation, h'(Phi)(D-Phi)=-4pi(rho_tot+p_tot).
    rho_d=3*c/(8*math.pi*a**3)
    angular=(1-y*y)**(-1.5)*(D-y)+4*math.pi*(rho_d+rho+p)
    return a,v,rho,p,angular
fd_ward=[];angular_checks=[]
eps=T*1e-4
for t in np.linspace(.01*T,.99*T,31):
    a,v,rho,p,angular=values(t)
    drho=(-values(t+2*eps)[2]+8*values(t+eps)[2]-8*values(t-eps)[2]+values(t-2*eps)[2])/(12*eps)
    fd_ward.append(abs(drho+3*v/a*(rho+p))/(1+abs(rho)+abs(p)))
    angular_checks.append(abs(angular)/(1+abs(rho)+abs(p)))

# Independent published covariant vacuum tensor, Herzog-Huang Eq. (23).
# Use mixed Ricci eigenvalues and compare every diagonal stress component.
tensor_errors=[]
for y in [.001,.1,.5,.9]:
  for D in [-.8,0.,.3,1.]:
    ric=np.array([3*D,D+2*y,D+2*y,D+2*y]);R=float(np.sum(ric))
    mixed=-(kappa/6)*((R*R/2-float(np.sum(ric*ric)))+2*ric*ric-4*R*ric/3)
    rho,p=quantum_stress(y,D,kappa)
    tensor_errors.append(float(np.max(abs(mixed-np.array([-rho,p,p,p])))))
# The static Einstein-universe check fixes the zero-radiation conformal
# vacuum in the same zero-Box-R prescription: E=3*a_Euler/(4*radius).
radius=2.;a_euler=neff/360
rho_esu,p_esu=quantum_stress(1/radius**2,0.,kappa)
casimir=dict(radius=radius,energy=rho_esu*2*math.pi**2*radius**3,
             expected=3*a_euler/(4*radius),pressure=p_esu)

allrows=rows+state_rows+fixed_mass
checks=dict(response_threshold=abs(opt.fun/lambda_star-1)<1e-11,
  bounce_evolution=all(r['bounce_radius_error']<1e-7 and r['velocity_return_error']<1e-7 for r in allrows),
  constraint=all(r['max_constraint_error']<1e-8 for r in allrows),
  conservation=all(r['max_ward_residual']<1e-12 for r in allrows),
  anomaly=all(r['max_anomaly_residual']<1e-11 for r in allrows),
  sufficient_curvature_bound=all(r['max_kretschmann']<=24*(1+1e-7) for r in rows if r['lam']<curv.fun),
  excessive_coupling_rejected=bool(approach[-1]['K']>1000*approach[0]['K'] and abs(Fp(yfold,lam))<1e-10),
  low_curvature_limit=bool(low[-1]['relative_difference']<1e-9),
  discrete_bounce=all(r['bounced'] for r in discrete),
  discrete_convergence=all(3.5<discrete[i]['final_error']/discrete[i+1]['final_error']<4.5 for i in range(2)),
  de_sitter_anomaly=all(abs(r['trace']-r['expected_trace'])<1e-14 for r in limits))
checks.update(ward_independent_finite_difference=max(fd_ward)<1e-8,
              full_angular_equation=max(angular_checks)<1e-10,
              analytic_SM_uniform_bound=m<1,
              covariant_quantum_tensor=max(tensor_errors)<1e-13,
              conformal_casimir_vacuum=abs(casimir['energy']-casimir['expected'])<1e-13)
checks={k:bool(v) for k,v in checks.items()}
out=dict(config=config,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
         free_SM_proxy=dict(N_eff=neff,lambda_at_alpha_planck_area=sm_lam,conditional_alpha_min_over_planck_area=sm_lam/lambda_star),
         bounds=dict(response_lambda_exact=lambda_star,response_lambda_numeric=opt.fun,response_y=opt.x,
                     sufficient_dust_curvature_lambda=curv.fun,sufficient_dust_curvature_y=curv.x),
         cosmologies=rows,state_sensitivity=state_rows,low_curvature=low,failure_control=failure,
         fixed_dust_comparison=fixed_mass,analytic_SM_uniform_bound=uniform_bound,
         independent_checks=dict(max_fd_ward=max(fd_ward),max_angular_equation=max(angular_checks),max_covariant_tensor_error=max(tensor_errors)),
         conformal_vacuum_casimir=casimir,
         discrete_steps=discrete,external_de_sitter_limits=limits,checks=checks)
(ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
np.savez_compressed(ROOT/'trajectories.npz',**{f'{name}__{k}':v for name,curve in allcurves.items() for k,v in curve.items()})
print(json.dumps(dict(bounds=out['bounds'],sm=out['free_SM_proxy'],checks=checks,
     sm_cases=[r for r in rows if r['lam']==sm_lam],failure=failure,discrete=discrete),indent=2))
assert all(checks.values()),checks
