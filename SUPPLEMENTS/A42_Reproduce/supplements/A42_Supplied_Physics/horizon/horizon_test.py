"""Radial 1+1 massless transport across a supplied Schwarzschild horizon.
Painleve-Gullstrand metric ds²=-dt²+(dr+sqrt(rs/r)dt)².
v_s=s-sqrt(rs/r); H_s=-i(v_s partial_r+v_s'/2).
No angular potential, Hawking-vacuum calculation or central completion.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp

RS=2.;LOW=.25;HIGH=6.25;T=.7
def vel(r,s):return s-np.sqrt(RS/r)
def deriv(r):return .5*np.sqrt(RS)*r**(-1.5)
def initial(r,center):
    z=(r-center)/.30
    out=np.zeros_like(r,dtype=complex);inside=abs(z)<1
    out[inside]=np.exp(-1/(1-z[inside]**2))*np.exp(6j*(r[inside]-center))
    return out
def reference(r,s,center):
    # Backward characteristic and its Jacobian, independent of lattice gates.
    n=len(r)
    def rhs(t,z):
        q=z[:n];j=z[n:];return np.r_[-vel(q,s),-deriv(q)*j]
    sol=solve_ivp(rhs,(0,T),np.r_[r,np.ones(n)],rtol=3e-12,atol=3e-13)
    q,j=sol.y[:,-1].reshape(2,n)
    return initial(q,center)*np.sqrt(j)
def run(n,s,center):
    a=(HIGH-LOW)/n;r=LOW+(np.arange(n)+.5)*a
    mid=r[:-1]+a/2
    steps=int(np.ceil(T/(.4*a)));steps+=steps%2;dt=T/steps
    theta=np.arcsin(dt*vel(mid,s)/(2*a))
    assert np.isfinite(theta).all()
    psi=np.sqrt(a)*initial(r,center);norm=np.linalg.norm(psi);psi/=norm
    layers={}
    for parity in [0,1]:
        j=np.arange(parity,n-1,2)
        for fraction in [.5,1.]:
            ang=fraction*theta[j];layers[parity,fraction]=(j,np.cos(ang),np.sin(ang))
    for t in range(steps):
        p=t%2
        for par,frac in [(p,.5),(1-p,1.),(p,.5)]:
            j,c,b=layers[par,frac];u=psi[j].copy();v=psi[j+1].copy()
            psi[j]=c*u-b*v;psi[j+1]=b*u+c*v
    exact=np.sqrt(a)*reference(r,s,center)/norm
    ray=solve_ivp(lambda t,y:vel(y,s),(0,T),[center],rtol=1e-12,atol=1e-13).y[0,-1]
    row={'sites':n,'chirality':s,'initial_center':center,'steps':steps,
         'state_error':float(np.linalg.norm(psi-exact)),
         'norm_error':float(abs(np.linalg.norm(psi)-1)),
         'reference_norm_error':float(abs(np.linalg.norm(exact)-1)),
         'numerical_inside_probability':float(np.sum(abs(psi[r<RS])**2)),
         'reference_inside_probability':float(np.sum(abs(exact[r<RS])**2)),
         'mean_r':float(np.sum(r*abs(psi)**2)),'reference_mean_r':float(np.sum(r*abs(exact)**2)),
         'central_ray_r':float(ray),'boundary_probability':float(np.sum(abs(psi[(r<.5)|(r>6.)])**2))}
    return row
def main():
    rows=[]
    for s,c in [(-1,3.),(1,1.5),(1,3.)]:
        for n in [256,512,1024,2048,4096,8192]:
            row=run(n,s,c);rows.append(row);print(json.dumps(row),flush=True)
    # Purely omitting r=0 leaves curvature cut-off dependence K~a^-6.
    cutoff=[{'a':a,'first_radius':a/2,'K_at_first_cell':12*RS**2/(a/2)**6} for a in [.2,.1,.05,.025]]
    ratios=[cutoff[i+1]['K_at_first_cell']/cutoff[i]['K_at_first_cell'] for i in range(3)]
    assert np.allclose(ratios,64)
    for s,c in [(-1,3.),(1,1.5),(1,3.)]:
        z=[r for r in rows if r['chirality']==s and r['initial_center']==c]
        assert z[-1]['state_error']<z[-2]['state_error']/3.5
    assert max(r['norm_error'] for r in rows)<1e-11
    result={'scope':'Horizon-regular radial massless 1+1 supplied-background propagation. Open finite interval with packets away from boundaries. No central-singularity resolution or Hawking spectrum.',
            'metric':'ds²=-dt²+(dr+sqrt(rs/r)dt)², rs=2; v_plus/minus=plus/minus1-sqrt(rs/r)',
            'time':T,'rows':rows,'curvature_cutoff':cutoff,'curvature_refinement_ratios':ratios,'all_assertions_passed':True}
    Path(__file__).with_name('horizon_test_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
