"""Independent algebra/metric-variation audit for pass 17; no main-code imports."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent

def check():
    mass=.4;M2=1.7;V=2.3;Lambda=.05
    K=np.array([0.,.4,3.,8.])
    qq=np.array([.6,.2,.4,.8]);pp=np.array([.7,1.3,2.,.9]);qp=np.array([.1,-.2,.1,-.3])
    a=1.23;pa=-3.14
    def energy(a):
        return .5*np.sum(pp/a**3+(a*K+a**3*mass**2)*qq)+V*Lambda*a**3
    def constraint(a,pa):
        return -pa**2/(12*M2*V*a)+energy(a)
    rho=energy(a)/(V*a**3)
    pressure=.5*np.sum(pp/a**6-K*qq/(3*a**2)-mass**2*qq)/V-Lambda
    dh=1e-4
    Ea=(energy(a-2*dh)-8*energy(a-dh)+8*energy(a+dh)-energy(a+2*dh))/(12*dh)
    pressure_fd=-Ea/(3*V*a*a)
    rows=[]
    for s in [0,1]:
        b=(3-s)/2;X=a**b;P=pa/(b*a**(b-1));gamma=(s+1)/b;delta=(s+3)/b
        H=-b*b*P*P/(12*M2*V)+.5*np.sum(pp/X**2+(K*X**gamma+mass**2*X**delta)*qq)+V*Lambda*X**delta
        Xdot=-b*b*P/(6*M2*V)
        Pdot=np.sum(pp)/X**3-.5*np.sum((gamma*K*X**(gamma-1)+delta*mass**2*X**(delta-1))*qq)-delta*V*Lambda*X**(delta-1)
        adot=Xdot/(b*a**(b-1))
        padot=b*((b-1)*a**(b-2)*adot*P+a**(b-1)*Pdot)
        Ca=pa**2/(12*M2*V*a*a)+.5*np.sum(-3*pp/a**4+(K+3*a*a*mass**2)*qq)+3*V*Lambda*a*a
        padot_expected=-a**s*Ca-(s*a**(s-1)*constraint(a,pa) if s else 0.)
        adot_expected=-a**s*pa/(6*M2*V*a)
        rows.append(dict(lapse_power=s,hamiltonian_error=float(abs(H-a**s*constraint(a,pa))),adot_error=float(abs(adot-adot_expected)),padot_error=float(abs(padot-padot_expected))))
    # Continuity identity holds off the gravitational constraint as well.
    Hubble=-pa/(6*M2*V*a*a)
    adot=Hubble*a
    qqdot=2*qp/a**3
    ppdot=-2*(a*K+a**3*mass**2)*qp
    rhodot=.5*np.sum(ppdot/a**6-6*pp*adot/a**7+K*qqdot/a**2-2*K*qq*adot/a**3+mass**2*qqdot)/V
    continuity=float(abs(rhodot+3*Hubble*(rho+pressure)))
    out=dict(gauge_checks=rows,pressure_metric_derivative_error=float(abs(pressure-pressure_fd)),continuity_error=continuity)
    out['all_checks_passed']=all(max(row[k] for k in ['hamiltonian_error','adot_error','padot_error'])<1e-12 for row in rows) and abs(pressure-pressure_fd)<1e-10 and continuity<1e-12
    return out

if __name__=='__main__':
    out=check();(ROOT/'gravity_formula_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));assert out['all_checks_passed']
