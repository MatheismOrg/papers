"""Pass18: RST endpoint candidates with explicitly paired asymptotic charges.

This is an analytical reduced-model screening program, not a new BH solver or
SM/quantum-gravity simulator. Units choose A=kappa*lambda/4=1 and lambda=1.
The dimensionless incoming mass is mu=M/A (NOT the mu of arXiv2510.18298).
U=log q is shifted affine retarded time, r=q/(1+q)=expit(U).
Traditional: B/A=mu-log(1+q)-r, F/(A lambda)=r(2-r).
Quantum-corrected: B/A=mu-log(1+q)+r, F/(A lambda)=r^2.
Each charge has its OWN flux. Their difference is a boundary/quantum term;
negative traditional charge is not treated as a physical instability theorem.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent

def charge(U,mu,corrected=True):
 r=expit(U);return mu-np.logaddexp(0.,U)+(r if corrected else -r)

def flux(U,corrected=True):
 r=expit(U);return r*r if corrected else r*(2-r)

def weak_z(gap):
 if gap<0:raise ValueError('No real RST dilaton at W<1')
 if gap==0:return 1.
 return brentq(lambda z:z-np.log(z)-1-gap,1.,max(4.,2+2*gap),xtol=5e-15)

def grid(mu,h):
 end=float(np.log(np.expm1(mu)));start=-24.
 t=start+np.arange(int(np.ceil((end-start)/h))+1)*h
 t[-1]=end
 return t

def ledger(mu,h,definition):
 t=grid(mu,h);g,w=leggauss(12)
 correct_B=definition!='traditional';correct_F=definition=='corrected'
 mid=(t[1:]+t[:-1])/2;half=np.diff(t)/2
 increments=(flux(mid[:,None]+half[:,None]*g,correct_F)*w).sum(axis=1)*half
 # Include analytic early tail; avoid a hidden finite start-time budget error.
 r0=expit(t[0]);emitted0=np.logaddexp(0.,t[0])+(-r0 if correct_F else r0)
 emitted=np.r_[emitted0,emitted0+np.cumsum(increments)]
 B=charge(t,mu,correct_B);residual=B+emitted-mu
 return dict(mu=mu,coordinate_step=h,definition=definition,steps=len(t)-1,
  endpoint_charge=float(B[-1]),emitted=float(emitted[-1]),
  max_abs_budget_residual=float(np.max(abs(residual))),
  endpoint_relative_budget_residual=float(residual[-1]/mu)),dict(U=t,charge=B,emitted=emitted,residual=residual)

def mass_clip(mu):
 # The deliberately naive rule cuts off flux when TRADITIONAL B first hits0.
 y0=brentq(lambda y:y+1-np.exp(-y)-mu,0.,mu,xtol=1e-14)
 r0=-np.expm1(-y0);q0=np.expm1(y0);z=weak_z(r0)
 return dict(mu=mu,y_clip=y0,q_clip=q0,U_clip=float(np.log(q0)),
  U_endpoint=float(np.log(np.expm1(mu))),
  horizon_W_gap=r0,horizon_z=z,X_horizon=q0/mu,postshock_horizon=bool(q0/mu>=1),
  corrected_charge_at_old_zero=2*r0,
  interpretation='Only the traditional-charge clipping rule. The formal horizon value is physical only if postshock_horizon=true; then apparent horizon is above critical z=1 and clipping cannot declare complete evaporation. No claim all Bondi definitions vanish here.')

def main():
 rows=[];clips=[];endpoints=[]
 for mu in [1.,2.,4.,8.,20.,40.]:
  d=-np.expm1(-mu)
  endpoints.append(dict(mu=mu,traditional_endpoint=-d,corrected_endpoint=d,
   traditional_radiation=mu+d,corrected_radiation=mu-d,
   mixing_budget_error=2*d,mixing_relative_error=2*d/mu))
  clips.append(mass_clip(mu))
  for h in [1.,.5,.25,.125]:
   for definition in ['traditional','corrected','mixed']:
    row,hist=ledger(mu,h,definition);rows.append(row)
    np.savez_compressed(ROOT/f'ledger_mu{mu:g}_h{h:g}_{definition}.npz',**hist)
 max_matched=max(x['max_abs_budget_residual'] for x in rows if x['definition']!='mixed')
 max_mixed_error=max(abs(x['endpoint_relative_budget_residual']-2*(-np.expm1(-x['mu']))/x['mu']) for x in rows if x['definition']=='mixed')
 result=dict(description=__doc__,rows=rows,endpoints=endpoints,naive_clips=clips,
  matched_max_budget_error=max_matched,mixed_prediction_max_error=max_mixed_error,
  all_checks_passed=max_matched<1e-11 and max_mixed_error<1e-12 and all(x['horizon_z']>1 for x in clips))
 (ROOT/'ledger_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['rows','description']},indent=2))

if __name__=='__main__':main()
