"""Independent real-space check of the Bloch-reduced gate computation."""
import json
from pathlib import Path
import numpy as np
import expanding_dirac as e

def full_chain(n,k,tau,coordinate):
 dx=2*np.pi/n;start=-e.BOUND*tau;end=e.BOUND*tau
 duration=end-start if coordinate=='eta' else e.cosmic(end,tau)
 steps=int(np.ceil(duration/(.5*dx)));steps+=steps%2;dt=duration/steps
 seed=np.array([.3+.1j,.2-.4j,.7+.2j,-.1+.2j]);seed/=np.linalg.norm(seed)
 phase=np.exp(2j*k*dx*np.arange(n//2))/np.sqrt(n//2)
 psi=np.empty((n,2),complex)
 psi[::2]=phase[:,None]*seed[:2];psi[1::2]=phase[:,None]*seed[2:]
 for step in range(steps):
  eta=start+(step+.5)*dt if coordinate=='eta' else e.inverse((step+.5)*dt,tau)
  A=e.scale(eta,tau);mu=e.MASS*A if coordinate=='eta' else e.MASS
  speed=1 if coordinate=='eta' else 1/A
  theta=np.arcsin(dt*speed/(2*dx))
  def mass():
   u=psi[:,0].copy();v=psi[:,1].copy();c=np.cos(dt*mu/2);s=np.sin(dt*mu/2)
   psi[:,0]=c*u-1j*s*v;psi[:,1]=c*v-1j*s*u
  mass();par=step%2
  for p,f in [(par,.5),(1-par,1),(par,.5)]:
   j=np.arange(p,n,2);j2=(j+1)%n;c=np.cos(f*theta);s=np.sin(f*theta)
   u=psi[j].copy();v=psi[j2].copy();sign=np.array([1,-1])
   psi[j]=c*u-s*v*sign;psi[j2]=s*u*sign+c*v
  mass()
 u,_=e.circuit(n,k,e.MASS,tau,coordinate);target=u@seed
 expected=np.empty_like(psi);expected[::2]=phase[:,None]*target[:2];expected[1::2]=phase[:,None]*target[2:]
 return float(np.linalg.norm(psi-expected))

def main():
 chain=[{'clock':clock,'state_error':full_chain(64,2.,.15,clock)} for clock in ['eta','t']]
 r=json.load(open(Path(__file__).with_name('expanding_dirac_results.json')))
 assert max(x['state_error'] for x in chain)<1e-11
 assert max(x['two_clock_reference_error'] for x in r['reference_audits'])<1e-10
 assert max(x['tight_tolerance_error'] for x in r['reference_audits'])<1e-10
 assert max(x['unitarity_error'] for x in r['rows'])<1e-10
 assert all(3.8<x['ratios'][-1]<4.2 for x in r['error_ratios'])
 assert max(x.get('norm_error',0) for x in r['fock_results'])<1e-10
 assert r['controls']['massless']['pair_mean_two_branches']<1e-10
 assert r['controls']['static']['pair_mean_two_branches']<1e-10
 assert r['controls']['spin_connection']['rescaled_state_error']<1e-10
 assert abs(r['controls']['spin_connection']['physical_norm']-1)<1e-10
 assert abs(r['controls']['sudden_approach'][-1]['pair_probability']-r['controls']['sudden_reference_probability'])<2e-6
 assert r['controls']['adiabatic_approach'][-1]['pair_probability']<2e-9
 assert all(x['propagator_error']>2.4 for x in r['controls']['wrong_mass'])
 checks={'full_chain_vs_Bloch':chain,'all_assertions_passed':True,
         'scope':'Numerical implementation consistency, not new physics or proof of a ToE.'}
 Path(__file__).with_name('audit_results.json').write_text(json.dumps(checks,indent=2)+'\n')
 print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
