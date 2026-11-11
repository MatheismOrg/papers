"""Pass 10: supplied 1+1 FLRW expansion, fermion modes, two time coordinates.

ds^2=A(eta)^2(-deta^2+dx^2)=-dt^2+A(t)^2 dx^2; dt=A deta.
chi=sqrt(A)*psi: H_eta=k Z+m A X; H_t=(k/A) Z+m X.
The exact 4x4 Bloch reduction of the previous nearest-neighbor circuit is used.
Both parity branches are retained; no removal of doublers or quantum gravity.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
from scipy.optimize import brentq

I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex);Z=np.diag([1,-1]).astype(complex)
AI=.5;AF=3.;MASS=.7;MOMENTA=(1.,2.);BOUND=10.

def scale(eta,tau,static=False):
 return AI if static else (AI+AF)/2+(AF-AI)/2*np.tanh(eta/tau)
def derivative(eta,tau):return (AF-AI)/(2*tau)/np.cosh(eta/tau)**2
def cosmic(eta,tau,static=False):
 if static:return AI*(eta+BOUND*tau)
 return ((AI+AF)/2*(eta+BOUND*tau)+(AF-AI)/2*tau*
         (np.logaddexp(eta/tau,-eta/tau)-np.logaddexp(BOUND,-BOUND)))
def inverse(t,tau,static=False):
 return brentq(lambda e:cosmic(e,tau,static)-t,-BOUND*tau,BOUND*tau,xtol=1e-14)

def ref(k,m,tau,coordinate='eta',static=False,rtol=2e-12):
 start=-BOUND*tau;end=BOUND*tau
 if coordinate=='eta':
  def rhs(e,y):return (-1j*(k*Z+m*scale(e,tau,static)*X)@y.reshape(2,2)).ravel()
  sol=solve_ivp(rhs,(start,end),I.ravel(),method='DOP853',rtol=rtol,atol=rtol*.03)
  return sol.y[:,-1].reshape(2,2)
 # Integrate eta(t) jointly, independently of the analytic clock inversion.
 def rhs(t,y):
  e=y[0].real;A=scale(e,tau,static)
  return np.r_[1/A,(-1j*(k/A*Z+m*X)@y[1:].reshape(2,2)).ravel()]
 sol=solve_ivp(rhs,(0,cosmic(end,tau,static)),np.r_[complex(start),I.ravel()],method='DOP853',rtol=rtol,atol=rtol*.03)
 assert abs(sol.y[0,-1].real-end)<1e-8
 return sol.y[1:,-1].reshape(2,2)

def layers(theta,q):
 c=np.cos(theta);s=np.sin(theta)
 return (np.block([[c*I,-s*Z],[s*Z,c*I]]),
         np.block([[c*I,s*np.exp(-1j*q)*Z],[-s*np.exp(1j*q)*Z,c*I]]))

def circuit(n,k,m,tau,coordinate='eta',static=False,wrong_mass=False):
 dx=2*np.pi/n;start=-BOUND*tau;end=BOUND*tau
 duration=(end-start) if coordinate=='eta' else cosmic(end,tau,static)
 steps=int(np.ceil(duration/(.5*dx)));steps+=steps%2;dt=duration/steps
 U=np.eye(4,dtype=complex);q=2*dx*k
 for j in range(steps):
  e=(start+(j+.5)*dt) if coordinate=='eta' else inverse((j+.5)*dt,tau,static)
  A=scale(e,tau,static)
  v=1. if coordinate=='eta' else 1/A
  mu=m*A if coordinate=='eta' else m
  if wrong_mass:mu=m # negative control applies only in conformal time
  theta=np.arcsin(dt*v/(2*dx))
  eh,oh=layers(theta/2,q);ef,of=layers(theta,q)
  mass=np.kron(I,np.cos(dt*mu/2)*I-1j*np.sin(dt*mu/2)*X)
  kin=(eh@of@eh) if j%2==0 else (oh@ef@oh)
  U=mass@kin@mass@U
 return U,{'sites':n,'steps':steps,'time_coordinate':coordinate,'spacing':dx}

def basis(n,k):
 dx=2*np.pi/n;d=np.diag([1,1,np.exp(1j*k*dx),np.exp(1j*k*dx)])
 return np.column_stack([d@np.kron(np.array([1,p])/np.sqrt(2),I[:,j]) for p in [1.,-1.] for j in range(2)])

def energy_basis(k,m,A):
 return block_diag(*[np.linalg.eigh(p*k*Z+m*A*X)[1] for p in [1,-1]])

def observables(U,k,m,tau,static=False):
 bi=energy_basis(k,m,scale(-BOUND*tau,tau,static));bf=energy_basis(k,m,scale(BOUND*tau,tau,static))
 orbitals=bf.conj().T@U@bi[:,[0,2]]
 C=orbitals@orbitals.conj().T
 return {'pair_mean_two_branches':float(np.trace(C[np.ix_([1,3],[1,3])]).real),
         'pair_probability_per_branch':np.diag(C).real[[1,3]].tolist(),
         'unitarity_error':float(np.linalg.norm(U.conj().T@U-np.eye(4)))},orbitals

def fock(orbitals):
 # 8 one-body orbitals, four occupied: all 70 amplitudes of that Fock sector.
 return np.array([np.linalg.det(orbitals[list(rows),:]) for rows in combinations(range(8),4)])

def unrescaled_check(k,m,tau):
 start=-BOUND*tau;end=BOUND*tau
 spin=np.linalg.eigh(k*Z+m*scale(start,tau)*X)[1][:,0]
 psi0=spin/np.sqrt(scale(start,tau))
 def rhs(e,y):
  A=scale(e,tau)
  return -1j*(k*Z+m*A*X)@y-.5*derivative(e,tau)/A*y
 s=solve_ivp(rhs,(start,end),psi0,method='DOP853',rtol=2e-12,atol=5e-14)
 chi=np.sqrt(scale(end,tau))*s.y[:,-1]
 target=ref(k,m,tau)@spin
 return {'rescaled_state_error':float(np.linalg.norm(chi-target)),
         'physical_norm':float(scale(end,tau)*np.vdot(s.y[:,-1],s.y[:,-1]).real),
         'physical_norm_if_connection_omitted':float(scale(end,tau)/scale(start,tau))}

def main():
 rows=[];audits=[];all_refs={};all_fock={}
 for tau in [.15,.6]:
  reference_blocks=[]
  for k in MOMENTA:
   re=block_diag(ref(k,MASS,tau),ref(-k,MASS,tau))
   rt=block_diag(ref(k,MASS,tau,'t'),ref(-k,MASS,tau,'t'))
   rr=block_diag(ref(k,MASS,tau,rtol=3e-13),ref(-k,MASS,tau,rtol=3e-13))
   audit={'tau':tau,'k':k,'two_clock_reference_error':float(np.linalg.norm(re-rt)),
          'tight_tolerance_error':float(np.linalg.norm(re-rr))}
   stats,orb=observables(re,k,MASS,tau);audit.update(stats);audits.append(audit)
   reference_blocks.append(orb);all_refs[tau,k]=re
  reference_fock=fock(block_diag(*reference_blocks))
  for n in [64,128,256,512,1024]:
   clock_outputs={}
   for coordinate in ['eta','t']:
    fockblocks=[]
    for k in MOMENTA:
     u,meta=circuit(n,k,MASS,tau,coordinate);b=basis(n,k);u=b.conj().T@u@b
     stats,orb=observables(u,k,MASS,tau);fockblocks.append(orb)
     row={**meta,'tau':tau,'k':k,**stats,
          'propagator_error':float(np.linalg.norm(u-all_refs[tau,k])),
          'reference_pair_mean_two_branches':observables(all_refs[tau,k],k,MASS,tau)[0]['pair_mean_two_branches']}
     rows.append(row)
    state=fock(block_diag(*fockblocks));phase=np.angle(np.vdot(reference_fock,state))
    clock_outputs[coordinate]=state
    all_fock[str((tau,n,coordinate))]={'tau':tau,'sites':n,'coordinate':coordinate,
       'phase_aligned_state_error':float(np.linalg.norm(state*np.exp(-1j*phase)-reference_fock)),
       'norm_error':float(abs(np.vdot(state,state).real-1))}
   a,b=clock_outputs['eta'],clock_outputs['t'];phase=np.angle(np.vdot(a,b))
   all_fock[str((tau,n,'comparison'))]={'tau':tau,'sites':n,'coordinate':'comparison',
       'two_clock_state_difference':float(np.linalg.norm(b*np.exp(-1j*phase)-a))}
  print(json.dumps({'tau':tau,'finest':[r for r in rows if r['tau']==tau and r['sites']==1024]}),flush=True)
 # Analytic limiting controls and independent geometric rescaling.
 controls={}
 for label,m,static in [('massless',0.,False),('static',MASS,True)]:
  u,_=circuit(1024,1.,m,.15,static=static);b=basis(1024,1.)
  stats,_=observables(b.conj().T@u@b,1.,m,.15,static);controls[label]=stats
 sudden=.5*(1-(1+MASS*MASS*AI*AF)/(np.hypot(1,MASS*AI)*np.hypot(1,MASS*AF)))
 controls['sudden_reference_probability']=sudden
 controls['sudden_approach']=[]
 controls['adiabatic_approach']=[]
 for tau in [.15,.03,.006,.0012,1.,2.,4.]:
  u=block_diag(ref(1.,MASS,tau),ref(-1.,MASS,tau))
  stat,_=observables(u,1.,MASS,tau)
  row={'tau':tau,'pair_probability':stat['pair_probability_per_branch'][0]}
  controls['sudden_approach' if tau<1 else 'adiabatic_approach'].append(row)
 controls['spin_connection']=unrescaled_check(1.,MASS,.15)
 controls['wrong_mass']=[]
 for n in [128,512,1024]:
  u,_=circuit(n,1.,MASS,.15,wrong_mass=True);b=basis(n,1.);u=b.conj().T@u@b
  stat,_=observables(u,1.,MASS,.15)
  controls['wrong_mass'].append({'sites':n,'propagator_error':float(np.linalg.norm(u-all_refs[.15,1.])),**stat})
 ratios=[]
 for tau in [.15,.6]:
  for k in MOMENTA:
   for coordinate in ['eta','t']:
    z=[r for r in rows if r['tau']==tau and r['k']==k and r['time_coordinate']==coordinate]
    ratios.append({'tau':tau,'k':k,'clock':coordinate,'ratios':[z[j]['propagator_error']/z[j+1]['propagator_error'] for j in range(len(z)-1)]})
 result={'scope':'Free Dirac modes in a supplied expanding 1+1 FLRW geometry. Two selected momenta and both parity branches. Quantum/Fock structure imposed; no interactions, backreaction, full vacuum sum or quantum gravity.',
         'scale_initial':AI,'scale_final':AF,'mass':MASS,'momenta':MOMENTA,
         'scale':'A(eta)=(Ai+Af)/2+(Af-Ai)/2*tanh(eta/tau); eta endpoints +/-10 tau',
         'reference_audits':audits,'rows':rows,'fock_results':list(all_fock.values()),'controls':controls,'error_ratios':ratios}
 Path(__file__).with_name('expanding_dirac_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'reference_audits':audits,'controls':controls,'error_ratios':ratios},indent=2))

if __name__=='__main__':main()
