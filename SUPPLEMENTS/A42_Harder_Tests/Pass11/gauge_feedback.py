"""Pass 11: finite open staggered U(1) chain with dynamical electric links.

Established lattice Hamiltonian benchmark, not a derivation from the previous
Dirac walk or a curved-spacetime/quantum-gravity model. hbar=1, lattice units.
H = -w sum(c_i^dag U_i c_{i+1}+h.c.) + m sum((-1)^i n_i)
    + (g^2/2) sum(E_i^2). [E,U]=U with hard electric cutoff.
Gauss G_i=E_i-E_{i-1}-(n_i-background_i), background_i=i mod 2.
Open boundary fields E_-1=E_(N-1)=1. Initial bare staggered vacuum, E=1.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh

N=6;W=1.;M=.4;G=1.2;BOUNDARY=1;T=1.5

def basis(cutoff,full=False):
 dims=tuple(x for i in range(N) for x in ((2,2*cutoff+1) if i<N-1 else (2,)))
 if full:
  digits=np.array(np.unravel_index(np.arange(np.prod(dims)),dims)).T
  ns=digits[:,::2];es=digits[:,1::2]-cutoff
 else:
  ns=[];es=[]
  for matter in product([0,1],repeat=N):
   q=np.array(matter)-np.arange(N)%2
   e=BOUNDARY+np.cumsum(q)[:-1]
   if sum(q)==0 and np.max(abs(e))<=cutoff:ns.append(matter);es.append(e)
  ns=np.array(ns);es=np.array(es)
  digits=np.empty((len(ns),2*N-1),int);digits[:,::2]=ns;digits[:,1::2]=es+cutoff
 codes=np.ravel_multi_index(digits.T,dims)
 index={int(c):j for j,c in enumerate(codes)}
 gauss=np.diff(np.column_stack([np.full(len(ns),BOUNDARY),es,np.full(len(ns),BOUNDARY)]),axis=1)-(ns-np.arange(N)%2)
 energy=M*(ns*((-1.)**np.arange(N))).sum(axis=1)+G*G/2*(es*es).sum(axis=1)
 initial=np.flatnonzero(np.all(ns==np.arange(N)%2,axis=1)&np.all(es==BOUNDARY,axis=1))[0]
 return {'n':ns,'e':es,'digits':digits,'dims':dims,'codes':codes,'index':index,'gauss':gauss,'energy':energy,'initial':initial,'cutoff':cutoff}

def bonds(b,mode='correct'):
 pairs=[];L=b['cutoff'];digits=b['digits'];ns=b['n'];es=b['e']
 for i in range(N-1):
  allowed=(ns[:,i]==0)&(ns[:,i+1]==1)
  if mode=='correct':allowed&=es[:,i]<L
  ids=np.flatnonzero(allowed);target=digits[ids].copy()
  target[:,2*i]=1;target[:,2*i+2]=0
  if mode=='correct':target[:,2*i+1]+=1
  if mode=='wrap':target[:,2*i+1]=(target[:,2*i+1]+1)%(2*L+1)
  codes=np.ravel_multi_index(target.T,b['dims'])
  pairs.append((ids,np.array([b['index'][int(c)] for c in codes],dtype=int)))
 return pairs

def circuit(b,steps,mode='correct',time=T):
 dt=time/steps;psi=np.zeros(len(b['n']),complex);psi[b['initial']]=1
 phase=np.exp(-.5j*dt*b['energy']);pairs=bonds(b,mode)
 def layer(parity,fraction):
  c=np.cos(W*dt*fraction);s=1j*np.sin(W*dt*fraction)
  for i in range(parity,N-1,2):
   left,right=pairs[i];u=psi[left].copy();v=psi[right].copy()
   psi[left]=c*u+s*v;psi[right]=s*u+c*v
 for _ in range(steps):
  psi*=phase;layer(0,.5);layer(1,1);layer(0,.5);psi*=phase
 return psi

def reference(cutoff,time=T,coupling=G):
 # Independent basis and Hamiltonian construction from matter bit strings.
 # Electric fields are reconstructed algebraically by Gauss's law.
 configs=[];fields=[]
 for word in range(2**N):
  n=np.array([(word>>i)&1 for i in range(N)])
  q=n-np.arange(N)%2;e=BOUNDARY+np.cumsum(q)[:-1]
  if q.sum()==0 and max(abs(e))<=cutoff:configs.append(n);fields.append(e)
 configs=np.array(configs);fields=np.array(fields);lookup={tuple(n):i for i,n in enumerate(configs)}
 H=np.diag(M*(configs*((-1.)**np.arange(N))).sum(axis=1)+coupling**2/2*(fields**2).sum(axis=1)).astype(float)
 for j,n in enumerate(configs):
  for i in range(N-1):
   if n[i]!=n[i+1]:
    changed=n.copy();changed[i],changed[i+1]=changed[i+1],changed[i]
    dest=lookup.get(tuple(changed))
    if dest is not None:H[dest,j]-=W
 initial=np.flatnonzero(np.all(configs==np.arange(N)%2,axis=1))[0]
 values,vectors=eigh(H);psi=vectors@(np.exp(-1j*time*values)*vectors[initial].conj())
 return configs,fields,psi,H

def mapped_reference(b,cutoff,time=T,coupling=G):
 ns,es,psi,H=reference(cutoff,time,coupling)
 out=np.zeros(len(b['n']),complex)
 for n,e,amp in zip(ns,es,psi):
  d=np.empty(2*N-1,int);d[::2]=n;d[1::2]=e+b['cutoff']
  if min(d[1::2])<0 or max(d[1::2])>=2*b['cutoff']+1:continue
  code=int(np.ravel_multi_index(d,b['dims']));out[b['index'][code]]=amp
 return out

def stats(b,psi):
 p=abs(psi)**2;n=b['n'];e=b['e'];q=n-np.arange(N)%2
 return {'norm_error':float(abs(p.sum()-1)),
         'gauss_squared':float(p@np.sum(b['gauss']**2,axis=1)),
         'physical_sector_probability':float(p[np.all(b['gauss']==0,axis=1)].sum()),
         'net_charge_squared':float(p@(q.sum(axis=1)**2)),
         'pair_mean':float(p@n[:,::2].sum(axis=1)),
         'mean_electric_field':(p@e).tolist(),
         'electric_energy':float(G*G/2*(p@(e*e).sum(axis=1))),
         'matter_mass_energy':float(M*(p@(n*((-1.)**np.arange(N))).sum(axis=1))),
         'vacuum_survival':float(p[b['initial']])}

def phase_error(a,b):
 phase=np.angle(np.vdot(b,a));return float(np.linalg.norm(a*np.exp(-1j*phase)-b))

def main():
 full=basis(1,True);phys=basis(1);rows=[]
 target=mapped_reference(full,1)
 for steps in [8,16,32,64,128,256,512,1024]:
  state=circuit(full,steps);rows.append({'steps':steps,'state_error':phase_error(state,target),**stats(full,state)})
 print(json.dumps({'full_dimension':len(full['n']),'physical_dimension':len(phys['n']),'finest':rows[-1]}),flush=True)
 # Local gate support itself preserves all G_i; this checks the whole Hilbert space.
 local_gauss_change=max(float(np.max(abs(full['gauss'][j]-full['gauss'][k]))) for j,k in bonds(full))
 reduced_state=circuit(phys,128);full_state=circuit(full,128)
 lifted=np.zeros(len(full_state),complex);lifted[phys['codes']]=reduced_state
 negative=[]
 for mode in ['frozen','wrap']:
  v=circuit(full,256,mode);negative.append({'mode':mode,**stats(full,v)})
 # Cutoff test: fixed six-site system, independent exact physical Hamiltonians.
 large=basis(4);exact_large=mapped_reference(large,4)
 cutoffs=[]
 for L in [1,2,3,4]:
  v=mapped_reference(large,L)
  cutoffs.append({'cutoff':L,'physical_dimension':len(basis(L)['n']),
     'state_error_vs_cutoff4':phase_error(v,exact_large),**stats(large,v)})
 resolved=basis(3);resolved_target=mapped_reference(resolved,3);resolved_rows=[]
 for steps in [16,32,64,128,256,512,1024,2048,4096]:
  v=circuit(resolved,steps)
  resolved_rows.append({'steps':steps,'state_error':phase_error(v,resolved_target),**stats(resolved,v)})
 feedback=[]
 for g in [0.,.6,1.2,2.]:
  v=mapped_reference(large,4,coupling=g)
  feedback.append({'coupling':g,'pair_mean':stats(large,v)['pair_mean'],'mean_electric_field':stats(large,v)['mean_electric_field']})
 timeline=[]
 for time in np.linspace(0,T,16):
  v=mapped_reference(large,4,time=time);timeline.append({'time':float(time),**stats(large,v)})
 ratios=[rows[i]['state_error']/rows[i+1]['state_error'] for i in range(len(rows)-1)]
 checks={'local_gate_max_gauss_change':local_gauss_change,'full_vs_reduced_error':float(np.linalg.norm(full_state-lifted))}
 assert local_gauss_change==0
 assert checks['full_vs_reduced_error']<1e-11
 assert max(x['norm_error'] for x in rows)<1e-10
 assert max(x['gauss_squared'] for x in rows)<1e-20
 assert all(3.8<x<4.2 for x in ratios[-3:])
 assert min(x['gauss_squared'] for x in negative)>.1
 assert max(x['norm_error'] for x in negative)<1e-10
 assert max(x['net_charge_squared'] for x in negative)<1e-20
 assert cutoffs[-2]['state_error_vs_cutoff4']<1e-12
 output={'scope':'Established finite 1+1 staggered U(1) lattice model with dynamical links. New benchmark branch, not derived from the prior free Dirac circuit. No continuum, curved-spacetime, full-QED or quantum-gravity claim.',
 'parameters':{'sites':N,'hopping':W,'mass':M,'coupling':G,'time':T,'boundary_field':BOUNDARY},
 'full_hilbert_dimension':len(full['n']),'physical_dimension_at_cutoff1':len(phys['n']),
 'rows':rows,'resolved_cutoff3_rows':resolved_rows,'error_ratios':ratios,'reference_cutoff1':stats(full,target),'negative_controls':negative,
 'cutoff_study':cutoffs,'coupling_study':feedback,'timeline_cutoff4':timeline,'audit':checks,'all_assertions_passed':True}
 Path(__file__).with_name('gauge_feedback_results.json').write_text(json.dumps(output,indent=2)+'\n')
 print(json.dumps({'cutoff_study':cutoffs,'negative_controls':negative,'coupling_study':feedback,'audit':checks,'ratios':ratios},indent=2))
if __name__=='__main__':main()
