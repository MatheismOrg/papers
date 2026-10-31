"""Checks against dense diagonalization, free-chain spectrum, and local gates."""
from pathlib import Path
import json,importlib.util
import numpy as np
from scipy.linalg import eigh
import lattice_scaling as s

root=Path(__file__).parent;dense=[];free=[]
for N in [6,8,10]:
    b=s.make_basis(N)
    for improved in [False,True]:
        ham=s.H(b,4.,improved=improved)
        sparse,_,res=s.low_states(ham,k=3)
        exact=eigh(ham.toarray(),eigvals_only=True,subset_by_index=[0,2])
        dense.append(dict(N=N,improved=improved,max_eigenvalue_error=float(np.max(abs(exact-sparse)))))
    ham=s.H(b,4.,coupling=0.)
    vals,_,_=s.low_states(ham)
    a=4/(N+1);prediction=2/a*np.sin(np.pi/(2*(N+1)))
    free.append(dict(N=N,calculated_gap=float(vals[1]-vals[0]),analytic_gap=float(prediction),
                     error=abs(float(vals[1]-vals[0])-prediction)))
assert max(z['max_eigenvalue_error'] for z in dense)<1e-10
assert max(z['error'] for z in free)<1e-10

# Independent six-site previous implementation uses a different basis builder.
prior=root.parent/'toe_pass12/curved_gauge.py'
if not prior.exists():prior=root/'curved_gauge_pass12.py'
spec=importlib.util.spec_from_file_location('pass12',prior);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
b=s.make_basis(6,boundary=1);reference=old.reference_operators()
order=[np.flatnonzero(np.all(b['n']==n,axis=1))[0] for n in reference['n']]
errors=[]
for omega in [.5,1.,1.5]:
    ham=s.H(b,3.5,.4,np.sqrt(2.88),omega).toarray()
    errors.append(float(np.max(abs(ham[np.ix_(order,order)]-old.hamiltonian(reference,omega)))))
assert max(errors)<1e-12
_,oldstate,_=old.reference_evolution(reference,.6)
z,newstate=s.evolve(b)
old_error=float(np.linalg.norm(newstate[order]*np.exp(-1j*np.angle(np.vdot(oldstate,newstate[order])))-oldstate))
assert old_error<2e-8

def local_circuit(b,steps,improved):
    length=3.5;mass=.4;coupling=np.sqrt(2.88);tau=.6;lo=-10*tau;hi=10*tau
    _,vec,_=s.low_states(s.H(b,length,mass,coupling,s.profile(lo,tau),improved))
    psi=vec[:,0].astype(complex);dt=(hi-lo)/steps
    for step in range(steps):
        omega=s.profile(lo+(step+.5)*dt,tau)
        w,diag=s.components(b,length,mass,coupling,omega,improved)
        phase=np.exp(-.5j*dt*diag);psi*=phase
        for parity,fraction in [(0,.5),(1,1.),(0,.5)]:
            co=np.cos(w*dt*fraction);si=1j*np.sin(w*dt*fraction)
            for i in range(parity,b['N']-1,2):
                left,right=b['pairs'][i];u=psi[left].copy();v=psi[right].copy()
                psi[left]=co*u+si*v;psi[right]=si*u+co*v
        psi*=phase
    return psi

gates=[];b=s.make_basis(10,boundary=1.)
for improved in [False,True]:
    _,refstate=s.evolve(b,improved=improved,tol=2e-11)
    rows=[]
    for steps in [1024,2048,4096]:
        state=local_circuit(b,steps,improved)
        error=float(np.linalg.norm(state*np.exp(-1j*np.angle(np.vdot(refstate,state)))-refstate))
        rows.append(dict(steps=steps,state_error=error,norm_error=float(abs(np.vdot(state,state).real-1))))
    ratios=[rows[i]['state_error']/rows[i+1]['state_error'] for i in range(2)]
    assert all(3.9<v<4.1 for v in ratios),ratios
    gates.append(dict(N=10,improved=improved,rows=rows,ratios=ratios))

result=dict(dense_diagonalization=dense,free_chain=free,pass12_generator_errors=errors,
    pass12_reference_state_error=old_error,local_gate_checks=gates,all_assertions_passed=True)
(root/'audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
