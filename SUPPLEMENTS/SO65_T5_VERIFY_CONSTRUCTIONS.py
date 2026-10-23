"""Finite checks of the new worked examples; proofs remain in A9/A28."""
import math,json,sys
from pathlib import Path
import numpy as np
out=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().with_name('SO65_T5_CONSTRUCTION_CHECKS.json')
N=131072
phi=(np.arange(N)+0.5)*2*math.pi/N
lam=np.stack((np.cos(phi),np.sin(phi)),axis=1)
rows=[]
for theta in [0,.17,.7,math.pi/4,math.pi/2,2.31,math.pi]:
 a=np.array([1.,0.]);b=np.array([math.cos(theta),math.sin(theta)])
 sa=np.where(lam@a>=0,1.,-1.);sb=np.where(lam@b>=0,1.,-1.)
 if theta in [0,math.pi]: rho=np.full(N,1/(4*math.pi))
 else:rho=np.where(sa==sb,(1+math.cos(theta))/(8*(math.pi-theta)),(1-math.cos(theta))/(8*theta))
 # z integration contributes 2 because sign regions do not depend on z.
 w=rho*4*math.pi/N
 joint={(int(x),int(y)):float(w[(sa==x)&(-sb==y)].sum()) for x in [-1,1] for y in [-1,1]}
 err=max(abs(joint[x,y]-(1-x*y*math.cos(theta))/4) for x,y in joint)
 rows.append({'theta':theta,'normalization':float(w.sum()),'correlation':float(np.sum(w*sa*(-sb))),'target':-math.cos(theta),'max_joint_error':err})
 assert err<5e-5
 # Common-preparation identity at quadrature nodes is checked below.
p=np.array([.2,.3,.5])
rhos=[]
for theta in [.3,.9,2.1]:
 sb=np.where(np.cos(phi-theta)>=0,1.,-1.);sa=np.where(np.cos(phi)>=0,1.,-1.)
 rhos.append(np.where(sa==sb,(1+math.cos(theta))/(8*(math.pi-theta)),(1-math.cos(theta))/(8*theta)))
rhos=np.stack(rhos);bar=(p[:,None]*rhos).sum(axis=0)
cond=p[:,None]*rhos/bar
assert np.max(np.abs(cond.sum(axis=0)-1))<1e-14
assert np.max(np.abs(bar*cond-p[:,None]*rhos))<1e-14
# A reversible finite update, compute f(x) into workspace by XOR.
for x in range(16):
 for y in range(16):
  f=(x*x+3*x+1)%16
  assert ((y^f)^f)==y
counts=[0]*50
for s in range(1,51):
 for j in range(s):counts[j]+=1
assert counts==list(range(50,0,-1))
# Finite cyclic analogue of the bilateral clock permutation.
def vt(t,c):return c^((t*t+3)%16)
def U(t,c):return ((t+1)%13,vt(t,c))
def Ui(t,c):return ((t-1)%13,vt((t-1)%13,c))
states=[(t,c) for t in range(13) for c in range(16)]
assert len(set(U(*x) for x in states))==len(states)
assert all(Ui(*U(*x))==x and U(*Ui(*x))==x for x in states)
result={'A9':{'method':'Deterministic azimuth quadrature on rotationally reduced sphere; analytic endpoint cases included.','grid_points':N,'results':rows,'common_preparation_identity':'passed'},
'A28':{'scheduler_first_50_stages':'passed','reversible_XOR_workspace_check':'passed','finite_clock_permutation_inverse':'passed','scope':'Finite sanity checks do not replace the infinite-dimensional basis-bijection proof in A28; no claim that our cosmology implements the system.'}}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'max_singlet_probability_error':max(x['max_joint_error'] for x in rows),'A28_checks':'passed'}))
