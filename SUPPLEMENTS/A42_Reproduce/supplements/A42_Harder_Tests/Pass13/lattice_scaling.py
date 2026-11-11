"""Pass 13: spatial/volume tests of the established staggered Schwinger model.

No electric truncation: enumerate the complete neutral matter sector and solve
Gauss's law exactly, with fixed equal boundary fluxes. No new-physics claim.
Hamiltonian uses real nearest-neighbor hopping, phase-equivalent to standard KS.
The known mass improvement m_lat=m-g^2*a/8 is an independent prescription, not a
fit. Its homogeneous conformal extension substitutes g->g*Omega, m->m*Omega.
"""
from pathlib import Path
from itertools import combinations
import json,time
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp

def make_basis(N,boundary=0.):
    words=np.sort(np.fromiter((sum(1<<i for i in comb) for comb in combinations(range(N),N//2)),dtype=np.int64))
    n=((words[:,None]>>np.arange(N))&1).astype(np.int8)
    q=n-np.arange(N)%2;e=boundary+np.cumsum(q,axis=1)[:,:-1]
    rows=[];cols=[];pairs=[]
    for i in range(N-1):
        left=np.flatnonzero((n[:,i]==0)&(n[:,i+1]==1))
        right=np.searchsorted(words,words[left]^((1<<i)|(1<<(i+1))))
        pairs.append((left,right))
        rows.extend([left,right]);cols.extend([right,left])
    rows=np.concatenate(rows);cols=np.concatenate(cols)
    hopping=coo_matrix((-np.ones(len(rows)),(rows,cols)),shape=(len(words),len(words))).tocsr()
    assert np.max(abs(np.diff(np.column_stack([np.full(len(n),boundary),e,np.full(len(n),boundary)]),axis=1)-q))==0
    return dict(N=N,words=words,n=n,e=e,Q=(n*((-1)**np.arange(N))).sum(axis=1),
                V=(e*e).sum(axis=1),K=hopping,pairs=pairs,boundary=boundary)

def components(b,length,mass=0.,coupling=1.,omega=1.,improved=False,wall_offset=1):
    a=length/(b['N']+wall_offset);w=1/(2*a);J=coupling**2*a/2
    ml=mass*omega-(coupling**2*a*omega**2/8 if improved else 0.)
    return w,ml*b['Q']+J*omega**2*b['V']

def H(b,length,mass=0.,coupling=1.,omega=1.,improved=False,wall_offset=1):
    w,d=components(b,length,mass,coupling,omega,improved,wall_offset)
    return w*b['K']+diags(d)

def low_states(ham,k=2):
    dim=ham.shape[0];start=np.random.default_rng(4213).normal(size=dim);start/=np.linalg.norm(start)
    vals,vec=eigsh(ham,k=k,which='SA',tol=2e-11,v0=start,ncv=max(24,2*k+1))
    order=np.argsort(vals);vals=vals[order];vec=vec[:,order]
    residuals=[float(np.linalg.norm(ham@v-val*v)) for val,v in zip(vals,vec.T)]
    assert max(residuals)<2e-8,residuals
    return vals,vec,residuals

def spectrum(b,length,improved,wall_offset=1):
    start=time.time();ham=H(b,length,improved=improved,wall_offset=wall_offset);vals,vec,res=low_states(ham,k=3)
    target=float(np.sqrt(1/np.pi+(np.pi/length)**2));gap=float(vals[1]-vals[0])
    return dict(N=b['N'],length=length,a=length/(b['N']+wall_offset),wall_offset=wall_offset,dimension=len(b['words']),
       improved=improved,ground_energy=float(vals[0]),gap=gap,second_gap=float(vals[2]-vals[0]),
       continuum_box_gap=target,relative_gap_error=(gap-target)/target,
       infinite_volume_gap=float(1/np.sqrt(np.pi)),max_eigen_residual=max(res),seconds=time.time()-start)

def profile(eta,tau):return 1+.5*np.tanh(eta/tau)

def evolve(b,length=3.5,mass=.4,coupling=np.sqrt(2.88),tau=.6,improved=False,tol=2e-9,wall_offset=1):
    start=time.time();lo=-10*tau;hi=10*tau;om0=profile(lo,tau);omf=profile(hi,tau)
    h0=H(b,length,mass,coupling,om0,improved,wall_offset);hf=H(b,length,mass,coupling,omf,improved,wall_offset)
    val0,vec0,res0=low_states(h0);valf,vecf,resf=low_states(hf)
    initial=vec0[:,0].astype(complex)
    # Subtract a real scalar reference energy to reduce a rapidly rotating phase.
    # It cannot change any observable or state error after global-phase alignment.
    a=length/(b['N']+wall_offset);w=1/(2*a);J=coupling**2*a/2
    K=w*b['K'];q=b['Q'];v=b['V'];p0=abs(initial)**2
    kinetic0=float(np.vdot(initial,K@initial).real)
    def rhs(eta,state):
        om=profile(eta,tau);ml=mass*om-(coupling**2*a*om**2/8 if improved else 0.)
        diag=ml*q+J*om**2*v;scalar=kinetic0+float(p0@diag)
        return -1j*(K@state+(diag-scalar)*state)
    sol=solve_ivp(rhs,(lo,hi),initial,method='DOP853',rtol=tol,atol=tol/100,t_eval=[hi])
    assert sol.success
    psi=sol.y[:,-1];norm=float(np.vdot(psi,psi).real)
    assert abs(norm-1)<2e-6,norm
    survival=float(abs(np.vdot(vecf[:,0],psi))**2)
    excess=float(np.vdot(psi,hf@psi).real-valf[0]*norm)
    out=dict(N=b['N'],length=length,a=a,dimension=len(initial),improved=improved,tau=tau,
       ground_survival=survival,minus_log_survival_per_comoving_length=-float(np.log(survival))/length,
       final_excess_energy=excess,proper_excess_energy_density=excess/(length*omf**2),
       initial_gap=float(val0[1]-val0[0]),final_gap=float(valf[1]-valf[0]),
       initial_bare_even_occupation=float(p0@b['n'][:,::2].sum(axis=1)),
       final_bare_even_occupation=float(abs(psi)**2@b['n'][:,::2].sum(axis=1)),
       mean_flux_change=float((abs(psi)**2-p0)@b['e'].mean(axis=1)),
       norm_error=abs(norm-1),rtol=tol,evaluations=sol.nfev,
       max_eigen_residual=max(res0+resf),seconds=time.time()-start)
    return out,psi

def main():
    root=Path(__file__).parent;rows=[]
    # Fixed volume and varying a; several volumes separate boundary/momentum shift.
    for N in [6,8,10,12,14,16,18,20]:
        b=make_basis(N)
        for length in [4.,8.]:
            for improved in [False,True]:
                z=spectrum(b,length,improved);rows.append(z);print(json.dumps(z),flush=True)
        if N in [10,18]:
            for length in [4.,8.]:
                for improved in [False,True]:
                    z=spectrum(b,length,improved,wall_offset=0);rows.append(z);print(json.dumps(z),flush=True)
        (root/'spectral_results.json').write_text(json.dumps(dict(rows=rows,complete=False),indent=2)+'\n')
    (root/'spectral_results.json').write_text(json.dumps(dict(rows=rows,complete=True),indent=2)+'\n')

if __name__=='__main__':main()
