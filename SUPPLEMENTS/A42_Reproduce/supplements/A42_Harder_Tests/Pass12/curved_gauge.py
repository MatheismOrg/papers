"""Pass 12: supplied conformal geometry and interacting staggered U(1).

Known-model benchmark, not a quantum-gravity model or continuum extrapolation.
ds^2 = Omega(eta,x)^2 (-deta^2 + dx^2).
H_eta = K + m sum((-1)^i Omega_i n_i) + J sum(Omega_link^2 E_i^2).
For homogeneous expansion dt=Omega deta, H_t=K/Omega+M+Omega V.
See Budd et al., arXiv:2601.02331v2, Eqs. (80), (89), (98).
The real hopping convention is related to imaginary staggered hopping by a
site phase rotation. Metric at link midpoints defines our finite regulator.
All reference calculations concern this SAME finite spatial lattice.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, expm

N=6; W=1.; M=.4; J=.72; BOUNDARY=1; A=.5

def local_basis(cutoff=3, full=False):
    dims=tuple(q for i in range(N) for q in ((2,2*cutoff+1) if i<N-1 else (2,)))
    if full:
        digits=np.array(np.unravel_index(np.arange(np.prod(dims)),dims)).T
        ns=digits[:,::2]; es=digits[:,1::2]-cutoff
    else:
        ns=[]; es=[]
        for n in product([0,1],repeat=N):
            charge=np.array(n)-np.arange(N)%2
            field=BOUNDARY+np.cumsum(charge)[:-1]
            if charge.sum()==0 and max(abs(field))<=cutoff:
                ns.append(n); es.append(field)
        ns=np.array(ns); es=np.array(es)
        digits=np.empty((len(ns),2*N-1),int)
        digits[:,::2]=ns; digits[:,1::2]=es+cutoff
    codes=np.ravel_multi_index(digits.T,dims)
    index={int(code):i for i,code in enumerate(codes)}
    pairs=[]
    for i in range(N-1):
        ids=np.flatnonzero((ns[:,i]==0)&(ns[:,i+1]==1)&(es[:,i]<cutoff))
        dest=digits[ids].copy();dest[:,2*i]=1;dest[:,2*i+2]=0;dest[:,2*i+1]+=1
        targets=np.ravel_multi_index(dest.T,dims)
        pairs.append((ids,np.array([index[int(c)] for c in targets],int)))
    gauss=np.diff(np.column_stack([np.full(len(ns),BOUNDARY),es,np.full(len(ns),BOUNDARY)]),axis=1)-(ns-np.arange(N)%2)
    return dict(n=ns,e=es,pairs=pairs,codes=codes,dims=dims,gauss=gauss,cutoff=cutoff)

def reference_operators(cutoff=3):
    # Separate construction: bit words + Gauss-eliminated matter Hamiltonian.
    configs=[]; fields=[]
    for word in range(2**N):
        n=np.array([(word>>i)&1 for i in range(N)])
        electric=[]; E=BOUNDARY
        for i in range(N-1):
            E+=int(n[i])-i%2;electric.append(E)
        if sum(n)==N//2 and max(abs(np.array(electric)))<=cutoff:
            configs.append(n);fields.append(electric)
    configs=np.array(configs);fields=np.array(fields)
    lookup={tuple(n):i for i,n in enumerate(configs)};dim=len(configs)
    K=np.zeros((dim,dim))
    for j,n in enumerate(configs):
        for i in range(N-1):
            if n[i]!=n[i+1]:
                q=n.copy();q[i],q[i+1]=q[i+1],q[i]
                dest=lookup.get(tuple(q))
                if dest is not None:K[dest,j]=-W
    return dict(n=configs,e=fields,K=K,
                mass=M*(configs*(-1.)**np.arange(N)),electric=J*fields**2)

def to_local(v, ref, local):
    lookup={tuple(n):z for n,z in zip(ref['n'],v)}
    return np.array([lookup.get(tuple(n),0) if np.all(G==0) else 0 for n,G in zip(local['n'],local['gauss'])],complex)

def profile(eta,tau):
    z=np.tanh(eta/tau)
    return 1+.5*z,.5/tau*(1-z*z)

def cosmic_time(eta,tau):
    logcosh=lambda z:np.logaddexp(z,-z)-np.log(2.)
    return eta+10*tau+.5*tau*(logcosh(eta/tau)-logcosh(10.))

def inverse_clock(t,tau):
    eta=np.array(t,dtype=float)-10*tau
    for _ in range(16):eta-=(cosmic_time(eta,tau)-t)/profile(eta,tau)[0]
    assert np.max(abs(cosmic_time(eta,tau)-t))<2e-13*max(1,tau)
    return eta

def hamiltonian(ref,omega_site=1.,omega_link=None,clock='eta',power=2):
    if omega_link is None:omega_link=omega_site
    diagonal=(ref['mass']*omega_site).sum(axis=1)+(ref['electric']*np.asarray(omega_link)**power).sum(axis=1)
    H=ref['K']+np.diag(diagonal)
    return H if clock=='eta' else H/omega_site

def derivative(ref,omega,prime,clock='eta',power=2):
    mass=ref['mass'].sum(axis=1);electric=ref['electric'].sum(axis=1)
    if clock=='eta':return np.diag(prime*(mass+power*omega**(power-1)*electric))
    return -prime/omega**3*ref['K']+np.diag(prime*(power-1)*omega**(power-3)*electric)

def reference_evolution(ref,tau,clock='eta',power=2,tol=2e-12):
    lo=-10*tau;hi=10*tau;omega0=profile(lo,tau)[0]
    # Every candidate starts from the same correct interacting initial vacuum.
    values,vectors=eigh(hamiltonian(ref,omega0))
    initial=vectors[:,0].astype(complex);dim=len(initial)
    y0=np.r_[initial,0j,complex(lo)]
    def rhs(time,y):
        eta=time if clock=='eta' else y[-1].real
        omega,prime=profile(eta,tau);H=hamiltonian(ref,omega,clock=clock,power=power)
        dH=derivative(ref,omega,prime,clock,power)
        return np.r_[-1j*(H@y[:dim]),complex(np.vdot(y[:dim],dH@y[:dim]).real),complex(1 if clock=='eta' else 1/omega)]
    span=(lo,hi) if clock=='eta' else (0.,float(cosmic_time(hi,tau)))
    sol=solve_ivp(rhs,span,y0,method='DOP853',rtol=tol,atol=tol/20)
    assert sol.success
    state=sol.y[:dim,-1];work=sol.y[-2,-1].real
    H0=hamiltonian(ref,omega0,clock=clock,power=power)
    Hf=hamiltonian(ref,profile(hi,tau)[0],clock=clock,power=power)
    energy_change=np.vdot(state,Hf@state).real-np.vdot(initial,H0@initial).real
    return initial,state,dict(work=float(work),energy_change=float(energy_change),
         work_balance_error=float(abs(energy_change-work)),norm_error=float(abs(np.vdot(state,state).real-1)),
         final_eta_error=float(abs(sol.y[-1,-1].real-hi)) if clock=='t' else 0.,evaluations=sol.nfev)

def split_step(psi,local,dt,diagonal,hop):
    state=psi*np.exp(-.5j*dt*diagonal)
    for parity,fraction in [(0,.5),(1,1.),(0,.5)]:
        c=np.cos(hop*dt*fraction);s=1j*np.sin(hop*dt*fraction)
        for i in range(parity,N-1,2):
            left,right=local['pairs'][i];u=state[left].copy();v=state[right].copy()
            state[left]=c*u+s*v;state[right]=s*u+c*v
    return state*np.exp(-.5j*dt*diagonal)

def expectation(psi,local,diagonal,hop):
    z=float(abs(psi)**2@diagonal)
    for left,right in local['pairs']:z-=2*hop*np.vdot(psi[left],psi[right]).real
    return float(z)

def circuit(local,initial,tau,steps,clock='eta',power=2,measure_work=True):
    lo=-10*tau;hi=10*tau
    start,end=(lo,hi) if clock=='eta' else (0.,float(cosmic_time(hi,tau)))
    dt=(end-start)/steps;mid=start+(np.arange(steps)+.5)*dt
    eta=mid if clock=='eta' else inverse_clock(mid,tau)
    omegas,primes=profile(eta,tau)
    md=M*(local['n']*(-1.)**np.arange(N)).sum(axis=1)
    ed=J*(local['e']**2).sum(axis=1)
    psi=initial.copy();work=0.
    for omega,prime in zip(omegas,primes):
        diag=omega*md+omega**power*ed;hop=W
        if clock=='t':diag=diag/omega;hop/=omega
        if measure_work:
            half=split_step(psi,local,.5*dt,diag,hop)
            if clock=='eta':d_diag=prime*(md+power*omega**(power-1)*ed);d_hop=0.
            else:d_diag=prime*(power-1)*omega**(power-3)*ed;d_hop=-W*prime/omega**3
            work+=dt*expectation(half,local,d_diag,d_hop)
        psi=split_step(psi,local,dt,diag,hop)
    def energy(state,eta):
        omega=profile(eta,tau)[0];diag=omega*md+omega**power*ed;hop=W
        if clock=='t':diag/=omega;hop/=omega
        return expectation(state,local,diag,hop)
    change=energy(psi,hi)-energy(initial,lo)
    return psi,dict(work=float(work),energy_change=float(change),work_balance_error=float(abs(work-change)))

def phase_error(a,b):
    return float(np.linalg.norm(a*np.exp(-1j*np.angle(np.vdot(b,a)))-b))

def observables(psi,ref,H):
    values,vectors=eigh(H)
    return dict(ground_state_survival=float(abs(np.vdot(vectors[:,0],psi))**2),
                excess_energy=float(np.vdot(psi,H@psi).real-values[0]),
                bare_even_site_occupation=float(abs(psi)**2@ref['n'][:,::2].sum(axis=1)),
                mean_electric_field=(abs(psi)**2@ref['e']).tolist())

def main():
    local=local_basis();ref=reference_operators();allrows=[];summaries=[]
    for tau in [.15,.6]:
        initial,exact,check=reference_evolution(ref,tau)
        _,cosmic,ccheck=reference_evolution(ref,tau,'t')
        init=to_local(initial,ref,local);target=to_local(exact,ref,local)
        Hf=hamiltonian(ref,profile(10*tau,tau)[0])
        summary=dict(tau=tau,reference_eta=check,reference_t=ccheck,reference_clock_error=phase_error(exact,cosmic),
            initial_bare_occupation=float(abs(initial)**2@ref['n'][:,::2].sum(axis=1)),
            final=observables(exact,ref,Hf))
        for steps in [128,256,512,1024,2048,4096]:
            states=[]
            for clock in ['eta','t']:
                state,work=circuit(local,init,tau,steps,clock)
                states.append(state)
                allrows.append(dict(tau=tau,clock=clock,steps=steps,state_error=phase_error(state,target),
                     norm_error=float(abs(np.vdot(state,state).real-1)),gauss_squared=float(abs(state)**2@np.sum(local['gauss']**2,axis=1)),**work))
            allrows[-1]['circuit_clock_error']=phase_error(*states)
        summaries.append(summary)
        print(json.dumps({'reference':summary,'finest':allrows[-2:]}),flush=True)
    # Deliberately incorrect electric metric exponents: still gauge invariant/unitary.
    wrong=[];tau=.6
    initial,correct,_=reference_evolution(ref,tau)
    for power in [0,1,2]:
        _,state,check=reference_evolution(ref,tau,power=power)
        _,cosmic,_=reference_evolution(ref,tau,'t',power=power)
        wrong.append(dict(electric_exponent=power,state_error_vs_correct=phase_error(state,correct),
             own_clock_agreement=phase_error(state,cosmic),reference=check,
             final=observables(state,ref,hamiltonian(ref,profile(10*tau,tau)[0]))))
    # Static curved geometry, flat limit and constant conformal-factor clock relation.
    x=(np.arange(N)-(N-1)/2)*A;xl=(x[:-1]+x[1:])/2
    om=1+.25*np.cos(2*np.pi*x/(N*A));ol=1+.25*np.cos(2*np.pi*xl/(N*A))
    bare=np.array([float(np.all(n==np.arange(N)%2)) for n in ref['n']],complex)
    init=to_local(bare,ref,local);H=hamiltonian(ref,om,ol)
    target=to_local(expm(-1.5j*H)@bare,ref,local);static=[]
    diagonal=M*(local['n']*((-1.)**np.arange(N))*om).sum(axis=1)+J*(local['e']**2*ol**2).sum(axis=1)
    for steps in [128,256,512,1024,2048]:
        psi=init.copy()
        for _ in range(steps):psi=split_step(psi,local,1.5/steps,diagonal,W)
        static.append(dict(steps=steps,state_error=phase_error(psi,target),
           energy_error=abs(expectation(psi,local,diagonal,W)-expectation(init,local,diagonal,W))))
    omega=1.3;H=hamiltonian(ref,omega)
    time_rescaling=phase_error(expm(-1j*.8*H)@bare,expm(-1j*(.8*omega)*(H/omega))@bare)
    # Static interacting vacuum: it must not be reported as newly created particles.
    H=hamiltonian(ref);vals,vecs=eigh(H);vac=vecs[:,0]
    static_vac=observables(expm(-1.5j*H)@vac,ref,H)
    # Full tensor-product check at L=1: no projection is used by any gate.
    full=local_basis(1,True);small=local_basis(1);smallref=reference_operators(1)
    iv,sv,_=reference_evolution(smallref,.15)
    fullstate,_=circuit(full,to_local(iv,smallref,full),.15,128,measure_work=False)
    smallstate,_=circuit(small,to_local(iv,smallref,small),.15,128,measure_work=False)
    lift=np.zeros(len(fullstate),complex);lift[small['codes']]=smallstate
    checks=dict(full_dimension=len(fullstate),complete_physical_dimension=len(local['n']),
        full_reduced_agreement=float(np.linalg.norm(fullstate-lift)),
        full_gauss_squared=float(abs(fullstate)**2@np.sum(full['gauss']**2,axis=1)),
        all_local_gate_gauss_change=max(float(np.max(abs(full['gauss'][i]-full['gauss'][j]))) for i,j in full['pairs']),
        constant_scale_clock_error=time_rescaling,static_vacuum=static_vac,
        max_electric_field_complete_sector=int(np.max(abs(local['e']))))
    assert checks['all_local_gate_gauss_change']==0 and checks['full_gauss_squared']==0
    assert checks['full_reduced_agreement']<1e-12
    assert time_rescaling<1e-12
    assert abs(static_vac['ground_state_survival']-1)<1e-12
    assert max(z['reference_clock_error'] for z in summaries)<2e-9
    assert max(z['norm_error'] for z in allrows)<2e-10
    assert max(z['gauss_squared'] for z in allrows)==0
    for tau in [.15,.6]:
        for clock in ['eta','t']:
            rows=[z for z in allrows if z['tau']==tau and z['clock']==clock]
            ratios=[rows[i]['state_error']/rows[i+1]['state_error'] for i in range(len(rows)-1)]
            assert all(3.8<r<4.2 for r in ratios[-3:]),(tau,clock,ratios)
    assert min(z['state_error_vs_correct'] for z in wrong[:2])>.1
    result=dict(scope=__doc__,parameters=dict(sites=N,w=W,m=M,J=J,lattice_spacing=A,
           physical_gauge_coupling=float(np.sqrt(2*J/A)),boundary_field=BOUNDARY,
           cutoff=3,initial_state='interacting ground state for expansion, bare state for static-curved transport'),
           expansion_summaries=summaries,rows=allrows,wrong_electric_exponents=wrong,
           static_curved=dict(omega_sites=om.tolist(),omega_links=ol.tolist(),rows=static),
           audit=checks,all_assertions_passed=True)
    Path(__file__).with_name('curved_gauge_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'wrong_models':wrong,'static':static[-1],'audit':checks},indent=2))

if __name__=='__main__':main()
