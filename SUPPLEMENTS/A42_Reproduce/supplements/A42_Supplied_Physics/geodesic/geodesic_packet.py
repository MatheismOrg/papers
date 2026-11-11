"""Supplied static geometry: massive Dirac packets compared with timelike rays.

The quantum reference is an FFT Hermitian spectral Hamiltonian, independent of
the candidate's local bond circuit. A simultaneous narrow/semiclassical family
sigma=0.8/sqrt(m) is used; finite-width quantum packets need not follow exact
point-particle paths at finite m. NumPy and SciPy, no quantum computer.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import LinearOperator, expm_multiply


def initial_state(n,m):
    x=np.arange(n)*2*np.pi/n;dx=2*np.pi/n
    k=np.fft.fftfreq(n,d=1/n)
    sigma=.8/np.sqrt(m);center=np.pi/2
    distance=(x-center+np.pi)%(2*np.pi)-np.pi
    envelope=np.exp(-distance**2/(4*sigma**2)).astype(complex)
    # Positive-energy projection in local constant-coefficient Dirac basis.
    # p0=0, beta+ seed; this fixes the preparation, without fitting the result.
    seed=np.stack([envelope,envelope])/np.sqrt(2)
    sp=np.fft.fft(seed,axis=1,norm='ortho')
    e=np.sqrt(k*k+m*m)
    hp=np.stack([k*sp[0]+m*sp[1],m*sp[0]-k*sp[1]])
    sp=.5*(sp+hp/e)
    initial=np.fft.ifft(sp,axis=1,norm='ortho')
    initial/=np.linalg.norm(initial)
    return initial


def quantum(n,m,mass_lapse=True,times=(0.,.25,.5,.75,1.)):
    x=np.arange(n)*2*np.pi/n;dx=2*np.pi/n
    lapse=1+.25*np.cos(x);k=np.fft.fftfreq(n,d=1/n)
    sigma=.8/np.sqrt(m);center=np.pi/2
    distance=(x-center+np.pi)%(2*np.pi)-np.pi
    initial=initial_state(n,m)
    coefficient=m*lapse if mass_lapse else np.full(n,m)

    def h(v):
        arr=np.asarray(v).reshape(2,n,-1)
        def p(z):return np.fft.ifft(k[:,None]*np.fft.fft(z,axis=0),axis=0)
        kinetic=.5*(lapse[:,None]*p(arr[0])+p(lapse[:,None]*arr[0]))
        kinetic1=-.5*(lapse[:,None]*p(arr[1])+p(lapse[:,None]*arr[1]))
        out=np.stack([kinetic+coefficient[:,None]*arr[1],kinetic1+coefficient[:,None]*arr[0]])
        return out.reshape(np.asarray(v).shape)
    A=LinearOperator((2*n,2*n),matvec=lambda v:-1j*h(v),rmatvec=lambda v:1j*h(v),
                     matmat=lambda v:-1j*h(v),rmatmat=lambda v:1j*h(v),dtype=np.complex128)
    assert np.allclose(np.diff(times),times[1]-times[0])
    states=expm_multiply(A,initial.ravel(),start=times[0],stop=times[-1],num=len(times),traceA=0.)
    rows=[]
    for t,v in zip(times,states):
        psi=v.reshape(2,n);rho=np.sum(abs(psi)**2,axis=0)
        # Packets are remote from the periodic seam; local unwrapping explicit.
        xpos=center+distance
        mean=float(np.sum(xpos*rho));width=float(np.sqrt(np.sum((xpos-mean)**2*rho)))
        rows.append({'time':float(t),'mean_x':mean,'rms_width':width,
                     'norm_error':float(abs(np.linalg.norm(v)-1)),
                     'probability_near_periodic_seam':float(np.sum(rho[(x<.2)|(x>2*np.pi-.2)]))})
    return {'sites':n,'mass':m,'mass_uses_lapse':mass_lapse,'initial_sigma_parameter':sigma,
            'initial_state':'local constant-coefficient positive-energy projection of zero-momentum beta+ Gaussian',
            'trajectory':rows},states[-1]


def ray(times):
    # eta=asinh(p/m); v_local=tanh(eta). Mass cancels from these equations.
    def rhs(t,y):
        x,eta=y;return [(1+.25*np.cos(x))*np.tanh(eta),.25*np.sin(x)]
    sol=solve_ivp(rhs,(0,float(max(times))),[np.pi/2,0.],t_eval=times,rtol=2e-12,atol=2e-13)
    energy=(1+.25*np.cos(sol.y[0]))*np.cosh(sol.y[1])
    return [{'time':float(t),'x':float(x),'rapidity':float(e),'energy_over_mass':float(en)}
            for t,x,e,en in zip(sol.t,sol.y[0],sol.y[1],energy)]


def main():
    times=np.linspace(0,1,5);rays=ray(times);rows=[]
    for mass in [8.,16.,32.,64.]:
        row,_=quantum(1024,mass,True,times)
        row['mean_x_minus_point_ray']=[q['mean_x']-r['x'] for q,r in zip(row['trajectory'],rays)]
        rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='initial_state'}),flush=True)
    bad,_=quantum(1024,64.,False,times)
    coarse,cstate=quantum(512,64.,True,times)
    fine,fstate=quantum(2048,64.,True,times)
    # Compare sampled fine spectral reference, with grid-normalization factor.
    restricted=fstate.reshape(2,2048)[:,::4]*2
    refinement=float(np.linalg.norm(cstate.reshape(2,512)-restricted))
    spectral_mean_difference=abs(coarse['trajectory'][-1]['mean_x']-fine['trajectory'][-1]['mean_x'])
    assert refinement<1e-9
    assert max(r['norm_error'] for z in rows+[bad,coarse,fine] for r in z['trajectory'])<1e-9
    out={'scope':'Quantum reference and semiclassical timelike-ray control on supplied static curved background; no local circuit result in this file.',
         'H':'-i sigma_z (N partial_x+Nprime/2)+m b(x) sigma_x; N=1+.25cos(x), b=N correct, b=1 negativecontrol',
         'ray_equations':'xdot=N tanh(eta), etadot=-Nprime, eta=asinh(p/m)',
         'ray':rays,'quantum_packets':rows,'wrong_mass_coupling':bad,
         'spectral_refinement':{'mass':64.,'coarse_n':512,'fine_n':2048,'full_state_error':refinement,'mean_x_difference':spectral_mean_difference},
         'all_assertions_passed':True}
    Path(__file__).with_name('geodesic_packet_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'spectral_refinement':out['spectral_refinement'],'ray_final':rays[-1],'wrong_mass_final':bad['trajectory'][-1]},indent=2),flush=True)


if __name__=='__main__':main()
