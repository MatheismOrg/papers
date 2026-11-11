"""Massive two-spinor transport on a supplied static curved background.

Independent reference: Hermitian Fourier pseudospectral generator evaluated
by scipy.sparse.linalg.expm_multiply, with reference refinement and exact
flat-space mode checks. Candidate: calibrated alternating bond rotations,
each microstep sandwiched by the supplied m*N(x) mass coin.
No fit, dynamical gravity, QFT, or backreaction is claimed.
"""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def geometry(n,amplitude):
    dx=2*np.pi/n;x=np.arange(n)*dx
    N=1+amplitude*np.cos(x);Nmid=1+amplitude*np.cos(x+dx/2)
    root=np.sqrt(1-amplitude**2);Y=2*np.pi/root
    y=2/root*np.arctan2(np.sqrt((1-amplitude)/(1+amplitude))*np.sin(x/2),np.cos(x/2))
    return dx,x,N,Nmid,y,Y


def initial(n,amplitude,center):
    dx,x,N,Nmid,y,Y=geometry(n,amplitude)
    modes=np.arange(1,10);freq=2*np.pi*modes/Y
    coeff=np.exp(-.5*((modes-5)/1.5)**2)
    coeff/=np.sqrt(Y*np.sum(coeff**2))
    f=np.sum(coeff[:,None]*np.exp(1j*freq[:,None]*(y[None,:]-center)),axis=0)/np.sqrt(N)
    spin=np.array([1.,.35+.2j]);spin/=np.linalg.norm(spin)
    psi=np.sqrt(dx)*f[:,None]*spin[None,:]
    return psi/np.linalg.norm(psi)


def spectral(n,amplitude,mass,center,T=1.,wrong_mass=False):
    k=np.rint(np.fft.fftfreq(n,d=1/n)).astype(int)
    ii=np.arange(n)
    # Multiplication by 1+a*cos(x) in the finite Fourier basis.
    N=sparse.eye(n,format='csr')+sparse.csr_matrix((np.full(2*n,amplitude/2),
        (np.r_[ii,ii],np.r_[(ii+1)%n,(ii-1)%n])),shape=(n,n))
    P=sparse.diags(k,format='csr')
    A=(N@P+P@N)/2
    M=mass*(sparse.eye(n,format='csr') if wrong_mass else N)
    h=sparse.bmat([[A,M],[M,-A]],format='csr')
    defect=h-h.conj().T
    assert defect.nnz==0 or np.max(abs(defect.data))<1e-14
    psi0=initial(n,amplitude,center)
    v0=np.fft.fft(psi0,axis=0,norm='ortho').T.ravel()
    v=expm_multiply(-1j*T*h,v0,traceA=0.)
    coeff=v.reshape(2,n).T
    psi=np.fft.ifft(coeff,axis=0,norm='ortho')
    flat_error=None
    if amplitude==0:
        init=v0.reshape(2,n).T;E=np.sqrt(k*k+mass*mass)
        exact=np.cos(T*E)[:,None]*init
        action=np.column_stack([k*init[:,0]+mass*init[:,1],mass*init[:,0]-k*init[:,1]])
        exact-=1j*(np.sin(T*E)/E)[:,None]*action
        flat_error=float(np.linalg.norm(coeff-exact))
        assert flat_error<1e-10
    return psi,{'sites':n,'mass':mass,'lapse_amplitude':amplitude,'center':center,
                'wrong_mass':wrong_mass,'norm_error':float(abs(np.linalg.norm(psi)-1)),
                'high_frequency_probability':float(np.sum(abs(coeff[abs(k)>.4*n])**2)),
                'flat_exact_mode_error':flat_error}


def mass_coin(psi,c,s):
    u=psi[:,0].copy();v=psi[:,1].copy()
    psi[:,0]=c*u-1j*s*v;psi[:,1]=c*v-1j*s*u


def bond(psi,c,s,j,k):
    u=psi[j].copy();v=psi[k].copy()
    z=np.array([1.,-1.])
    psi[j]=c[:,None]*u-s[:,None]*v*z
    psi[k]=s[:,None]*u*z+c[:,None]*v


def local_circuit(n,amplitude,mass,center,T=1.,wrong_mass=False,reverse_check=False,initial_state=None):
    dx,x,N,Nmid,y,Y=geometry(n,amplitude)
    steps=int(np.ceil(T/dx));steps+=steps%2;dt=T/steps
    psi=(initial(n,amplitude,center) if initial_state is None
         else np.asarray(initial_state,dtype=complex).copy())
    assert psi.shape==(n,2) and np.all(np.isfinite(psi))
    assert abs(np.linalg.norm(psi)-1)<1e-10
    start=psi.copy()
    even=np.arange(0,n,2);odd=np.arange(1,n,2)
    theta=np.arcsin(dt*Nmid/(2*dx))
    cache={}
    for parity,indices in enumerate([even,odd]):
        for fraction in [.5,1.]:
            angle=fraction*theta[indices]
            cache[parity,fraction]=(np.cos(angle),np.sin(angle),indices,(indices+1)%n)
    mass_angle=.5*dt*mass*(np.ones(n) if wrong_mass else N)
    cm,sm=np.cos(mass_angle),np.sin(mass_angle)
    def microstep(index,backward=False):
        sign=-1 if backward else 1
        mass_coin(psi,cm,sign*sm)
        parity=index%2
        for p,f in [(parity,.5),(1-parity,1.),(parity,.5)]:
            c,s,j,k=cache[p,f];bond(psi,c,sign*s,j,k)
        mass_coin(psi,cm,sign*sm)
    for step in range(steps):microstep(step)
    result=psi.copy();inverse_error=None
    if reverse_check:
        for step in reversed(range(steps)):microstep(step,True)
        inverse_error=float(np.linalg.norm(psi-start))
        assert inverse_error<1e-10
    return result,{'sites':n,'steps':steps,'dt':dt,'dt_over_dx':dt/dx,
                   'norm_error':float(abs(np.linalg.norm(result)-1)),
                   'inverse_error':inverse_error}


def downsample(psi,n):
    assert len(psi)%n==0
    factor=len(psi)//n
    return psi[::factor]*np.sqrt(factor)


def compare(state,target):
    overlap=np.vdot(target,state)
    return {'state_error':float(np.linalg.norm(state-target)),
            'phase_aligned_state_error':float(np.linalg.norm(state*np.exp(-1j*np.angle(overlap))-target)),
            'density_L1_error':float(np.sum(abs(np.sum(abs(state)**2,axis=1)-np.sum(abs(target)**2,axis=1)))),
            'chirality_probabilities':np.sum(abs(state)**2,axis=0).tolist(),
            'reference_chirality_probabilities':np.sum(abs(target)**2,axis=0).tolist()}


def main():
    cases=[{'id':'curved_m0p5_c0p9','a':.25,'m':.5,'center':.9},
           {'id':'curved_m2_c0p9','a':.25,'m':2.,'center':.9},
           {'id':'curved_m2_c3p5','a':.25,'m':2.,'center':3.5},
           {'id':'flat_m2_c0p9','a':0.,'m':2.,'center':.9}]
    results=[];reference_audits=[];plots={}
    for case in cases:
        refs={};audits=[]
        for n in [512,1024,2048]:
            refs[n],audit=spectral(n,case['a'],case['m'],case['center']);audits.append(audit)
        audit_errors=[float(np.linalg.norm(refs[n]-downsample(refs[2048],n))) for n in [512,1024]]
        assert max(audit_errors)<1e-9
        reference_audits.append({'case':case['id'],'runs':audits,'reference_errors_against_2048':audit_errors})
        rows=[]
        for wrong in ([False,True] if case['a'] else [False]):
            for n in [64,128,256,512,1024]:
                psi,record=local_circuit(n,case['a'],case['m'],case['center'],wrong_mass=wrong,reverse_check=(n==128 and not wrong))
                target=downsample(refs[2048],n)
                record.update(compare(psi,target));record['wrong_mass']=wrong
                rows.append(record)
                if case['id']=='curved_m2_c0p9' and n==512:
                    key='wrong' if wrong else 'correct'
                    plots[key]={'x':(np.arange(n)*2*np.pi/n).tolist(),
                                'density':(np.sum(abs(psi)**2,axis=1)*n/(2*np.pi)).tolist(),
                                'reference_density':(np.sum(abs(target)**2,axis=1)*n/(2*np.pi)).tolist()}
        correct=[r for r in rows if not r['wrong_mass']]
        ratios=[correct[i]['state_error']/correct[i+1]['state_error'] for i in range(len(correct)-1)]
        assert all(3.5<ratio<4.6 for ratio in ratios[1:])
        floor=None
        if case['a']:
            wrong_ref,wrong_audit=spectral(2048,case['a'],case['m'],case['center'],wrong_mass=True)
            floor=compare(wrong_ref,refs[2048]);floor['reference_audit']=wrong_audit
            assert floor['state_error']>.01
        item={'case':case,'reference_sites':2048,'results':rows,'correct_error_ratios':ratios,
              'wrong_mass_continuum_difference':floor}
        results.append(item)
        print(json.dumps({'case':case['id'],'correct_errors':[(r['sites'],r['state_error']) for r in correct],
                          'ratios':ratios,'wrong_mass_floor':None if floor is None else floor['state_error']}),flush=True)
    assert max(r['norm_error'] for case in results for r in case['results'])<1e-10
    out={'scope':'Supplied 1+1 static metric ds²=-N(x)²dt²+dx². Massive first-quantized Dirac two-spinor; no field quantization, source backreaction, dynamical metric or ToE.',
         'target_H':'-i sigma_z [N partial_x + Nprime/2] + m N sigma_x',
         'local_rule':'Pass7 arcsine-calibrated alternating even/odd bond microsteps, sandwiched at every microstep by exp(-i dt m N sigma_x/2).',
         'reference':'Independent Hermitian Fourier pseudospectral H=(N P+P N) sigma_z/2 + mN sigma_x; scipy expm_multiply. Both numerical methods have a cutoff and are tested on smooth low-frequency packets.',
         'reference_audits':reference_audits,'cases':results,'plot_samples':plots,'all_assertions_passed':True}
    Path(__file__).with_name('massive_curved_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
