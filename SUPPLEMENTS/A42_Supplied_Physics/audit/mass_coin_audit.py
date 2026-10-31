"""Independent homogeneous 4x4 cell audit of the exterior mass-coin sandwich."""
from pathlib import Path
import json
import numpy as np

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)


def unitary(a,nu,N,m,p):
    dt=nu*a;theta=np.arcsin(nu*N/2);k=2*a*p
    c=np.cos(theta/2);s=np.sin(theta/2)
    even=np.block([[c*I,-s*SZ],[s*SZ,c*I]])
    c=np.cos(theta);s=np.sin(theta)
    odd=np.block([[c*I,s*np.exp(-1j*k)*SZ],[-s*np.exp(1j*k)*SZ,c*I]])
    angle=m*N*dt/2
    mass=np.kron(I,np.cos(angle)*I-1j*np.sin(angle)*SX)
    return mass@even@odd@even@mass


def main():
    nu=.8;N=1.2;m=.7;p=.4
    alpha=np.kron(SX,SZ);beta=np.kron(I,SX)
    assert np.linalg.norm(alpha@beta+beta@alpha)==0
    target=N*(p*alpha+m*beta)
    projector=np.vstack([I,I])/np.sqrt(2)
    rows=[]
    for a in [.1,.05,.025,.0125,.00625]:
        dt=nu*a;u=unitary(a,nu,N,m,p)
        generator=(u.conj().T-u)/(2j*dt)
        rest=unitary(a,nu,N,m,0)
        expected=np.kron(I,np.cos(m*N*dt)*I-1j*np.sin(m*N*dt)*SX)
        resterr=float(np.linalg.norm(rest-expected,2))
        assert resterr<1e-14
        rows.append({'a':a,'dt':dt,'full_centered_generator_error':float(np.linalg.norm(generator-target,2)),
            'smooth_branch_centered_generator_error':float(np.linalg.norm(projector.conj().T@(generator-target)@projector,2)),
            'exact_rest_mass_gate_error':resterr})
    assert rows[-1]['full_centered_generator_error']<rows[0]['full_centered_generator_error']/10
    out={'scope':'Constant lapse, two-site cell and two internal spin components. No import of production code.',
         'placement':'Mass half-coins outside the entire kinetic palindrome.',
         'expected_H':'N[p sigma_x(sublattice) tensor sigma_z(spin) + m I tensor sigma_x(spin)]',
         'rows':rows,'all_assertions_passed':True}
    Path(__file__).with_name('mass_coin_audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
