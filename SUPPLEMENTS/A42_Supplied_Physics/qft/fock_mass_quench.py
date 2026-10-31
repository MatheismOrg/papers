"""Explicit eight-mode fermionic Fock test of the uniform local Dirac circuit.

Two fixed low cell momenta, both regulator parity branches, four occupied
negative-energy orbitals. No interactions; no projection during evolution.
The complete N=4 sector has 70 states inside the 256-dimensional CAR space.
Python 3, NumPy, SciPy. Runtime is a few seconds.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
from scipy.linalg import block_diag, expm
from scipy.sparse import csr_matrix, eye

I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Z=np.diag([1,-1]).astype(complex)
MOMENTA=(1.,2.)
MI=.4
MF=3.
TIME=1.3
SUBSETS=list(combinations(range(8),4))
MASKS=np.array([sum(1<<i for i in rows) for rows in SUBSETS])
POSITIVE=(1,3,5,7)


def layers(theta,k):
    c=np.cos(theta);s=np.sin(theta)
    even=np.block([[c*I2,-s*Z],[s*Z,c*I2]])
    odd=np.block([[c*I2,s*np.exp(-1j*k)*Z],[-s*np.exp(1j*k)*Z,c*I2]])
    return even,odd


def circuit_block(n,p,m,t=TIME,calibrated=True):
    a=2*np.pi/n
    steps=int(np.ceil(t/a));steps+=steps%2
    dt=t/steps;nu=dt/a
    theta=np.arcsin(nu/2) if calibrated else nu/2
    eh,oh=layers(theta/2,2*a*p)
    ef,of=layers(theta,2*a*p)
    mass=np.kron(I2,np.cos(m*dt/2)*I2-1j*np.sin(m*dt/2)*X)
    first=mass@eh@of@eh@mass
    second=mass@oh@ef@oh@mass
    u=np.linalg.matrix_power(second@first,steps//2)
    return u,{'sites':n,'steps':steps,'dt':dt,'dt_over_dx':nu}


def continuum_basis(n,p,m):
    # Rephase odd-site amplitudes to their physical positions within the cell.
    a=2*np.pi/n
    d=np.diag([1,1,np.exp(1j*p*a),np.exp(1j*p*a)])
    columns=[]
    energies=[]
    for parity in [1.,-1.]:
        vals,vec=np.linalg.eigh(parity*p*Z+m*X)
        cell=np.array([1,parity],complex)/np.sqrt(2)
        for j in range(2):
            columns.append(d@np.kron(cell,vec[:,j]));energies.append(vals[j])
    return np.column_stack(columns),np.array(energies)


def exterior_state(orbitals):
    return np.array([np.linalg.det(orbitals[list(rows),:]) for rows in SUBSETS])


def exterior_unitary(u):
    # Independent full N=4-sector lift for one selected check.
    return np.array([[np.linalg.det(u[np.ix_(rows,cols)]) for cols in SUBSETS] for rows in SUBSETS])


def statistics(state):
    probabilities=abs(state)**2
    counts=np.array([sum((int(mask)>>j)&1 for j in POSITIVE) for mask in MASKS])
    holes=np.array([4-sum((int(mask)>>j)&1 for j in [0,2,4,6]) for mask in MASKS])
    distribution=np.array([probabilities[counts==j].sum() for j in range(5)])
    mean=float(counts@probabilities)
    return {'norm_error':float(abs(probabilities.sum()-1)),
            'mean_pairs':mean,'pair_variance':float((counts*counts)@probabilities-mean*mean),
            'vacuum_persistence':float(distribution[0]),
            'pair_distribution':distribution.tolist(),
            'particle_minus_hole_squared':float(((counts-holes)**2)@probabilities)}


def analytic_predictions(mi=MI,mf=MF):
    probabilities=[]
    for p in MOMENTA:
        probability=.5*(1-(p*p+mi*mf)/(np.hypot(p,mi)*np.hypot(p,mf)))
        probabilities.extend([probability,probability])
    distribution=np.array([1.])
    for b in probabilities:distribution=np.convolve(distribution,[1-b,b])
    return {'single_channel_pair_probabilities':probabilities,
            'mean_pairs':float(sum(probabilities)),
            'pair_variance':float(sum(b*(1-b) for b in probabilities)),
            'vacuum_persistence':float(distribution[0]),
            'pair_distribution':distribution.tolist()}


def run(n,calibrated=True,no_quench_control=False):
    init_blocks=[];final_bases=[];u_blocks=[];exact_blocks=[]
    for p in MOMENTA:
        bi,_=continuum_basis(n,p,MI);bf,ef=continuum_basis(n,p,MF)
        init_blocks.append(bi[:,[0,2]])
        final_bases.append(bf)
        u,meta=circuit_block(n,p,MI if no_quench_control else MF,calibrated=calibrated)
        u_blocks.append(u)
        exact_blocks.append((bf*np.exp(-1j*ef*TIME))@bf.conj().T)
    initial=block_diag(*init_blocks)
    out_basis=block_diag(*final_bases)
    u=block_diag(*u_blocks);exact=block_diag(*exact_blocks)
    numerical_orbitals=out_basis.conj().T@u@initial
    reference_orbitals=out_basis.conj().T@exact@initial
    state=exterior_state(numerical_orbitals)
    reference=exterior_state(reference_orbitals)
    phase=np.angle(np.vdot(reference,state))
    stats=statistics(state)
    analytic=analytic_predictions()
    correlation=numerical_orbitals@numerical_orbitals.conj().T
    correlation_reference=reference_orbitals@reference_orbitals.conj().T
    row={**meta,'calibrated':calibrated,'no_quench_control':no_quench_control,
         **stats,'phase_aligned_many_fermion_state_error':float(np.linalg.norm(state*np.exp(-1j*phase)-reference)),
         'one_body_correlation_error':float(np.linalg.norm(correlation-correlation_reference)),
         'mean_pair_error':float(abs(stats['mean_pairs']-analytic['mean_pairs'])),
         'vacuum_persistence_error':float(abs(stats['vacuum_persistence']-analytic['vacuum_persistence'])),
         'distribution_l1_error':float(np.linalg.norm(np.array(stats['pair_distribution'])-analytic['pair_distribution'],ord=1)),
         'one_body_unitarity_error':float(np.linalg.norm(u.conj().T@u-np.eye(8))),
         'projector_idempotence_error':float(np.linalg.norm(correlation@correlation-correlation))}
    if n==64 and calibrated and not no_quench_control:
        uout=out_basis.conj().T@u@out_basis
        lift=exterior_unitary(uout)
        initial_state=exterior_state(out_basis.conj().T@initial)
        row['explicit_fock_lift_agreement']=float(np.linalg.norm(lift@initial_state-state))
        row['explicit_fock_lift_unitarity']=float(np.linalg.norm(lift.conj().T@lift-np.eye(70)))
    assert row['norm_error']<1e-11
    return row


def car_check():
    operators=[]
    for j in range(8):
        rows=[];cols=[];data=[]
        for s in range(256):
            if s&(1<<j):
                rows.append(s^(1<<j));cols.append(s)
                data.append((-1)**((s&((1<<j)-1)).bit_count()))
        operators.append(csr_matrix((data,(rows,cols)),shape=(256,256),dtype=complex))
    maximum=0.
    for i,c in enumerate(operators):
        for j,d in enumerate(operators):
            for residual in [c@d+d@c,c@d.conj().T+d.conj().T@c-(eye(256) if i==j else 0)]:
                if residual.nnz:maximum=max(maximum,float(max(abs(residual.data))))
    assert maximum==0
    return {'fermionic_modes':8,'full_fock_dimension':256,'occupied_fermions':4,
            'occupied_sector_dimension':70,'maximum_CAR_residual':maximum}


def main():
    rows=[run(n) for n in [32,64,128,256,512,1024,2048]]
    bad=[run(n,calibrated=False) for n in [128,512,2048]]
    no_quench=run(2048,no_quench_control=True)
    ratios=[rows[j]['phase_aligned_many_fermion_state_error']/rows[j+1]['phase_aligned_many_fermion_state_error'] for j in range(len(rows)-1)]
    assert all(3.5<r<4.5 for r in ratios[-3:])
    output={'scope':'Free fermionic mass-quench control, two selected momentum sectors, both parity branches retained. CAR/Fock quantization supplied. No interactions, no electromagnetic Schwinger effect, no complete vacuum/UV regulator limit.',
            'mi':MI,'mf':MF,'time':TIME,'momenta':MOMENTA,
            'continuum_occupation_formula':'b_p=1/2 [1-(p^2+mi mf)/(sqrt(p^2+mi^2) sqrt(p^2+mf^2))]',
            'analytic':analytic_predictions(),'CAR':car_check(),'convergent_results':rows,
            'successive_many_fermion_state_error_ratios':ratios,
            'uncalibrated_control':bad,'no_quench_dynamics_control':no_quench,
            'all_assertions_passed':True}
    Path(__file__).with_name('fock_mass_quench_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
