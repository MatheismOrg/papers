"""Fixed-box refinement of the pass-12 expansion, with/without mass improvement."""
from pathlib import Path
import json
import numpy as np
import lattice_scaling as s

def main():
    root=Path(__file__).parent;rows=[];checks=[]
    for N in [6,8,10,12,14,16]:
        b=s.make_basis(N,boundary=1.)
        for improved in [False,True]:
            z,state=s.evolve(b,improved=improved);rows.append(z);print(json.dumps(z),flush=True)
            if N in [6,16]:
                tighter,other=s.evolve(b,improved=improved,tol=1e-11)
                diff=float(np.linalg.norm(state*np.exp(-1j*np.angle(np.vdot(other,state)))-other))
                checks.append(dict(N=N,improved=improved,tight_reference_state_difference=diff,
                    survival_difference=abs(z['ground_survival']-tighter['ground_survival']),
                    excess_density_difference=abs(z['proper_excess_energy_density']-tighter['proper_excess_energy_density']),
                    tighter=tighter))
                assert diff<2e-6
        (root/'dynamics_results.json').write_text(json.dumps(dict(rows=rows,checks=checks,complete=False),indent=2)+'\n')
    # Increase the box at fixed a=0.5, separating volume from spacing changes.
    volume=[]
    for N in [10,14,16]:
        b=s.make_basis(N,boundary=1.)
        for improved in [False,True]:
            z,_=s.evolve(b,length=.5*(N+1),improved=improved);volume.append(z)
            print(json.dumps(dict(volume_test=z)),flush=True)
    result=dict(parameters=dict(fixed_wall_separation=3.5,boundary_flux=1,mass=.4,gauge_coupling=float(np.sqrt(2.88)),tau=.6,
            prescription='a=L/(N+1), m_eta=m*Omega-g^2*a*Omega^2/8 when improved'),
            rows=rows,checks=checks,volume_rows=volume,complete=True)
    (root/'dynamics_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
