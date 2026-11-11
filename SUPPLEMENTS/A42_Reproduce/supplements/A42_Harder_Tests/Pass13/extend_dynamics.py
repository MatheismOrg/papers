"""Extend the saved spatial sweep without repeating completed calculations."""
from pathlib import Path
import json
import numpy as np
import lattice_scaling as s
root=Path(__file__).parent;path=root/'dynamics_results.json';r=json.loads(path.read_text())
for N in [18,20]:
    b=s.make_basis(N,boundary=1.)
    for improved in [False,True]:
        if any(z['N']==N and z['improved']==improved for z in r['rows']):continue
        z,state=s.evolve(b,improved=improved);r['rows'].append(z)
        print(json.dumps(z),flush=True)
        if N==20 and improved:
            tight,other=s.evolve(b,improved=True,tol=2e-10)
            diff=float(np.linalg.norm(state*np.exp(-1j*np.angle(np.vdot(other,state)))-other))
            r['checks'].append(dict(N=N,improved=True,tight_reference_state_difference=diff,
                survival_difference=abs(z['ground_survival']-tight['ground_survival']),
                excess_density_difference=abs(z['proper_excess_energy_density']-tight['proper_excess_energy_density']),tighter=tight))
            assert diff<2e-6
        r['largest_requested_N_completed']=N;path.write_text(json.dumps(r,indent=2)+'\n')
print('EXTENSION COMPLETE',flush=True)
