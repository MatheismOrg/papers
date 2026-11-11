"""Independent eigenvector-overlap check of the quench occupation formula."""
from pathlib import Path
import json
import numpy as np

X=np.array([[0,1],[1,0]],complex)
Z=np.diag([1,-1]).astype(complex)
errors=[]
for p in [1.,2.]:
    for sign in [-1,1]:
        _,vi=np.linalg.eigh(sign*p*Z+.4*X)
        _,vf=np.linalg.eigh(sign*p*Z+3.*X)
        measured=abs(np.vdot(vf[:,1],vi[:,0]))**2
        target=.5*(1-(p*p+.4*3.)/(np.hypot(p,.4)*np.hypot(p,3.)))
        errors.append(float(abs(measured-target)))
assert max(errors)<1e-14
r={'analytic_quench_overlap_max_error':max(errors),'both_parity_branches_checked':True,'all_assertions_passed':True}
Path(__file__).with_name('quench_overlap_audit_results.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
