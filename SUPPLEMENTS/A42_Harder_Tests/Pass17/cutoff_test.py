"""Fixed-volume bare-vacuum cutoff diagnostic; not renormalized predictions."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
rows=[]
m=.4;L=4.
for n in [4,8,16,32,64]:
 eps=L/n
 one=4*np.sin(np.pi*np.arange(n)/n)**2/eps**2
 K=(one[:,None,None]+one[None,:,None]+one[None,None,:]).ravel()
 omega=np.sqrt(K+m*m)
 rho=float(omega.sum()/(2*L**3));p=float((K/omega).sum()/(6*L**3))
 rows.append(dict(n=n,spacing=eps,rho_vacuum=rho,p_vacuum=p,
                  H_initial=np.sqrt(rho/3),w=p/rho,scaled_rho=rho*eps**4,
                  Lambda_cancel_initial_rho=-rho,pressure_after_constant_cancel=p+rho))
r=dict(mass=m,length=L,volume=L**3,rows=rows,
       conclusion='Bare vacuum grows approximately epsilon^-4. A constant Lambda chosen to cancel rho at a=1 leaves p+rho, not zero pressure. No renormalized continuum claim.')
(ROOT/'cutoff_results.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
