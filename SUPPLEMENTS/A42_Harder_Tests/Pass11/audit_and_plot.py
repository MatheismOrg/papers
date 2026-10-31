from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gauge_feedback as g

root=Path(__file__).parent;r=json.load(open(root/'gauge_feedback_results.json'))
b=g.basis(3);state=g.circuit(b,4096)
n,e,ref,H=g.reference(3)
ordered=[]
for nn,ee in zip(n,e):
 d=np.empty(2*g.N-1,int);d[::2]=nn;d[1::2]=ee+3
 ordered.append(state[b['index'][int(np.ravel_multi_index(d,b['dims']))]])
ordered=np.array(ordered)
i=np.flatnonzero(np.all(n==np.arange(g.N)%2,axis=1))[0]
direct=expm(-1j*g.T*H)[:,i]
initial_energy=float(H[i,i]);final_energy=float(np.vdot(ordered,H@ordered).real)
checks={'eigen_vs_matrix_exponential_error':float(np.linalg.norm(ref-direct)),
 'initial_total_energy':initial_energy,'circuit_total_energy':final_energy,
 'energy_absolute_error':abs(final_energy-initial_energy),
 'resolved_state_error':g.phase_error(ordered,ref),
 'pair_mean_relative_cutoff1_error':abs(r['cutoff_study'][0]['pair_mean']/r['cutoff_study'][-1]['pair_mean']-1)}
assert checks['eigen_vs_matrix_exponential_error']<1e-12
assert checks['energy_absolute_error']<1e-6
assert max(x['norm_error'] for x in r['resolved_cutoff3_rows'])<1e-10
assert max(x['gauss_squared'] for x in r['resolved_cutoff3_rows'])==0
z=r['resolved_cutoff3_rows'];ratios=[z[j]['state_error']/z[j+1]['state_error'] for j in range(len(z)-1)]
assert all(3.9<x<4.1 for x in ratios[-3:])
checks['resolved_refinement_ratios']=ratios;checks['all_assertions_passed']=True
(root/'audit_results.json').write_text(json.dumps(checks,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(7.2,3.4),layout='constrained')
ax.loglog([x['steps'] for x in r['rows']],[x['state_error'] for x in r['rows']],'o-',label='Time-step error, cutoff 1')
ax.loglog([x['steps'] for x in z],[x['state_error'] for x in z],'s-',label='Time-step error, cutoff 3')
ax.axhline(r['cutoff_study'][0]['state_error_vs_cutoff4'],color='#a43b31',ls='--',label='Cutoff-1 discrepancy from complete sector')
ax.set_xlabel('Number of time steps over T = 1.5');ax.set_ylabel('Phase-aligned state error')
ax.grid(alpha=.18,which='both');ax.legend(loc='lower left',fontsize=8.7)
fig.savefig(root/'error_budget.png',dpi=190)
print(json.dumps(checks,indent=2))
