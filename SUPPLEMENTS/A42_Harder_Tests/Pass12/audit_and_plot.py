"""Independent audits and figures for pass 12. No parameter fitting."""
from pathlib import Path
import importlib.util
import json
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import curved_gauge as c

root=Path(__file__).parent;r=json.loads((root/'curved_gauge_results.json').read_text())
ref=c.reference_operators();local=c.local_basis();tau=.6
initial,normal,_=c.reference_evolution(ref,tau)
_,tight,tightcheck=c.reference_evolution(ref,tau,tol=8e-14)
tolerance_error=c.phase_error(normal,tight)
# Compare independently assembled local generator with the Gauss-eliminated one.
omega_site=1+.2*np.sin(np.arange(c.N));omega_link=1+.15*np.cos(np.arange(c.N-1))
diag=c.M*(local['n']*((-1.)**np.arange(c.N))*omega_site).sum(axis=1)+c.J*(local['e']**2*omega_link**2).sum(axis=1)
H=np.diag(diag)
for left,right in local['pairs']:H[left,right]=-c.W;H[right,left]=-c.W
order=[np.flatnonzero(np.all(local['n']==n,axis=1))[0] for n in ref['n']]
generator_error=float(np.max(abs(H[np.ix_(order,order)]-c.hamiltonian(ref,omega_site,omega_link))))
# Check analytic work derivatives independently using centered differences.
derivatives=[]
for eta in [-.7,0.,.5]:
    omega,prime=c.profile(eta,tau)
    for clock in ['eta','t']:
        eps=2e-6
        if clock=='eta':lo,hi=eta-eps,eta+eps
        else:
            time=c.cosmic_time(eta,tau)
            lo,hi=c.inverse_clock(np.array([time-eps,time+eps]),tau)
        numeric=(c.hamiltonian(ref,c.profile(hi,tau)[0],clock=clock)-c.hamiltonian(ref,c.profile(lo,tau)[0],clock=clock))/(2*eps)
        error=float(np.max(abs(numeric-c.derivative(ref,omega,prime,clock))))
        derivatives.append(dict(eta=eta,clock=clock,max_entry_error=error))
# Flat limit must reproduce the preceding pass, including all hopping signs.
previous_path=root.parent/'toe_pass11/gauge_feedback.py'
if not previous_path.exists():
    previous_path=root/'gauge_feedback_pass11.py'
if not previous_path.exists():
    raise FileNotFoundError('The cross-pass audit needs the included gauge_feedback_pass11.py.')
spec=importlib.util.spec_from_file_location('pass11',previous_path)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
old_n,_,_,old_H=old.reference(3)
reorder=[np.flatnonzero(np.all(old_n==n,axis=1))[0] for n in ref['n']]
flat_error=float(np.max(abs(old_H[np.ix_(reorder,reorder)]-c.hamiltonian(ref))))
# Sudden and slow expansion: actual initial interacting ground state throughout.
rates=[]
for width in [.001,.15,.6,2.,6.]:
    iv,sv,checks=c.reference_evolution(ref,width)
    Hf=c.hamiltonian(ref,c.profile(10*width,width)[0])
    vals,vecs=eigh(Hf)
    sudden_survival=float(abs(np.vdot(vecs[:,0],iv))**2)
    rates.append(dict(tau=width,**c.observables(sv,ref,Hf),reference_work_error=checks['work_balance_error']))

ratios=[]
for width in [.15,.6]:
    for clock in ['eta','t']:
        rows=[z for z in r['rows'] if z['tau']==width and z['clock']==clock]
        ratios.append(dict(tau=width,clock=clock,
            state_error_ratios=[rows[i]['state_error']/rows[i+1]['state_error'] for i in range(len(rows)-1)],
            work_error_ratios=[rows[i]['work_balance_error']/rows[i+1]['work_balance_error'] for i in range(len(rows)-1)]))
assert tolerance_error<1e-10
assert generator_error<1e-13 and flat_error<1e-13
assert max(z['max_entry_error'] for z in derivatives)<1e-7
assert abs(rates[0]['ground_state_survival']-sudden_survival)<1e-4
assert rates[-1]['ground_state_survival']>rates[1]['ground_state_survival']
assert max(z['reference_work_error'] for z in rates)<1e-8
audit=dict(tight_reference_state_difference=tolerance_error,tight_reference=tightcheck,
    local_vs_independent_generator_max_error=generator_error,flat_limit_vs_pass11_max_error=flat_error,
    energy_derivative_checks=derivatives,refinement_ratios=ratios,
    sudden_limit_survival=sudden_survival,expansion_rate_study=rates,all_assertions_passed=True)
(root/'audit_results.json').write_text(json.dumps(audit,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
for clock,label,marker in [('eta','Conformal clock','o'),('t','Physical clock','s')]:
    rows=[z for z in r['rows'] if z['tau']==.6 and z['clock']==clock]
    axes[0].loglog([z['steps'] for z in rows],[z['state_error'] for z in rows],marker+'-',label=label)
axes[0].set_xlabel('Number of time steps');axes[0].set_ylabel('Phase-aligned state error')
axes[0].set_title('Correct coupling converges');axes[0].legend(fontsize=9);axes[0].grid(alpha=.2,which='both')
wrong=r['wrong_electric_exponents'];colors=['#bd5748','#c98f45','#367e96']
axes[1].bar(['No factor','Factor Ω','Factor Ω²'],[z['final']['ground_state_survival'] for z in wrong],color=colors,width=.65)
axes[1].set_ylabel('Final ground-state probability');axes[1].set_ylim(0,.6)
axes[1].set_title('Charge and clock tests pass all three')
for i,z in enumerate(wrong):axes[1].text(i,z['final']['ground_state_survival']+.015,f"{z['final']['ground_state_survival']:.4f}",ha='center',fontsize=10)
axes[1].text(.5,.94,'Same initial state; τ = 0.6',transform=axes[1].transAxes,ha='center',fontsize=9,color='#5a6670')
fig.savefig(root/'curved_gauge_diagnostics.png',dpi=190)
print(json.dumps(audit,indent=2))
