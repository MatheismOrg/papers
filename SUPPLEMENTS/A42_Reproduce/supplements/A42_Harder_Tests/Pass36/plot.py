import json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import Fp

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results.json').read_text());data=np.load(ROOT/'trajectories.npz')
fig,axs=plt.subplots(1,3,figsize=(13.4,4.3),layout='constrained')
for row,color,label in zip(r['fixed_dust_comparison'],['#666666','#126e91'],['Classical gravity model','With 4D quantum stress']):
    key=row['label'];t=data[key+'__t']-row['bounce_time'];a=data[key+'__a']
    axs[0].plot(t,a,color=color,label=label,lw=2)
axs[0].set(xlabel=r'$(t-t_{bounce})/\sqrt{\alpha}$',ylabel=r'$a/\sqrt{\alpha}$',title='Same dust constant; both evolve through a bounce')
axs[0].legend(fontsize=8,frameon=False)
key=r['fixed_dust_comparison'][1]['label'];row=r['fixed_dust_comparison'][1]
t=data[key+'__t']-row['bounce_time']
for name,color,label in [('rho_q','#126e91',r'$\rho_q$'),('p_q','#b34a29',r'$p_q$, including angular pressure')]:
    axs[1].plot(t,data[key+'__'+name]*8*math.pi/3,color=color,label=label,lw=2)
axs[1].set(xlabel=r'$(t-t_{bounce})/\sqrt{\alpha}$',ylabel=r'$(8\pi G\alpha/3)\,T$',title='Quantum density and pressure remain finite')
axs[1].legend(fontsize=8,frameon=False)
y=np.linspace(0,.82,500)
for lam,color,label in [(r['free_SM_proxy']['lambda_at_alpha_planck_area'],'#126e91','Free SM-count proxy'),(1.6,'#b34a29','Excessive-coupling control')]:
    axs[2].plot(y,Fp(y,lam),label=label,color=color,lw=2)
axs[2].axhline(0,color='#333',lw=.8)
axs[2].set(xlabel=r'$y=\alpha\Phi$',ylabel=r'$F_\lambda^\prime(y)$',ylim=(-.25,3.5),title='Loss of response occurs above the bound')
axs[2].legend(fontsize=8,frameon=False)
for ax in axs:ax.grid(alpha=.18)
fig.suptitle('A42 Pass36: a homogeneous quantum-backreaction benchmark',fontsize=15)
fig.text(.01,-.035,'Free massless conformal fields; specified state and zero total Box R coefficient. No inhomogeneous black-hole claim.',fontsize=9)
fig.savefig(ROOT/'quantum_backreaction.png',dpi=170,bbox_inches='tight')
