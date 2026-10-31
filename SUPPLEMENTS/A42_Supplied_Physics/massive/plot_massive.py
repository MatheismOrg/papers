"""Static figure from completed massive-packet results."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
d=json.loads((root/'massive_curved_results.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(11.5,4.3),constrained_layout=True)
c=d['plot_samples']['correct'];w=d['plot_samples']['wrong']
axs[0].plot(c['x'],c['reference_density'],color='#172b41',lw=2.5,label='Fourier reference, mN coupling')
axs[0].plot(c['x'],c['density'],'--',color='#20799e',lw=1.7,label='Local circuit, mN coupling')
axs[0].plot(w['x'],w['density'],':',color='#d77029',lw=2,label='Same transport, wrong m coupling')
axs[0].set(xlabel='Position x',ylabel='Probability density',title='Massive packet, m=2, time 1\nSupplied curved metric, 512 sites')
axs[0].legend(frameon=False,fontsize=8.3)
for case,color,label in [(d['cases'][0],'#754197','Correct mN, mass 0.5'),(d['cases'][1],'#20799e','Correct mN, mass 2')]:
    rows=[x for x in case['results'] if not x['wrong_mass']]
    axs[1].loglog([x['sites'] for x in rows],[x['state_error'] for x in rows],'o-',color=color,label=label)
case=d['cases'][1];rows=[x for x in case['results'] if x['wrong_mass']]
axs[1].loglog([x['sites'] for x in rows],[x['state_error'] for x in rows],'o-',color='#d77029',label='Wrong m coupling, mass 2')
axs[1].axhline(case['wrong_mass_continuum_difference']['state_error'],ls='--',color='#a97043',lw=1,label='Wrong-coupling continuum difference')
axs[1].set(xlabel='Number of spatial sites',ylabel='Full spinor state norm error',title='Correct mass coupling converges\nThe same light cone is insufficient')
axs[1].legend(frameon=False,fontsize=8.1)
for ax in axs:
    ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
fig.savefig(root/'massive_curved.png',dpi=180)
print(root/'massive_curved.png')
