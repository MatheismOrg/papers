from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;B=R
j=json.loads((B/'dm_unified_repair/results/results.json').read_text())['models_by_dataset']['augmented']['models']
k=json.loads((B/'dm_unified_regular/results/results.json').read_text())['models_by_dataset']['augmented']['models']['mixed_regular']
j=dict(j,mixed_regular=k)
keys=['RAR','power_sum','mixed_gate','mixed_gate_range','mixed_regular']
labels=['Common RAR','Earlier power sum','Mass + potential, gated','Squared outer cutoff','Regularized shared form']
fig,axes=plt.subplots(1,2,figsize=(10.5,4.3),gridspec_kw={'width_ratios':[1.35,1]})
ys=np.arange(len(keys));h=.34
for off,section,color,label in [(-h/2,'cv','#24627a','Whole systems excluded'),(h/2,'population','#b67843','Whole populations excluded')]:
 vals=[j[z]['CV']['20261027']['metrics']['B_dex'] if section=='cv' else j[z]['population_holdout']['metrics']['B_dex'] for z in keys]
 axes[0].barh(ys+off,vals,h,color=color,label=label)
 for yy,val in zip(ys+off,vals):axes[0].text(val+.004,yy,f'{val:.3f}',va='center',fontsize=8)
axes[0].set_yticks(ys,labels);axes[0].set_ylim(len(keys)+.3,-.6);axes[0].set_xlim(0,.41);axes[0].set_xlabel('Balanced log acceleration RMS (dex)');axes[0].set_title('Shared fit across 286 systems',loc='left',fontsize=11,fontweight='bold');axes[0].legend(loc='lower left',fontsize=7,frameon=False)
for i,(tag,label) in enumerate([('power_sum','Earlier power sum'),('mixed_regular','Regularized shared form')]):
 vals=[j[tag]['CV']['20261027']['metrics']['populations'][p]['group_RMS_dex'] for p in ['SPARC','dSph','ETG','MaNGA_BCG','clusters']]
 axes[1].plot(np.arange(5),vals,'o-',label=label,color=['#89949b','#24627a'][i],lw=1.7,ms=5)
axes[1].set_xticks(np.arange(5),['Disks','Dwarfs','Early-type','Cluster\ncenters','Clusters'],rotation=25,ha='right');axes[1].set_ylabel('Whole-system-held-out RMS (dex)');axes[1].set_ylim(0,.47);axes[1].set_title('Where the remaining error lies',loc='left',fontsize=11,fontweight='bold');axes[1].legend(fontsize=7,frameon=False,loc='upper right')
for ax in axes:
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='x' if ax==axes[0] else 'y',alpha=.17);ax.set_axisbelow(True);ax.tick_params(labelsize=8)
fig.suptitle('One formula can improve a broad fit without winning every transfer test',fontsize=12,y=1.01)
fig.tight_layout();fig.savefig(R/'unified_regular_comparison.png',dpi=200,bbox_inches='tight');plt.close(fig)
print('Saved comparison figure')
