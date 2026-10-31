from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
B=R
m=json.loads((B/'dm_mond_repairs/results/results.json').read_text())
h=json.loads((B/'dm_halo_benchmark/results/results.json').read_text())
old=json.loads((B/'dm_unified_regular/results/results.json').read_text())['models_by_dataset']['augmented']['models']['mixed_regular']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(11.2,5.1))
names=['Common RAR','Earlier A43','Transient repair','Potential taper','Two windows']
data=[m['models']['RAR'],old,m['models']['mixed_transient'],m['models']['mixed_potential_taper'],m['models']['two_log_bumps']]
y=np.arange(len(data));dh=.33
for d,lab,col in [(-dh/2,'Excluded systems','#23669b'),(dh/2,'Excluded population','#cf7c35')]:
 vals=[v['CV']['20261027']['metrics']['B_dex'] if d<0 else v['population_holdout']['metrics']['B_dex'] for v in data]
 ax[0].barh(y+d,vals,height=dh,color=col,label=lab)
ax[0].set_yticks(y,names);ax[0].invert_yaxis();ax[0].set_xlim(0,.43)
ax[0].set_title('Shared formula: 286 systems',weight='bold',pad=39)
ax[0].set_xlabel('Population-balanced RMS (dex)')
ax[0].legend(loc='lower center',bbox_to_anchor=(.5,1.005),ncol=2,fontsize=9,frameon=False,handlelength=1.2,columnspacing=1.0)
names=['Shared RAR','A43, disk trained','NFW halo','Burkert halo','RAR + stellar scale']
keys=['RAR_shared','mixed_regular_disk','NFW_fixed','Burkert_fixed','RAR_stellar']
for d,lab,sk,col in [(-dh/2,'Outer third','outer_third','#23669b'),(dh/2,'Middle third','middle_third','#cf7c35')]:
 vals=[h['splits'][sk]['models'][k]['test']['galaxy_RMS_dex'] for k in keys]
 ax[1].barh(y+d,vals,height=dh,color=col,label=lab)
ax[1].set_yticks(y,names);ax[1].invert_yaxis();ax[1].set_xlim(0,.155)
ax[1].set_title('Withheld radii: 104 disk galaxies',weight='bold',pad=39)
ax[1].set_xlabel('Galaxy-balanced RMS (dex)')
ax[1].legend(loc='lower center',bbox_to_anchor=(.5,1.005),ncol=2,fontsize=9,frameon=False,handlelength=1.2,columnspacing=1.0)
for a in ax:a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
fig.subplots_adjust(left=.16,right=.985,bottom=.20,top=.78,wspace=.72)
fig.text(.5,.035,'Lower is better. The panels test different tasks; halos use two fitted coefficients per galaxy.\nThese conditional comparisons do not establish superiority over full published theories.',ha='center',fontsize=9.4)
fig.savefig(R/'best_models_comparison.png',dpi=220,facecolor='white')
