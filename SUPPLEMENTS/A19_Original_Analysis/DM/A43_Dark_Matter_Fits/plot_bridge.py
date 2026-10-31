"""Reproduce the bridge figure from frozen numeric outputs; no fitting.

Run from any directory: python /path/to/tmp/dm_bridge_integration/plot_bridge.py
Requires NumPy and Matplotlib, plus the four sibling result JSON files.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

ROOT = Path(__file__).resolve().parent
FILES = {
    'shared_halo': ROOT/'dm_shared_halo/results/results.json',
    'cosmological_halo': ROOT/'dm_cosmological_halo/results/results.json',
    'gas_bridge': ROOT/'dm_gas_bridge/results.json',
    'field_bridge': ROOT/'dm_field_bridge/results.json',
}
data = {name: json.loads(path.read_text()) for name,path in FILES.items()}
seed='20261027'
shared=data['shared_halo']; cosmo=data['cosmological_halo']
gas=data['gas_bridge']; field=data['field_bridge']
reference=shared['models']['RAR']['CV'][seed]['metrics']['galaxy_RMS_dex']
panel_a=[]
for label,key,source in [
    ('A43 disk\n7 coefficients','mixed_regular',shared),
    ('NFW: mass + size\n6 coefficients','NFW_size',shared),
    ('NFW + c(M)\n3 coefficients','NFW_Planck_size',cosmo),
]:
    score=source['models'][key]['CV'][seed]['metrics']['galaxy_RMS_dex']
    uncertainty=source['bootstrap']['cv_'+seed][key+' minus RAR']
    assert abs(score-reference-uncertainty['difference_dex'])<1e-12
    panel_a.append({'label':label,'model':key,'RMS_dex':score,**uncertainty})

panel_b=[]
for label,fixed,fitted in [('RAR','RAR_fixed','RAR_fit'),('NLG','NLG_frozen','NLG_fit'),('MOG','MOG_frozen','MOG_fit')]:
    panel_b.append({'label':label,
        'fixed_RMS_dex':field['models_by_ML']['2'][fixed]['CV'][seed]['metrics']['RMS_dex'],
        'refit_RMS_dex':field['models_by_ML']['2'][fitted]['CV'][seed]['metrics']['RMS_dex']})

gas_panels={}
for domain in ['disk','cluster']:
    separate=gas['scopes'][domain]['simple']['fitted']['CV'][seed]['metrics']
    gas_panels[domain]={
        'rows':separate['rows'],'systems':separate['systems'],
        'values':[
            {'label':'No extra gas','B_dex':gas['scopes'][domain]['simple']['zero']['CV'][seed]['metrics']['B_dex']},
            {'label':'Shared\ncoupling','B_dex':gas['scopes']['joint']['simple']['fitted']['CV'][seed]['metrics']['domains'][domain]['B_dex']},
            {'label':'Separate\ncalibrations','B_dex':separate['B_dex']},
        ]}

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':12.5,'axes.titlesize':14,
    'axes.labelsize':12.5,'xtick.labelsize':11.5,'ytick.labelsize':12,
    'axes.edgecolor':'#a6afb7','axes.linewidth':.8,'savefig.facecolor':'white',
})
fig,axes=plt.subplots(2,2,figsize=(8.8,7.6))
fig.subplots_adjust(left=.215,right=.98,top=.885,bottom=.115,wspace=.74,hspace=.91)
blue='#276c9e';teal='#167d82';orange='#cb7835';gray='#727b84';purple='#85629a'

def base(ax,title,subtitle):
    ax.set_title(title,loc='left',fontweight='bold',pad=42)
    ax.text(0,1.075,subtitle,transform=ax.transAxes,ha='left',va='bottom',fontsize=11.5)
    ax.set_axisbelow(True);ax.grid(axis='x',color='#dbe1e6',linewidth=.7)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=8)
    ax.tick_params(axis='x',length=3,color='#a6afb7')

ax=axes[0,0]
base(ax,'A  Unseen galaxies',f'149 galaxies · RAR {reference:.4f} dex')
ax.axvline(0,color='#7d858d',linewidth=1,linestyle='--',zorder=1)
for y,(entry,color) in enumerate(zip(panel_a,[purple,blue,teal])):
    mid=entry['difference_dex'];lo,hi=entry['percentile_95_dex']
    ax.errorbar(mid,y,xerr=np.array([[mid-lo],[hi-mid]]),fmt='o',markersize=7,
        color=color,elinewidth=2,capsize=4,zorder=3)
ax.set_yticks(range(len(panel_a)),[v['label'] for v in panel_a]);ax.set_ylim(2.65,-.65)
ax.set_xlim(-.015,.018);ax.xaxis.set_major_locator(MultipleLocator(.01))
ax.xaxis.set_major_formatter(FormatStrFormatter('%+.02f'))
ax.set_xlabel('RMS − RAR RMS (dex)\nLower is better',labelpad=9)

ax=axes[0,1]
base(ax,'B  Dwarf kernels','44 dwarfs · stellar M/L = 2')
yy=np.arange(len(panel_b));height=.28
for offset,column,label,color in [(-.18,'fixed_RMS_dex','Fixed',gray),(.18,'refit_RMS_dex','Dwarf refit',teal)]:
    vals=[v[column] for v in panel_b]
    ax.barh(yy+offset,vals,height=height,label=label,color=color)
    for y,value in zip(yy+offset,vals):ax.text(value+.045,y,f'{value:.3f}',va='center',fontsize=11)
ax.set_yticks(yy,[v['label'] for v in panel_b]);ax.set_ylim(2.55,-.60)
ax.set_xlim(0,2.42);ax.xaxis.set_major_locator(MultipleLocator(.5))
ax.set_xlabel('Acceleration RMS (dex)\nLower is better',labelpad=9)
ax.legend(loc='upper right',bbox_to_anchor=(1.04,1.005),frameon=False,fontsize=10.8,
          borderpad=0,labelspacing=.35,handlelength=1.1,handletextpad=.45)

for ax,domain,title,subtitle in [
    (axes[1,0],'disk','C  Cold-gas disks','61 galaxies · 402 selected radii'),
    (axes[1,1],'cluster','D  Hot-gas clusters','23 clusters · 86 estimates'),
]:
    base(ax,title,subtitle)
    entries=gas_panels[domain]['values'];vals=[v['B_dex'] for v in entries]
    ax.barh(np.arange(3),vals,height=.48,color=[gray,orange,teal])
    for y,value in enumerate(vals):ax.text(value+.012,y,f'{value:.3f}',va='center',fontsize=11.5)
    ax.set_yticks(range(3),[v['label'] for v in entries]);ax.set_ylim(2.65,-.65)
    ax.set_xlim(0,.56);ax.xaxis.set_major_locator(MultipleLocator(.1))
    ax.set_xlabel('Balanced RMS (dex)\nLower is better',labelpad=9)

fig.savefig(ROOT/'bridge_comparison.png',dpi=300)
fig.savefig(ROOT/'bridge_comparison.pdf')
plt.close(fig)

provenance={
    'seed':int(seed),'purpose':'Display-only comparison; no fits or model selection performed.',
    'source_sha256':{str(path.relative_to(ROOT.parent)):hashlib.sha256(path.read_bytes()).hexdigest() for path in FILES.values()},
    'A_reference_RAR_RMS_dex':reference,'A':panel_a,'B':panel_b,'C':gas_panels['disk'],'D':gas_panels['cluster'],
    'intervals':'A: paired galaxy-Poisson bootstrap 95% intervals of frozen primary CV residuals; no refitting or selection/systematics uncertainty.',
}
(ROOT/'bridge_figure_values.json').write_text(json.dumps(provenance,indent=2))
print('Wrote bridge_comparison.png, bridge_comparison.pdf and bridge_figure_values.json')
