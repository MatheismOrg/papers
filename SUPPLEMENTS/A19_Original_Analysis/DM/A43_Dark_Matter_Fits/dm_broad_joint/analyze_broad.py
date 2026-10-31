from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
j=json.loads((OUT/'broad_results.json').read_text())
rows=list(csv.DictReader((ROOT/'combined_input.csv').open()))
x=np.log10([float(r['gbar_m_s2']) for r in rows]);y=np.log10([float(r['gobs_m_s2']) for r in rows]);r=np.array([float(a['radius_kpc']) for a in rows]);p=np.array([a['population'] for a in rows]);sp=np.array([a['subpopulation'] for a in rows]);group=np.array([a['group'] for a in rows])
order=['SPARC','ETG','dSph','MaNGA_BCG','CLASH']
labels={'SPARC':'Rotating galaxies (SPARC)','ETG':'Early-type galaxies','dSph':'Dwarf spheroidals','MaNGA_BCG':'Central cluster galaxies','CLASH':'Cluster lensing (CLASH)'}
colors={'SPARC':'#4c78a8','ETG':'#e69f00','dSph':'#b04c9a','MaNGA_BCG':'#14856f','CLASH':'#b53535'}
def rar(g):return g/(-np.expm1(-np.sqrt(g/1.1603274962269242e-10)))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
fig,axs=plt.subplots(1,2,figsize=(12,5.2),layout='constrained')
for pop in order:
 m=p==pop;col=colors[pop];sz=9 if pop=='SPARC' else 26
 axs[0].scatter(x[m],y[m],s=sz,c=col,alpha=.3 if pop=='SPARC' else .78,label=labels[pop],edgecolors='none')
 axs[1].scatter(np.log10(r[m]),y[m]-np.log10(rar(10**x[m])),s=sz,c=col,alpha=.3 if pop=='SPARC' else .78,edgecolors='none')
xx=np.linspace(x.min()-.12,x.max()+.1,300)
axs[0].plot(xx,np.log10(rar(10**xx)),color='black',lw=1.3,label='Frozen galaxy RAR')
axs[0].plot(xx,xx,color='grey',ls=':',lw=1,label='Newtonian baryons')
axs[0].set_xlabel('log10 baryonic acceleration [m/s²]');axs[0].set_ylabel('log10 inferred total acceleration [m/s²]')
axs[0].set_title('Shared regularity and population differences')
axs[0].legend(loc='upper left',fontsize=8,frameon=False)
axs[1].axhline(0,color='black',lw=1);axs[1].set_xlabel('log10 physical radius [kpc]');axs[1].set_ylabel('log10(observed / galaxy-RAR prediction)')
axs[1].set_title('Radius alone does not organize all residuals')
for ax in axs:ax.grid(alpha=.13)
fig.suptitle('2,910 measurements or inferred estimates in 279 objects; five data populations',fontsize=13)
fig.savefig(OUT/'broad_acceleration_planes.png',dpi=180,bbox_inches='tight');plt.close(fig)

kinds=['frozen','RAR','RAR_radius_quad','RAR_potential_quad','RAR_density_quad','quadratic_2D','nested','RAR_piecewise']
names=['Frozen SPARC RAR*','Common fitted RAR','RAR + radius curvature','RAR + potential-proxy curvature','RAR + density-proxy curvature','Common 2D quadratic','Nested shared-model selection','Piecewise RAR by population']
values=[]
for k in kinds:
 z=j['frozen_SPARC_RAR'] if k=='frozen' else j['nested_shared_selection']['metrics'] if k=='nested' else j['models'][k]['cross_validation']
 values.append([100*z['populations'][q]['median_absolute_acceleration_fraction'] for q in order])
values=np.array(values)
fig,ax=plt.subplots(figsize=(10.7,5.6),layout='constrained')
im=ax.imshow(values,cmap='YlOrRd',vmin=0,vmax=100,aspect='auto')
ax.set_xticks(range(len(order)),['SPARC','Early type','Dwarf sph.','MaNGA BCG','CLASH'])
ax.set_yticks(range(len(names)),names)
for (i,k),v in np.ndenumerate(values):ax.text(k,i,f'{v:.1f}%',ha='center',va='center',color='white' if v>66 else 'black',fontsize=11)
ax.set_title('Median absolute acceleration error by population\nAll percentages here use acceleration, not rotation speed',fontsize=13,pad=13)
fig.colorbar(im,ax=ax,label='Median absolute acceleration error (%)',shrink=.85)
fig.text(.5,-.025,'*Frozen calibration: SPARC in-sample; other populations external. Other rows: grouped internal validation.',ha='center',fontsize=9)
fig.savefig(OUT/'broad_model_errors.png',dpi=180,bbox_inches='tight');plt.close(fig)

# Common coefficient models are compared across matched internal and harder population holdouts.
models=['RAR','RAR_radius_quad','RAR_potential_quad','RAR_density_quad','quadratic_2D']
names2=['Common RAR','Radius curvature','Potential-proxy curvature','Density-proxy curvature','2D quadratic']
a=[j['models'][k]['cross_validation']['balanced_population_system_rms_dex'] for k in models]
b=[float(np.sqrt(np.mean([v['models'][k]['metrics']['equal_system_rms_dex']**2 for v in j['leave_population_out'].values()]))) for k in models]
fig,ax=plt.subplots(figsize=(10,4.7),layout='constrained');yy=np.arange(len(models))
ax.barh(yy-.18,a,height=.32,color='#235d75',label='Whole-system folds within populations')
ax.barh(yy+.18,b,height=.32,color='#b87943',label='Entire population withheld')
ax.set_yticks(yy,names2);ax.invert_yaxis();ax.set_xlabel('Equal-population, equal-group RMS error [dex]')
ax.set_title('A good shared fit can still fail on an unseen population',fontsize=13)
ax.legend(frameon=False,fontsize=9);ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
fig.savefig(OUT/'population_holdout_comparison.png',dpi=180,bbox_inches='tight');plt.close(fig)

# Verify every main prediction CSV against the declared primary metrics independently.
checked=[]
for kind in list(j['models'])+['nested_shared']:
 rr=list(csv.DictReader((OUT/(kind+'_predictions.csv')).open()))
 expected=j['nested_shared_selection']['metrics'] if kind=='nested_shared' else j['models'][kind]['cross_validation']
 assert len(rr)==len(rows)
 bygroup={}
 for z in rr:bygroup.setdefault(z['group'],set()).add(z['fold'])
 assert all(len(v)==1 for v in bygroup.values())
 for pp in order:
  ss=[z for z in rr if z['population']==pp];rat=np.array([float(z['pred_gobs_m_s2'])/float(z['gobs_m_s2']) for z in ss])
  med=float(np.median(abs(rat-1)));rms=float(np.sqrt(np.mean(np.log10(rat)**2)))
  assert abs(med-expected['populations'][pp]['median_absolute_acceleration_fraction'])<1e-11
  assert abs(rms-expected['populations'][pp]['rms_log10_acceleration_dex'])<1e-11
 checked.append(kind)
out={'verified_models':checked,'rows':len(rows),'objects':len(set(z['system'] for z in rows)),'validation_groups':len(set(group)),
 'population_holdout_equal_weight_RMS_dex':dict(zip(models,b)),
 'caution':'Frozen SPARC RAR is calibrated on SPARC, so its SPARC row is descriptive training-sample performance. Its other population rows are frozen transfer. New model rows use matched whole-system internal validation. Error bars/covariances are not used in fitting.'}
(OUT/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
