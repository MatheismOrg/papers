"""Conditional paired-system bootstrap and figures for search_formulas.py outputs."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
s=json.loads((OUT/'search_results.json').read_text())
rows=list(csv.DictReader((OUT/'galaxy_only_equal_point_RAR_predictions.csv').open()))
g=np.array([r['system'] for r in rows]);ug=np.unique(g)
base=np.array([float(r['log10_residual']) for r in rows])
count=np.array([np.sum(g==a) for a in ug]);ssb=np.array([np.sum(base[g==a]**2) for a in ug])
rng=np.random.default_rng(123456);draw=rng.integers(0,len(ug),(10000,len(ug)))
results={}
for kind in ['quadratic','quadratic_radius','RAR_radius','RAR_exponent','RAR_radius_exponent']:
    rr=list(csv.DictReader((OUT/f'galaxy_only_equal_point_{kind}_predictions.csv').open()))
    assert [r['system'] for r in rr]==list(g)
    e=np.array([float(r['log10_residual']) for r in rr]);ssa=np.array([np.sum(e[g==a]**2) for a in ug])
    delta=np.sqrt(np.sum(ssa[draw],axis=1)/np.sum(count[draw],axis=1))-np.sqrt(np.sum(ssb[draw],axis=1)/np.sum(count[draw],axis=1))
    results[kind]={'point_delta_rms_dex':float(np.sqrt(ssa.sum()/count.sum())-np.sqrt(ssb.sum()/count.sum())),
     'percentile_95_interval_dex':np.quantile(delta,[.025,.975]).tolist()}
result={'status':'Descriptive paired whole-galaxy bootstrap of existing five-fold predictions,10000 draws,seed123456. Conditional on those trained models; omits retraining uncertainty, nuisance-parameter uncertainty, and candidate-selection uncertainty. Not a confirmatory p-value.',
 'quantity':'Point-weighted RMS(log10 acceleration) of candidate minus standard RAR; negative favors candidate.', 'comparisons':results}
(OUT/'paired_galaxy_bootstrap.json').write_text(json.dumps(result,indent=2))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
kinds=['RAR','quadratic','quadratic_radius','RAR_radius','RAR_exponent','RAR_radius_exponent']
labels=['Standard RAR','Quadratic','Quadratic + radius','RAR + radius','RAR + exponent','RAR + radius + exponent']
fig,ax=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
for j,(key,factor,title) in enumerate([('galaxy_point_rmse_dex',1,'Acceleration RMS error (dex)'),('galaxy_median_abs_speed_error',100,'Median absolute speed error (%)')]):
    vals=[s['galaxy_only_equal_point_training'][q]['galaxy_cv'][key]*factor for q in kinds]
    ax[j].barh(labels,vals,color=['#174a67']+['#6c9fae']*5)
    ax[j].invert_yaxis();ax[j].set_xlabel(title)
    ax[j].set_xlim(0,max(vals)*1.18)
    for i,v in enumerate(vals):ax[j].text(v+max(vals)*.018,i,f'{v:.4f}' if j==0 else f'{v:.2f}',va='center',fontsize=10)
    ax[j].grid(axis='x',alpha=.18);ax[j].set_axisbelow(True)
ax[1].set_yticklabels([])
fig.suptitle('Galaxy-only validation: extra flexibility gives no consistent improvement',fontsize=13)
fig.savefig(OUT/'galaxy_candidate_errors.png',bbox_inches='tight');plt.close(fig)

models=['constant_acceleration','radius_definition_g_equals_C_r']
fig,ax=plt.subplots(figsize=(8.5,4.2),layout='constrained')
for i,(name,label,col) in enumerate(zip(models,['Constant acceleration null','Radius-definition null: g = C r'],['#b87943','#174a67'])):
    rr=list(csv.DictReader((OUT/f'cluster_null_{name}_predictions.csv').open()))
    ratio=np.array([10**float(r['log10_residual']) for r in rr]);ids=[r['system'] for r in rr]
    ax.plot(np.arange(len(ids)),ratio,'o-',label=label,color=col)
rr=list(csv.DictReader((OUT/'equal_system_quadratic_radius_predictions.csv').open()))
rr=[r for r in rr if r['cluster']=='True'];ratio=np.array([10**float(r['log10_residual']) for r in rr])
ax.plot(np.arange(len(ids)),ratio,'s--',label='Radius-augmented quadratic',color='#80924a')
ax.axhline(1,color='black',lw=.8);ax.set_xticks(np.arange(len(ids)),ids);ax.set_ylabel('Held-out predicted / catalog acceleration')
ax.set_title('R500 carries mass information: a definitional null beats the fitted relation',fontsize=12)
ax.legend(loc='best',fontsize=9);ax.grid(axis='y',alpha=.18)
fig.savefig(OUT/'cluster_definition_null.png',bbox_inches='tight');plt.close(fig)
print(json.dumps(result,indent=2))
