from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'backreaction_results.json').read_text())
rows=r['rows'];positive=[x for x in rows if x['feedback']]
ratios=[]
for case in r['references']:
 for s in [0,1]:
  a=next(x for x in positive if x['case']==case and x['lapse_power']==s and x['h']==.01)
  b=next(x for x in positive if x['case']==case and x['lapse_power']==s and x['h']==.005)
  for k in ['a_relative_sup_error','H_relative_sup_error','rho_relative_sup_error','pressure_relative_sup_error','max_constraint']:
   ratios.append({'case':case,'s':s,'quantity':k,'ratio':a[k]/b[k]})
summary=dict(positive_runs=len(positive),negative_runs=len(rows)-len(positive),
 max_purity_error=max(x['max_purity_error'] for x in rows),
 min_last_refinement_ratio=min(x['ratio'] for x in ratios),
 max_last_refinement_ratio=max(x['ratio'] for x in ratios),ratios=ratios,
 finest_max_clock_a=max(x['a_clock_relative_sup_difference'] for x in r['clock_rows'] if x['h']==.005),
 finest_max_constraint=max(x['max_constraint'] for x in positive if x['h']==.005))
summary['all_checks_passed']=(summary['positive_runs']==30 and summary['negative_runs']==10 and summary['max_purity_error']<1e-12 and 3.8<summary['min_last_refinement_ratio']<=summary['max_last_refinement_ratio']<4.2 and summary['finest_max_constraint']<3e-5)
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(9,3.1),layout='constrained')
for case,label,color in [('vacuum','Vacuum','#17658c'),('squeezed','Squeezed','#b26422'),('vacuum_Lambda','Vacuum + Lambda','#528559')]:
 d=np.load(ROOT/f'reference_{case}.npz');axs[0].plot(d['time'],d['a'],label=label,color=color,lw=2)
axs[0].set(xlabel='Proper time',ylabel='Scale factor a',title='Geometry responds to quantum matter')
axs[0].legend(frameon=False,fontsize=8)
for feedback,s,label,style in [(True,0,'Full source, proper','-o'),(True,1,'Full source, conformal','-s'),(False,0,'Feedback omitted, proper','--^')]:
 rr=[x for x in rows if x['case']=='vacuum' and x['feedback']==feedback and x['lapse_power']==s]
 axs[1].loglog([x['h'] for x in rr],[x['max_constraint'] for x in rr],style,label=label,markersize=4)
axs[1].set(xlabel='Coordinate step h',ylabel='Max. relative Friedmann residual',title='Purity alone misses gravity failure')
axs[1].legend(frameon=False,fontsize=8)
fig.savefig(ROOT/'backreaction.png',dpi=190)
print(json.dumps({k:v for k,v in summary.items() if k!='ratios'},indent=2))
