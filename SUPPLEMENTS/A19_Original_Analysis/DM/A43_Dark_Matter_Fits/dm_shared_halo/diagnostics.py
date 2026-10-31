"""Post-fit descriptive diagnostics; does not refit or select any model."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
from scipy.spatial import Delaunay

R = Path(__file__).resolve().parent
D = R / 'diagnostics'
D.mkdir(exist_ok=True)
res = json.loads((R/'results/results.json').read_text())
features = {v['system']: {k: float(x) for k,x in v.items() if k != 'system'}
            for v in csv.DictReader((R/'galaxy_features.csv').open())}
rows = list(csv.DictReader((R/'input.csv').open()))
names = sorted(features)

def write(name, records):
    with (D/name).open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)

summary=[]; galaxy_scores=[]; boot=[]; bounds=[]; states=[]; extrap=[]
for model, val in res['models'].items():
    tests={'all': val['all'], **{'cv_'+k:v for k,v in val['CV'].items()},
           **val['tail_holdout']}
    for test, data in tests.items():
        summary.append({'model':model, 'test':test, **data['metrics']})
        for name, score in data['by_galaxy'].items():
            galaxy_scores.append({'model':model,'test':test,'system':name,**score})
    fits=[('all',val['selected_all_fit'])]
    for seed,data in val['CV'].items():
        fits.extend((f'cv_{seed}_{j}',fid) for j,fid in enumerate(data['selected_fit_ids']))
    fits.extend((test,data['selected_fit_id']) for test,data in val['tail_holdout'].items())
    if model=='RAR_fixed': continue
    spec=res['spec'][model]; lo=np.array(spec['lo']); hi=np.array(spec['hi'])
    for test,fid in fits:
        fit=res['fit_log'][fid]; c=np.array(fit['coefficients'])
        near=np.minimum((c-lo)/(hi-lo),(hi-c)/(hi-lo))<=1e-4
        bounds.append({'model':model,'test':test,'selected_fit_id':fid,
            'success':fit['success'],'near_bound_indices':';'.join(map(str,np.flatnonzero(near))),
            'solver_active_indices':';'.join(map(str,np.flatnonzero(fit['active_bounds']))),
            'coefficients':json.dumps(c.tolist()),'cost':fit['cost']})
        if not model.startswith(('NFW','Burkert')): continue
        train={rows[i]['system'] for i in fit['train_rows']}
        targets=names if test=='all' else [g for g in names if g not in train]
        for name in targets:
            v=features[name]; m=np.log10(v['mass_proxy_msun']/1e10); h=np.log10(v['reff_kpc']/3)
            if model.endswith('mass'): la=c[0]+c[1]*m; lr=c[2]+c[3]*m
            else: la=c[0]+c[1]*m+c[2]*h; lr=c[3]+c[4]*m+c[5]*h
            states.append({'model':model,'test':test,'system':name,
                           'A_m_s2':float(10**la*1e-10),'rs_kpc':float(10**lr)})

for test,models in res['bootstrap'].items():
    for model,v in models.items():
        boot.append({'test':test,'comparison':model,'difference_dex':v['difference_dex'],
                     'lower_95_dex':v['percentile_95_dex'][0],'upper_95_dex':v['percentile_95_dex'][1]})

for test, tail in res['tails'].items():
    te=tail['galaxies']; tr=[g for g in names if g not in te]
    fields=['mass_proxy_msun','reff_kpc','sigma_proxy_msun_kpc2']
    entry={'test':test,'train_galaxies':len(tr),'test_galaxies':len(te),'test_rows':len(tail['rows'])}
    for field in fields:
        tv=np.array([features[g][field] for g in tr]); vv=np.array([features[g][field] for g in te])
        entry.update({field+'_train_min':float(tv.min()),field+'_train_max':float(tv.max()),
                      field+'_test_min':float(vv.min()),field+'_test_max':float(vv.max()),
                      field+'_test_outside_train_range':int(np.sum((vv<tv.min())|(vv>tv.max())))})
    def points(gg):return np.log10([[features[g]['mass_proxy_msun']/1e10,features[g]['reff_kpc']/3] for g in gg])
    entry['test_outside_train_log_mass_log_size_convex_hull']=int(np.sum(Delaunay(points(tr)).find_simplex(points(te))<0))
    extrap.append(entry)

write('score_summary.csv',summary);write('per_galaxy_scores.csv',galaxy_scores)
write('paired_bootstrap.csv',boot);write('selected_fit_bounds.csv',bounds)
write('predicted_halo_states.csv',states);write('tail_extrapolation.csv',extrap)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
models=['RAR_fixed','mixed_regular','NFW_mass','Burkert_mass','NFW_size','Burkert_size']
tests=['cv_20261027','cv_7331','cv_48105','mass_low','mass_high','surface_low','surface_high']
array=np.array([[res['bootstrap'][test][model+' minus RAR']['difference_dex'] for test in tests] for model in models])
fig,ax=plt.subplots(figsize=(10.3,4.6),layout='constrained')
image=ax.imshow(array,cmap='RdBu_r',vmin=-.055,vmax=.055,aspect='auto')
ax.set_yticks(range(len(models)),['Fixed RAR','A43 disk refit','NFW: mass','Burkert: mass','NFW: mass + size','Burkert: mass + size'])
ax.set_xticks(range(len(tests)),['CV seed 1','CV seed 2','CV seed 3','Low mass','High mass','Low surface\ndensity','High surface\ndensity'])
for i in range(len(models)):
    for j in range(len(tests)):
        ax.text(j,i,f'{array[i,j]:+.4f}',ha='center',va='center',fontsize=9,color='white' if abs(array[i,j])>.035 else 'black')
ax.set_title('Whole-galaxy prediction: difference from refitted RAR\nEqual-galaxy RMS log acceleration (dex); lower is better',fontsize=12)
fig.colorbar(image,ax=ax,label='RMS difference (dex)',shrink=.8)
fig.savefig(D/'comparison.png',dpi=190)
plt.close(fig)

boundcounts={model:sum(bool(v['near_bound_indices']) for v in bounds if v['model']==model) for model in res['spec']}
report={'method':'Post-fit descriptive diagnostics only; no model or fit changes.',
        'near_bound_rule':'within 1e-4 of coefficient range, in addition to recorded solver active flags',
        'selected_fits_with_near_bounds_of_20_per_model':boundcounts,
        'input_feature_ranges':{k:[min(v[k] for v in features.values()),max(v[k] for v in features.values())] for k in ['mass_proxy_msun','reff_kpc','sigma_proxy_msun_kpc2']},
        'tail_extrapolation':extrap}
(D/'diagnostics.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'bounds':boundcounts,'tail_outside_hulls':{r['test']:r['test_outside_train_log_mass_log_size_convex_hull'] for r in extrap}},indent=2))
