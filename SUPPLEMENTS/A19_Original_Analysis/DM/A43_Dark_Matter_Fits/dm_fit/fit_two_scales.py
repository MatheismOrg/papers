"""Exploratory common polynomial for galaxies and seven cluster R500 estimates."""
from pathlib import Path
import json,csv
import numpy as np
ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'results/held_out_predictions.csv').open()))
c=json.loads((ROOT/'results/cluster_test.json').read_text())['clusters']
groups=np.array([r['galaxy'] for r in rows]+[r['cluster'] for r in c])
x=np.log10(np.array([float(r['gbar_m_s2']) for r in rows]+[r['gbar'] for r in c])/1e-10)
y=np.log10(np.array([float(r['gobs_m_s2']) for r in rows]+[r['gobs'] for r in c])/1e-10)
z=np.log10(np.array([float(r['radius_kpc']) for r in rows]+[r['r500_kpc'] for r in c]))
fold=np.array([int(r['fold']) for r in rows]+[i%5 for i in np.random.default_rng(20261027).permutation(len(c))])
cluster=np.arange(len(x))>=len(rows)
w=np.array([1/np.sum(groups==g) for g in groups])
res={}
for kind in ['quadratic_acceleration','quadratic_acceleration_radius']:
    X=np.column_stack([np.ones(len(x)),x,x*x] if kind=='quadratic_acceleration' else [np.ones(len(x)),x,z,x*x,x*z,z*z])
    pred=np.zeros(len(y))
    for k in range(5):
        train=fold!=k;test=~train
        assert not set(groups[train])&set(groups[test])
        coef=np.linalg.lstsq(X[train]*np.sqrt(w[train,None]),y[train]*np.sqrt(w[train]),rcond=None)[0]
        pred[test]=X[test]@coef
    coef=np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0]
    e=pred-y
    res[kind]={'coefficients':coef.tolist(),
        'galaxy_point_cv_rmse_dex':float(np.sqrt(np.mean(e[~cluster]**2))),
        'galaxy_median_abs_velocity_error':float(np.median(abs(10**(e[~cluster]/2)-1))),
        'cluster_cv_rmse_dex':float(np.sqrt(np.mean(e[cluster]**2))),
        'cluster_median_abs_acceleration_error':float(np.median(abs(10**e[cluster]-1))),
        'cluster_cv_observed_over_predicted':(10**(-e[cluster])).tolist()}
    with (ROOT/'results'/f'{kind}_combined_predictions.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['system','cluster','fold','x','z','y','held_out_prediction_y'])
        writer.writerows(zip(groups,cluster,fold,x,z,y,pred))
result={'status':'Exploratory, selected after examining galaxy-to-cluster transfer failure; requires new independent data.',
        'fit_weights':'Every galaxy and cluster has total weight 1. Cluster subset is seven R500 estimates, one per cluster.',
        'validation':'Five whole-system folds; some clusters are now in training. This is not the galaxy-only transfer test.',
        'basis':'x=log10(gbar/1e-10 m/s^2); z=log10(radius/kpc); 2D basis [1,x,z,x^2,x*z,z^2]',
        'models':res}
(ROOT/'results/two_scale_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
