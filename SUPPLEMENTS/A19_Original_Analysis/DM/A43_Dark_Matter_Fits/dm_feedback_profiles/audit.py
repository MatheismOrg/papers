"""Independent checks of saved predictions, exclusions and fitted mass profiles.

Run after run.py. Metric arithmetic and DC14 mass quadrature are independent of
the production implementation. The imported module supplies parameters/inputs
and predictions to compare, not the numerical quadrature or scoring formula.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.integrate import quad
import run as model

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
s=json.loads((OUT/'summary.json').read_text())
records=[json.loads(p.read_text()) for p in sorted(OUT.glob('*_record.json'))]
tab=pd.read_csv(ROOT/'input.csv')
assert len(tab)==2700 and tab.system.nunique()==149
checks={}
for name in ['input.csv','galaxy_features.csv']:
 a=(ROOT/name).read_bytes();b=(ROOT.parent/'dm_cosmological_halo'/name).read_bytes()
 assert a==b
 checks[name+'_sha256']=hashlib.sha256(a).hexdigest()

largest_metric_difference=0.
for seed in model.SEEDS:
 old=pd.read_csv(ROOT.parent/f'dm_cosmological_halo/results/NFW_Planck_size_cv_{seed}.csv')
 rr=[r for r in records if r['name'].startswith(f'cv_{seed}_')]
 assert len(rr)==5
 visits=np.zeros(len(tab),dtype=int)
 for r in rr:
  te=np.array(r['test']);tr=np.array(r['train']);fold=int(r['name'].split('_')[-1])
  assert not set(tab.system.iloc[te])&set(tab.system.iloc[tr])
  assert sorted(np.r_[te,tr])==list(range(len(tab)))
  assert np.array_equal(te,np.where(old.fold.to_numpy()==fold)[0])
  visits[te]+=1
 assert np.all(visits==1)
 for m in s['means']:
  t=pd.read_csv(OUT/f'{m}_cv_{seed}.csv')
  assert np.array_equal(t.system,tab.system)
  e=np.log10(t.gpred_m_s2/t.gobs_m_s2)
  mse=pd.DataFrame({'system':t.system,'e2':e**2}).groupby('system').e2.mean()
  rms=float(np.sqrt(mse.mean()));speed=float(100*np.median(np.abs(np.sqrt(t.gpred_m_s2/t.gobs_m_s2)-1)))
  largest_metric_difference=max(largest_metric_difference,abs(rms-s['repeats'][str(seed)][m]['rms']),abs(speed-s['repeats'][str(seed)][m]['median_speed']))
assert largest_metric_difference<1e-11
checks['maximum_recomputed_metric_difference']=largest_metric_difference

selected_failures=[];failures=0;attempts=0
for r in records:
 for m,x in r['models'].items():
  a=x['fit'].get('attempts',[]);attempts+=len(a)
  failures+=sum(not z.get('success',True) for z in a)
  if a and 'loss' in a[0] and not min(a,key=lambda z:z['loss'])['success']:
   selected_failures.append([r['name'],m])
 assert r['models']['DC14_response']['fit']['loss']<=min(r['models'][m]['fit']['loss'] for m in ['NFW','DC14'])+1e-13
assert not selected_failures
checks.update(stored_optimizer_attempts=attempts,failed_attempts=failures,selected_failures=selected_failures)

# Integrate actual all-data fitted DC14 densities at observed radii and compare
# mass normalization to direct quadrature. No hypergeometric functions here.
relative_errors=[];min_increment=1.;finite_errors=[]
for m in ['DC14','DC14_response']:
 p=s['all_coefficients'][m]['p'];mn,md,h=model.components(p)
 for j,name in enumerate(model.NAMES):
  a,b,g=h['a'][j],h['b'][j],h['g'][j];rs=h['dc_rs'][j]
  radii=model.DATA['r'][model.DATA['gi']==j]
  xmax=h['Rv'][j]/rs
  def integ(x):
   return quad(lambda t:np.exp((3-g)*t-(b-g)/a*np.logaddexp(0,a*t)),np.log(x)-50,np.log(x),epsabs=1e-12,epsrel=2e-11)[0]
  norm=integ(xmax)
  for row in np.flatnonzero(model.DATA['gi']==j)[[0,-1]]:
   val=h['mh'][j]*integ(min(model.DATA['r'][row],h['Rv'][j])/rs)/norm
   relative_errors.append(abs(val/md[row]-1))
  x=np.geomspace(1e-7,xmax,150)
  mass=model.mass_integral(x,a,b,g)/model.mass_integral(xmax,a,b,g)
  assert np.all(np.diff(mass)>0) and np.all(mass>0)
  min_increment=min(min_increment,float(np.min(np.diff(mass))))
  finite_errors.append(abs(float(mass[-1])-1))
checks['actual_profile_quadrature_cases']=len(relative_errors)
checks['actual_profile_max_mass_relative_error']=max(relative_errors)
checks['actual_profile_min_mass_increment']=min_increment
checks['actual_profile_boundary_mass_error']=max(finite_errors)
assert max(relative_errors)<2e-7

# Independent speed-to-acceleration identity and cross-population scores.
d=pd.read_csv(OUT/'frozen_dwarf_transfer.csv')
assert len(d)==44 and d.system.nunique()==44
for m in s['means']:
 v=float(np.sqrt(np.mean(np.log10(d[m]/d.gobs_m_s2)**2)))
 assert abs(v-s['dwarfs'][m]['rms'])<1e-12
checks['dwarf_rows']=len(d)
checks['baseline_max_log_prediction_difference']=s['baseline_reproduction']['maximum_absolute_log_prediction_difference']
assert checks['baseline_max_log_prediction_difference']<1e-6
checks['passed']=True
(OUT/'audit.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
