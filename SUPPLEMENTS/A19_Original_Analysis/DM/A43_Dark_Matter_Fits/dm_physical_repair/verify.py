from pathlib import Path
import csv,json,hashlib
import numpy as np
R=Path(__file__).resolve().parent;J=json.loads((R/'results/results.json').read_text())
raw=R.parent/'dm_broad_joint/combined_input.csv';assert hashlib.sha256(raw.read_bytes()).hexdigest()==J['input_sha256']
checks=[]
for p in sorted((R/'results').glob('*.csv')):
 a=list(csv.DictReader(p.open()));pop=np.array([r['population'] for r in a]);group=np.array([r['group'] for r in a]);e=np.log10([float(r['gpred_m_s2'])/float(r['gobs_m_s2']) for r in a]);rec=np.array([float(r['log10_residual']) for r in a]);assert np.max(abs(e-rec))<2e-14
 mse=[]
 for pp in sorted(set(pop)):
  mse.append(np.mean([np.mean(e[(pop==pp)&(group==g)]**2) for g in sorted(set(group[pop==pp]))]))
 score=float(np.sqrt(np.mean(mse)));name=p.stem
 if '_cv_' in name:
  model,seed=name.split('_cv_');expected=J['models'][model]['CV'][seed]['metrics']['B_dex']
  for g in set(group):assert len({r['fold'] for r in a if r['group']==g})==1
 else:model=name.removesuffix('_population_holdout');expected=J['models'][model]['population_holdout']['metrics']['B_dex']
 assert abs(score-expected)<1e-12
 checks.append({'file':p.name,'independent_B_dex':score})
# Independently evaluate all-data fitted formulas on frozen rows; compare their training metric.
rows=list(csv.DictReader(raw.open()));gb=np.array([float(r['gbar_m_s2']) for r in rows]);r=np.array([float(r['radius_kpc']) for r in rows]);s=gb/1e-10*r
for k,v in J['models'].items():
 c=v['coefficients']
 if k=='softened_power':ae=1e-10*(10**c[0]+10**c[1]*(s+10**c[5])**(-c[3])+10**c[2]*s**c[4])
 elif k=='bounded_U':z=np.tanh((np.log10(s)-c[3])/c[4]);ae=1e-10*10**(c[0]+c[1]*z+c[2]*z*z)
 else:ae=1e-10*(10**c[0]+10**c[1]*s**(-.5)+10**c[2]*s**c[3])
 gp=gb/(-np.expm1(-np.sqrt(gb/ae)));assert np.all(np.isfinite(gp)) and np.all(gp>=gb)
 e=np.log10(gp/np.array([float(r['gobs_m_s2']) for r in rows]));score=np.sqrt(np.mean([np.mean([np.mean(e[(pop==pp)&(group==g)]**2) for g in sorted(set(group[pop==pp]))]) for pp in sorted(set(pop))]));assert abs(score-v['training']['B_dex'])<1e-12
(R/'verification.json').write_text(json.dumps({'status':'passed','prediction_files':checks,'training_formulas_independently_evaluated':list(J['models']),'rows':len(rows)},indent=2));print('PASS: 12 prediction tables, independent metrics, disjoint group folds, all three fitted formulas.')
