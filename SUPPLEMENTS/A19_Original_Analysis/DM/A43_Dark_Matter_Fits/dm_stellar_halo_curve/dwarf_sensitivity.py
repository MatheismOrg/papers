"""Reuse the earlier high-precision dwarf cut; no tuning of its threshold."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent;d=pd.read_csv(R.parent/'dm_feedback_profiles/dwarfs_benchmark.csv');p=pd.read_csv(R/'results/dwarf_comparison.csv');assert np.array_equal(d.system,p.system)
mask=d[['sigma_err_plus_km_s','sigma_err_minus_km_s']].max(axis=1)/d.sigma_los_km_s<=.25
mask&=~d.possible_tidal_or_classification_caution
assert int(mask.sum())==25
models=['old_RAR','old_NFW','old_RAR_NFW_blend','stellar_power','stellar_curved','stellar_blend','Moster_fixed']
out={}
for label,selection in [('all44',np.ones(len(d),bool)),('archived_high_precision25',mask.to_numpy())]:
 z={}
 for model in models:
  residual=np.log10(p[model]/p.gobs_m_s2).to_numpy()[selection];rel=(p[model]/p.gobs_m_s2-1).to_numpy()[selection]
  z[model]=dict(rms=float(np.sqrt(np.mean(residual*residual))),median_absolute_log_error=float(np.median(abs(residual))),median_acceleration_percent=float(100*np.median(abs(rel))))
 out[label]=z
base=np.log10(p.old_NFW/p.gobs_m_s2).to_numpy()
out['individual_dwarfs_improved_log_error_vs_old_NFW']={m:int(np.sum(abs(np.log10(p[m]/p.gobs_m_s2))<abs(base))) for m in ['stellar_power','stellar_blend']}
out['subset_names']=d.system[mask].tolist();out['interpretation']='This reuses a previously declared source-quality sensitivity. It is not an independent sample and does not account for host-correlated systematics or all mass-estimator uncertainty.'
(R/'results/dwarf_sensitivity.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
