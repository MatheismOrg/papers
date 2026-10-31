"""Frozen X-COP transfer. Run from the complete A43 supplement, no network."""
from pathlib import Path
import csv, json, hashlib, sys
import numpy as np
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R.parent/'dm_unified_limits'))
import run as model
src=R.parent/'dm_fixed_radius/fixed_radius_inputs.csv'
fitfile=R.parent/'dm_unified_limits/results/results.json'
fit=json.loads(fitfile.read_text())
assert set(fit['models'])==set(model.KINDS) and 'optimizer_calls' in fit
rows=list(csv.DictReader(src.open()))
group=np.array([r['cluster'] for r in rows])
u=np.log10([float(r['radius_kpc']) for r in rows])
assert len(rows)==21 and len(np.unique(group))==7
out={'protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'fit_result_sha256':hashlib.sha256(fitfile.read_bytes()).hexdigest(),'cases':{}}
records=[]
for target in ['forward','nfw','einasto']:
 for baryon in ['gbar','gbar_no_stars','gbar_full_projected_stars']:
  xx=np.log10([float(r[baryon])/1e-10 for r in rows]);yy=np.log10([float(r['gobs_'+target])/1e-10 for r in rows]);case=target+'__'+baryon
  out['cases'][case]={}
  for k in model.KINDS:
   cc=fit['models'][k]['coefficients'];pr=model.rawpredict(k,cc,xx,u);err=pr-yy
   out['cases'][case][k]={'RMS_dex':float(np.sqrt(np.mean([np.mean(err[group==g]**2) for g in np.unique(group)]))),'median_acceleration_error_percent':float(100*np.median(abs(10**err-1))),'per_cluster_RMS_dex':{g:float(np.sqrt(np.mean(err[group==g]**2))) for g in np.unique(group)}}
   for i,r in enumerate(rows):records.append({'model':k,'case':case,'cluster':r['cluster'],'radius_kpc':float(r['radius_kpc']),'gbar_m_s2':float(r[baryon]),'gobs_m_s2':float(r['gobs_'+target]),'gpred_m_s2':float(1e-10*10**pr[i]),'log10_residual':float(err[i])})
(R/'results.json').write_text(json.dumps(out,indent=2))
with (R/'predictions.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
print(json.dumps(out['cases']['forward__gbar'],indent=2))
