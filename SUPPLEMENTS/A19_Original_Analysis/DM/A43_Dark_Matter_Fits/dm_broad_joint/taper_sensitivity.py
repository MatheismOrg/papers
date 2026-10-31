"""Post-search high-acceleration taper sensitivity; no parameters refitted.

The fixed g_cut=1e-8 m/s^2 is a chosen regulator, not a learned physical scale.
It is evaluated after the main search and is not an independently selected model.
This checks a central-force Newtonian limit only, not a relativistic solar test.
"""
from pathlib import Path
import csv,json
import numpy as np
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
j=json.loads((OUT/'broad_results.json').read_text())
rr=list(csv.DictReader((OUT/'RAR_potential_quad_predictions.csv').open()))
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rr]);u=np.log10([float(r['radius_kpc']) for r in rr]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rr])
pop=np.array([r['population'] for r in rr]);grp=np.array([r['group'] for r in rr]);fold=np.array([int(r['fold']) for r in rr]);I=np.arange(len(rr))
cut=1e-8
def pred(x,u,c):
 q=x+u;A=(c[0]+c[1]*q+c[2]*q*q)/(1+(1e-10*10**x/cut)**2)
 return x-np.log10(-np.expm1(-10**((x-A)/2)))
def metrics(pr,idx):
 e=pr-y[idx]
 return {'rows':len(idx),'median_absolute_acceleration_fraction':float(np.median(abs(10**e-1))),
  'rms_log10_acceleration_dex':float(np.sqrt(np.mean(e**2))),
  'equal_group_rms_dex':float(np.sqrt(np.mean([np.mean(e[grp[idx]==g]**2) for g in set(grp[idx])])))}
pr=np.zeros(len(rr))
for k,c in enumerate(j['models']['RAR_potential_quad']['fold_coefficients']):
 m=fold==k;pr[m]=pred(x[m],u[m],c)
res={'status':'Post-search regulator sensitivity with existing fitted coefficients; no refitting, independent data, or confirmatory model selection.',
 'formula':'g=gbar/[1-exp(-sqrt(gbar/a_eff))]; log10(a_eff/1e-10)=(a+b*q+c*q*q)/(1+(gbar/1e-8)^2); q=x+u',
 'chosen_g_cut_m_s2':cut,'cross_validation':{},'leave_population_out':{},'solar_central_force_smoke_checks':[]}
for p in sorted(set(pop)):
 idx=I[pop==p];res['cross_validation'][p]=metrics(pr[idx],idx)
res['balanced_group_cv_rms_dex']=float(np.sqrt(np.mean([v['equal_group_rms_dex']**2 for v in res['cross_validation'].values()])))
for p,v in j['leave_population_out'].items():
 idx=I[pop==p];c=v['models']['RAR_potential_quad']['coefficients'];pp=pred(x[idx],u[idx],c)
 res['leave_population_out'][p]=metrics(pp,idx)
res['balanced_population_holdout_rms_dex']=float(np.sqrt(np.mean([v['equal_group_rms_dex']**2 for v in res['leave_population_out'].values()])))
c=j['models']['RAR_potential_quad']['all_data_coefficients'];AU=149597870700.;KPC=3.085677581491367e19;GM=1.32712440018e20
for au in [1.,9.58,30.07]:
 r=au*AU;gb=GM/r**2;xx=np.log10(gb/1e-10);uu=np.log10(r/KPC);q=xx+uu
 A=(c[0]+c[1]*q+c[2]*q*q)/(1+(gb/cut)**2);ae=1e-10*10**A;t=np.sqrt(gb/ae)
 # Fractional correction=e^-t/(1-e^-t); for t>700 log10 correction=-t/ln(10).
 logcorrection=float(-t/np.log(10)-np.log10(-np.expm1(-t)))
 res['solar_central_force_smoke_checks'].append({'radius_AU':au,'gbar_m_s2':float(gb),'a_eff_m_s2':float(ae),'log10_fractional_acceleration_correction':logcorrection})
res['solar_scope']='Only spherical central Sun acceleration is evaluated. Does not calculate external-field quadrupole, planetary ephemerides, lensing, relativistic effects or CMB. Radius is far outside the empirical training range.'
(OUT/'taper_sensitivity.json').write_text(json.dumps(res,indent=2))
with (OUT/'tapered_potential_predictions.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['population','group','fold','gbar_m_s2','gobs_m_s2','pred_gobs_m_s2'])
 for i,r in enumerate(rr):w.writerow([pop[i],grp[i],fold[i],r['gbar_m_s2'],r['gobs_m_s2'],1e-10*10**pr[i]])
print(json.dumps(res,indent=2))
