"""Independent direct-algebra diagnostic audit, no analysis module imports."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
B=Path(__file__).resolve().parent.parent;D=B/'dm_mond_repairs';O=Path(__file__).resolve().parent
J=json.loads((D/'results/results.json').read_text());Q=json.loads((D/'diagnostics.json').read_text());LN=np.log(10.)

def ae(k,c,b,r):
 if k in ['RAR','simple','standard','nu_n']:return np.zeros_like(b)+10**c[0]
 s=b*r;U=b*r*r
 if k=='two_log_bumps':
  q=np.log10(s);return 10**c[0]+10**c[1]*np.exp(-((q-c[2])/10**c[3])**2/2)+10**c[4]*np.exp(-((q-c[5])/10**c[6])**2/2)
 S=1/(1+(10**c[5]/r)**4);W=1/(1+(b/100)**2)
 if k=='mixed_transient':return 10**c[0]+(10**c[1]*U**-c[3]*S*W+10**c[2]*s**c[4])/(1+r/10**c[6])
 return 10**c[0]+S*W*(10**c[1]*U**-c[3]+10**c[2]*s**c[4])*s/(s+10**c[6])

def lognu(k,c,b,r):
 a=ae(k,c,b,r)
 if k in ['simple','standard','nu_n']:
  n=1 if k=='simple' else 2 if k=='standard' else c[1];z=(a/b)**n
  return np.log1p(2*z/(1+np.sqrt(1+4*z)))/n
 return -np.log(-np.expm1(-np.sqrt(b/a)))

def pred(k,c,b,r):return np.log10(b)+lognu(k,c,b,r)/LN

def logextra(k,c,b,r):
 if k in ['simple','standard','nu_n']:return np.log(np.expm1(lognu(k,c,b,r)))/LN
 t=np.sqrt(b/ae(k,c,b,r));return -(t+np.log(-np.expm1(-t)))/LN
out={};G=6.67430e-11;MS=1.98847e30;KPC=3.085677581491367e19
xx=np.linspace(-20,20,801);uu=np.linspace(-9,9,181);X,U=np.meshgrid(xx,uu)
for k,model in J['models'].items():
 c=model['coefficients'];yy=pred(k,c,10**X,10**U);sl=np.diff(yy,axis=1)/np.diff(xx);num=int(np.sum(sl < -1e-7));low=float(sl.min())
 d=Q['models'][k]['fixed_radius_monotonicity'];assert num==d['negative_intervals'];assert abs(low-d['minimum_logarithmic_slope'])<1e-9
 ur=np.linspace(-9,9,2001);radius=10**ur;bad=0;tested=0;minimum=np.inf;badm=[]
 for mass in 10.**np.arange(-6,17,2):
  b=G*MS*mass/(KPC*radius)**2/1e-10;v=np.diff(logextra(k,c,b,radius))/np.diff(ur);valid=np.isfinite(v);tested+=sum(valid);bad+=sum(v[valid]<-1e-6);minimum=min(minimum,min(v[valid]));
  if np.any(v[valid]<-1e-6):badm.append(mass)
 d=Q['models'][k]['point_source_inferred_extra_mass'];assert int(bad)==d['negative_intervals'];assert tested==d['finite_slope_intervals'];assert np.isclose(minimum,d['minimum_logarithmic_slope'],rtol=1e-8,atol=1e-8)
 observed=list(csv.DictReader((D/'input.csv').open()));ob=np.array([float(v['gbar_m_s2'])/1e-10 for v in observed]);orr=np.array([float(v['radius_kpc']) for v in observed]);eps=1e-5;oslope=(pred(k,c,ob*10**eps,orr)-pred(k,c,ob*10**(-eps),orr))/(2*eps);onum=int(np.sum(oslope< -1e-7))
 out[k]={'observed_coordinate_negative_response_intervals':onum,'minimum_observed_coordinate_response_slope':float(oslope.min()),'fixed_radius_negative_response_intervals':num,'minimum_fixed_radius_response_slope':low,'point_source_negative_extra_density_intervals':int(bad),'tested_point_source_intervals':int(tested)}
# Check each delivered solar scalar probe, including robust tiny excess logs.
maxsolar=0.
for row in csv.DictReader((D/'solar_stress.csv').open()):
 k=row['model'];c=J['models'][k]['coefficients'];r=float(row['r_AU'])*149597870700./KPC;b=G*MS/(r*KPC)**2/1e-10
 aa=float(ae(k,c,np.array(b),np.array(r)))*1e-10;logex=float(logextra(k,c,np.array(b),np.array(r)));err=abs(logex-float(row['log10_relative_excess']));maxsolar=max(maxsolar,err)
 assert np.isclose(aa,float(row['a_eff_m_s2']),rtol=1e-10);assert np.isclose(logex,float(row['log10_relative_excess']),rtol=1e-10,atol=1e-10)
report={'status':'passed','models':out,'solar_probe_rows':42,'max_solar_log_excess_absolute_difference':maxsolar,'interpretation':'Transient repair passes these grids and has an analytic nonnegative point-source extra-density proof for all declared positive coefficients and alpha,beta<1. This is not an arbitrary-profile, field-theory, lensing, binary or external-field solar validation. Taper and bumps fail fixed-radius response monotonicity; bumps also fail point-source density positivity. Formal far-field BTFR normalization does not establish observed finite-radius BTFR agreement.','diagnostic_result_sha256':hashlib.sha256((D/'diagnostics.json').read_bytes()).hexdigest()}
(O/'physics_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
