from pathlib import Path
import csv,json,hashlib
import numpy as np
import run as m
R=Path(__file__).resolve().parent;D=json.loads((R/'results/results.json').read_text());LN=np.log(10.)
G=6.67430e-11;MSUN=1.98847e30;KPC=3.085677581491367e19;AU=149597870700.
def exact_pred(k,c,g,r):
 xx=np.log10(g/1e-10);uu=np.log10(r);ae=1e-10*10**m.logae(k,c,np.array([xx]),np.array([uu]))[0]
 if k in ['simple','standard','nu_n']:
  nn=1 if k=='simple' else 2 if k=='standard' else c[1]
  z=(g/ae)**(-nn);lnnu=np.log1p(2*z/(np.sqrt(1+4*z)+1))/nn
  rel=np.expm1(lnnu);lr=np.log10(rel)
 else:
  t=np.sqrt(g/ae);lr=-(t+np.log(-np.expm1(-t)))/LN
  rel=10**lr if lr>-320 else 0.
 return g*(1+rel),ae,rel,lr

def log_extra_ratio(k,c,xx,uu):
 if k not in ['simple','standard','nu_n']:
  ae=m.logae(k,c,xx,uu);t=10**np.clip(.5*(xx-ae),-290,290)
  # log(expm1(t))=t+log1p(-exp(-t)); log(-expm1(-t)) stable at low t.
  return -(t+np.log(-np.expm1(-t)))/LN
 # Stable identity avoids subtracting nearly equal force values at high g.
 nn=1 if k=='simple' else 2 if k=='standard' else c[1]
 z=10**np.clip(-nn*(xx-c[0]),-290,290)
 lnu=np.log1p(2*z/(np.sqrt(1+4*z)+1))/nn
 return np.log(np.expm1(lnu))/LN

out={'status':'Illustrative isolated scalar force diagnostics; no external field or observation likelihood.','models':{}}
sol=[];cores=[];far=[]
for k in m.KINDS:
 c=D['models'][k]['coefficients'];s={};a_inf=10**c[0]*1e-10
 s['asymptotic_a_inf_m_s2']=a_inf;s['a_inf_over_canonical_1p2e10']=a_inf/1.2e-10
 s['velocity_ratio_to_canonical_at_fixed_mass']=(a_inf/1.2e-10)**.25
 # fixed positive radius, g sweeping 40 decades; compute local slope finite differences.
 xx=np.linspace(-20,20,801);uu=np.linspace(-9,9,181);X,U=np.meshgrid(xx,uu)
 yy=m.rawpredict(k,c,X.ravel(),U.ravel()).reshape(X.shape);sl=np.diff(yy,axis=1)/np.diff(xx)[None,:]
 s['fixed_radius_monotonicity']={'grid_points':int(X.size),'slope_intervals':int(sl.size),'negative_intervals':int(np.sum(sl < -1e-7)),'minimum_logarithmic_slope':float(sl.min()),'x_range':[-20,20],'u_range':[-9,9]}
 # Derivative at observed input coordinates, to distinguish interpolation-domain problems from extreme extrapolation.
 eps=1e-5;obs_slope=(m.rawpredict(k,c,m.x+eps,m.u)-m.rawpredict(k,c,m.x-eps,m.u))/(2*eps)
 s['observed_input_monotonicity']={'points':len(obs_slope),'negative_slopes':int(np.sum(obs_slope < -1e-7)),'minimum_logarithmic_slope':float(obs_slope.min()),'difference_step_dex':eps}
 j,i=np.unravel_index(np.argmin(sl),sl.shape)
 s['fixed_radius_monotonicity']['minimum_location']={'gbar_m_s2':float(10**((xx[i]+xx[i+1])/2-10)),'radius_kpc':float(10**uu[j])}
 # Solar masses including tiny and massive point sources. Algebraic density test only.
 ugrid=np.linspace(-9,9,2001);masses=10.**np.arange(-6,17,2);bad=0;tested=0;min_slope=np.inf;neg_cases=[]
 for mass in masses:
  gg=G*MSUN*mass/(KPC*10**ugrid)**2;xxg=np.log10(gg/1e-10);lmer=log_extra_ratio(k,c,xxg,ugrid);delta=np.diff(lmer)
  finite=np.isfinite(delta);v=delta[finite]/np.diff(ugrid)[finite];tested+=int(finite.sum());bad+=int(np.sum(v < -1e-6))
  if len(v):min_slope=min(min_slope,float(v.min()))
  if np.any(v < -1e-6):neg_cases.append(float(mass))
 s['point_source_inferred_extra_mass']={'masses_solar':masses.tolist(),'radius_kpc_range':[1e-9,1e9],'radii_per_mass':len(ugrid),'finite_slope_intervals':tested,'negative_intervals':bad,'minimum_logarithmic_slope':min_slope,'masses_with_negative_intervals':neg_cases,'caveat':'Stable log extra mass used for every family; inferred density is a point-source diagnostic, not a full field equation.'}
 for au in [.4,1,10,30,1000,10000]:
  rr=au*AU/KPC;gg=G*MSUN/(au*AU)**2;pr,ae,rel,lr=exact_pred(k,c,gg,rr)
  sol.append({'model':k,'r_AU':au,'gbar_m_s2':gg,'gpred_m_s2':float(pr),'a_eff_m_s2':float(ae),'relative_excess':float(rel),'log10_relative_excess':float(lr),'extra_acceleration_m_s2':float(gg*rel)})
 # Smooth core path g_b=a_ref*r/kpc. This is a shape test, no density fitted.
 for uu0 in np.linspace(-12,0,49):
  xx0=uu0;pr=m.rawpredict(k,c,np.array([xx0]),np.array([uu0]))[0]
  cores.append({'model':k,'radius_kpc':10**uu0,'gbar_m_s2':1e-10*10**xx0,'gpred_m_s2':1e-10*10**pr})
 # Point-mass finite-r normalized V^4, distinct from formal asymptote.
 for mass in [1,1e6,1e8,1e10,1e12]:
  for rr in [1,10,100,1000,10000,1e6,1e9]:
   gg=G*MSUN*mass/(KPC*rr)**2;pr,ae,rel,lr=exact_pred(k,c,gg,rr)
   far.append({'model':k,'mass_solar':mass,'radius_kpc':rr,'a_eff_over_a_inf':float(ae/a_inf),'V4_over_GM_a_inf':float((KPC*rr)**2*pr**2/(G*MSUN*mass*a_inf))})
 s['boundary_contacts_all']=[i for i,(z,lo,hi) in enumerate(zip(c,m.SPEC[k]['lo'],m.SPEC[k]['hi'])) if min(abs(z-lo),abs(z-hi))<1e-6]
 out['models'][k]=s
for fn,rr in [('solar_stress.csv',sol),('smooth_core.csv',cores),('finite_radius_btfr.csv',far)]:
 with (R/fn).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
(R/'diagnostics.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
