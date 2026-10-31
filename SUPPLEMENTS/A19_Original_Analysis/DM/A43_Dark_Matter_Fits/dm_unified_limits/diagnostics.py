"""Mathematical spot checks of frozen fits; no observational Solar likelihood."""
from pathlib import Path
import json,csv
import numpy as np
import run as b
R=Path(__file__).resolve().parent;J=json.loads((R/'results/results.json').read_text());GM=1.32712440018e20;AU=149597870700.;KPC=3.085677581491367e19;G=6.67430e-11;MSUN=GM/G;MREF=1e-10*KPC*KPC/G/MSUN
solar=[];diag={}
for k,d in J['models'].items():
 c=np.array(d['coefficients']);a=np.array([.4,1.,10.,30.]);r=a*AU;gb=GM/r**2;x=np.log10(gb/1e-10);u=np.log10(r/KPC);pr=b.rawpredict(k,c,x,u);delta=np.expm1((pr-x)*np.log(10.))
 for i in range(len(a)):solar.append({'model':k,'radius_AU':a[i],'gbar_m_s2':gb[i],'gpred_m_s2':10**pr[i]*1e-10,'fractional_excess':delta[i],'anomalous_acceleration_m_s2':gb[i]*delta[i]})
 # Point-source total inferred mass derivative sign, broad but declared finite grid.
 logsM=np.linspace(0,15,31);uu=np.linspace(-9,4,161);MM,UU=np.meshgrid(logsM,uu,indexing='ij');XX=MM-np.log10(MREF)-2*UU;h=1e-4
 predhi=b.rawpredict(k,c,XX.ravel()-2*h,UU.ravel()+h);predlo=b.rawpredict(k,c,XX.ravel()+2*h,UU.ravel()-h);slope=(predhi-predlo)/(2*h);md=2+slope
 # Fixed-radius response monotonicity across declared extended two-input grid.
 xx=np.linspace(-8,8,181);uu2=np.linspace(-9,4,101);XX2,UU2=np.meshgrid(xx,uu2,indexing='ij');der=(b.rawpredict(k,c,XX2.ravel()+h,UU2.ravel())-b.rawpredict(k,c,XX2.ravel()-h,UU2.ravel()))/(2*h)
 ddata=(b.rawpredict(k,c,b.x+h,b.u)-b.rawpredict(k,c,b.x-h,b.u))/(2*h)
 z={'point_source_grid':{'mass_Msun_log10':[0,15,31],'radius_kpc_log10':[-9,4,161],'min_dlog_inferred_mass_dlog_radius':float(min(md)),'negative_below_minus_1e_5_count':int(sum(md<-1e-5)),'grid_points':len(md)},'fixed_radius_monotonicity':{'log10_gbar_over_aref_range':[-8,8,181],'log10_radius_kpc_range':[-9,4,101],'min_dlogg_dloggbar':float(min(der)),'negative_count':int(sum(der<-1e-5)),'grid_points':len(der),'observed_input_min_derivative':float(min(ddata)),'observed_input_negative_count':int(sum(ddata<-1e-5))}}
 coreu=np.array([-8.,-10.,-12.,-14.]);corepr=b.rawpredict(k,c,coreu,coreu)
 z['uniform_density_core_path']={'definition':'gbar/a_ref=r/kpc; finite central baryonic density, r->0','log10_radius_kpc':coreu.tolist(),'log10_gpred_over_aref':corepr.tolist(),'last_interval_dlogg_dlogr':float((corepr[-1]-corepr[-2])/(coreu[-1]-coreu[-2]))}
 if k in b.NEW:
  mass=np.array([1.,1e5,1e8,1e10,1e12,1e15]);mm=np.log10(mass/MREF)
  def asym(logm):
   if k=='mass_bounded':zz=np.tanh((logm-c[3])/c[4]);return 10**(c[0]+c[1]*zz+c[2]*zz*zz)
   val=10**c[0]+10**(c[1]-c[3]*logm)
   if k in ['mass_power','mass_additive','mass_screened'] or k=='index_power' and abs(c[5]-2)<1e-9:val+=10**(c[2]+c[4]*logm)
   return val
  aa=asym(mm);ad=(np.log10(asym(mm+h))-np.log10(asym(mm-h)))/(2*h);z['point_mass_asymptotic_BTFR']={'mass_Msun':mass.tolist(),'a_infinity_m_s2':(aa*1e-10).tolist(),'dlog_vflat_dlog_M':((1+ad)/4).tolist(),'formal_dlog_M_dlog_vflat':(4/(1+ad)).tolist(),'meaning':'v_flat^4=G*M*a_infinity(M); algebraic asymptote only, not a fitted BTFR likelihood'}
 diag[k]=z
with (R/'results/solar_spot_checks.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(solar[0]));w.writeheader();w.writerows(solar)
out={'status':'frozen-fit mathematical diagnostics, not Solar-System or BTFR likelihood fits','Mref_Msun':MREF,'solar_constants':{'GM_sun_SI':GM,'AU_m':AU,'kpc_m':KPC},'models':diag}
(R/'diagnostics.json').write_text(json.dumps(out,indent=2));print('Saved solar and mathematical diagnostics')
for k in b.KINDS:
 vals=[r for r in solar if r['model']==k];print(k,'Solar excess@10AU',next(r['fractional_excess'] for r in vals if r['radius_AU']==10),'density neg',diag[k]['point_source_grid']['negative_below_minus_1e_5_count'],'monotonicity observed neg',diag[k]['fixed_radius_monotonicity']['observed_input_negative_count'])
