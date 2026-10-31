from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.special import logsumexp, expit

R=Path(__file__).resolve().parent;J=json.loads((R/'results/results.json').read_text());LN=np.log(10.)
G=6.67430e-11;MSUN=1.98847e30;KPC=3.0856775814913673e19;AU=149597870700.

def evaluate(c,x,u):
    x,u=np.broadcast_arrays(np.asarray(x,float),np.asarray(u,float));m=x+2*u;q=x+u
    # Independently written dimensional-log components; run.py is not imported.
    h=2*LN*(x-2);j=4*LN*(c[5]-u);o=LN*(u-c[6])
    logs=np.stack([np.zeros_like(x)+LN*c[0],LN*(c[1]-c[3]*m)-np.logaddexp(0,h)-np.logaddexp(0,j),LN*(c[2]+c[4]*q)-np.logaddexp(0,o)])
    la=logsumexp(logs,axis=0);w=np.exp(logs-la)
    a_slope=w[1]*(-c[3]-2*expit(h))+w[2]*c[4]
    # Fixed-source radius slope: m constant, q decreases at rate1.
    a_r_slope=w[1]*(4*expit(j)+4*expit(h))+w[2]*(-c[4]-expit(o))
    t=np.exp(np.clip(.5*(LN*x-la),-700,700))
    den=-np.expm1(-t);logg=LN*x-np.log(den)
    excess=np.zeros_like(t);mask=t<700;excess[mask]=1/np.expm1(t[mask])
    dlogg_dloggb=1-.5*(1-a_slope)*(t*excess)
    return {'log10_gpred_over_aref':logg/LN,'a_eff':np.exp(la)*1e-10,'relative_excess':excess,'a_r_slope':a_r_slope,'dlogg_dloggb':dlogg_dloggb}

def write_csv(path,rows):
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

solar=[];core=[];report={'constants':{'G':G,'solar_mass_kg':MSUN,'kpc_m':KPC,'AU_m':AU},'datasets':{}}
for name,z in J['models_by_dataset'].items():
    m=z['models']['mixed_regular'];c=m['coefficients'];dd={}
    for r in [.4,1.,10.,30.,1000.,10000.]:
        rm=r*AU;gb=G*MSUN/rm**2;v=evaluate(c,np.log10(gb/1e-10),np.log10(rm/KPC))
        solar.append({'dataset':name,'radius_AU':r,'gbar_m_s2':gb,'a_eff_m_s2':float(v['a_eff']),'relative_acceleration_excess':float(v['relative_excess']),'relative_circular_speed_excess':float(np.sqrt(1+v['relative_excess'])-1),'extra_acceleration_m_s2':float(gb*v['relative_excess'])})
    # Smooth central density1 Msun/pc^3. At the center g_b~r, M_b~r^3.
    rho=MSUN/(KPC/1000)**3;rr=np.logspace(-12,-1,111);gb=4*np.pi*G*rho*(rr*KPC)/3
    vv=evaluate(c,np.log10(gb/1e-10),np.log10(rr));pg=1e-10*10**vv['log10_gpred_over_aref']
    for i in range(len(rr)):core.append({'dataset':name,'radius_kpc':float(rr[i]),'gbar_m_s2':float(gb[i]),'gpred_m_s2':float(pg[i]),'a_eff_m_s2':float(vv['a_eff'][i])})
    dd['smooth_core']={'central_density_Msun_pc3':1.,'smallest_radius_kpc':float(rr[0]),'smallest_gpred_m_s2':float(pg[0]),'small_radius_log_slope':float(np.polyfit(np.log10(rr[:10]),np.log10(pg[:10]),1)[0]),'analytic_limit':'g_pred~sqrt(a_ref*C0*gbar)~sqrt(r), except beta=0 adds finiteC2 to scale; g_pred->0.'}
    # Exact analytic density-sign factor checked on a wide numerical domain.
    point=[]
    for mass in [1.,1e3,1e6,1e9,1e12,1e15]:
        rr=np.logspace(-12,7,1501);gb=G*mass*MSUN/(rr*KPC)**2
        vv=evaluate(c,np.log10(gb/1e-10),np.log10(rr));factor=1+.5*vv['a_r_slope']
        # M_extra=M/(exp(t)-1). Its derivative has the sign of factor.
        mextra=mass*vv['relative_excess'];decreases=np.diff(mextra)<-1e-10*np.maximum(mextra[1:],1)
        point.append({'mass_Msun':mass,'minimum_extra_mass_derivative_sign_factor':float(np.min(factor)),'numerically_decreasing_extra_mass_steps':int(np.sum(decreases))})
    dd['point_mass_density']={'radius_grid_kpc':[1e-12,1e7],'radii_per_mass':1501,'results':point,'analytic_guarantee':'Each positive a_eff component has radius-log-slope>-2 for0<=beta<=.95. Thus dM_extra/dr>=0 for every fixed point-source mass and positive radius.'}
    xx,uu=np.meshgrid(np.linspace(-30,15,251),np.linspace(-12,7,241),indexing='ij');vv=evaluate(c,xx,uu)
    dd['input_monotonicity_grid']={'log10_gbar_over_aref':[-30,15],'log10_radius_kpc':[-12,7],'points':int(xx.size),'minimum_dlogg_dloggb':float(np.min(vv['dlogg_dloggb'])),'negative_count':int(np.sum(vv['dlogg_dloggb']<0))}
    # Isolated fixed-M asymptote, with m=0 as a numerical illustration.
    ur=np.array([8.,9.,12.,13.]);fr=evaluate(c,-2*ur,ur)['log10_gpred_over_aref']
    dd['far_field_slopes']=[float(fr[1]-fr[0]),float(fr[3]-fr[2])]
    dd['coefficients']=c;dd['length_l_pc']=1000*10**c[5];dd['range_R_kpc']=10**c[6]
    dd['all_data_bound_contacts']=[i for i,(v,lo,hi) in enumerate(zip(c,J['spec']['mixed_regular']['lo'],J['spec']['mixed_regular']['hi'])) if min(abs(v-lo),abs(v-hi))<1e-6]
    selected=[]
    for split in sorted({r['split'] for r in J['fit_log'] if r['dataset']==name}):
        rr=[r for r in J['fit_log'] if r['dataset']==name and r['split']==split];best=min(rr,key=lambda z:z['cost'])
        selected.append({'split':split,'selected_start':best['start_index'],'success':best['success'],'active_bounds':best['active_bounds'],'bound_contacts_absolute_1e-6':[i for i,(v,lo,hi) in enumerate(zip(best['coefficients'],J['spec']['mixed_regular']['lo'],J['spec']['mixed_regular']['hi'])) if min(abs(v-lo),abs(v-hi))<1e-6]})
    dd['selected_fit_bound_diagnostics']=selected
    report['datasets'][name]=dd

write_csv(R/'solar_stress.csv',solar);write_csv(R/'smooth_core.csv',core)
report['limitations']='No planetary or binary likelihood; no Galactic external-field quadrupole. Large wide-binary MOND-like response from C0 remains. Positive inferred density is proved only for point-source radial profiles, not arbitrary extended systems or a covariant theory.'
report['results_sha256']=hashlib.sha256((R/'results/results.json').read_bytes()).hexdigest()
(R/'diagnostics.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:{z:d[z] for z in ['length_l_pc','range_R_kpc','far_field_slopes','input_monotonicity_grid','smooth_core']} for k,d in report['datasets'].items()},indent=2))
print('Solar stress rows:',solar)
