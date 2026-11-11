"""SA73 T50: how much of the dwarfs' error and of its trend with size is measurement noise? Take the formula's prediction
as the truth, perturb each dwarf's dispersion, half-light radius, distance and luminosity by its published errors (LVDB),
recompute the 'observed' acceleration (sigma^2 / R_e) and the prediction at the perturbed inputs, 2000 times."""
import numpy as np, json
import harness_T48 as M8
M=M8.P; H=M.H; B=M.B; B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
rows=H.dwarf_rows(); a=H.arr(rows)
rr,_=H.V4.load_rows(str(H.PKG/'data/LVDB/dwarf_all.csv'),str(H.PKG/'blind_test')); R={x['key']:x for x in rr}
def f(x,k,d=np.nan):
    try: return float(x[k])
    except: return d
es=[];er=[];ed=[];em=[]
for n in a['name']:
    x=R[n]; s=f(x,'vlos_sigma'); es.append(.5*(f(x,'vlos_sigma_em',0.)+f(x,'vlos_sigma_ep',0.))/s)
    er.append(.5*(f(x,'rhalf_em',0.)+f(x,'rhalf_ep',0.))/f(x,'rhalf')); ed.append(.5*(f(x,'distance_modulus_em',.1)+f(x,'distance_modulus_ep',.1)))
    em.append(.5*(f(x,'apparent_magnitude_v_em',.2)+f(x,'apparent_magnitude_v_ep',.2)))
es,er,ed,em=map(lambda v:np.nan_to_num(np.array(v),nan=0.),(es,er,ed,em))
print('median frac err sigma',np.median(es).round(2),'Re',np.median(er).round(2),'dm',np.median(ed).round(2),'mV',np.median(em).round(2))
g0=H.gpred(M.MODEL,th,a); e_obs=np.log10(a['go']/g0)
def stats(e,lre): return np.sqrt(np.mean(e**2)), np.corrcoef(lre,e)[0,1], np.polyfit(lre,e,1)[0]
print('observed rms, corr, slope vs logRe', np.round(stats(e_obs,np.log10(a['re'])),3))
rng=np.random.default_rng(50); out=[]
for it in range(2000):
    fs=np.clip(1+es*rng.standard_normal(len(es)),.05,None); dd=10**(0.2*ed*rng.standard_normal(len(es)))   # distance factor
    fr=np.clip(1+er*rng.standard_normal(len(es)),.1,None)*dd; fm=10**(-0.4*em*rng.standard_normal(len(es)))*dd**2
    b=dict(a); b['re']=a['re']*fr; b['r']=a['r']*fr; b['ms']=a['ms']*fm; b['mp']=a['mp']*fm; b['gb']=a['gb']*fm/fr**2
    go=g0*fs**2/fr; e=np.log10(go/H.gpred(M.MODEL,th,b)); out.append(stats(e,np.log10(b['re'])))
out=np.array(out); print('noise only: rms median',np.median(out[:,0]).round(3),'[16,84]',np.percentile(out[:,0],[16,84]).round(3),
  '| corr',np.median(out[:,1]).round(3),np.percentile(out[:,1],[2.5,97.5]).round(3),'| slope',np.median(out[:,2]).round(3),np.percentile(out[:,2],[2.5,97.5]).round(3))
r_obs=stats(e_obs,np.log10(a['re']))[0]; rn=np.median(out[:,0]); print('intrinsic rms (quadrature)',np.sqrt(max(r_obs**2-rn**2,0)).round(3))
json.dump(dict(observed=list(map(float,stats(e_obs,np.log10(a['re'])))),noise_rms=float(rn),noise_corr=float(np.median(out[:,1])),noise_slope=float(np.median(out[:,2]))),open('res50_noise_dwarfs.json','w'))
