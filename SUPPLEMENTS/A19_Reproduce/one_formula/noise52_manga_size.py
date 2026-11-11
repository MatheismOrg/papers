"""SA73 T52: would errors in the MaNGA BCGs' half-light radii alone produce the residual trends? Truth = the formula;
each galaxy's R_e redrawn with a log-normal error s (the stellar g_bar at the fixed measured radius scaled as R_e^-2,
the size entering the formula too), plus the catalog's g_obs error; 1000 draws per s."""
import numpy as np, csv, json
import harness_T48 as M8
M=M8.P; H=M.H; B=M.B; B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
rows=[r for r in H.cluster_rows() if r['kind']=='MaNGA_BCG']; a=H.arr(rows)
E={r['system']:float(r['loggobs_err_dex']) for r in csv.DictReader(open(H.PKG/'data/Clusters/manga_bcg_standardized.csv'))}
eo=np.array([E.get(n,E.get('MaNGA_'+n,np.nan)) for n in a['name']]); eo=np.nan_to_num(eo,nan=np.nanmedian(eo))
g0=H.gpred(M.MODEL,th,a); e=np.log10(a['go']/g0)
obs=dict(rms=float(np.sqrt(np.mean(e**2))),corr_Re=float(np.corrcoef(np.log10(a['re']),e)[0,1]),corr_gbar=float(np.corrcoef(np.log10(a['gb']),e)[0,1]),slope_Re=float(np.polyfit(np.log10(a['re']),e,1)[0]))
print('observed',{k:round(v,3) for k,v in obs.items()}); rng=np.random.default_rng(52); out={'observed':obs}
for s in (0.,.05,.09,.12):
    R=[]
    for it in range(1000):
        d=s*rng.standard_normal(50); b=dict(a); b['re']=a['re']*10**d; b['gb']=a['gb']*10**(-2*d)
        ee=np.log10(g0*10**(eo*rng.standard_normal(50))/H.gpred(M.MODEL,th,b))
        R.append((np.sqrt(np.mean(ee**2)),np.corrcoef(np.log10(b['re']),ee)[0,1],np.corrcoef(np.log10(b['gb']),ee)[0,1],np.polyfit(np.log10(b['re']),ee,1)[0]))
    R=np.array(R); out[str(s)]=dict(rms=float(np.median(R[:,0])),corr_Re=np.percentile(R[:,1],[2.5,50,97.5]).round(3).tolist(),corr_gbar=np.percentile(R[:,2],[2.5,50,97.5]).round(3).tolist(),slope_Re=float(np.median(R[:,3])))
    print('size error',s,json.dumps(out[str(s)]),flush=True)
json.dump(out,open('res52_manga_size.json','w'),indent=1)
