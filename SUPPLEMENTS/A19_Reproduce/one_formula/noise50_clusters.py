"""SA73 T50: noise floors for MaNGA BCGs and CLASH from the catalogs' own errors (loggbar_err, loggobs_err): truth = the
formula; perturb log g_obs and log g_bar (and M* with g_bar for MaNGA, one point each), 2000 times. Also dwarfs with
errors symmetric in the log (robustness)."""
import numpy as np, csv, json
import harness_T48 as M8
M=M8.P; H=M.H; B=M.B; B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
D=H.PKG/'data/Clusters'
E={}
for fn in ('manga_bcg_standardized.csv','clash_standardized.csv'):
    for r in csv.DictReader(open(D/fn)): E[(r['canonical_cluster_name'] if 'clash' in fn else r['system'], round(float(r['radius_kpc']),1))]=(float(r['loggbar_err_dex']),float(r['loggobs_err_dex']))
cl=H.cluster_rows(); rng=np.random.default_rng(51); out={}
for pop,sel in (('MaNGA',lambda r:r['kind']=='MaNGA_BCG'),('CLASH_centre',lambda r:r['kind']=='CLASH_centre'),('CLASH_outer',lambda r:r['kind']=='CLASH_outer')):
    rows=[r for r in cl if sel(r)]; a=H.arr(rows); keys=list(E)
    def err(r):
        k=(r['name'],round(float(r['r']),1))
        if k in E: return E[k]
        c=[v for (n,rad),v in E.items() if (n==r['name'] or n.endswith(r['name']) or r['name'].endswith(n)) and abs(rad-r['r'])<.06*r['r']+.2]
        return c[0] if c else (np.nan,np.nan)
    eb,eo=np.array([err(r) for r in rows]).T; miss=np.isnan(eb).sum(); eb=np.nan_to_num(eb,nan=np.nanmedian(eb)); eo=np.nan_to_num(eo,nan=np.nanmedian(eo))
    g0=H.gpred(M.MODEL,th,a); e_obs=np.log10(a['go']/g0); hs=np.unique(a['name'])
    def rms(e): return np.sqrt(np.mean([np.mean(e[a['name']==h]**2) for h in hs]))
    sims=[]; cors=[]
    for it in range(1000):
        db=eb*rng.standard_normal(len(eb)); do=eo*rng.standard_normal(len(eo)); b=dict(a); b['gb']=a['gb']*10**db
        if pop=='MaNGA': b['ms']=a['ms']*10**db; b['mp']=a['mp']*10**db
        e=np.log10(g0*10**do/H.gpred(M.MODEL,th,b)); sims.append(rms(e)); cors.append(np.corrcoef(np.log10(b['gb']),e)[0,1])
    out[pop]=dict(observed=float(rms(e_obs)),noise=float(np.median(sims)),p16_84=np.percentile(sims,[16,84]).round(4).tolist(),missing_err=int(miss),
                  med_err_gbar=float(np.median(eb)),med_err_gobs=float(np.median(eo)),corr_gbar_obs=float(np.corrcoef(np.log10(a['gb']),e_obs)[0,1]),corr_gbar_noise=np.percentile(cors,[2.5,50,97.5]).round(3).tolist())
    print(pop,json.dumps(out[pop]),flush=True)
json.dump(out,open('res50_noise_clusters.json','w'),indent=1)
