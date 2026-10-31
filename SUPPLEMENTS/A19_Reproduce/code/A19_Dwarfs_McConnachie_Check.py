"""A19 supplement: independent check of Sections 3.3 and 4 on the original 44 dwarf spheroidals.
Input: dwarfs_benchmark.csv built by the original pipeline from McConnachie's January 2021 FITS
catalog (McConnachie 2012, AJ 144, 4), circularized radii, mass-to-light ratio 2; and the
results file of A19_Disk_Fit_SPARC.py. Written independently of the original code.
Usage: python3 A19_Dwarfs_McConnachie_Check.py dwarfs_benchmark.csv A19_Disk_Fit_SPARC_Results.json out.json
"""
import sys, os, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import g_nfw, rar, moster_mref, eff_frac, FB, G, ACC, NLG, MOG
from scipy.optimize import minimize, minimize_scalar
rms=lambda v: float(np.sqrt(np.mean(np.square(v))))
def main(bench, diskjson, out):
    D=json.load(open(diskjson))['frozen_on_all_disks']
    r=list(csv.DictReader(open(bench))); f=lambda k: np.array([float(x[k]) for x in r]); n=len(r)
    R=f('radius_kpc'); gb=f('gbar_m_s2'); go=f('gobs_m_s2'); a=f('projected_Re_circularized_kpc'); L=f('luminosity_V_Lsun'); Ms=f('stellar_mass_Msun'); lo=np.log10(go)
    H=lambda p: gb+g_nfw(R,Ms/FB+10**p[0]*(Ms/1e10)**p[1]*(a/3)**p[2])
    res={'stellar_rule':lo-np.log10(H(D['stellar_ABC'])),'baryonic_rule':lo-np.log10(H(D['baryonic_ABC'])),
         'rar_disk_a0':lo-np.log10(rar(gb,D['rar_a0'])),'moster':lo-np.log10(np.array([g+g_nfw(rr,max(moster_mref(m),m/FB)) for g,rr,m in zip(gb,R,Ms)]))}
    cl=L>=1e5; out_={'n':n,'section3_rmse_dex':{k:round(rms(v),5) for k,v in res.items()},
        'classical':{k:round(rms(v[cl]),4) for k,v in res.items()},'ultrafaint':{k:round(rms(v[~cl]),4) for k,v in res.items()}}
    ce={}
    for ml in (2.0,3.0):
        M=ml*L; fr=R**3/(R*R+a*a)**1.5; gp=G*fr*M/R**2*ACC
        ce['ML_%g'%ml]={'exceed_NLG':int(np.sum(go/gp>1+NLG[0]/fr)),'exceed_MOG':int(np.sum(go/gp>1+MOG[0]/fr))}
    out_['section4_ceiling_counts']=ce
    M=2*L; fenc=R**3/(R*R+a*a)**1.5; gN=G*M*fenc/R**2*ACC; MU=np.logspace(-3,3,121)
    F={k:np.array([[eff_frac(R[i],a[i],mu,k) for mu in MU] for i in range(n)]) for k in ('NLG','MOG')}
    # gN (1 + alpha F / fenc): F is per unit of the total baryonic mass, fenc = M_b(<R) / M_b (turn 13 correction, as in A19_Dwarfs_LVDB.py)
    kp=lambda kind,al,lmu,idx: gN[idx]*(1+al*np.array([np.interp(lmu,np.log10(MU),F[kind][i]) for i in idx])/fenc[idx])
    s4={'rar_published_a0':round(rms(lo-np.log10(rar(gN,1.2e-10))),4)}
    for kind,(al,mu) in (('NLG',NLG),('MOG',MOG)): s4[kind+'_published']=round(rms(lo-np.log10(kp(kind,al,np.log10(mu),np.arange(n)))),4)
    cv={'rar_refit':[],'NLG_refit':[],'MOG_refit':[]}
    for seed in (1,2,3):
        idx=np.random.default_rng(seed).permutation(n); pr={k:np.zeros(n) for k in cv}
        for fo in np.array_split(idx,5):
            tr=np.array([i for i in idx if i not in set(fo)])
            la=minimize_scalar(lambda x: np.mean((lo[tr]-np.log10(rar(gN[tr],10**x)))**2),bounds=(-12,-7),method='bounded').x
            pr['rar_refit'][fo]=lo[fo]-np.log10(rar(gN[fo],10**la))
            for kind in ('NLG','MOG'):
                best=min((minimize(lambda p: np.mean((lo[tr]-np.log10(kp(kind,np.exp(p[0]),p[1],tr)))**2),[np.log(10),l0],method='L-BFGS-B',bounds=[(np.log(1e-3),np.log(1e4)),(-3,3)]) for l0 in np.log10(MU)[::4]),key=lambda o:o.fun)
                pr[kind+'_refit'][fo]=lo[fo]-np.log10(kp(kind,np.exp(best.x[0]),best.x[1],fo))
        for k in cv: cv[k].append(rms(pr[k]))
    s4.update({k+'_heldout':round(float(np.mean(v)),4) for k,v in cv.items()}); out_['section4_plummer_ML2']=s4
    json.dump(out_,open(out,'w'),indent=1); print(json.dumps(out_,indent=1))
if __name__=='__main__': main(*sys.argv[1:4])
