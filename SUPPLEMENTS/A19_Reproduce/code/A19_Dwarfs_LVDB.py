"""A19 supplement: dwarf spheroidal tests on the Local Volume Database.

Inputs: dwarf_all.csv of the Local Volume Database (Pace 2025, Open Journal of Astrophysics 8,
doi:10.33232/001c.144859), github.com/apace7/local_volume_database commit
c9fb49c620adab3df32c7d206e8599f4bb624a3a, sha256
c5c029dce9970e7b70fc329cc46dbf2ed02b308901579a3f2aa756887d542b26; and the results file of
A19_Disk_Fit_SPARC.py, which supplies the rules frozen on the 149 SPARC disks.
Usage: python3 A19_Dwarfs_LVDB.py dwarf_all.csv A19_Disk_Fit_SPARC_Results.json out_prefix
Nothing here is fitted to the dwarfs except where a row says "refitted on dwarfs", and those
rows are scored on dwarfs held out of the fit.
"""
import csv, json, sys, numpy as np
from scipy.optimize import minimize, minimize_scalar
G=4.30091e-6; ACC=3.24078e-14; h=0.671; H0=100*h; FB=0.0490/0.3175
KEYS=['fornax_1','sculptor_1','leo_1','leo_2','sextans_1','carina_1','draco_1','ursa_minor_1',
 'canes_venatici_1','crater_2','antlia_2','andromeda_01','andromeda_02','andromeda_03','andromeda_05',
 'andromeda_06','andromeda_07','andromeda_09','andromeda_14','andromeda_15','andromeda_16','andromeda_18',
 'andromeda_21','andromeda_23','andromeda_25','bootes_1','coma_berenices_1','canes_venatici_2',
 'hercules_1','leo_4','leo_5','ursa_major_1','ursa_major_2','segue_1','segue_2','willman_1',
 'reticulum_2','eridanus_2','horologium_1','hydrus_1','andromeda_10','andromeda_11','andromeda_12',
 'andromeda_13']
NLG=(10.94,0.059); MOG=(8.89,0.042)   # Rahvar & Mashhoon 2014; Moffat & Rahvar 2013 (alpha, mu per kpc)

def g_nfw(r,Mref):
    rho=3*(H0*1e-3)**2/(8*np.pi*G); r200=(3*Mref/(4*np.pi*200*rho))**(1/3)
    c=10**0.905*(Mref/(1e12/h))**-0.101; x=np.minimum(r/(r200/c),c); m=lambda y: np.log(1+y)-y/(1+y)
    return G*(1-FB)*Mref*m(x)/m(c)/r**2*ACC
def halo(p,r,gbar,MX,Mp,Re):
    A,B,C=p; return gbar+g_nfw(r,Mp/FB+10**A*(MX/1e10)**B*(Re/3.0)**C)
def rar(gbar,a0): return gbar/(1-np.exp(-np.sqrt(gbar/a0)))
def moster_mref(ms):
    M1,N,be,ga=10**11.590,0.0351,1.376,0.608
    f=lambda l:(np.log10(2*N*10**l/((10**l/M1)**-be+(10**l/M1)**ga))-np.log10(ms))**2
    return 10**minimize_scalar(f,bounds=(8,15),method='bounded').x
rms=lambda v: float(np.sqrt(np.mean(np.square(v))))

# --- kernel theories: enclosed effective mass of a Plummer sphere (scale a) inside radius r ---
XG,WG=np.polynomial.legendre.leggauss(48)
MF=(np.arange(4000)+0.5)/4000                    # equal-mass shells of the Plummer sphere
def Kcum(D,mu,kind):                            # fraction of the kernel within distance D
    x=mu*D; return 1-np.exp(-x)*(1+x/2) if kind=='NLG' else 1-np.exp(-x)*(1+x)
def kdens(d,mu,kind):
    x=mu*d; return mu*(1+x)*np.exp(-x)/2 if kind=='NLG' else mu*mu*d*np.exp(-x)
def eff_frac(r,a,mu,kind):
    """M_eff(<r) / (alpha M_b) for a Plummer sphere of scale a."""
    s=a/np.sqrt(MF**(-2/3)-1)
    inner=np.where(s<r,Kcum(np.abs(r-s),mu,kind),0.0)
    lo=np.abs(r-s); hi=r+s; d=(hi-lo)[:,None]/2*XG[None,:]+(hi+lo)[:,None]/2
    P=(r*r-(s[:,None]-d)**2)/(4*s[:,None]*d)
    mid=((hi-lo)/2)*np.sum(WG[None,:]*kdens(d,mu,kind)*P,axis=1)
    return float(np.mean(inner+mid))

def load(path):
    db={r['key']:r for r in csv.DictReader(open(path))}
    rows=[db[k] for k in KEYS if db[k]['vlos_sigma'] and db[k]['M_V'] and db[k]['rhalf_physical']]
    return rows,[k for k in KEYS if db[k] not in rows]

def main(lvdb,diskjson,out):
    D=json.load(open(diskjson))['frozen_on_all_disks']; STAR,BAR,A0D=D['stellar_ABC'],D['baryonic_ABC'],D['rar_a0']
    rows,dropped=load(lvdb); n=len(rows)
    L=10**(-0.4*(np.array([float(x['M_V']) for x in rows])-4.83))
    ReMaj=np.array([float(x['rhalf_physical']) for x in rows])/1000; ReSph=np.array([float(x['rhalf_sph_physical']) for x in rows])/1000
    s=np.array([float(x['vlos_sigma']) for x in rows]); host=np.array([x['host'] for x in rows])
    def models(ml=2.0,Re=ReMaj):
        Ms=ml*L; r=4/3*Re; gbar=G*(Ms/2)/r**2*ACC; lo=np.log10(3*s**2/r*ACC)
        mo=np.array([gb+g_nfw(rr,max(moster_mref(m),m/FB)) for gb,rr,m in zip(gbar,r,Ms)])
        return lo,{'stellar_rule':lo-np.log10(halo(STAR,r,gbar,Ms,Ms,Re)),'baryonic_rule':lo-np.log10(halo(BAR,r,gbar,Ms,Ms,Re)),
                   'rar_disk_a0':lo-np.log10(rar(gbar,A0D)),'rar_published_a0':lo-np.log10(rar(gbar,1.20e-10)),'moster':lo-np.log10(mo)}
    lo,R=models(); res={'n':n,'dropped_upper_limit_dispersion':dropped,'frozen_from_disks':D,'rmse_dex':{k:round(rms(v),4) for k,v in R.items()}}
    def boot(a,b,seed=20261030):
        rng=np.random.default_rng(seed); d=[]
        for _ in range(20000):
            i=rng.integers(0,n,n); d.append(rms(a[i])-rms(b[i]))
        return [round(rms(a)-rms(b),4)]+np.percentile(d,[2.5,97.5]).round(4).tolist()
    st=R['stellar_rule']
    res['paired_bootstrap']={'stellar_minus_baryonic':boot(st,R['baryonic_rule']),'stellar_minus_rar_disk':boot(st,R['rar_disk_a0'])}
    res['closer_under_stellar']={'vs_baryonic':int(np.sum(np.abs(st)<np.abs(R['baryonic_rule']))),'vs_rar_disk':int(np.sum(np.abs(st)<np.abs(R['rar_disk_a0'])))}
    sub={'LV_ge_1e5':L>=1e5,'LV_lt_1e5':L<1e5,'host_milky_way_incl_lmc':np.isin(host,['mw','lmc']),'host_andromeda':host=='m_031'}
    res['subsets']={k:dict(n=int(m.sum()),**{kk:round(rms(v[m]),4) for kk,v in R.items()}) for k,m in sub.items()}
    res['sensitivity']={lab:{k:round(rms(v),4) for k,v in models(**kw)[1].items()} for lab,kw in [('ML_1',dict(ml=1.0)),('ML_3',dict(ml=3.0)),('circularized_Re',dict(Re=ReSph))]}
    loo=[rms(np.delete(st,k))/rms(np.delete(R['rar_disk_a0'],k)) for k in range(n)]; res['leave_one_out_ratio_stellar_over_rar']=[round(min(loo),3),round(max(loo),3)]
    # --- Section 4: kernel ceilings and full integration on Plummer spheres ---
    res['ceilings']={}
    for ml in (2.0,3.0):
        Ms=ml*L; r=4/3*ReMaj; f=r**3/(r*r+ReMaj**2)**1.5; ratio=(3*s**2/r)/(G*f*Ms/r**2)
        res['ceilings']['ML_%g'%ml]={'plummer_fraction_at_4/3Re':round(float(f[0]),4),
            'exceed_NLG_1+alpha/f':int(np.sum(ratio>1+NLG[0]/f)),'exceed_MOG_1+alpha/f':int(np.sum(ratio>1+MOG[0]/f))}
    Ms=2.0*L; r=4/3*ReMaj; fenc=r**3/(r*r+ReMaj**2)**1.5; Mb=Ms*fenc; gN=G*Mb/r**2*ACC
    MU=np.logspace(-3,3,121); F={k:np.array([[eff_frac(r[i],ReMaj[i],mu,k) for mu in MU] for i in range(n)]) for k in ('NLG','MOG')}
    # g = G [M_b(<r) + alpha M_b F] / r^2 = gN (1 + alpha F / fenc), with F = M_eff(<r) / (alpha M_b) and fenc = M_b(<r) / M_b.
    # Turn 13 correction (found by an independent review): the turn 10 version omitted "/ fenc" here and in the multiplier below.
    def kern_pred(kind,al,lmu,idx): return gN[idx]*(1+al*np.array([np.interp(lmu,np.log10(MU),F[kind][i]) for i in idx])/fenc[idx])
    allidx=np.arange(n); K={}
    K['rar_published_a0']=lo-np.log10(rar(gN,1.20e-10))
    for kind,(al,mu) in (('NLG',NLG),('MOG',MOG)): K[kind+'_published']=lo-np.log10(kern_pred(kind,al,np.log10(mu),allidx))
    cvs={'rar_refit':[], 'NLG_refit':[], 'MOG_refit':[]}; bounds_hit=[]; mult=[]; nfit=0
    for seed in (1,2,3):
        idx=np.random.default_rng(seed).permutation(n); folds=np.array_split(idx,5); pr={k:np.zeros(n) for k in cvs}
        for fo in folds:
            tr=np.array([i for i in idx if i not in set(fo)])
            la=minimize_scalar(lambda x: np.mean((lo[tr]-np.log10(rar(gN[tr],10**x)))**2),bounds=(-12,-7),method='bounded').x
            pr['rar_refit'][fo]=lo[fo]-np.log10(rar(gN[fo],10**la))
            for kind in ('NLG','MOG'):
                best=None
                for lmu in np.log10(MU)[::4]:
                    o=minimize(lambda p: np.mean((lo[tr]-np.log10(kern_pred(kind,np.exp(p[0]),p[1],tr)))**2),[np.log(10),lmu],method='L-BFGS-B',bounds=[(np.log(1e-3),np.log(1e4)),(-3,3)])
                    if best is None or o.fun<best.fun: best=o
                nfit+=1; bounds_hit.append(best.x[1]>2.99); mult.append(1+np.exp(best.x[0])*np.interp(best.x[1],np.log10(MU),np.mean(F[kind]/fenc[:,None],axis=0)))
                pr[kind+'_refit'][fo]=lo[fo]-np.log10(kern_pred(kind,np.exp(best.x[0]),best.x[1],fo))
        for k in cvs: cvs[k].append(rms(pr[k]))
    res['section4_plummer_ML2']={**{k:round(rms(v),4) for k,v in K.items()},**{k+'_heldout':round(float(np.mean(v)),4) for k,v in cvs.items()},
        'kernel_refits':nfit,'refits_with_mu_at_upper_bound_1000_per_kpc':int(np.sum(bounds_hit)),'refit_multiplier_median':round(float(np.median(mult)),1),
        'published_kernel_median_miss_factor':{k:round(float(10**np.median(K[k+'_published'])),1) for k in ('NLG','MOG')}}
    # dwarf-refitted RAR on the point estimates of Section 3 (held out)
    cv=[]
    for seed in (1,2,3):
        idx=np.random.default_rng(seed).permutation(n); pr=np.zeros(n); gb=G*(Ms/2)/r**2*ACC
        for fo in np.array_split(idx,5):
            tr=np.array([i for i in idx if i not in set(fo)])
            la=minimize_scalar(lambda x: np.mean((lo[tr]-np.log10(rar(gb[tr],10**x)))**2),bounds=(-12,-7),method='bounded').x
            pr[fo]=lo[fo]-np.log10(rar(gb[fo],10**la))
        cv.append(rms(pr))
    res['rmse_dex']['rar_refit_on_dwarfs_heldout']=round(float(np.mean(cv)),4)
    json.dump(res,open(out+'_Results.json','w'),indent=1)
    with open(out+'.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['lvdb_key','name','host','L_V_Lsun','R_e_kpc','sigma_kms','log_g_obs']+['res_'+k for k in R]+['ref_structure','ref_vlos','ref_m_v'])
        for i,x in enumerate(rows): w.writerow([x['key'],x['name'],x['host'],f'{L[i]:.4g}',f'{ReMaj[i]:.4f}',x['vlos_sigma'],f'{lo[i]:.4f}']+[f'{R[k][i]:.4f}' for k in R]+[x['ref_structure'],x['ref_vlos'],x['ref_m_v']])
    print(json.dumps(res,indent=1))
if __name__=='__main__': main(*sys.argv[1:4])
