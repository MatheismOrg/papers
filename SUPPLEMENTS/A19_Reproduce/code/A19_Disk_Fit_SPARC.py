"""A19 supplement: disk fit and held-out scores from public SPARC files.
Inputs: SPARC_Lelli2016c.mrt (galaxy table) and Rotmod_LTG/*_rotmod.dat (Lelli, McGaugh and
Schombert 2016, AJ 152, 157). Usage: python3 A19_Disk_Fit_SPARC.py table1.mrt rotmod_dir out_prefix [folds.csv]
"""
import os, sys, json, csv, numpy as np
from scipy.optimize import least_squares, minimize_scalar, minimize
G=4.30091e-6; ACC=3.24078e-14; h=0.671; H0=100*h; FB=0.0490/0.3175

def table1(path):
    rows={}; lines=open(path).read().split('\n')
    last=max(i for i,l in enumerate(lines) if l.startswith('-----'))
    for l in lines[last+1:]:
        p=l.split()
        if len(p)<18: continue
        rows[p[0]]=dict(D=float(p[2]),Inc=float(p[5]),L36=float(p[7])*1e9,Reff=float(p[9]),MHI=float(p[13])*1e9,Q=int(p[17]))
    return rows

def load(t1, rdir):
    gal=[]; 
    for g,r in sorted(t1.items()):
        if not (r['Q']<3 and r['Inc']>=30): continue
        d=np.loadtxt(os.path.join(rdir,g+'_rotmod.dat'),comments='#',ndmin=2)
        R,V,eV,Vg,Vd,Vb=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]
        vb2=Vg*np.abs(Vg)+0.5*Vd*np.abs(Vd)+0.7*Vb*np.abs(Vb)
        k=(R>0)&(V>0)&(eV/np.where(V>0,V,1)<=0.1)&(vb2>0)
        if k.sum()==0: continue
        Ms=0.5*r['L36']; Mp=Ms+1.33*r['MHI']
        gal.append(dict(name=g,r=R[k],gobs=V[k]**2/R[k]*ACC,gbar=vb2[k]/R[k]*ACC,Ms=Ms,Mp=Mp,Re=r['Reff']))
    return gal

def g_nfw(r,Mref):
    rho=3*(H0*1e-3)**2/(8*np.pi*G); r200=(3*Mref/(4*np.pi*200*rho))**(1/3)
    c=10**0.905*(Mref/(1e12/h))**-0.101; x=np.minimum(r/(r200/c),c)
    m=lambda y: np.log(1+y)-y/(1+y)
    return G*(1-FB)*Mref*m(x)/m(c)/r**2*ACC

def halo(p, r, gbar, MX, Mp, Re):
    A,B,C=p; Mref=Mp/FB+10**A*(MX/1e10)**B*(Re/3.0)**C
    return gbar+g_nfw(r,Mref)

def rar(gbar,a0): return gbar/(1-np.exp(-np.sqrt(gbar/a0)))

def moster_mref(ms):
    M1,N,be,ga=10**11.590,0.0351,1.376,0.608
    f=lambda l:(np.log10(2*N*10**l/((10**l/M1)**-be+(10**l/M1)**ga))-np.log10(ms))**2
    return 10**minimize_scalar(f,bounds=(8,15),method='bounded').x

def stack(gs, mass):
    r=np.concatenate([g['r'] for g in gs]); gb=np.concatenate([g['gbar'] for g in gs]); go=np.concatenate([g['gobs'] for g in gs])
    MX=np.concatenate([np.full(len(g['r']),g['Ms' if mass=='star' else 'Mp']) for g in gs])
    Mp=np.concatenate([np.full(len(g['r']),g['Mp']) for g in gs]); Re=np.concatenate([np.full(len(g['r']),g['Re']) for g in gs])
    return r,gb,go,MX,Mp,Re

def fit_halo(gs, mass, p0=(11.7,0.8,-0.7)):
    r,gb,go,MX,Mp,Re=stack(gs,mass)
    w=np.concatenate([np.full(len(g['r']),1/np.sqrt(len(g['r']))) for g in gs])
    res=least_squares(lambda p: w*(np.log10(go)-np.log10(halo(p,r,gb,MX,Mp,Re))), p0, x_scale=[1,0.3,0.3])
    return res.x

def fit_a0(gs):
    f=lambda la: np.mean([np.mean((np.log10(g['gobs'])-np.log10(rar(g['gbar'],10**la)))**2) for g in gs])
    return 10**minimize_scalar(f,bounds=(-10.5,-9.3),method='bounded').x

def resid(model, gs, par):
    out=[]
    for g in gs:
        if model=='rar': p=rar(g['gbar'],par)
        elif model=='moster': p=g['gbar']+g_nfw(g['r'],max(moster_mref(g['Ms']),g['Mp']/FB))
        else: p=halo(par,g['r'],g['gbar'],g['Ms' if model=='star' else 'Mp'],g['Mp'],g['Re'])
        out.append(np.log10(g['gobs'])-np.log10(p))
    return out

if __name__=='__main__':
    t1=table1(sys.argv[1]); gal=load(t1,sys.argv[2])
    npts=sum(len(g['r']) for g in gal); print('galaxies',len(gal),'points',npts)
    full={'star':fit_halo(gal,'star'),'baryon':fit_halo(gal,'baryon'),'a0':fit_a0(gal)}
    print('fit on all disks: stellar',full['star'].round(3),' baryonic',full['baryon'].round(3),' a0 %.4e'%full['a0'])
    # score: each held-out galaxy weighs equally (root of the mean over galaxies of the mean squared residual)
    rmsg=lambda lst: float(np.sqrt(np.mean([np.mean(v**2) for v in lst])))
    cv={k:[] for k in ('rar','baryon','moster','star')}; spd={k:[] for k in cv}
    if len(sys.argv)>4:   # the original fold assignments (three complete assignments, five folds each)
        rows=list(csv.DictReader(open(sys.argv[4]))); pos={g['name']:i for i,g in enumerate(gal)}
        assignments=[[np.array([pos[r['galaxy']] for r in rows if r[c]==f]) for f in sorted(set(r[c] for r in rows))] for c in rows[0] if c.startswith('fold_')]
    else:
        assignments=[np.array_split(np.random.default_rng(seed).permutation(len(gal)),5) for seed in (1,2,3)]
    for folds in assignments:
        idx=np.concatenate(folds)
        res={k:[] for k in cv}
        for f in folds:
            te=[gal[i] for i in f]; tr=[gal[i] for i in idx if i not in set(f)]
            res['rar']+=resid('rar',te,fit_a0(tr)); res['baryon']+=resid('baryon',te,fit_halo(tr,'baryon'))
            res['star']+=resid('star',te,fit_halo(tr,'star')); res['moster']+=resid('moster',te,None)
        for k in cv:
            cv[k].append(rmsg(res[k])); v=np.concatenate(res[k]); spd[k].append(float(np.median(np.abs(10**(-v/2)-1)))*100)
    cvm={k:round(float(np.mean(v)),5) for k,v in cv.items()}; spm={k:round(float(np.mean(v)),2) for k,v in spd.items()}
    print('held-out disk RMSE, dex (5 folds x 3 assignments):',cvm); print('median speed error, percent:',spm)
    json.dump(dict(selection='Q<3, Inc>=30, R>0, Vobs>0, errV/Vobs<=0.1, gbar>0; Ups_disk 0.5, Ups_bulge 0.7',galaxies=len(gal),points=npts,
        frozen_on_all_disks=dict(stellar_ABC=[round(float(x),4) for x in full['star']],baryonic_ABC=[round(float(x),4) for x in full['baryon']],rar_a0=float('%.4e'%full['a0'])),
        heldout_rmse_dex=cvm,heldout_runs={k:[round(x,5) for x in v] for k,v in cv.items()},median_speed_error_percent=spm),open(sys.argv[3]+'_Results.json','w'),indent=1)
