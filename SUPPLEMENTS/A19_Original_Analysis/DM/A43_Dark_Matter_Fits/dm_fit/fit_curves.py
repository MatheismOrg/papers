"""Empirical SPARC fits; no galaxy-specific refitting in the validation set."""
from pathlib import Path
from zipfile import ZipFile
import io,json,csv,hashlib
import numpy as np
from scipy.optimize import minimize_scalar
from numpy.polynomial import Polynomial

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
meta={}
for line in (DATA/'SPARC_Lelli2016c.mrt').read_text().splitlines()[98:]:
    if not line.strip():continue
    parts=line.split()
    assert len(parts)==19
    meta[parts[0]]={'inc':float(parts[5]),'quality':int(parts[17])}
assert len(meta)==175
rows=[];total=0;filtered_galaxies=0
with ZipFile(DATA/'Rotmod_LTG.zip') as z:
    for name in sorted(z.namelist()):
        g=Path(name).name.replace('_rotmod.dat','');m=meta[g]
        a=np.loadtxt(io.BytesIO(z.read(name)),ndmin=2);total+=len(a)
        if m['quality']==3 or m['inc']<30:continue
        filtered_galaxies+=1
        for r,v,ev,gas,disk,bulge,*_ in a:
            vb2=gas*abs(gas)+0.5*disk*abs(disk)+0.7*bulge*abs(bulge)
            if r<=0 or v<=0 or ev/v>0.1 or vb2<=0:continue
            gb=vb2/r*1e6/3.085677581491367e19;go=v*v/r*1e6/3.085677581491367e19
            rows.append([g,r,v,ev,gb,go])
gal=np.array([r[0] for r in rows]);a=np.array([r[1:] for r in rows],float)
x=np.log10(a[:,3]/1e-10);y=np.log10(a[:,4]/1e-10)
groups=np.unique(gal);rng=np.random.default_rng(20261027)
shuffled=rng.permutation(groups);folds={str(g):i%5 for i,g in enumerate(shuffled)}
fold=np.array([folds[g] for g in gal])

def rar(x,loga):
    # a0 in units 1e-10 m/s^2
    return x-np.log10(-np.expm1(-np.sqrt(10.**(x-loga))))

def fit_predict(kind,xt,yt,xp,weights=None):
    if kind=='baryons_only':return xp.copy(),[]
    if kind.startswith('poly'):
        d=int(kind[4:]);p=Polynomial.fit(xt,yt,d,w=None if weights is None else np.sqrt(weights)).convert()
        return p(xp),p.coef.tolist()
    f=lambda la:np.average((rar(xt,la)-yt)**2,weights=weights)
    res=minimize_scalar(f,bounds=(-3,2),method='bounded',options={'xatol':1e-12})
    assert res.success
    return rar(xp,res.x),[float(10.**res.x*1e-10)]

kinds=['baryons_only','poly1','poly2','poly3','poly5','poly9','RAR']
summary={};preds={}
for kind in kinds:
    pred=np.zeros(len(x))
    for k in range(5):
        train=fold!=k;test=~train
        assert not set(gal[train])&set(gal[test])
        pred[test],_=fit_predict(kind,x[train],y[train],x[test])
    fitted,coef=fit_predict(kind,x,y,x)
    e=pred-y
    summary[kind]={'coefficients':coef,'train_rmse_dex':float(np.sqrt(np.mean((fitted-y)**2))),
        'cv_rmse_dex':float(np.sqrt(np.mean(e**2))),
        'cv_equal_galaxy_rmse_dex':float(np.sqrt(np.mean([np.mean(e[gal==g]**2) for g in groups]))),
        'cv_median_abs_velocity_fraction':float(np.median(abs(10**(e/2)-1))),
        'cv_multiplicative_acceleration_scatter':float(10**np.sqrt(np.mean(e**2))),
        'cv_fraction_predicted_below_baryonic':float(np.mean(pred<x))}
    preds[kind]=pred

# Equal-galaxy weighting checks whether densely sampled galaxies drive the trend.
weights=np.array([1/np.sum(gal==g) for g in gal])
for kind in ['poly3','RAR']:
    pred=np.zeros(len(x))
    for k in range(5):
        train=fold!=k;test=~train
        pred[test],_=fit_predict(kind,x[train],y[train],x[test],weights[train])
    summary[kind]['equal_galaxy_weight_fit_cv_rmse_dex']=float(np.sqrt(np.mean((pred-y)**2)))

res={'raw_points':total,'quality_inclination_selected_galaxies':filtered_galaxies,'galaxies':len(groups),'points':len(x),
    'selection':'Q < 3, inclination >= 30 deg, errV/Vobs <= 0.1, positive radius/observed speed/net inward baryonic acceleration',
    'mass_to_light':'disk 0.5, bulge 0.7 in solar units at 3.6 micron; published distances/inclinations fixed',
    'objective':'Unweighted least squares in log10 acceleration; fixed observed x values; no error-in-variables or nuisance marginalization',
    'cross_validation':'5 whole-galaxy folds, numpy default_rng seed 20261027, sorted galaxy IDs then permutation',
    'gbar_range_m_s2':[float(a[:,3].min()),float(a[:,3].max())], 'radius_range_kpc':[float(a[:,0].min()),float(a[:,0].max())],
    'models':summary,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.suffix in ['.zip','.mrt']}}
(OUT/'fit_results.json').write_text(json.dumps(res,indent=2))
with (OUT/'held_out_predictions.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['galaxy','radius_kpc','vobs_km_s','err_v_km_s','gbar_m_s2','gobs_m_s2','fold']+[f'{k}_pred_gobs_m_s2' for k in kinds])
    for i,row in enumerate(rows):w.writerow(row+[int(fold[i])]+[float(1e-10*10**preds[k][i]) for k in kinds])
print(json.dumps(res,indent=2))
