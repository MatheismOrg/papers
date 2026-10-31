"""External ETG and stacked weak-lensing checks of the frozen SPARC RAR.

No model parameter is fitted to these samples. The input files are official public
tables/releases. Lensing conversion follows Brouwer2021's SIS approximation.
"""
from pathlib import Path
from zipfile import ZipFile
import json,csv,io,hashlib,re
import numpy as np
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
BASE=json.loads((ROOT.parent/'dm_fit/results/fit_results.json').read_text())
A0=BASE['models']['RAR']['coefficients'][0]
KPC=3.085677581491367e19
def rar(g):return g/(-np.expm1(-np.sqrt(g/A0)))
def metrics(g,p):
    g=np.array(g);p=np.array(p);ok=(g>0)&(p>0)
    e=np.log10(p[ok]/g[ok]);return {'n_total':len(g),'n_positive_for_fraction_and_log_metrics':int(ok.sum()),
      'median_abs_acceleration_error':float(np.median(abs(p[ok]/g[ok]-1))),
      'rms_log10_acceleration_dex':float(np.sqrt(np.mean(e*e))),
      'median_signed_log10_predicted_over_observed':float(np.median(e)),
      'median_abs_circular_speed_equivalent_error':float(np.median(abs(np.sqrt(p[ok]/g[ok])-1)))}

meta={}
for line in (DATA/'ETG_Lelli2017.mrt').read_text().splitlines():
    if not line or line.startswith('#'):continue
    t=line.split();assert len(t)==20
    meta[t[0]]=np.array(t[1:],float)
etg=[]
with ZipFile(DATA/'Rotmod_ETG.zip') as arc:
    for name in sorted(n for n in arc.namelist() if n.endswith('_rotmod.dat')):
        system=name.replace('_rotmod.dat','');a=np.loadtxt(io.BytesIO(arc.read(name)));assert a.shape==(2,6)
        for j,(r,v,ev,gas,disk,bulge) in enumerate(a):
            go=v*v/r*1e6/KPC;gb=(gas*abs(gas)+.8*disk*abs(disk)+.8*bulge*abs(bulge))/r*1e6/KPC
            m=meta[system];tgo=m[11+2*j];tgb=m[15+2*j]
            assert abs(np.log10(go)-tgo)<.03,(system,j,np.log10(go),tgo)
            assert abs(np.log10(gb)-tgb)<.03,(system,j,np.log10(gb),tgb)
            # Use the published acceleration pairs; radii in the two-point
            # rotation files are rounded and can shift reconstructed logs ~0.01dex.
            gb=10.**tgb;go=10.**tgo
            etg.append({'population':'ETG','system':system,'tracer':'inner_JAM' if j==0 else 'outer_HI',
              'radius_kpc':r,'gbar_m_s2':gb,'gobs_m_s2':go,'gbar_error_dex':m[16+2*j],
              'gobs_error_dex':m[12+2*j],'vobs_km_s':v,'vobs_error_km_s':ev,
              'gpred_frozen_RAR_m_s2':float(rar(gb)),
              'log10_baryonic_potential_proxy_gbar_r_over_a0_kpc':float(np.log10(gb*r/A0)),
              'log10_baryonic_density_proxy_gbar_over_r_over_a0_per_kpc':float(np.log10(gb/(r*A0))),
              'interpretation':'inner JAM model velocity, not direct circular motion; outer HI velocity; fixed stellar M/L=0.8 at3.6micron',
              'source':'Lelli2017 official ETG tables'})
assert len(etg)==32 and len(set(r['system'] for r in etg))==16
with (OUT/'etg_external_pairs.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(etg[0]));w.writeheader();w.writerows(etg)
summary={'status':'External sample transfer with frozen SPARC RAR a0; no refit. Published model-derived observables retain source assumptions.',
 'frozen_a0_m_s2':A0,'source_fit':'149SPARCgalaxies,2700radii,unweightedlogaccelerationfit',
 'etg':{},'lensing':{},'group_dispersions':{}}
for tracer in ['all','inner_JAM','outer_HI']:
    r=etg if tracer=='all' else [q for q in etg if q['tracer']==tracer]
    d=metrics([q['gobs_m_s2'] for q in r],[q['gpred_frozen_RAR_m_s2'] for q in r])
    d['radius_kpc_range']=[min(q['radius_kpc'] for q in r),max(q['radius_kpc'] for q in r)]
    d['gbar_m_s2_range']=[min(q['gbar_m_s2'] for q in r),max(q['gbar_m_s2'] for q in r)]
    e=np.array([np.log10(q['gobs_m_s2']/q['gpred_frozen_RAR_m_s2']) for q in r])
    d['descriptive_spearman_residual_correlations']={}
    for key in ['radius_kpc','gbar_m_s2','log10_baryonic_potential_proxy_gbar_r_over_a0_kpc','log10_baryonic_density_proxy_gbar_over_r_over_a0_per_kpc']:
        d['descriptive_spearman_residual_correlations'][key]=float(spearmanr([q[key] for q in r],e).statistic)
    summary['etg'][tracer]=d

# Official README Eq7 convention, with h70=1. Constants intentionally reproduce release.
conv=4*4.52e-30*3.086e16
profiles={'KiDS_isolated':'Fig-4-5-C1_RAR-KiDS-isolated','GAMA_isolated':'Fig-4-C1_RAR-GAMA-isolated','KiDS_hotgas_sensitivity':'Fig-4_RAR-KiDS-isolated_hotgas'}
lens=[]
for sample,stem in profiles.items():
    p=DATA/'kids'/f'{stem}_Nobins.txt';a=np.loadtxt(p);rawcov=np.loadtxt(DATA/'kids'/f'{stem}_covmatrix.txt')
    n=len(a);assert rawcov.shape[0]==n*n
    gb=a[:,0];go=conv*a[:,1]/a[:,4];sigma=conv*a[:,3]/a[:,4]
    cov=np.zeros((n,n))
    for row in rawcov:
        i=np.argmin(abs(np.log(gb/row[2])));j=np.argmin(abs(np.log(gb/row[3])))
        assert abs(np.log(gb[i]/row[2]))<1e-4 and abs(np.log(gb[j]/row[3]))<1e-4
        cov[i,j]=conv*conv*row[4]/row[6]
    assert np.allclose(cov,cov.T,rtol=1e-9,atol=1e-36)
    assert np.linalg.eigvalsh(cov).min()>0
    assert np.allclose(np.sqrt(np.diag(cov)),sigma,rtol=1e-4)
    np.savetxt(OUT/f'{sample}_gobs_covariance_m2_s4.csv',cov,delimiter=',')
    for i in range(n):
        lens.append({'population':'weak_lensing_stack','system':sample,'stack_bin':i,
           'radius_kpc':'','gbar_m_s2':gb[i],'gobs_m_s2':go[i],'gobs_error_m_s2':sigma[i],
           'gbar_error_dex':'','gpred_frozen_RAR_m_s2':rar(gb[i]),
           'isolation_reliable_by_brouwer2021':bool(sample=='GAMA_isolated' or gb[i]>=1e-13),
           'interpretation':'stacked bin; SIS deprojection; covariance supplied; no uniquephysicalradius; KiDS/GAMA are overlapping analyses notindependentpopulations',
           'source':p.name})
    sets={}
    for selection,mask in [('all_bins',np.ones(n,bool)),('gbar_ge_1e-13',gb>=1e-13)]:
        d=metrics(go[mask],rar(gb[mask]));c=cov[np.ix_(mask,mask)];res=go[mask]-rar(gb[mask])
        d['chi2_statistical_covariance_only']=float(res@np.linalg.solve(c,res));d['dof_no_refit']=int(mask.sum())
        c_sis=c+np.diag((np.log(10)*.1*abs(go[mask]))**2)
        d['chi2_with_independent_diagonal_0p1dex_SIS_uncertainty']=float(res@np.linalg.solve(c_sis,res))
        d['mass_scale_sensitivity_stat_chi2']={}
        for shift in [-.2,0,.2]:
            rr=go[mask]-rar(gb[mask]*10**shift)
            d['mass_scale_sensitivity_stat_chi2'][str(shift)]=float(rr@np.linalg.solve(c,rr))
        d['gbar_m_s2_range']=[float(gb[mask].min()),float(gb[mask].max())]
        sets[selection]=d
    summary['lensing'][sample]=sets
with (OUT/'lensing_external_pairs.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(lens[0]));w.writeheader();w.writerows(lens)

# Group table is a genuinely different observable. Do not manufacture local RAR
# pairs from global dispersions and an arbitrary virial coefficient.
text=(DATA/'Milgrom2019.txt').read_text();part=text[text.index('TABLE I: The 67 groups'):]
group=[];G=6.67430e-11;MS=1.988409870698051e30
for line in part.splitlines():
    t=line.split()
    if len(t)!=13 or not re.match(r'^(NGC|ESO)',t[0]):continue
    vals=np.array(t[1:],float);name=t[0].replace('–','-');excluded='*' in name
    nv,sig,r,loglk=vals[:4];mass=10**loglk
    sigpred=(4/81*G*(mass*MS)*A0)**.25/1000
    # Reconstruct the source table's published prediction at a0=1.2e-10.
    sigsource=(4/81*G*(mass*MS)*1.2e-10)**.25/1000
    assert abs(sigsource-vals[6])<1.6,(name,sigsource,vals[6])
    group.append({'population':'galaxy_group_dispersion','system':name.replace('*',''),'source_excluded':excluded,
      'n_member_velocities':int(nv),'harmonic_radius_kpc':r,'log10_LK_solar':loglk,
      'assumed_baryonic_mass_msun':mass,'sigma_observed_km_s':sig,'sigma_predicted_frozen_a0_km_s':sigpred,
      'sigma_source_table_at_a0_1p2e10_km_s':vals[6],
      'interpretation':'global deepMONDvirial estimate; M/L_K=1; no gas inventory; equilibrium/isotropy assumptions; notlocalRARacceleration'})
assert len(group)==66 and sum(not q['source_excluded'] for q in group)==55
summary['group_table_count_note']='The paper text/caption say67total/56selected, but the supplied v2PDFTableI contains66tabulatedrows/55unstarredrows; we useonlythoseactualrows and do not invent the missingrow. Tablepageswerevisuallychecked.'
with (OUT/'group_velocity_dispersion_test.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(group[0]));w.writeheader();w.writerows(group)
for selection in ['all66_tabulated','source_selected55_tabulated']:
    rr=group if selection=='all66_tabulated' else [q for q in group if not q['source_excluded']]
    obs=np.array([q['sigma_observed_km_s'] for q in rr]);pred=np.array([q['sigma_predicted_frozen_a0_km_s'] for q in rr])
    summary['group_dispersions'][selection]={'n':len(rr),'median_abs_dispersion_error':float(np.median(abs(pred/obs-1))),
      'rms_log10_dispersion_dex':float(np.sqrt(np.mean(np.log10(pred/obs)**2))),
      'median_predicted_over_observed':float(np.median(pred/obs)),
      'harmonic_radius_kpc_range':[min(q['harmonic_radius_kpc'] for q in rr),max(q['harmonic_radius_kpc'] for q in rr)]}
summary['source_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.is_file() and p.suffix in ['.pdf','.mrt','.zip','.tar']}
summary['limitations']=['The samples are external to fitting, but all RAR choices were informed by existing literature; not blindconfirmatorytesting.',
 'ETGinnerJAM gobs is modeldependent and totalmass follows light in sourceinference; outerHI providesdifferenttracer.',
 'ETG stellar M/L=0.8 differsfromSPARC disk0.5/bulge0.7 bysourcepopulation prescription.',
 'Lensing gobs=4GDeltaSigma is the published singular-isothermal-sphere approximation, not directmeasurementof3Dacceleration.',
 'Lensing statisticalcovariance omitsfullbaryonicmass, gas, isolation anddeprojection systematics. Diagonal0.1dex is a sensitivitymodel, not a validatedfullcovariance.',
 'KiDS isolation deteriorates belowgbar~1e-13; GAMA spectroscopicselectionextendslower. KiDSandGAMAprofilesoverlap and arenotindependenttests.',
 'Group sourceexcludes10richgroups/clusters andonebad-distanceobject; stellar-onlybaryonicproxy undercountsgas especiallyrichsystems. Sourceindividualsigmaerrorsnotprovided.',
 'Globalpressure-supportedvirial relation is a separateMONDprediction, not proofthatlocaldiskRARworksunchanged.',
 'Noadditionalfittohistory,potential,ordensityproxiesperformed; proxycorrelationsaredescriptiveandsharegbarwithresidual.']
(OUT/'external_validation_results.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
