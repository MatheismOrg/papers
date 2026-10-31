"""Frozen external dSph benchmark from McConnachie's January 2021 catalog.

No dwarf is used to train any formula here. Requires numpy and astropy.
Projected major-axis radii are circularized before the Wolf mass estimator.
"""
from pathlib import Path
import csv, json, re, hashlib, warnings
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
G=6.67430e-11; MSUN=1.988409870698051e30; KPC=3.085677581491367e19
MSUN_V=4.83

def roman(s):
    vals={'I':1,'V':5,'X':10};out=0
    for i,c in enumerate(s):out+=-vals[c] if i+1<len(s) and vals[c]<vals[s[i+1]] else vals[c]
    return str(out)

def norm(s):
    s=s.replace('(','').replace(')','').replace('*','').strip()
    s=re.sub(r'([IVX]+)$',lambda m:roman(m[1]),s)
    return re.sub(r'[^a-z0-9]','',s.lower())

# Names and classifications are read from the author's actual updated Table 1.
morph={}
for line in (DATA/'morphology.txt').read_text().splitlines():
    if not re.search(r'\bdSph\??\b',line):continue
    parts=re.split(r'\s{2,}',line.strip())
    if len(parts)<3:continue
    name=parts[0]
    mt=next((s for s in parts[1:] if re.fullmatch(r'(?:dE/|dIrr/)?dSph\??',s)),None)
    if mt is None:continue
    host=next((s for s in parts[1:] if s in ['G','A','L','N','A/L','G/L']), 'unknown')
    morph[norm(name)]={'table1_name':name,'morphology':mt,'host':host,'source_line':line.strip()}

def finite(v):return np.isfinite(v) and abs(v)<900
def num(v):return float(v) if finite(v) else None
def errmean(r,k):
    v=[float(r[k+c]) for c in ['+','-']]
    return np.mean(v) if all(finite(x) and x>=0 for x in v) else None

selected=[]; excluded=[]; candidates=[]
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    with fits.open(DATA/'NearbyGalaxies_Jan2021_PUBLIC.fits') as hdus:
        (OUT/'fits_field_descriptions.xml').write_bytes(bytes(hdus[0].data))
        for row in hdus[1].data:
            name=row['GalaxyName'].strip();meta=morph.get(norm(name)); reasons=[]
            if meta is None or 'dSph' not in meta['morphology']:continue
            # Restrict to unambiguous dSph classification; archive all exclusions.
            if meta['morphology']!='dSph':reasons.append('not_unambiguous_pure_dSph')
            if name.startswith('*'):reasons.append('catalog_starred_candidate')
            if name in ['SagittariusdSph','AndromedaXXVII']:
                reasons.append('table1_explicit_tidal_disruption')
            if name=='Willman1':reasons.append('table1_cluster_classification_uncertain')
            if name=='Segue2':reasons.append('catalog_comments_newer_sigma_upper_limit')
            for k in ['dmod','Vmag','rh','e=1-b/a','sigma_s']:
                if not finite(row[k]):reasons.append('missing_'+k)
            if not (row['rh']>0 and row['sigma_s']>0 and 0<=row['e=1-b/a']<1):
                reasons.append('invalid_radius_dispersion_or_ellipticity')
            if not all(finite(row['sigma_s'+c]) and row['sigma_s'+c]>0 for c in ['+','-']):
                reasons.append('missing_dispersion_error_or_possible_limit')
            if row['sigma_s']-row['sigma_s-']<=0:
                reasons.append('dispersion_not_resolved_above_zero_at_quoted_lower_error')
            if finite(row['vrot_s']) and row['vrot_s']>0.5*row['sigma_s']:
                reasons.append('reported_rotation_over_half_dispersion')
            if reasons:
                excluded.append({'system':name,'reasons':';'.join(reasons),**meta});continue
            dm=float(row['dmod']); d=10**((dm-10)/5)
            # McConnachie2012 section4/Table3: apparent V is foreground-corrected.
            mv=float(row['Vmag'])-dm; lum=10**((MSUN_V-mv)/2.5)
            major=float(row['rh'])*np.pi/(180*60)*d
            Re=major*np.sqrt(1-float(row['e=1-b/a']))
            rhalf=4*Re/3
            sig=float(row['sigma_s']); go=3*sig**2*1e6/(rhalf*KPC)
            star=2*lum; gb=G*(star/2)*MSUN/(rhalf*KPC)**2
            hi=num(row['MHI']); hiratio=None if hi is None else hi*1e6/star
            if hiratio is not None and hiratio>0.1:
                excluded.append({'system':name,'reasons':'catalog_HI_mass_over_0.1_nominal_stellar_mass',**meta});continue
            sd=errmean(row,'dmod'); sr=errmean(row,'rh'); se=errmean(row,'e'); sv=errmean(row,'Vmag'); ss=errmean(row,'sigma_s')
            elog_r=None if None in (sd,sr,se) else np.sqrt((np.log(10)/5*sd)**2+(sr/row['rh'])**2+(se/(2*(1-row['e=1-b/a'])))**2)/np.log(10)
            elog_go=None if elog_r is None or ss is None else np.sqrt((2*ss/(sig*np.log(10)))**2+elog_r**2)
            # Distance cancels from inferred stellar acceleration; do not count it twice.
            elog_gb=None if None in (sv,sr,se) else np.sqrt((0.4*sv)**2+(2*sr/(row['rh']*np.log(10)))**2+(se/((1-row['e=1-b/a'])*np.log(10)))**2)
            out={'population':'ultrafaint_dSph' if lum<1e5 else 'classical_dSph','system':name,'host':meta['host'],
                 'radius_kpc':rhalf,'gbar_m_s2':gb,'gobs_m_s2':go,'projected_Re_circularized_kpc':Re,
                 'projected_major_axis_rh_arcmin':float(row['rh']),'ellipticity':float(row['e=1-b/a']),
                 'distance_kpc':d,'distance_modulus':dm,'Vmag_foreground_corrected':float(row['Vmag']),
                 'absolute_MV':mv,'luminosity_V_Lsun':lum,'stellar_ML_V_assumed':2.0,'stellar_mass_Msun':star,
                 'sigma_los_km_s':sig,'sigma_err_plus_km_s':float(row['sigma_s+']),'sigma_err_minus_km_s':float(row['sigma_s-']),
                 'rh_err_plus_arcmin':num(row['rh+']),'rh_err_minus_arcmin':num(row['rh-']),
                 'distance_modulus_err_plus':num(row['dmod+']),'distance_modulus_err_minus':num(row['dmod-']),
                 'Vmag_err_plus':num(row['Vmag+']),'Vmag_err_minus':num(row['Vmag-']),
                 'ellipticity_err_plus':num(row['e+']),'ellipticity_err_minus':num(row['e-']),
                 'err_log10_gobs_dex_stat_approx':elog_go,'err_log10_gbar_dex_stat_approx':elog_gb,
                 'gbar_ML1_m_s2':gb/2,'gbar_ML3_m_s2':1.5*gb,'Mhalf_Wolf_Msun':3*sig**2*1e6*rhalf*KPC/(G*MSUN),
                 'catalog_HI_1e6Msun':hi,'HI_over_assumed_stellar_mass':hiratio,
                 'catalog_stellar_rotation_km_s':num(row['vrot_s']),
                 'possible_tidal_or_classification_caution':bool(re.search(r'tidal|remnant|cluster\?',meta['source_line'],re.I)),
                 'morphology':meta['morphology'],'catalog_reference_numbers':row['References'].strip()}
            selected.append(out)

assert len(selected)>=20
def writecsv(name, rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
writecsv('dwarfs_benchmark.csv',selected);writecsv('excluded_dwarfs.csv',excluded)
(OUT/'morphology_lookup.json').write_text(json.dumps(morph,indent=2))

def rar(g,a0):return g/(-np.expm1(-np.sqrt(g/a0)))
def metrics(pred,obs):
    e=np.log10(pred/obs)
    return {'n':len(e),'rmse_log10_acceleration_dex':float(np.sqrt(np.mean(e*e))),
            'median_abs_acceleration_fraction':float(np.median(abs(pred/obs-1))),
            'median_abs_dispersion_fraction_under_same_estimator':float(np.median(abs(np.sqrt(pred/obs)-1))),
            'median_observed_over_predicted':float(np.median(obs/pred))}
fit=json.loads((ROOT.parent/'dm_fit/results/fit_results.json').read_text())
a0=fit['models']['RAR']['coefficients'][0]
summary={}
for pop in ['all','classical_dSph','ultrafaint_dSph','higher_precision']:
    subset=[r for r in selected if pop=='all' or r['population']==pop or
            (pop=='higher_precision' and max(r['sigma_err_plus_km_s'],r['sigma_err_minus_km_s'])/r['sigma_los_km_s']<=0.25
             and not r['possible_tidal_or_classification_caution'])]
    gb=np.array([r['gbar_m_s2'] for r in subset]);go=np.array([r['gobs_m_s2'] for r in subset]);r=np.array([r['radius_kpc'] for r in subset])
    summary[pop]={'radius_range_kpc':[float(r.min()),float(r.max())],
                  'gbar_range_m_s2':[float(gb.min()),float(gb.max())],
                  'gobs_range_m_s2':[float(go.min()),float(go.max())],
                  'RAR_frozen_SPARC':metrics(rar(gb,a0),go),
                  'RAR_literature_a0_1p2e_10':metrics(rar(gb,1.2e-10),go),
                  'ML1_RAR_frozen':metrics(rar(gb/2,a0),go),
                  'ML3_RAR_frozen':metrics(rar(gb*1.5,a0),go)}
for d in selected:
    d['RAR_frozen_SPARC_pred_m_s2']=float(rar(d['gbar_m_s2'],a0))
    d['RAR_observed_over_predicted']=d['gobs_m_s2']/d['RAR_frozen_SPARC_pred_m_s2']
writecsv('dwarfs_predictions.csv',selected)
metadata={
 'status':'External population transfer of frozen galaxy RAR; no dwarf fitting.',
 'catalog':'McConnachie2012 author maintained January2021 FITS update,144catalog rows.',
 'catalog_url':'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/PANDAS/NearbyGals/NearbyGalaxies_Jan2021_PUBLIC.fits',
 'catalog_page':'https://www.cadc-ccda.hia.nrc-cnrc.gc.ca/en/community/nearby/',
 'morphology_url':'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/PANDAS/NearbyGals/table1_OCT2019.pdf',
 'catalog_reference':'McConnachie2012 AJ144,4; DOI10.1088/0004-6256/144/1/4; arXiv1204.1562',
 'estimator_reference':'Wolf etal2010 MNRAS406,1220; DOI10.1111/j.1365-2966.2010.16753.x; arXiv0908.2995',
 'selection':'Pure dSph morphology, unstarred catalog name; measured positive radius/sigma/ellipticity, positive finite sigma errors and positive lower1sigma endpoint; excludes reported rotation>0.5sigma, explicit Sagittarius/AndXXVII tidal disruption, ambiguousWillman1, upperlimitSegue2; catalog HI/nominalMstar<=0.1 when available. Missing HI permitted for dSph morphology.',
 'formula':'Re=D*rh_arcmin*pi/(180*60)*sqrt(1-ellipticity); rhalf=4Re/3; L_V=10^((4.83-(Vmag-dmod))/2.5); Mstar=2L_V; gbar=G*Mstar/(2*rhalf^2); gobs=3*sigma_los^2/rhalf.',
 'photometry':'Vmag foreground extinction corrected per McConnachie2012 section4; no additional reddening correction.',
 'population_definition':'ultrafaint_dSph iff L_V<100000Lsun; otherwise classical_dSph, an operational luminosity split independent of residuals.',
 'higher_precision_selection':'max(quoted sigma plus/minus error)/sigma <=0.25 and no possible tidal/classification caution flag.',
 'assumptions':'Spherical approximately equilibrium dispersion-supported tracer; circularized half-light radius; 4/3deprojection; approximately flat luminosity-weighted dispersion nearhalf-light; constantstellarM/L; gas ignored; missingrotationnotproofofzero. Catalogcompiledheterogeneousobservationsnotuniformreanalysis.',
 'uncertainty':'Input asymmetricerrors retained. Derived dexerrors arefirstorder propagation withmean+/−errors, independenceassumed; covarianceandWolfmodelerrornotincluded. Distancecancelsfromgbar. ML1and3sensitivitynotconfidencebounds. Binarycontamination,tides,equilibrium,externalfields,contaminants remainunmodeled.',
 'RAR_a0_from_SPARC_m_s2':a0,'summary':summary,
 'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.suffix in ['.fits','.pdf','.dat','.log']},
 'excluded_count':len(excluded),'selected_count':len(selected)}
(OUT/'dwarfs_summary.json').write_text(json.dumps(metadata,indent=2))
print(json.dumps({'selected':len(selected),'excluded':len(excluded),'summary':summary},indent=2))
