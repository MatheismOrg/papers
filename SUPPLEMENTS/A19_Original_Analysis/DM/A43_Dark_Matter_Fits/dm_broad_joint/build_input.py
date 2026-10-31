from pathlib import Path
import csv,json,re,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;TMP=ROOT.parent
files={
'SPARC':TMP/'dm_fit/results/held_out_predictions.csv',
'dSph':TMP/'dm_broad_dwarfs/results/dwarfs_benchmark.csv',
'ETG':TMP/'dm_broad_groups/results/etg_external_pairs.csv',
'MaNGA_BCG':TMP/'dm_broad_lensing/manga_bcg_standardized.csv',
'CLASH':TMP/'dm_broad_lensing/clash_standardized.csv'}

def canon(name):
    name=re.sub('[^A-Z0-9]','',name.upper())
    name=re.sub(r'^(NGC|IC|UGC)0+(\d)',r'\1\2',name)
    return name

out=[]
for p,path in files.items():
    for r in csv.DictReader(path.open()):
        name=r.get('system',r.get('galaxy'));gr=canon(name)
        if name in ['MaNGA_8989-12703','MaNGA_8991-6102']:
            gr='MANGA_NEARBY_BCG_PAIR_8989_8991'
        sub=r.get('population',p)
        if p=='dSph':env={'A':'M31','G':'Milky_Way','L':'Local_Group'}.get(r.get('host',''),r.get('host',''))
        else:env=r.get('environment','')
        row=dict(population=p,subpopulation=sub,system=name,group=gr,
          radius_kpc=float(r['radius_kpc']),gbar_m_s2=float(r['gbar_m_s2']),gobs_m_s2=float(r['gobs_m_s2']),
          environment=env,source=r.get('source',p),
          loggbar_err_dex=r.get('loggbar_err_dex',r.get('gbar_error_dex',r.get('err_log10_gbar_dex_stat_approx',''))),
          loggobs_err_dex=r.get('loggobs_err_dex',r.get('gobs_error_dex',r.get('err_log10_gobs_dex_stat_approx',''))),
          measurement_method=r.get('kinematic_method',r.get('tracer','Wolf_half_light_estimator' if p=='dSph' else 'rotation_curve')))
        assert all(np.isfinite(row[k]) and row[k]>0 for k in ['radius_kpc','gbar_m_s2','gobs_m_s2'])
        out.append(row)

with (ROOT/'combined_input.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)

summary={}
for p in files:
    rs=[r for r in out if r['population']==p]
    summary[p]={'rows':len(rs),'systems':len({r['system'] for r in rs}),'validation_groups':len({r['group'] for r in rs}),
      **{k:[min(r[k] for r in rs),max(r[k] for r in rs)] for k in ['radius_kpc','gbar_m_s2','gobs_m_s2']}}
seen={}
for r in out:seen.setdefault(r['group'],set()).add(r['population'])
metadata={'source_files':{p:{'path':str(f.relative_to(TMP)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for p,f in files.items()},
 'counts':summary,'total_rows':len(out),'total_physical_systems':len({r['system'] for r in out}),'total_validation_groups':len(seen),
 'cross_population_same_systems':{k:sorted(v) for k,v in seen.items() if len(v)>1},
 'population_definition':'CLASH inner BCG and outer cluster points share one training population and one group per cluster. Classical/ultrafaint dwarfs share dSph training population and are reported separately. MaNGA and CLASH redshift ranges do not overlap. Known names normalized across SPARC/ETG. Two nearby MaNGA BCGs (8989-12703 and8991-6102;2.694arcmin,delta_z.002943) conservatively share a validation group because a common hostenvironment is possible; host identities are not fully documented.',
 'reserved_inputs':'X-COP fixed-radius profiles remain a separate prior-inspected stress test; no X-COP R500 points. KiDS/GAMA stacked lensing have no unique physical radius; separate external one-input tests. Group dispersions are separate global virial observables. Globular central-dispersion data are controls, not pooled local RAR points.',
 'limitations':'Fixed-input heterogeneous estimators; no errors-in-variables or covariance inference. Physical system folds do not hold out shared MilkyWay/M31 environments, which receive a separate sensitivity check. Some shared calibration/survey uncertainties remain across folds.'}
(ROOT/'input_metadata.json').write_text(json.dumps(metadata,indent=2))
print(json.dumps(metadata,indent=2))
