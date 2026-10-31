"""Parse published numerical tables, never digitize plotted points.

Inputs: cached VizieR TSV files, and pdftotext -layout output of arXiv2402.12016.
Retains published log10 accelerations and one-sigma errors, with no invented covariance.
"""
from pathlib import Path
import re,csv,json,hashlib,sys
import numpy as np
P=Path(__file__).resolve().parent

def read_tsv(path):
    lines=[l for l in path.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return [{k.strip():v.strip() for k,v in zip(lines[0].split('\t'), l.split('\t'))} for l in lines[3:]]

meta=read_tsv(P/'clash_table1_idia.tsv')
byname={r['AName']:r for r in meta}
clash=[]
for row in read_tsv(P/'clash_idia.tsv'):
    host=row['AName']; m=byname[host];r=float(row['Rad']);x=float(row['log(gbar)']);y=float(row['log(gtot)'])
    pop='CLASH_BCG' if r<100 else 'CLASH_cluster'
    clash.append(dict(population=pop,system='CLASH_'+host,parent_system='CLASH_'+host,radius_kpc=r,
        gbar_m_s2=10**x,gobs_m_s2=10**y,loggbar_err_dex=float(row['e_log(gbar)']),loggobs_err_dex=float(row['e_log(gtot)']),
        loggbar=x,loggobs=y,redshift=float(m['z']),source='Tian2020_CDS_fig2',kinematic_method='lensing_NFW_spherical_mass',
        environment='cluster_central_BCG' if r<100 else 'cluster',canonical_cluster_name=m['Name'],
        ra_j2000=m['RAJ2000'],dec_j2000=m['DEJ2000']))
assert len(clash)==84 and len(set(r['parent_system'] for r in clash))==20
assert sum(r['population']=='CLASH_BCG' for r in clash)==20

if (P/'manga_bcg_table1_extracted.txt').exists():
    t=(P/'manga_bcg_table1_extracted.txt').read_text()
else:
    t=(P/'bcg.txt').read_text().split('Table 1.     Properties and Results of 50 MaNGA BCGs')[1].split('Note. —')[0]
lines=[l for l in t.splitlines() if re.match(r'^\s*\d{4}-\d{4,5}',l)]
assert len(lines)==50
(P/'manga_bcg_table1_extracted.txt').write_text('\n'.join(lines)+'\n')
bcg=[]
for l in lines:
    ident,rest=l.strip().split(maxsplit=1)
    vals=[float(v) for v in re.findall(r'[-+]?\d+(?:\.\d+)?',rest)]
    assert len(vals)==14,(ident,len(vals),vals)
    z,n,re_,rlast,m,b,lm,elm,s,es,x,ex,y,ey=vals
    name='MaNGA_'+ident.replace('∗','')
    bcg.append(dict(population='MaNGA_BCG',system=name,parent_system=name,radius_kpc=re_*rlast,
        gbar_m_s2=10**x,gobs_m_s2=10**y,loggbar_err_dex=ex,loggobs_err_dex=ey,loggbar=x,loggobs=y,
        redshift=z,effective_radius_kpc=re_,last_radius_over_Re=rlast,source='Tian2024_Table1',
        kinematic_method='Jeans_Abel_isotropic',environment='cluster_central_BCG'))

# Check identifiers against actual sky positions before calling samples disjoint.
position_lines=[l for l in (P/'manga_positions_sdss.csv').read_text().splitlines() if not l.startswith('#')]
positions={r['plateifu']:r for r in csv.DictReader(position_lines)}
sys.path.insert(0,str(P.parent/'dm_fit/deps'))
from astropy.io import fits
xcop=[]
for path in (P.parent/'dm_fit/data/xcop').glob('*_mstar.fits'):
    with fits.open(path) as f:
        h=f['MSTAR_RAW'].header
        xcop.append((path.stem.replace('_mstar',''),h['RA'],h['DEC'],h['REDSHIFT']))
crossmatch=[]
for r in bcg:
    pos=positions[r['system'].replace('MaNGA_','')]
    ra,dec=np.radians([float(pos['objra']),float(pos['objdec'])])
    candidates=[]
    for name,ra0,dec0,z0 in xcop:
        a,d=np.radians([ra0,dec0])
        angle=float(np.degrees(np.arccos(np.clip(np.sin(d)*np.sin(dec)+np.cos(d)*np.cos(dec)*np.cos(a-ra),-1,1))))
        candidates.append((angle,name,z0))
    angle,name,z0=min(candidates)
    crossmatch.append(dict(system=r['system'],ra_deg=float(pos['objra']),dec_deg=float(pos['objdec']),redshift=float(pos['z']),
        nearest_xcop=name,separation_arcmin=angle*60,redshift_difference=float(pos['z'])-z0))

def write(name,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (P/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
write('clash_standardized.csv',clash);write('manga_bcg_standardized.csv',bcg)
write('external_standardized.csv',clash+bcg)
write('manga_xcop_crossmatch.csv',crossmatch)

root=P.parent/'dm_fit'
fit=json.loads((root/'results/fit_results.json').read_text())
prev=json.loads((root/'results/two_scale_results.json').read_text())
fixed=json.loads((P.parent/'dm_fixed_radius/stress_test_results.json').read_text())
a0g=fit['models']['RAR']['coefficients'][0]
a0c=fixed['cases']['forward__gbar']['full_cluster_RAR_a0_m_s2']
pc=prev['models']['quadratic_acceleration_radius']['coefficients']
rr=[]
for r in clash+bcg:
    gb=r['gbar_m_s2']; x=r['loggbar']+10;z=np.log10(r['radius_kpc'])
    pred=dict(baryons_only=gb,galaxy_RAR=gb/(-np.expm1(-np.sqrt(gb/a0g))),
        xcop_fixed_radius_RAR=gb/(-np.expm1(-np.sqrt(gb/a0c))),
        published_CLASH_sqrt= float(np.sqrt(2.02e-9*gb)),
        galaxy_quadratic=1e-10*10**np.polynomial.polynomial.polyval(x,fit['models']['poly2']['coefficients']),
        old_R500_radius_quadratic=1e-10*10**(np.array([1,x,z,x*x,x*z,z*z])@pc))
    rr.append(dict(**r,**pred))
write('external_frozen_predictions.csv',rr)
results={}
for pop in ['CLASH_BCG','CLASH_cluster','CLASH_all','MaNGA_BCG']:
    subset=[r for r in rr if r['population']==pop or (pop=='CLASH_all' and r['population'].startswith('CLASH'))]
    out={}
    for key in pred:
        ratios=np.array([r[key]/r['gobs_m_s2'] for r in subset]);e=np.log10(ratios)
        out[key]=dict(median_absolute_acceleration_fraction=float(np.median(abs(ratios-1))),
            rmse_dex=float(np.sqrt(np.mean(e*e))),median_observed_over_predicted=float(np.median(1/ratios)))
    results[pop]=dict(points=len(subset),systems=len(set(r['parent_system'] for r in subset)),models=out)
result=dict(status='Frozen-model transfer to newly acquired published acceleration estimates; no fitting to these targets in this script.',
    notes=['CLASH central BCG and outer points share 20 parent clusters and must stay together in splits.',
        'CLASH mass is spherical NFW mass inferred from combined strong/weak lensing and magnification under assumed cosmology/gravity, not a theory-independent raw force measurement.',
        'CLASH 64 outer points are at fixed 100/200/400/600 kpc, avoiding the R500 target-definition identity.',
        'MaNGA 50 targets are published isotropic spherical Jeans/Abel estimates with photometric Sersic baryons. Radii are Re times rlast/Re from Table1.',
        'All stored logarithms and errors are log10 SI/dex. The CLASH paper equation8 uses natural logarithms; its intercept must not be applied to these log10 numbers.',
        'No full cross-radius covariance matrices are supplied in these tables; uncertainty is not marginalized here.',
        'CLASH names do not match any seven X-COP names, and MaNGA redshifts .027-.148 do not overlap CLASH .187-.686.',
        'SDSS DR17 mangaDrpAll sky positions matched to X-COP FITS cluster centers: no MaNGA/X-COP center overlap; nearest separation is 71.27 arcmin.',
        'MaNGA table lacks host-cluster IDs. Its system and parent_system identify each observed BCG. MaNGA8989-12703 and8991-6102 lie 2.694 arcmin apart with delta-z .002943; a conservative host-based split may keep them together.',
        'published_CLASH_sqrt uses the published g_dagger=2.02e-9 m/s^2 low-acceleration relation. Its CLASH score is in-sample to that publication, not independent validation; MaNGA was also previously compared with this relation by Tian2024.'],
    galaxy_RAR_a0=a0g,xcop_fixed_radius_RAR_a0=a0c,results=results,
    sources={'CLASH_catalog':'https://doi.org/10.26093/cds/vizier.18960070',
        'CLASH_fig2_download':'https://vizier.idia.ac.za/viz-bin/asu-tsv?-source=J%2FApJ%2F896%2F70%2Ffig2&-out.all&-out.max=100',
        'CLASH_table1_download':'https://vizier.idia.ac.za/viz-bin/asu-tsv?-source=J%2FApJ%2F896%2F70%2Ftable1&-out.all&-out.max=100',
        'CLASH_paper':'https://arxiv.org/abs/2001.08340','MaNGA_paper':'https://arxiv.org/abs/2402.12016',
        'MaNGA_positions':'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=select%20plateifu,objra,objdec,z%20from%20mangaDrpAll&format=csv'},
    hashes={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['clash_idia.tsv','clash_table1_idia.tsv','manga_bcg_table1_extracted.txt','manga_positions_sdss.csv']})
(P/'external_results.json').write_text(json.dumps(result,indent=2))
for pop,v in results.items():print(pop,v['points'],{k:round(100*m['median_absolute_acceleration_fraction'],2) for k,m in v['models'].items()})
