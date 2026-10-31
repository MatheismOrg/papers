"""Reproduce a descriptive globular-cluster structural control (not a RAR fit).

Primary author catalog snapshot downloaded 2026-10-27. No velocity-dispersion
mass estimator is applied: the tabulated sigma0 is CENTRAL, not global.
"""
from pathlib import Path
import csv, hashlib, html, json, re, urllib.request
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'
DATA.mkdir(exist_ok=True)
URLS = {
    'combined_table.txt': 'https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/combined_table.txt',
    'current_parameters.html': 'https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/parameter.html',
}
for name,url in URLS.items():
    if not (DATA/name).exists():
        (DATA/name).write_bytes(urllib.request.urlopen(url,timeout=45).read())

html_rows = {}
for tr in re.findall(r'<tr\b[^>]*>(.*?)</tr>', (DATA/'current_parameters.html').read_text(), re.S|re.I):
    cells = re.findall(r'<td\b[^>]*>(.*?)</td>', tr, re.S|re.I)
    if len(cells) != 27: continue
    anchor = re.search(r'<a\b([^>]*)>(.*?)</a>', cells[0], re.S|re.I)
    if not anchor: continue
    name = html.unescape(re.sub('<[^>]*>', '', anchor.group(2))).strip().replace(' ','_')
    # Author uses blue for fitted systems and green for a luminosity-estimated
    # mass. Inspect the object link itself; italic aliases/MF values are unrelated.
    fitted_link = bool(re.search(r'color\s*:\s*blue', anchor.group(1), re.I))
    html_rows[name] = {'fitted_link':fitted_link,'mass_has_reported_error':'&plusmn;' in cells[7]}

columns = ['name','ra_deg','dec_deg','distance_kpc','distance_error_kpc',
           'galactocentric_distance_kpc','galactocentric_distance_error_kpc',
           'n_radial_velocities','n_proper_motions','mass_model_Msun',
           'mass_model_error_Msun','apparent_V_mag','apparent_V_error_mag',
           'mass_to_light_V_Msun_Lsun','mass_to_light_V_error',
           'core_radius_pc','projected_half_light_radius_pc','half_mass_radius_3d_pc',
           'tidal_radius_pc','log_central_density_Msun_pc3','log_half_mass_density_Msun_pc3',
           'log_central_surface_density_Msun_pc2','log_half_mass_surface_density_Msun_pc2',
           'log_half_mass_relaxation_time_yr','log_initial_mass_Msun','remaining_lifetime_Gyr',
           'mass_function_lower_Msun','mass_function_upper_Msun','mass_function_slope',
           'mass_function_slope_error','central_sigma_1d_km_s','central_escape_speed_km_s',
           'core_mass_segregation_parameter','half_mass_segregation_parameter']
rows = []
for line in (DATA/'combined_table.txt').read_text().splitlines():
    if not line or line.startswith('#'): continue
    t = line.split()
    assert len(t) in (34,37), (t[0],len(t))
    d = {columns[0]:t[0]}
    d.update({k:float(v) for k,v in zip(columns[1:],t[1:34])})
    assert d['name'] in html_rows, d['name']
    h = html_rows[d['name']]
    d['author_link_indicates_fitted_mass'] = int(h['fitted_link'])
    d['reported_mass_error_in_html'] = int(h['mass_has_reported_error'])
    d['descriptive_control_selected'] = int(h['fitted_link'] and h['mass_has_reported_error']
        and d['mass_model_error_Msun'] > 0
        and (d['n_radial_velocities'] >= 20 or d['n_proper_motions'] >= 20))
    d['RAR_fit_eligible'] = 0
    rows.append(d)
assert len(rows) == 171
with (ROOT/'globular_cluster_controls.csv').open('w',newline='') as f:
    w = csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
selected = [r for r in rows if r['descriptive_control_selected']]
def quantiles(key):
    a=np.array([r[key] for r in selected]);return {str(q):float(np.quantile(a,q)) for q in [0,.16,.5,.84,1]}
ml = np.array([r['mass_to_light_V_Msun_Lsun'] for r in selected])
summary = {
    'status':'Descriptive structural and dynamical M/L control only; NOT included in acceleration regression or any fit accuracy claim.',
    'source':'Baumgardt author database, updated N-body parameter table, accessed 2026-10-27; methodology Baumgardt & Hilker 2018 plus author-listed subsequent updates.',
    'source_urls':URLS,
    'methodology_reference':{'authors':'H. Baumgardt and M. Hilker','year':2018,'title':'A catalogue of masses, structural parameters, and velocity dispersion profiles of 112 Milky Way globular clusters','doi':'10.1093/mnras/sty1057','arxiv':'https://arxiv.org/abs/1804.08359'},
    'total_catalog_rows':len(rows),'selected_control_rows':len(selected),
    'selection':'Author blue fitted-system link; mass error reported in HTML and positive in ASCII; >=20 radial-velocity or >=20 proper-motion measurements.',
    'mass_to_light_V_quantiles':quantiles('mass_to_light_V_Msun_Lsun'),
    'projected_half_light_radius_pc_quantiles':quantiles('projected_half_light_radius_pc'),
    'mass_Msun_quantiles':quantiles('mass_model_Msun'),
    'fraction_ML_between_1_and_3':float(np.mean((ml>=1)&(ml<=3))),
    'limitations':['Masses are inferred under N-body stellar models, not theory-independent measurements.',
       'No independent stellar mass or radial luminosity-weighted global dispersion is supplied by this compact table.',
       'Central sigma0 cannot be inserted into the Wolf half-light mass formula as if global.',
       'Half-mass and half-light radii differ due to mass segregation; they are not interchangeable.',
       'M/L variation is not an estimated dark-matter fraction. No direct RAR residual or dark-matter exclusion is calculated.',
       'Radius uncertainties are not supplied in the compact table; formal parameter errors do not include every systematic.'],
    'source_sha256':{k:hashlib.sha256((DATA/k).read_bytes()).hexdigest() for k in URLS}
}
(ROOT/'globular_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
