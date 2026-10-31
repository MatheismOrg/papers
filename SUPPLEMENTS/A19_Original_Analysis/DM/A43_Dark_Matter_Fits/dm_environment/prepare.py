from pathlib import Path
import csv,json,hashlib,re
import numpy as np
R=Path(__file__).resolve().parent
rows=list(csv.DictReader((R.parent/'dm_broad_joint/combined_input.csv').open()));dw={r['system']:r for r in csv.DictReader((R.parent/'dm_broad_dwarfs/results/dwarfs_benchmark.csv').open())}
h=[]
for r in csv.DictReader((R/'data/Weisz2015_matched_history.csv').open()):
 d=dw[r['system']];r.update({k:d[k] for k in ['radius_kpc','gbar_m_s2','gobs_m_s2','stellar_mass_Msun','err_log10_gobs_dex_stat_approx','err_log10_gbar_dex_stat_approx']});h.append(r)
assert len(h)==15 and len(set(r['system'] for r in h))==15
with (R/'history_input.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(h[0]));w.writeheader();w.writerows(h)
def norm(s):return re.sub(r'^(NGC|UGC|IC)0+(\d)',r'\1\2',re.sub('[^A-Z0-9]','',s.upper()))
field={}
for l in (R/'data/Table3_Chae_etal_2021.txt').read_text().splitlines():
 a=l.split()
 if len(a)!=5:continue
 try:v=list(map(float,a[1:]))
 except ValueError:continue
 field[norm(a[0])]=dict(zip(['log_eN_MC','err_log_eN_MC','log_eN_NC','err_log_eN_NC'],v))
assert len(field)==109
# Parse whitespace-delimited data rows of the original SPARC Table1; verify name matches.
dist={}
for l in (R.parent/'dm_fit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
 try:
  a=l.split();n=a[0];tt=int(a[1]);dd=float(a[2]);err=float(a[3]);assert len(a)>=18
 except (ValueError,IndexError,AssertionError):continue
 if n and dd>0:dist[norm(n)]=dd
s=[]
for i,r in enumerate(rows):
 if r['population']!='SPARC' or norm(r['system']) not in field:continue
 s.append(dict(r,original_row=i,distance_Mpc=dist[norm(r['system'])],**field[norm(r['system'])]))
with (R/'sparc_environment_input.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(s[0]));w.writeheader();w.writerows(s)
meta={'history_matches':len(h),'history_host_counts':{host:sum(r['host']==host for r in h) for host in ['MW','M31']},'environment_catalog':109,'SPARC_matched_galaxies':len(set(r['system'] for r in s)),'SPARC_matched_rows':len(s),'distance_band_counts':{str(j):len(set(r['system'] for r in s if np.digitize(r['distance_Mpc'],[10,20,40,80])==j)) for j in range(5)},'input_hashes':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'history_input.csv',R/'sparc_environment_input.csv',R/'data/Weisz2015_matched_history.csv',R/'data/Table3_Chae_etal_2021.txt']},'history_source':'Weisz et al.2015 arXiv1503.05195 Table1, p3; tau90 total errors second entries, log10 year units. Manually transcribed and visually checked against PDF. Host distances are older published photometric distances, not dynamical masses.','environment_source':'Chae et al.2021 ApJ921:104 Table3 from author download. Table2 fitted external fields explicitly excluded. eN normalized by fixed a0=1.2e-10 SI. Published MC/NC scenarios not direct field measurements.'}
(R/'provenance.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
