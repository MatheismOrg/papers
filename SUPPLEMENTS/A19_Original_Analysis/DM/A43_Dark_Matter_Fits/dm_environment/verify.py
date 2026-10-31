from pathlib import Path
import csv,json,hashlib
import numpy as np
R=Path(__file__).resolve().parent;j=json.loads((R/'results/results.json').read_text());prov=json.loads((R/'provenance.json').read_text())
for p,h in prov['input_hashes'].items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==h
checks=[]
for p in sorted((R/'results').glob('*.csv')):
 rows=list(csv.DictReader(p.open()));e=np.array([np.log10(float(r['gpred_m_s2'])/float(r['gobs_m_s2'])) for r in rows]);rec=np.array([float(r['log10_residual']) for r in rows]);assert max(abs(e-rec))<5e-14
 for name in {r['system'] for r in rows}:assert len({r['fold'] for r in rows if r['system']==name})==1
 if p.name.startswith('environment_'):
  scenario=p.stem.split('_')[1];end=next(s for s in ['distance_bands','20261028','7331','48105'] if p.stem.endswith(s));model=p.stem[len('environment_'+scenario+'_'):-len('_'+end)];v=j['environment'][scenario][model]['validation'][end]
  mse=np.mean([np.mean(e[[r['system']==g for r in rows]]**2) for g in {r['system'] for r in rows}]);expected=[]
  for r in rows:
   c=v['coefficients'][r['fold']];gb=float(r['gbar_m_s2']);a0=1e-10*10**c[0]
   if model=='RAR':gp=gb/(1-np.exp(-np.sqrt(gb/a0)))
   else:
    ge=0 if model=='simple' else 1.2e-10*10**float(r['log_eN_'+scenario])*(10**c[1] if model=='EFE_scale' else 1)
    y=gb/a0;ee=ge/a0;dd=1+ee/y
    gp=gb*(.5+np.sqrt(dd*dd/4+dd/y)-np.sqrt(ee*(1+ee/4))/y)
   expected.append(gp)
  assert np.max(abs(np.log10(np.array(expected)/np.array([float(r['gpred_m_s2']) for r in rows]))))<1e-11
 elif p.name.startswith('history_'):
  end='host_holdout' if p.stem.endswith('host_holdout') else 'leave_one_dwarf_out';model=p.stem[len('history_'):-len('_'+end)];v=j['history'][model]['validation'][end]
  mse=np.mean([np.mean(e[[r['host']==h for r in rows]]**2) for h in ['MW','M31']])
  if end=='host_holdout':assert all(r['fold']==r['host'] for r in rows)
  expected=[]
  for r in rows:
   c=v['coefficients'][r['fold']];aa=c[0];tau=2 if model=='history_2Gyr' else 10 if model=='history_10Gyr' else 5
   if model.startswith('history'):aa+=c[1]*np.exp(-min(10**(float(r['log_tau90_yr'])-9),13.8)/tau)
   if model in ['mass','history_mass']:aa+=c[-1]*np.log10(float(r['stellar_mass_Msun'])/1e6)
   if model in ['distance','history_distance']:aa+=c[-1]*np.log10(float(r['Dhost_kpc'])/100)
   gb=float(r['gbar_m_s2']);expected.append(gb/(1-np.exp(-np.sqrt(gb/(1e-10*10**aa)))))
  assert max(abs(np.log10(np.array(expected)/np.array([float(r['gpred_m_s2']) for r in rows]))))<1e-12
 else:
  v={'metrics':j['broad_power_sum_host_excluded']['metrics_on_15_matches']};mse=np.mean([np.mean(e[[r['host']==h for r in rows]]**2) for h in ['MW','M31']])
 assert abs(np.sqrt(mse)-v['metrics']['RMS_dex'])<1e-12
 checks.append(p.name)
assert len(checks)==49
assert j['broad_power_sum_host_excluded']['fits']['MW']['excluded_rows']==26
assert j['broad_power_sum_host_excluded']['fits']['M31']['excluded_rows']==18
assert not j['optimizer_warnings']
(R/'verification.json').write_text(json.dumps({'status':'passed','prediction_files_checked':len(checks),'checks':['input hashes','direct Chae Eq2 versus stable implementation','direct history predictions','independent equal-galaxy or equal-host RMS','unsplit validation groups','all44 dwarf host exclusions','no optimizer warnings'],'files':checks},indent=2));print('Verified',len(checks),'prediction tables and independent equations.')
