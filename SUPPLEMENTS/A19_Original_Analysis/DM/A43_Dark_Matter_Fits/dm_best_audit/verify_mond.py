"""Independent prediction, score, split and optimization-record audit.
Does not import the analysis implementation or use its prediction functions.
"""
from pathlib import Path
import csv, json, hashlib, collections, re
import numpy as np
B=Path(__file__).resolve().parent.parent; D=B/'dm_mond_repairs'; O=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def readcsv(p):return list(csv.DictReader(p.open()))
rows=readcsv(D/'input.csv'); N=len(rows); idx=np.arange(N)
pop=np.array([r['population'] for r in rows]); group=np.array([r['group'] for r in rows]); env=np.array([r['environment'] for r in rows]); src=np.array(['CLASH' if r['source']=='Tian2020_CDS_fig2' else r['source'] for r in rows])
b=np.array([float(r['gbar_m_s2'])/1e-10 for r in rows]); rad=np.array([float(r['radius_kpc']) for r in rows]); y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]); sparc=idx[pop=='SPARC']
J=json.loads((D/'results/results.json').read_text())
assert J['input_sha256']==sha(D/'input.csv')==sha(B/'dm_unified_regular/augmented_input.csv')
assert J['script_sha256']==sha(D/'run.py');assert J['protocol_sha256']==sha(D/'PROTOCOL.txt')

def pred(k,c,ii):
 z=b[ii];r=rad[ii]
 if k in ['simple','standard','nu_n']:
  n=1 if k=='simple' else 2 if k=='standard' else c[1]
  nu=((1+np.sqrt(1+4*(z/10**c[0])**(-n)))/2)**(1/n)
  return np.log10(z*nu)
 if k=='RAR': a=10**c[0]
 elif k=='two_log_bumps':
  q=np.log10(z*r)
  a=10**c[0]+10**c[1]*np.exp(-.5*((q-c[2])/10**c[3])**2)+10**c[4]*np.exp(-.5*((q-c[5])/10**c[6])**2)
 else:
  m=z*r*r;s=z*r;S=1/(1+(10**c[5]/r)**4);W=1/(1+(z/100)**2)
  if k=='mixed_transient':a=10**c[0]+(10**c[1]*m**(-c[3])*S*W+10**c[2]*s**c[4])/(1+r/10**c[6])
  elif k=='mixed_potential_taper':a=10**c[0]+S*W*(10**c[1]*m**(-c[3])+10**c[2]*s**c[4])*s/(s+10**c[6])
  else:raise KeyError(k)
 return np.log10(z/(-np.expm1(-np.sqrt(z/a))))

def weights(ii):
 res=np.zeros(len(ii)); ps=set(pop[ii])
 for p in ps:
  gs=set(group[ii][pop[ii]==p])
  for g in gs:
   mask=(pop[ii]==p)&(group[ii]==g);res[mask]=1/(len(ps)*len(gs)*np.count_nonzero(mask))
 assert abs(res.sum()-1)<1e-12
 return res

def scores(pr,ii):
 e=pr-y[ii]
 return {'B_dex':float(np.sqrt(np.dot(weights(ii),e*e))), 'RMS_dex':float(np.sqrt(np.mean(e*e))), 'median_acceleration_error_percent':float(100*np.median(abs(np.power(10,e)-1))), 'median_speed_error_percent':float(100*np.median(abs(np.power(10,e/2)-1)))}

def metriccheck(pr,ii,m):
 s=scores(pr,ii);assert abs(s['B_dex']-m['B_dex'])<1e-10
 assert m['rows']==len(ii); assert m['groups']==len(set(group[ii]))
 for label, field, labs in [('populations','population',pop),('sources','source',src)]:
  for v,d in m[label].items():
   use=labs[ii]==v;ss=scores(pr[use],ii[use]);assert d['rows']==sum(use);assert d['groups']==len(set(group[ii[use]]))
   for key in ['RMS_dex','median_acceleration_error_percent','median_speed_error_percent']:assert abs(ss[key]-d[key])<1e-8,(v,key)
   assert abs(ss['B_dex']-d['group_RMS_dex'])<1e-10

def folds(ii,seed):
 rng=np.random.default_rng(seed);assigned={}
 for p in sorted(set(pop[ii])):
  members=sorted(set(group[ii][pop[ii]==p]));count=[sum(assigned.get(g,-1)==f for g in members) for f in range(5)]
  new=[g for g in members if g not in assigned]
  for g in rng.permutation(new):
   f=int(np.argmin(count));assigned[g]=f;count[f]+=1
 return np.array([assigned[g] for g in group[ii]])

def trainset(tag):
 if tag=='all':return idx
 if tag=='sparc_only_all':return sparc
 if tag.startswith('sparc_only_cv_'):
  seed,fold=map(int,tag.removeprefix('sparc_only_cv_').split('_'));return sparc[folds(sparc,seed)!=fold]
 if tag.startswith('cv_'):
  seed,fold=map(int,tag.removeprefix('cv_').split('_'));return idx[folds(idx,seed)!=fold]
 if tag.startswith('population_'): test=idx[pop==tag.removeprefix('population_')]
 elif tag.startswith('host_'):test=idx[(pop=='dSph')&np.isin(env,{'MW':['Milky_Way','G/L'],'M31':['M31','A/L']}[tag.removeprefix('host_')])]
 elif tag.startswith('source_'):test=idx[src==tag.removeprefix('source_')]
 else:raise KeyError(tag)
 return idx[~np.isin(group,list(set(group[test])))]

maxerr=0.;tables=[];prediction_count=0;selected=[]
def checktable(name,k,c,ii,metric,foldids=None):
 global maxerr,prediction_count
 rr=readcsv(D/'results'/(name+'.csv'));assert len(rr)==len(ii)
 pp=np.array([float(r['gpred_m_s2']) for r in rr]);lp=np.log10(pp/1e-10)
 expected=np.empty(len(ii))
 for j,(record,i) in enumerate(zip(rr,ii)):
  assert int(record['row'])==i
  for field in ['population','system','group','environment']:assert record[field]==rows[i][field]
  assert record['source']==src[i];assert record['source_detail']==rows[i]['source']
  for field in ['gbar_m_s2','gobs_m_s2','radius_kpc']:assert np.isclose(float(record[field]),float(rows[i][field]),rtol=1e-13,atol=0)
  assert abs(float(record['log10_residual'])-(lp[j]-y[i]))<1e-11
  if foldids is None: assert record['fold']==''
  else: assert record['fold']==str(foldids[j])
 if isinstance(c,dict):
  for label,cc in c.items():
   mask=np.asarray(foldids).astype(str)==str(label);expected[mask]=pred(k,cc,ii[mask])
 else:expected=pred(k,c,ii)
 er=float(max(abs(expected-lp)));maxerr=max(maxerr,er);assert er<1e-10,(name,er)
 metriccheck(lp,ii,metric);tables.append(name);prediction_count+=len(ii)

for k,d in J['models'].items():
 checktable(k+'_all',k,d['coefficients'],idx,d['training']);selected.append((k,'all',d['coefficients']))
 for seed,cv in d['CV'].items():
  f=folds(idx,int(seed));cc={j:c for j,c in enumerate(cv['coefficients'])}
  checktable(k+'_cv_'+seed,k,cc,idx,cv['metrics'],f)
  for j,c in cc.items():
   tr=trainset('cv_'+seed+'_'+str(j));te=idx[f==j];assert not set(group[tr])&set(group[te]);selected.append((k,'cv_'+seed+'_'+str(j),c))
 cc=d['population_holdout']['coefficients'];checktable(k+'_population_holdout',k,cc,idx,d['population_holdout']['metrics'],pop)
 for p,c in cc.items():
  tr=trainset('population_'+p);assert not set(group[tr])&set(group[pop==p]);selected.append((k,'population_'+p,c))
 for prefix,section in [('host','host_holdout'),('source','source_holdout')]:
  for label,h in d[section].items():
   tr=trainset(prefix+'_'+label);ii=idx[~np.isin(idx,tr)]
   # Every excluded group belongs entirely to the target category in this dataset.
   if prefix=='host':assert len(ii)=={'MW':26,'M31':18}[label]
   else:assert len(ii)=={'XCOP':21,'CLASH':84}[label]
   assert len(ii)==h['excluded_rows'];assert len(tr)==h['training_rows'];assert not set(group[tr])&set(group[ii])
   checktable(k+'_'+prefix+'_'+label,k,h['coefficients'],ii,h['metrics'],[label]*len(ii));selected.append((k,prefix+'_'+label,h['coefficients']))
for k,d in J['sparc_only'].items():
 checktable(k+'_sparc_only_all',k,d['coefficients'],sparc,d['training']);selected.append((k,'sparc_only_all',d['coefficients']))
 for seed,cv in d['CV'].items():
  f=folds(sparc,int(seed));cc={j:c for j,c in enumerate(cv['coefficients'])};checktable(k+'_sparc_only_cv_'+seed,k,cc,sparc,cv['metrics'],f)
  for j,c in cc.items():
   tr=trainset('sparc_only_cv_'+seed+'_'+str(j));te=sparc[f==j];assert not set(group[tr])&set(group[te]);selected.append((k,'sparc_only_cv_'+seed+'_'+str(j),c))
assert len(tables)==79;assert len(selected)==239
attempts=collections.defaultdict(list);maxcosterr=0.
for entry in J['fit_log']:
 k=entry['model'];tag=entry['split'];ii=trainset(tag);assert len(ii)==entry['train_rows']
 c=np.array(entry['coefficients']);spec=J['spec'][k];assert np.all(c>=np.array(spec['lo'])-1e-10) and np.all(c<=np.array(spec['hi'])+1e-10)
 e=pred(k,c,ii)-y[ii];cost=.5*np.dot(weights(ii),e*e);ce=abs(cost-entry['cost']);maxcosterr=max(maxcosterr,ce);assert ce<1e-10
 attempts[k,tag].append(entry)
for k,tag,c in selected:
 aa=attempts[k,tag];assert len(aa)==len(J['spec'][k]['starts']);assert {a['start_index'] for a in aa}==set(range(len(aa)))
 best=min(aa,key=lambda a:a['cost']);assert np.max(abs(np.array(c)-best['coefficients']))<1e-12
assert len(J['fit_log'])==J['optimizer_calls']==471;assert J['fit_calls']==len(attempts)==239
# Check prior comparator files use exactly the same rows, folds, groups and targets.
for split in ['cv_20261027','population_holdout']:
 a=readcsv(D/'results'/('RAR_'+split+'.csv'));b0=readcsv(B/'dm_unified_regular/results/augmented'/('mixed_regular_'+split+'.csv'))
 for a0,b1 in zip(a,b0):
  for f in ['row','group','population','fold','gbar_m_s2','gobs_m_s2','radius_kpc']:assert a0[f]==b1[f]
report={'status':'passed','prediction_tables':len(tables),'prediction_rows':prediction_count,'fit_sets':len(attempts),'optimizer_attempts':len(J['fit_log']),'failed_optimizer_attempts':sum(not a['success'] for a in J['fit_log']),'max_independent_prediction_error_dex':maxerr,'max_independent_cost_error':maxcosterr,'input_hash':sha(D/'input.csv'),'result_hash':sha(D/'results/results.json'),'script_hash':sha(D/'run.py'),'audit_script_hash':sha(Path(__file__)),'all_groups_hosts_sources_and_bounds_checked':True,'limitations':'No refitting, likelihood/systematic validation or universal-law claim. Physical limit audit is separate.'}
(O/'mond_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
