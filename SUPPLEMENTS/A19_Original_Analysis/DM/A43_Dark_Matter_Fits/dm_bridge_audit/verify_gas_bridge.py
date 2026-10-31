"""Independent gas-source audit; no analysis prediction code imported."""
from pathlib import Path
import csv,json,hashlib,collections
import numpy as np
from scipy.integrate import quad
B=Path(__file__).resolve().parent.parent;D=B/'dm_gas_bridge';O=Path(__file__).resolve().parent
G=6.67430e-11;MSUN=1.988409870698051e30;KPC=3.085677581491367e19;A0=1.2e-10
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open()))
J=json.loads((D/'results.json').read_text());rows=read(D/'inputs.csv');full=read(D/'full_disk_components.csv');cl=read(D/'cluster_inputs.csv');quadr=json.loads((D/'cluster_quadrature.json').read_text())
for key,name in [('protocol','PROTOCOL.txt'),('script','run.py'),('input','inputs.csv')]:assert J[key+'_sha256']==sha(D/name)
for path,value in J['source_hashes'].items():assert sha(B/path)==value
for a,b in [('full_disk_components.csv','dm_halo_benchmark/matched_components.csv'),('cluster_inputs.csv','dm_cluster_repair/inputs.csv'),('cluster_quadrature.json','dm_cluster_repair/quadrature.json')]:assert sha(D/a)==sha(B/b)
selected=[]
for r in full:
 gg=float(r['gas_accel_over_aref'])*1e-10;gs=float(r['stellar_accel_over_aref'])*1e-10
 assert np.isclose(gg+gs,float(r['gbar_m_s2']),rtol=1e-10,atol=0)
 if gg>0 and gg/(gg+gs)>=.5 and float(r['radius_kpc'])<=30:selected.append(r)
assert len(selected)==402 and len(rows)==488
for i,r in enumerate(rows):
 assert int(r['row'])==i
 if r['domain']=='disk':
  s=selected[i];assert int(s['row'])==int(r['origin_row']) and s['galaxy']==r['system']
  for key in ['gbar_m_s2','gobs_m_s2','radius_kpc']:assert float(r[key])==float(s[key])
  gg=float(s['gas_accel_over_aref'])*1e-10;gs=float(s['stellar_accel_over_aref'])*1e-10;assert float(r['ggas_m_s2'])==gg and float(r['gstar_m_s2'])==gs and float(r['gas_force_fraction'])==gg/(gg+gs)
 else:
  s=cl[int(r['cluster_quadrature_row'])];assert r['system']==s['system'] and int(r['origin_row'])==int(s['broad_row'])
  for key in ['gbar_m_s2','gobs_m_s2','radius_kpc']:assert float(r[key])==float(s[key])
N=len(rows);I=np.arange(N);sys=np.array([r['system'] for r in rows]);source=np.array([r['source'] for r in rows]);dom=np.array([r['domain'] for r in rows]);radius=np.array([float(r['radius_kpc']) for r in rows]);bar=np.array([float(r['gbar_m_s2']) for r in rows]);gas=np.array([float(r['ggas_m_s2']) for r in rows]);target=np.log10([float(r['gobs_m_s2']) for r in rows]);DI=I[dom=='disk'];CI=I[dom=='cluster'];SC={'disk':DI,'cluster':CI,'joint':I}
# Integrate published density and piecewise mass derivatives directly; these
# independent adaptive integrals test the frozen quadrature, including core.
pars={r['system']:r for r in read(B/'dm_cluster_repair/clash_gas_parameters.csv')};prof=collections.defaultdict(list)
for r in read(B/'dm_cluster_repair/gas_profiles.csv'):prof[r['system']].append((float(r['radius_kpc']),float(r['mgas_msun'])))
def mass_integral(row,rc):
 r=float(row['radius_kpc'])
 if row['source']=='CLASH':
  p=pars[row['system'].replace('CLASH_','')];n0,r0,al,re,beta,n1,re1,beta1=[float(p[k]) for k in ['ne0','r0','alpha','re0','beta0','ne1','re1','beta1']]
  def fun(u):
   s=r*u;rho=n0*1e5*(s/r0)**(-al)*(1+(s/re)**2)**(-1.5*beta)
   if n1>0:rho+=n1*1e5*(1+(s/re1)**2)**(-1.5*beta1)
   return 4*np.pi*r*s*s*rho*np.exp(-s/rc)
  return quad(fun,0,1,epsabs=1e-3,epsrel=3e-11,limit=200)[0]
 points=prof[row['system']];r0,m0=points[0]
 mass=quad(lambda s:3*m0*s*s/r0**3*np.exp(-s/rc),0,min(r,r0),epsabs=1e-3,epsrel=2e-11)[0]
 for (lo,ml),(hi,mh) in zip(points[:-1],points[1:]):
  if lo>=r:break
  p=np.log(mh/ml)/np.log(hi/lo)
  mass+=quad(lambda s:p*ml*(s/lo)**p/s*np.exp(-s/rc),lo,min(hi,r),epsabs=1e-3,epsrel=2e-11)[0]
 return mass
maxquad=0.;masschecks=0
for row,q in zip(cl,quadr):
 assert int(row['row'])==q['row'];nodes=np.array(q['radius_nodes_kpc']);ww=np.array(q['mass_weights_msun']);assert np.all(nodes>=0) and np.all(nodes<=float(row['radius_kpc'])) and np.all(ww>=0)
 for rc in [10.,100.,400.,1000.,10000.,np.inf]:
  num=np.dot(ww,np.exp(-nodes/rc));direct=mass_integral(row,rc);err=abs(num/direct-1);maxquad=max(maxquad,err);masschecks+=1
  assert err<2e-6,(row['system'],row['radius_kpc'],rc,err)
def weight(ii):
 result=np.zeros(len(ii));domains=sorted(set(dom[ii]))
 for d in domains:
  sources=sorted(set(source[ii][dom[ii]==d]))
  for s in sources:
   hosts=sorted(set(sys[ii][(dom[ii]==d)&(source[ii]==s)]))
   for h in hosts:
    z=(dom[ii]==d)&(source[ii]==s)&(sys[ii]==h);result[z]=1/(len(domains)*len(sources)*len(hosts)*sum(z))
 assert abs(sum(result)-1)<1e-12
 return result
def pred(nu,c,ii):
 eta,logrc=c;K=gas[ii].copy()
 for j,i in enumerate(ii):
  if dom[i]=='cluster':
   q=quadr[int(rows[i]['cluster_quadrature_row'])];mass=np.dot(np.array(q['mass_weights_msun']),np.exp(-np.array(q['radius_nodes_kpc'])/10**logrc));K[j]=G*MSUN*mass/(radius[i]*KPC)**2
 src=bar[ii]+eta*K;assert np.all(src>0)
 val=(src+np.sqrt(src*src+4*A0*src))/2 if nu=='simple' else src/(1-np.exp(-np.sqrt(src/A0)))
 return np.log10(val)
def metrics(p,ii,nested=True):
 e=p-target[ii];out={'B_dex':float(np.sqrt(np.dot(weight(ii),e*e))),'rows':len(ii),'systems':len(set(sys[ii])),'median_acceleration_error_percent':float(100*np.median(abs(10**e-1))),'median_speed_error_percent':float(100*np.median(abs(10**(.5*e)-1)))}
 if nested:
  for key,arr in [('domains',dom),('sources',source)]:out[key]={s:metrics(p[arr[ii]==s],ii[arr[ii]==s],False) for s in sorted(set(arr[ii]))}
 return out
def comp(a,b):
 assert set(a)==set(b),(set(a),set(b))
 for k,v in a.items():
  if isinstance(v,dict):comp(v,b[k])
  else:assert np.isclose(v,b[k],rtol=2e-10,atol=1e-12),(k,v,b[k])
def folds(seed):
 rng=np.random.default_rng(seed);mp={}
 for s in sorted(set(source)):
  perm=rng.permutation(np.unique(sys[source==s]));mp.update({g:j%5 for j,g in enumerate(perm)})
 return np.array([mp[g] for g in sys])
fits=collections.defaultdict(list);maxcost=0.
for i,f in enumerate(J['fits']):
 assert f['fit_id']==i and f['success'];ii=np.array(f['train_rows']);p=pred(f['nu'],f['coefficients'],ii);value=np.dot(weight(ii),(p-target[ii])**2)/2;maxcost=max(maxcost,abs(value-f['cost']));assert abs(value-f['cost'])<1e-10
 lo,hi=map(np.array,f['bounds']);c=np.array(f['fitted_parameters']);assert np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10)
 assert f['disk_only']==bool(np.all(dom[ii]=='disk'))
 if f['disk_only']:assert len(c)==1 and abs(f['coefficients'][1]-np.log10(400))<1e-14
 fits[f['tag']].append(f)
tables=collections.defaultdict(list)
for r in read(D/'predictions.csv'):tables[r['scope'],r['nu'],r['model'],r['split']].append(r)
def checkfit(fid,tag,ii,cf):
 f=J['fits'][fid];assert f['tag']==tag and f['train_rows']==ii.tolist() and np.array_equal(f['coefficients'],cf)
 assert f['cost']==min(z['cost'] for z in fits[tag] if z['success'])
assert len(tables)==104;maxpred=0.;counts=0;pc={}
prior=json.loads((B/'dm_cluster_repair/results.json').read_text())['cases']['combined_forward_mid']
for (scope,nu,model,split),tab in tables.items():
 ii=np.array([int(r['row']) for r in tab]);p=np.log10([float(r['gpred_m_s2']) for r in tab]);expected=np.empty(len(ii));si=SC[scope]
 if split=='all':
  assert np.array_equal(ii,si);d=J['scopes'][scope][nu][model];co=d['coefficients'];score=d['all']
  if model=='fitted':checkfit(d['selected_fit_id'],scope+':'+nu+':all',ii,co)
  elif model=='zero':assert co==[0.,np.log10(400.)]
  elif model=='literature':assert co==[10.,np.log10(400.)]
  else:
   old=prior[nu+'_repair']['source_holdout']['XCOP']['coefficients'];assert co==[10**old[0],old[1]]
  expected=pred(nu,co,ii)
 elif split.startswith('cv_'):
  assert np.array_equal(ii,si);fd=folds(int(split[3:]));assert [int(r['fold']) for r in tab]==fd[ii].tolist();d=J['scopes'][scope][nu][model]['CV'][split[3:]];score=d['metrics']
  for f in range(5):
   test=ii[fd[ii]==f];train=si[fd[si]!=f];co=d['fold_coefficients'][str(f)]['coefficients'];fid=d['fold_coefficients'][str(f)]['selected_fit_id'];assert not set(sys[test])&set(sys[train])
   if model=='fitted':checkfit(fid,scope+':'+nu+':'+split+':'+str(f),train,co)
   expected[fd[ii]==f]=pred(nu,co,test)
 else:
  d=J['transfers'][nu+':'+split];co=d['coefficients'];score=d['metrics'];assert ii.tolist()==d['test_rows'];train=I[~np.isin(sys,sys[ii])];assert train.tolist()==d['train_rows'] and not set(sys[train])&set(sys[ii]);checkfit(d['selected_fit_id'],'joint:'+nu+':'+split,train,co)
  if split.startswith('source_'):assert set(source[ii])=={split[7:]} and split[7:] not in set(source[train])
  else:assert set(dom[ii])=={split[7:]} and split[7:] not in set(dom[train])
  assert d['Rc_externally_fixed_400kpc']==bool(split=='domain_cluster');expected=pred(nu,co,ii)
 err=float(np.max(abs(p-expected)));maxpred=max(maxpred,err);assert err<1e-10;comp(metrics(p,ii),score)
 for j,r in enumerate(tab):
  i=ii[j];assert r['system']==sys[i] and r['domain']==dom[i] and r['source']==source[i];assert abs(float(r['log10_residual'])-(p[j]-target[i]))<1e-12
  cf=json.loads(r['coefficients_json']);assert abs(pred(nu,cf,np.array([i]))[0]-p[j])<1e-10
  fid=int(r['selected_fit_id'])
  if fid>=0:assert cf==J['fits'][fid]['coefficients']
 counts+=len(ii)
 if split=='cv_20261027':pc[scope,nu,model]=p
# Reproduce paired host resampling, preserving original domain/source weights.
rng=np.random.default_rng(20261029);bschecks=0
def bs(p,ii,draws):
 ee=(p-target[ii])**2;ds=[]
 for d in sorted(set(dom[ii])):
  ss=[]
  for s in sorted(set(source[ii][dom[ii]==d])):
   gs=np.unique(sys[ii][source[ii]==s]);v=np.array([np.mean(ee[sys[ii]==g]) for g in gs]);ss.append(np.mean(v[draws[d,s]],axis=1))
  ds.append(np.mean(ss,axis=0))
 return np.sqrt(np.mean(ds,axis=0))
for nu in ['simple','exp']:
 cases=[('disk_CLASH_frozen_minus_zero',DI,pc['disk',nu,'CLASH_frozen'],pc['disk',nu,'zero']),('disk_literature_minus_zero',DI,pc['disk',nu,'literature'],pc['disk',nu,'zero']),('disk_joint_minus_separate',DI,pc['joint',nu,'fitted'][DI],pc['disk',nu,'fitted']),('cluster_joint_minus_separate',CI,pc['joint',nu,'fitted'][CI],pc['cluster',nu,'fitted']),('joint_fitted_minus_zero',I,pc['joint',nu,'fitted'],pc['joint',nu,'zero'])]
 for name,ii,a,b in cases:
  draws={}
  for d in sorted(set(dom[ii])):
   for s in sorted(set(source[ii][dom[ii]==d])):ng=len(set(sys[ii][source[ii]==s]));draws[d,s]=rng.integers(0,ng,(5000,ng))
  delta=bs(a,ii,draws)-bs(b,ii,draws);saved=J['descriptive_bootstrap'][nu+':'+name]
  assert np.allclose(np.percentile(delta,[2.5,97.5]),saved['percentile_95_dex'],rtol=1e-10,atol=1e-12);assert abs(metrics(a,ii)['B_dex']-metrics(b,ii)['B_dex']-saved['difference_dex'])<1e-12;bschecks+=1
# Independently preserve every nonpositive extrapolated gas-source force.
allfull=read(D/'full_disk_predictions.csv');ft=collections.defaultdict(list)
for r in allfull:ft[r['nu'],r['model']].append(r)
invalid={}
for (nu,model),tab in ft.items():
 assert len(tab)==2700;good=[];bad=[]
 for r,s in zip(tab,full):
  assert int(r['origin_row'])==int(s['row']);eta=float(r['eta']);gg=float(s['gas_accel_over_aref'])*1e-10;gs=float(s['stellar_accel_over_aref'])*1e-10;src=gs+(1+eta)*gg
  assert np.isclose(src,float(r['gsource_m_s2']),rtol=1e-14,atol=0);assert (r['valid_positive_source']=='True')==(src>0)
  if src>0:
   pp=(src+np.sqrt(src*src+4*A0*src))/2 if nu=='simple' else src/(1-np.exp(-np.sqrt(src/A0)));assert np.isclose(pp,float(r['gpred_m_s2']),rtol=1e-12,atol=0);ee=np.log10(pp/float(s['gobs_m_s2']));assert abs(ee-float(r['log10_residual']))<1e-10;good.append(ee)
  else:assert r['gpred_m_s2']==r['log10_residual']=='';bad.append(r)
 saved=J['full_disk_diagnostics'][nu+':'+model];assert saved['invalid_rows']==len(bad) and saved['total_rows']==2700 and saved['not_a_whole_sample_score_when_invalid']==bool(bad);assert saved['invalid_origin_rows']==[int(r['origin_row']) for r in bad];assert saved['invalid_systems']==sorted(set(r['system'] for r in bad));assert abs(saved['valid_conditional_row_RMS_dex']-np.sqrt(np.mean(np.array(good)**2)))<1e-10;invalid[nu+':'+model]=len(bad)
assert counts==J['prediction_rows']==32380 and len(fits)==J['fit_sets']==104 and len(J['fits'])==J['optimizer_attempts']==382 and J['optimizer_failures']==[]
report={'status':'PASS','independent_method':'Raw component and source joins; adaptive integration of cluster gas density and mass-profile derivatives; all forecasts, fits, objective minima, metrics, exclusions, paired host resampling and nonpositive-force retention. No analysis-code imports.','prediction_tables':len(tables),'prediction_rows':counts,'fit_sets':len(fits),'optimizer_attempts':len(J['fits']),'cluster_density_integrals':masschecks,'max_cluster_quadrature_relative_discrepancy':maxquad,'max_prediction_discrepancy_dex':maxpred,'max_objective_discrepancy':maxcost,'descriptive_bootstrap_comparisons':bschecks,'full_disk_diagnostic_rows':len(allfull),'full_disk_invalid_counts':invalid,'limits':['One eta is used across gas phases, but disk constant-weight source scaling and spherical finite-Rc cluster integral are different geometry calculations.','A disk evaluation radius cutoff cannot bound finite-Rc source approximation error: exterior sources contribute. No exact finite-Rc disk reconstruction is present.','The cold-HI versus hot-ICM extension is the tested hypothesis; the published hot-gas-only model is not falsified by its failure.','Cluster zero values in disk-only ggas/gstar columns are unused placeholders; measured cluster gas enters through nonzero source-density quadrature.','Disk-only fits do not identify Rc; cold-to-cluster transfer assumes external fixed Rc=400 kpc.','Frozen CLASH coefficients include CLASH targets from earlier work, so CLASH scores are descriptive. All newly refitted test targets are excluded as declared.','Cluster speed percentages are algebraic square-root transforms, not measured disk orbital-speed errors.','Resampling intervals are descriptive conditional on frozen forecasts and omit calibration/systematic/model-selection uncertainty.'],'result_sha256':sha(D/'results.json'),'audit_script_sha256':sha(Path(__file__))}
(O/'gas_bridge_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
