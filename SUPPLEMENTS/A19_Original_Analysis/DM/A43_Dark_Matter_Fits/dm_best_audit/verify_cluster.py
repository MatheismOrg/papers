"""Independent cluster density-integral, provenance, split and fit audit.
No imports from analysis scripts or external author code.
"""
from pathlib import Path
import json,csv,hashlib,ast,collections,sys
import numpy as np
from scipy.integrate import quad
B=Path(__file__).resolve().parent.parent;D=B/'dm_cluster_repair';O=Path(__file__).resolve().parent
sys.path.insert(0,str(B/'dm_fit/deps'))
from astropy.io import fits
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();readcsv=lambda p:list(csv.DictReader(p.open()))
J=json.loads((D/'results.json').read_text());Q=json.loads((D/'quadrature.json').read_text());rr=readcsv(D/'inputs.csv');P=readcsv(D/'predictions.csv');V=json.loads((D/'provenance.json').read_text());SP=json.loads((D/'source_provenance.json').read_text());DI=json.loads((D/'diagnostics.json').read_text());N=len(rr);idx=np.arange(N)
assert N==86;assert sha(D/'PROTOCOL.txt')==J['protocol_sha256'];assert sha(D/'inputs.csv')==J['input_sha256'];assert sha(D/'quadrature.json')==J['quadrature_sha256']
assert sha(D/'clash_gas_parameters.csv')==SP['parameter_csv_sha256']
for fn,h in V['source_hashes'].items():
 p=D/fn if fn=='clash_gas_parameters.csv' else B/fn;assert sha(p)==h
params={r['system']:np.array([float(r[k]) for k in ['ne0','r0','alpha','re0','beta0','ne1','re1','beta1']]) for r in readcsv(D/'clash_gas_parameters.csv')}
rawsource=D/'primary/ClusterInfo.py';primary_checked=False
if rawsource.exists():
 assert sha(rawsource)==SP['sha256']
 for node in ast.parse(rawsource.read_text()).body:
  if isinstance(node,ast.Assign) and isinstance(node.value,ast.List) and any(isinstance(t,ast.Name) and t.id=='cluster_data' for t in node.targets):literal=ast.literal_eval(node.value);break
 assert len(literal)==len(params)==16
 for row in literal:
  c=np.array(list(map(float,row[2:10])));
  if row[0]=='MACS1720':c[5]=0
  assert np.array_equal(c,params[row[0]])
 primary_checked=True
source=np.array([r['source'] for r in rr]);host=np.array([r['system'] for r in rr]);radius=np.array([float(r['radius_kpc']) for r in rr]);bar=np.array([float(r['gbar_m_s2']) for r in rr]);target=np.array([float(r['gobs_m_s2']) for r in rr])
assert sum(source=='CLASH')==65 and sum(source=='XCOP')==21;assert len(set(host))==23
original=readcsv(B/'dm_unified_regular/augmented_input.csv');oldxc={(r['cluster'],float(r['radius_kpc'])):r for r in readcsv(B/'dm_fixed_radius/fixed_radius_inputs.csv')}
profrows=readcsv(D/'gas_profiles.csv');profiles={}
G=6.67430e-11;MS=1.988409870698051e30;KPC=3.085677581491367e19
for h in sorted(set(host[source=='XCOP'])):
 name=h.removeprefix('XCOP_');path=B/'dm_fit/data/xcop'/(name+'_fgas_profile.fits')
 with fits.open(path) as f:
  rs=np.asarray(f['FGAS'].data['RADIUS'],float)*f['FGAS'].header['R500'];ms=np.asarray(f['FGAS'].data['MGAS'],float)
 delivered=[r for r in profrows if r['system']==h];assert len(delivered)==len(rs)
 assert np.array_equal(rs,np.array([float(r['radius_kpc']) for r in delivered]));assert np.array_equal(ms,np.array([float(r['mgas_msun']) for r in delivered]));assert np.all(np.diff(rs)>0) and np.all(np.diff(ms)>0);profiles[h]=(rs,ms)
for i,row in enumerate(rr):
 old=original[int(row['broad_row'])];assert row['system']==old['system'];assert int(row['row'])==i
 for f in ['radius_kpc','gbar_m_s2','gobs_m_s2']:assert float(row[f])==float(old[f])
 if source[i]=='XCOP':
  oo=oldxc[host[i].removeprefix('XCOP_'),radius[i]]
  for f in ['gbar_no_stars','gbar_full_projected_stars','gobs_forward','gobs_nfw','gobs_einasto']:assert float(row[f])==float(oo[f])
 else:
  for f in ['gbar_no_stars','gbar_full_projected_stars']:assert float(row[f])==bar[i]
  for f in ['gobs_forward','gobs_nfw','gobs_einasto']:assert float(row[f])==target[i]
 q=Q[i];assert q['row']==i;s=np.array(q['radius_nodes_kpc']);w=np.array(q['mass_weights_msun']);assert np.all((s>0)&(s<radius[i]));assert np.all(w>0)
 assert np.isclose(w.sum(),float(row['mgas_newprofile_msun']),rtol=1e-13)
 assert np.isclose(w.sum()/(bar[i]*(radius[i]*KPC)**2/G/MS),float(row['mgas_over_old_mb']),rtol=1e-13)

def adaptive_mass(i,rc):
 r=radius[i]
 if source[i]=='CLASH':
  n0,r0,al,re,beta,n1,re1,beta1=params[host[i].removeprefix('CLASH_')]
  def integr(t):
   if t==0:return 0.
   s=r*t;dens=n0*1e5*(s/r0)**(-al)*(1+(s/re)**2)**(-1.5*beta)
   if n1>0:dens+=n1*1e5*(1+(s/re1)**2)**(-1.5*beta1)
   return 4*np.pi*r**3*t*t*dens*(1 if np.isinf(rc) else np.exp(-s/rc))
  return quad(integr,0,1,epsabs=1e-4,epsrel=2e-11,limit=250)[0]
 rs,ms=profiles[host[i]];r0=rs[0];m0=ms[0]
 out=quad(lambda t:3*m0*t*t*(1 if np.isinf(rc) else np.exp(-r0*t/rc)),0,1,epsabs=1e-4,epsrel=2e-11)[0]
 for lo,hi,ml,mh in zip(rs[:-1],rs[1:],ms[:-1],ms[1:]):
  if lo>=r:break
  upper=min(r,hi);power=np.log(mh/ml)/np.log(hi/lo);loglength=np.log(upper/lo)
  out+=quad(lambda v:power*ml*np.exp(power*v)*(1 if np.isinf(rc) else np.exp(-lo*np.exp(v)/rc)),0,loglength,epsabs=1e-4,epsrel=2e-11)[0]
  if hi>=r:break
 return out
quaderr=0.;quadchecks=0
for i in idx:
 q=Q[i];s=np.array(q['radius_nodes_kpc']);w=np.array(q['mass_weights_msun'])
 for rc in [10.,400.,1e4,np.inf]:
  direct=adaptive_mass(i,rc);approx=np.sum(w*(1 if np.isinf(rc) else np.exp(-s/rc)));err=abs(approx/direct-1);quaderr=max(quaderr,err);assert err<2e-8,(i,rc,err);quadchecks+=1

def case(cname):
 if cname=='combined_forward_mid':return idx,bar,np.log10(target)
 if cname=='CLASH_forward_mid':return idx[source=='CLASH'],bar,np.log10(target)
 _,tar,key=cname.split('_',2);return idx[source=='XCOP'],np.array([float(r[key]) for r in rr]),np.log10([float(r['gobs_'+tar]) for r in rr])

def weights(ii):
 w=np.zeros(len(ii));srcs=set(source[ii])
 for s in srcs:
  hosts=set(host[ii][source[ii]==s])
  for h in hosts:
   z=host[ii]==h;w[z]=1/(len(srcs)*len(hosts)*sum(z))
 assert abs(w.sum()-1)<1e-12;return w

def model(k,c,ii,bb,refined=False,inner_bound=None):
 nu,typ=k.split('_',1);a=1.2e-10;g=bb[ii].copy()
 if typ=='refit_a0':a=1e-10*10**c[0]
 if typ in ['repair','frozen_gas']:
  eta,rc=(10.,400.) if typ=='frozen_gas' else 10**np.array(c)
  extra=[]
  for i in ii:
   q=Q[i];s=np.array(q['check_nodes_kpc' if refined else 'radius_nodes_kpc']);w=np.array(q['check_weights_msun' if refined else 'mass_weights_msun']);m=np.dot(w,np.exp(-s/rc))
   if inner_bound is not None and source[i]=='XCOP':
    r0=float(rr[i]['inner_radius_kpc']);m0=float(rr[i]['inner_mass_msun']);m-=np.dot(w[s<r0],np.exp(-s[s<r0]/rc));m+=m0*(np.exp(-r0/rc) if inner_bound=='lower' else 1)
   extra.append(m)
  g+=eta*G*MS*np.array(extra)/(radius[ii]*KPC)**2
 if nu=='simple':v=g/2+np.sqrt((g/2)**2+a*g)
 else:v=g/(-np.expm1(-np.sqrt(g/a)))
 return np.log10(v)

def metric(pr,yy,ii):
 e=pr-yy[ii];return {'B_dex':float(np.sqrt(np.dot(weights(ii),e*e))),'median_acceleration_error_percent':float(100*np.median(abs(10**e-1))),'rows':len(ii),'systems':len(set(host[ii])),'sources':{s:{'RMS_dex':float(np.sqrt(np.mean(e[source[ii]==s]**2))),'median_acceleration_error_percent':float(100*np.median(abs(10**e[source[ii]==s]-1)))} for s in set(source[ii])}}

def checkmetric(pr,yy,ii,d):
 m=metric(pr,yy,ii)
 for k in ['B_dex','median_acceleration_error_percent','rows','systems']:assert abs(m[k]-d[k])<1e-8
 for s,v in m['sources'].items():
  for k,a in v.items():assert abs(a-d['sources'][s][k])<1e-8
attempts=collections.defaultdict(list);maxcost=0.
for f in J['fits']:
 bits=f['tag'].split(':');cname=bits[0];ii,bb,yy=case(cname)
 if bits[1]=='host':ii=ii[host[ii]!=bits[2]]
 if bits[1]=='source':ii=ii[source[ii]!=bits[2]]
 assert ii.tolist()==f['train_rows'];c=np.array(f['coefficients']);lo=np.array([-3,1] if f['model'].endswith('repair') else [-3]);hi=np.array([3,4] if f['model'].endswith('repair') else [3]);assert np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10)
 e=model(f['model'],c,ii,bb)-yy[ii];cost=.5*np.dot(weights(ii),e*e);error=abs(cost-f['cost']);maxcost=max(maxcost,error);assert error<1e-10
 attempts[f['tag'],f['model']].append(f)
assert len(attempts)==J['fit_calls']==460;assert len(J['fits'])==J['optimizer_attempts']==1610
selected={}
for key,aa in attempts.items():
 count=4 if key[1].endswith('repair') else 3;assert len(aa)==count;assert {a['start'] for a in aa}==set(range(count));selected[key]=min([a for a in aa if a['success']],key=lambda a:a['cost'])['coefficients']
tables=collections.defaultdict(list)
for v in P:tables[v['case'],v['model'],v['split']].append(v)
assert len(tables)==J['prediction_tables']==192;assert len(P)==J['prediction_rows']==6128
maxpred=0.;maxrefined=0.;maxcore=0.;gasrows=0
for (cname,k,split),tab in tables.items():
 ii,bb,yy=case(cname);d=J['cases'][cname][k]
 if split.startswith('source_'):ii=ii[source[ii]==split.removeprefix('source_')]
 assert len(ii)==len(tab);pp=np.log10([float(r['gpred_m_s2']) for r in tab]);calc=[]
 for j,(row,i) in enumerate(zip(tab,ii)):
  assert int(row['row'])==i and row['system']==host[i] and row['source']==source[i]
  assert np.isclose(float(row['gbar_m_s2']),bb[i],rtol=1e-13);assert np.isclose(float(row['gobs_m_s2']),10**yy[i],rtol=1e-13);assert float(row['radius_kpc'])==radius[i]
  c=json.loads(row['coefficients_json']);tag=cname+(':all' if split=='all' else ':host:'+host[i] if split=='host_holdout' else ':source:'+row['source'])
  assert row['fold']==('all' if split=='all' else host[i] if split=='host_holdout' else row['source'])
  if k.endswith('fixed') or k.endswith('frozen_gas'):assert c==[]
  else:assert c==selected[tag,k]
  recorded=d['coefficients'] if split=='all' else d['host_coefficients'][host[i]] if split=='host_holdout' else d['source_holdout'][row['source']]['coefficients'];assert c==recorded
  p=model(k,c,np.array([i]),bb)[0];calc.append(p);assert abs(float(row['log10_residual'])-(pp[j]-yy[i]))<1e-10
  if k.endswith('repair') or k.endswith('frozen_gas'):
   gasrows+=1;ref=model(k,c,np.array([i]),bb,refined=True)[0];maxrefined=max(maxrefined,abs(ref-p))
   lower=model(k,c,np.array([i]),bb,inner_bound='lower')[0];upper=model(k,c,np.array([i]),bb,inner_bound='upper')[0];assert lower<=p+1e-12 and p<=upper+1e-12;maxcore=max(maxcore,upper-lower)
 error=max(abs(np.array(calc)-pp));maxpred=max(maxpred,error);assert error<1e-10
 met=d['all'] if split=='all' else d['host_holdout'] if split=='host_holdout' else d['source_holdout'][split.removeprefix('source_')]['metrics'];checkmetric(pp,yy,ii,met)
assert abs(maxrefined-DI['maximum_refined_quadrature_change_dex'])<1e-12
assert abs(maxcore-DI['maximum_unresolved_core_bound_span_dex'])<1e-12
assert gasrows==DI['checks']==3064
report={'status':'passed','prediction_tables':len(tables),'prediction_rows':len(P),'fit_sets':len(attempts),'optimizer_attempts':len(J['fits']),'failed_optimizer_attempts':sum(not f['success'] for f in J['fits']),'independent_adaptive_integral_checks':quadchecks,'max_adaptive_integral_relative_discrepancy':quaderr,'max_prediction_error_dex':maxpred,'max_cost_error':maxcost,'gas_prediction_rows_refined_and_inner_bounds_checked':gasrows,'maximum_refined_quad_change_dex':maxrefined,'maximum_unresolved_core_bound_span_dex':maxcore,'primary_parameter_literal_match_performed':primary_checked,'raw_XCOP_direct_MGAS_matched':True,'all_case_splits_and_target_exclusions_checked':True,'result_hash':sha(D/'results.json'),'audit_script_hash':sha(Path(__file__)),'limitations':['Frozen literature gas parameters are CLASH-informed; repeating them with a source_CLASH label is not independent CLASH-source exclusion.','Fits add a phenomenological positive source component before MOND; neither its material identity nor gravity/memory derivation is established.','Gas profiles and reconstructed total accelerations share observational/modeling systematics; algebraic nonleakage is not statistical independence.','Checks do not establish superiority to full published cluster likelihoods.']}
(O/'cluster_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
