"""Independent NFW/Burkert and radial-split audit; never imports run.py."""
from pathlib import Path
import csv,json,hashlib,io,zipfile,collections
import numpy as np
from scipy.integrate import quad
B=Path(__file__).resolve().parent.parent;D=B/'dm_halo_benchmark';O=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
readcsv=lambda p:list(csv.DictReader(p.open()))
J=json.loads((D/'results/results.json').read_text()); rr=readcsv(D/'augmented_input.csv');N=len(rr);idx=np.arange(N)
assert J['input_sha256']==sha(D/'augmented_input.csv')==sha(B/'dm_unified_regular/augmented_input.csv')
assert J['protocol_sha256']==sha(D/'PROTOCOL.txt');assert J['script_sha256']==sha(D/'run.py')
raw=B/'dm_fit/data/Rotmod_LTG.zip';assert J['raw_rotmod_sha256']==sha(raw)
assert J['matched_components_sha256']==sha(D/'matched_components.csv')
pop=np.array([v['population'] for v in rr]);gal=np.array([v['system'] for v in rr]);groups=np.array([v['group'] for v in rr]);r=np.array([float(v['radius_kpc']) for v in rr]);b=np.array([float(v['gbar_m_s2'])/1e-10 for v in rr]);y=np.log10([float(v['gobs_m_s2'])/1e-10 for v in rr]);sg=sorted(set(gal[pop=='SPARC']));eligible=[g for g in sg if sum((pop=='SPARC')&(gal==g))>=9]
assert eligible==J['eligibility']['galaxies'];assert len(eligible)==104
components=readcsv(D/'matched_components.csv');assert len(components)==2700
gas=np.full(N,np.nan);star=np.full(N,np.nan);lamlow={};rawerr=0.
with zipfile.ZipFile(raw) as z:
 for v in components:
  i=int(v['row']);assert v['galaxy']==gal[i];assert float(v['radius_kpc'])==r[i]
  dat=np.loadtxt(io.BytesIO(z.read(gal[i]+'_rotmod.dat')),ndmin=2);line=dat[np.isclose(dat[:,0],r[i],atol=1e-9,rtol=0)];assert len(line)==1
  a=line[0];f=1e16/(3.085677581491367e19*a[0]);gg=a[3]*abs(a[3])*f;ss=(.5*a[4]*abs(a[4])+.7*a[5]*abs(a[5]))*f
  gas[i]=gg;star[i]=ss;assert np.isclose(gg,float(v['gas_accel_over_aref']),rtol=1e-12,atol=1e-16);assert np.isclose(ss,float(v['stellar_accel_over_aref']),rtol=1e-12,atol=1e-16)
  assert np.isclose(gg+ss,b[i],rtol=1e-12);assert np.isclose(a[1]**2*f,10**y[i],rtol=1e-12);assert float(v['vobs_km_s'])==a[1];assert float(v['errv_km_s'])==a[2]
  lamlow[gal[i]]=float(v['stellar_multiplier_lower_bound'])
for g in sg:
 ii=idx[(pop=='SPARC')&(gal==g)];ss=star[ii];gg=gas[ii];valid=ss>0
 low=max(.2,max(-gg[valid]/ss[valid],default=-np.inf)+1e-6)
 assert abs(low-lamlow[g])<1e-12
# Independent halo mass integral: Gauss-Legendre quadrature for cancellation-prone x.
node,weight=np.polynomial.legendre.leggauss(12);node=(node+1)/2;weight=weight/2

def mass(kind,xx):
 z=np.atleast_1d(np.array(xx,dtype=np.longdouble));out=np.empty(len(z),dtype=np.longdouble);small=z<.03
 if np.any(small):
  zz=z[small,None];t=zz*node
  density=t/(1+t)**2 if kind=='NFW' else t*t/((1+t)*(1+t*t))
  out[small]=zz[:,0]*np.sum(density*weight,axis=1)
 q=z[~small]
 if kind=='NFW':out[~small]=np.log1p(q)-q/(1+q)
 else:out[~small]=(2*np.log1p(q)+np.log1p(q*q)-2*np.arctan(q))/4
 return np.array(out,dtype=float).reshape(np.shape(xx))
quaderr=0.
for kind in ['NFW','Burkert']:
 for xx in np.geomspace(1e-5,1e4,45):
  integ=quad(lambda t:(t/(1+t)**2 if kind=='NFW' else t*t/((1+t)*(1+t*t))),0,xx,epsabs=1e-28,epsrel=2e-12,limit=250)[0]
  rel=abs(float(mass(kind,xx))/integ-1);quaderr=max(quaderr,rel);assert rel<1e-10,(kind,xx,rel)
 assert abs(float(mass(kind,1))/(1*float(mass(kind,1)))-1)<1e-15

def rar(bb,a):return np.log10(bb/(-np.expm1(-np.sqrt(bb/a))))
def predict(k,c,ii):
 bb=b[ii];rad=r[ii]
 if k=='RAR_fixed':return rar(bb,1.2)
 if k in ['RAR_shared','RAR_individual']:return rar(bb,10**c[0])
 if k=='RAR_stellar':return rar(gas[ii]+10**c[0]*star[ii],1.2)
 if k.startswith('mixed_regular'):
  U=bb*rad**2;s=bb*rad;S=1/(1+(10**c[5]/rad)**4);W=1/(1+(bb/100)**2)
  a=10**c[0]+10**c[1]*U**(-c[3])*S*W+10**c[2]*s**c[4]/(1+rad/10**c[6]);return rar(bb,a)
 kind=k.split('_')[0];xx=rad/10**c[1];dm=10**c[0]*mass(kind,xx)/(xx*xx*mass(kind,1.))
 if k.endswith('_stellar'):bb=gas[ii]+10**c[2]*star[ii]
 assert np.all(bb>0);return np.log10(bb+dm)

def metric(pr,ii):
 e=pr-y[ii];gm=[np.mean(e[gal[ii]==g]**2) for g in sorted(set(gal[ii]))]
 return {'galaxy_RMS_dex':float(np.sqrt(np.mean(gm))),'row_RMS_dex':float(np.sqrt(np.mean(e*e))),'median_acceleration_error_percent':float(100*np.median(abs(10**e-1))),'median_speed_error_percent':float(100*np.median(abs(10**(e/2)-1))),'rows':len(ii),'galaxies':len(gm)}

def checkmetric(pr,ii,d):
 for k,v in metric(pr,ii).items():assert abs(v-d[k])<1e-8,(k,v,d[k])

def weights(ii,typ):
 if typ=='equal_point':return np.ones(len(ii))/len(ii)
 ps=sorted(set(pop[ii]));w=np.zeros(len(ii))
 for p in ps:
  gg=sorted(set(groups[ii][pop[ii]==p]))
  for g in gg:
   sel=(pop[ii]==p)&(groups[ii]==g);w[sel]=1/(len(ps)*len(gg)*sum(sel))
 assert abs(w.sum()-1)<1e-12
 return w
logs=J['fit_log'];sel=set(J['selected_ids']);assert len(logs)==3762;assert len(sel)==1254
maxcosterr=0.;allattempts=collections.defaultdict(list);expected_splits={}
for split,dd in J['splits'].items():
 test=[];gtrain={};gtest={}
 for g in eligible:
  ii=idx[(pop=='SPARC')&(gal==g)];ii=ii[np.argsort(r[ii],kind='stable')];k=len(ii)//3
  te=ii[-k:] if split=='outer_third' else ii[k:2*k];tr=ii[~np.isin(ii,te)];test.extend(te.tolist());gtrain[g]=tr;gtest[g]=te
  assert tr.tolist()==dd['galaxy_train_rows'][g];assert te.tolist()==dd['galaxy_test_rows'][g]
 te=np.array(sorted(test));tr=idx[~np.isin(idx,te)];disk=tr[pop[tr]=='SPARC'];assert len(te)==788
 assert te.tolist()==dd['test_rows'];assert tr.tolist()==dd['broad_train_rows'];assert disk.tolist()==dd['disk_train_rows']
 expected_splits[split]=(te,tr,disk,gtrain,gtest)
for entry in logs:
 i=entry['id'];assert logs[i] is entry;k=entry['model'];g=entry['galaxy'];te,tr,disk,gt,ge=expected_splits[entry['split']]
 ii=tr if k=='mixed_regular_broad' else disk if g=='shared' else gt[g]
 assert ii.tolist()==entry['train_rows'];assert not set(ii)&set(te)
 assert entry['weighting']==('equal_population_galaxy_point' if k=='mixed_regular_broad' else 'equal_galaxy_point' if g=='shared' else 'equal_point')
 c=np.array(entry['coefficients']);lo,hi=map(np.array,entry['bounds']);assert np.all(c>=lo-1e-11) and np.all(c<=hi+1e-11)
 error=predict(k,c,ii)-y[ii];cost=.5*np.dot(weights(ii,entry['weighting']),error*error);dif=abs(cost-entry['cost']);maxcosterr=max(maxcosterr,dif);assert dif<1e-9,(entry['id'],dif)
 allattemps=None;allattempts[entry['tag']].append(entry)
assert len(allattempts)==1254
for aa in allattempts.values():
 assert len(aa)==3;assert {a['start'] for a in aa}=={0,1,2}
 best=min(aa,key=lambda a:a['cost']);assert best['id'] in sel
maxerr=0.;tables=0;predrows=0;bootcount=0
for split,dd in J['splits'].items():
 te,tr,disk,gt,ge=expected_splits[split];gg=sorted(set(gal[te]));resamples=np.random.default_rng(20261028).integers(0,len(gg),size=(5000,len(gg)));boots={};points={}
 for k,d in dd['models'].items():
  pp=readcsv(D/'results'/f'{split}_{k}.csv');assert len(pp)==788;lp=np.log10([float(v['gpred_m_s2'])/1e-10 for v in pp]);ex=np.empty(len(te));seen=set();trainii=[];trainpr=[]
  for j,(v,i) in enumerate(zip(pp,te)):
   assert int(v['row'])==i;assert v['system']==gal[i];assert v['group']==groups[i];assert v['split']==split;assert v['model']==k
   for f,expected in [('radius_kpc',r[i]),('gbar_m_s2',b[i]*1e-10),('gobs_m_s2',10**y[i]*1e-10)]:assert np.isclose(float(v[f]),expected,rtol=1e-12,atol=0)
   fid=int(v['selected_fit_id']);seen.add(fid)
   if k=='RAR_fixed':assert fid==-1;c=[]
   else:
    assert fid in sel;a=logs[fid];assert a['split']==split and a['model']==k;assert a['galaxy'] in ['shared',gal[i]];c=a['coefficients'];label='shared' if a['galaxy']=='shared' else gal[i];assert c==d['coefficients'][label]
   ex[j]=predict(k,c,np.array([i]))[0]
   assert abs(float(v['log10_residual'])-(lp[j]-y[i]))<1e-10
  assert sorted(seen)==d['selected_fit_ids'];err=float(max(abs(ex-lp)));maxerr=max(maxerr,err);assert err<1e-9,(split,k,err)
  checkmetric(lp,te,d['test'])
  for g in eligible:
   ii=gt[g];c=[] if k=='RAR_fixed' else d['coefficients']['shared' if 'shared' in d['coefficients'] else g];trainii.extend(ii.tolist());trainpr.extend(predict(k,c,ii).tolist());mask=gal[te]==g;checkmetric(lp[mask],te[mask],d['by_galaxy'][g])
  checkmetric(np.array(trainpr),np.array(trainii),d['train_eligible_radii'])
  e=lp-y[te];eg=np.array([np.mean(e[gal[te]==g]**2) for g in gg]);boots[k]=np.sqrt(np.mean(eg[resamples],axis=1));points[k]=metric(lp,te)['galaxy_RMS_dex'];tables+=1;predrows+=len(te)
 for key,bd in dd['paired_galaxy_bootstrap'].items():
  k,base=key.split(' minus ');delta=boots[k]-boots[base];assert abs((points[k]-points[base])-bd['difference_dex'])<1e-12;assert max(abs(np.percentile(delta,[2.5,97.5])-bd['percentile_95_dex']))<1e-12;bootcount+=1
report={'status':'passed','prediction_tables':tables,'prediction_rows':predrows,'fit_sets':len(sel),'optimizer_attempts':len(logs),'failed_optimizer_attempts':sum(not a['success'] for a in logs),'max_independent_prediction_error_dex':maxerr,'max_independent_cost_error':maxcosterr,'raw_sparc_rows_independently_matched':len(components),'independent_halo_integral_relative_error':quaderr,'descriptive_bootstrap_comparisons_checked':bootcount,'all_contiguous_splits_and_training_targets_checked':True,'result_hash':sha(D/'results/results.json'),'audit_script_hash':sha(Path(__file__)),'limitations':'Conditional radial prediction within measured galaxies, not prediction of entirely new systems. Counts/metrics do not equal a cosmological likelihood or best published nuisance-marginalized analysis.'}
(O/'halo_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
