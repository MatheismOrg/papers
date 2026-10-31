"""Uniform gas-source extension compatibility test; see frozen PROTOCOL.txt."""
from pathlib import Path
import csv,json,hashlib,shutil
import numpy as np
from scipy.optimize import least_squares

R=Path(__file__).resolve().parent;B=R.parent
G=6.67430e-11;MSUN=1.988409870698051e30;KPC=3.085677581491367e19;A0=1.2e-10
SEEDS=[20261027,7331,48105]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def readcsv(p):return list(csv.DictReader(Path(p).open()))
def writecsv(p,rr):
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)

for src,dst in [(B/'dm_halo_benchmark/matched_components.csv',R/'full_disk_components.csv'),(B/'dm_cluster_repair/inputs.csv',R/'cluster_inputs.csv'),(B/'dm_cluster_repair/quadrature.json',R/'cluster_quadrature.json')]:shutil.copy2(src,dst)
full=readcsv(R/'full_disk_components.csv');rows=[]
for z in full:
 gg=float(z['gas_accel_over_aref'])*1e-10;gs=float(z['stellar_accel_over_aref'])*1e-10;frac=gg/(gg+gs)
 if gg>0 and frac>=.5 and float(z['radius_kpc'])<=30:
  rows.append({'row':len(rows),'origin_row':int(z['row']),'system':z['galaxy'],'domain':'disk','source':'SPARC','radius_kpc':float(z['radius_kpc']),'gbar_m_s2':float(z['gbar_m_s2']),'ggas_m_s2':gg,'gstar_m_s2':gs,'gas_force_fraction':frac,'gobs_m_s2':float(z['gobs_m_s2']),'cluster_quadrature_row':-1})
assert len(rows)==402 and len({r['system'] for r in rows})==61
cl=readcsv(R/'cluster_inputs.csv');quadr=json.loads((R/'cluster_quadrature.json').read_text())
for z in cl:rows.append({'row':len(rows),'origin_row':int(z['broad_row']),'system':z['system'],'domain':'cluster','source':z['source'],'radius_kpc':float(z['radius_kpc']),'gbar_m_s2':float(z['gbar_m_s2']),'ggas_m_s2':0.,'gstar_m_s2':0.,'gas_force_fraction':'','gobs_m_s2':float(z['gobs_m_s2']),'cluster_quadrature_row':int(z['row'])})
writecsv(R/'inputs.csv',rows)
n=len(rows);I=np.arange(n);domain=np.array([z['domain'] for z in rows]);source=np.array([z['source'] for z in rows]);system=np.array([z['system'] for z in rows]);rad=np.array([z['radius_kpc'] for z in rows]);bar=np.array([z['gbar_m_s2'] for z in rows]);gas=np.array([z['ggas_m_s2'] for z in rows]);target=np.log10([z['gobs_m_s2'] for z in rows]);DI=I[domain=='disk'];CI=I[domain=='cluster'];maxnodes=max(len(z['radius_nodes_kpc']) for z in quadr)
nodes=np.zeros((n,maxnodes));massweights=np.zeros_like(nodes)
for j,z in enumerate(rows):
 if z['domain']=='cluster':
  q=quadr[z['cluster_quadrature_row']];assert q['row']==z['cluster_quadrature_row'];L=len(q['radius_nodes_kpc']);nodes[j,:L]=q['radius_nodes_kpc'];massweights[j,:L]=q['mass_weights_msun']
prev=json.loads((B/'dm_cluster_repair/results.json').read_text())['cases']['combined_forward_mid']
FROZEN={nu:[10**prev[key]['source_holdout']['XCOP']['coefficients'][0],prev[key]['source_holdout']['XCOP']['coefficients'][1]] for nu,key in [('simple','simple_repair'),('exp','exp_repair')]}

def weights(idx):
 w=np.zeros(len(idx));dd=np.unique(domain[idx])
 for d in dd:
  ss=np.unique(source[idx][domain[idx]==d])
  for so in ss:
   gg=np.unique(system[idx][source[idx]==so])
   for g in gg:
    mask=system[idx]==g;w[mask]=1/(len(dd)*len(ss)*len(gg)*sum(mask))
 assert np.isclose(sum(w),1)
 return w

def logmond(src,nu):
 assert np.all(src>0)
 pred=.5*src+np.sqrt(.25*src*src+A0*src) if nu=='simple' else src/(-np.expm1(-np.sqrt(src/A0)))
 return np.log10(pred)

def predict(nu,c,idx):
 eta,logrc=c;rc=10**logrc;K=gas[idx].copy();iscl=domain[idx]=='cluster'
 if np.any(iscl):
  q=idx[iscl];K[iscl]=G*MSUN*np.sum(massweights[q]*np.exp(-nodes[q]/rc),axis=1)/(rad[q]*KPC)**2
 return logmond(bar[idx]+eta*K,nu)

def metrics(pr,idx):
 e=pr-target[idx];out={'B_dex':float(np.sqrt(np.sum(weights(idx)*e*e))),'rows':len(idx),'systems':len(set(system[idx])),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'median_speed_error_percent':float(np.median(abs(10**(e/2)-1))*100),'domains':{},'sources':{}}
 for name,vals in [('domains',domain),('sources',source)]:
  for key in np.unique(vals[idx]):
   z=vals[idx]==key;ii=idx[z];ee=e[z];out[name][key]={'B_dex':float(np.sqrt(np.sum(weights(ii)*ee*ee))),'median_acceleration_error_percent':float(np.median(abs(10**ee-1))*100),'median_speed_error_percent':float(np.median(abs(10**(ee/2)-1))*100),'rows':len(ii),'systems':len(set(system[ii]))}
 return out

FOLDS={}
for seed in SEEDS:
 rng=np.random.default_rng(seed);mp={}
 for so in sorted(set(source)):
  gs=np.unique(system[source==so])
  for j,g in enumerate(rng.permutation(gs)):mp[g]=j%5
 FOLDS[seed]=np.array([mp[g] for g in system])

FITLOG=[];TABLE=[]
def fit(nu,idx,tag):
 diskonly=not np.any(domain[idx]=='cluster');starts=[[0.],[1.],[10.]] if diskonly else [[.1,2.],[1.,2.6],[10.,2.6],[100.,1.5]];lo=[0.] if diskonly else [0.,1.];hi=[1000.] if diskonly else [1000.,4.];rr=[];ids=[]
 for j,start in enumerate(starts):
  rrj=least_squares(lambda p:np.sqrt(weights(idx))*(predict(nu,[p[0],np.log10(400.)] if diskonly else p,idx)-target[idx]),start,bounds=(lo,hi),max_nfev=2000,ftol=1e-11,xtol=1e-11,gtol=1e-11)
  c=[float(rrj.x[0]),float(np.log10(400.))] if diskonly else rrj.x.tolist();fid=len(FITLOG);ids.append(fid)
  FITLOG.append({'fit_id':fid,'tag':tag,'nu':nu,'start_index':j,'initial':start,'bounds':[lo,hi],'disk_only':diskonly,'parameter_names':['eta'] if diskonly else ['eta','log10_Rc_kpc'],'coefficients':c,'fitted_parameters':rrj.x.tolist(),'active_bounds':rrj.active_mask.tolist(),'train_rows':idx.tolist(),'success':bool(rrj.success),'message':rrj.message,'nfev':int(rrj.nfev),'cost':float(rrj.cost)})
  rr.append(rrj)
 ok=[j for j,v in enumerate(rr) if v.success];assert ok
 k=min(ok,key=lambda j:rr[j].cost);return FITLOG[ids[k]]['coefficients'],ids[k]

def save(scope,nu,model,split,pr,idx,fold,cfs,fids):
 for j,i in enumerate(idx):TABLE.append({'scope':scope,'nu':nu,'model':model,'split':split,'row':int(i),'system':system[i],'source':source[i],'domain':domain[i],'fold':str(fold[j]),'selected_fit_id':fids[j],'coefficients_json':json.dumps(cfs[j]),'radius_kpc':rad[i],'gbar_m_s2':bar[i],'gobs_m_s2':10**target[i],'gpred_m_s2':10**pr[j],'log10_residual':pr[j]-target[i]})

def run():
 out={'status':'Exploratory uniform-weight disk gas extension; not exact finite-cutoff disk field and not the published hot-only claim.','protocol_sha256':sha(R/'PROTOCOL.txt'),'script_sha256':sha(Path(__file__)),'source_hashes':{str(p.relative_to(B)):sha(p) for p in [B/'dm_halo_benchmark/matched_components.csv',B/'dm_cluster_repair/inputs.csv',B/'dm_cluster_repair/quadrature.json',B/'dm_cluster_repair/results.json']},'input_sha256':sha(R/'inputs.csv'),'frozen_CLASH_parameters':FROZEN,'counts':{'rows':n,'systems':len(set(system)),'disk_rows':len(DI),'disk_systems':len(set(system[DI])),'cluster_rows':len(CI),'cluster_systems':len(set(system[CI]))},'scopes':{},'fits':FITLOG,'transfers':{}}
 predcache={}
 for scope,idx in [('disk',DI),('cluster',CI),('joint',I)]:
  sd={};out['scopes'][scope]=sd
  for nu in ['simple','exp']:
   nd={};sd[nu]=nd
   for model,cf in [('zero',[0.,np.log10(400.)]),('literature',[10.,np.log10(400.)]),('CLASH_frozen',FROZEN[nu]),('fitted',None)]:
    cc,fid=fit(nu,idx,f'{scope}:{nu}:all') if cf is None else (cf,-1);pr=predict(nu,cc,idx);v={'coefficients':cc,'selected_fit_id':fid,'all':metrics(pr,idx),'CV':{}};nd[model]=v;save(scope,nu,model,'all',pr,idx,['all']*len(idx),[cc]*len(idx),[fid]*len(idx))
    for seed in SEEDS:
     fd=FOLDS[seed][idx];pp=np.full(len(idx),np.nan);cfs=[None]*len(idx);fids=[-1]*len(idx);fc={}
     for f in range(5):
      tr=idx[fd!=f];te=idx[fd==f];assert not set(system[tr])&set(system[te]);cc,fid=fit(nu,tr,f'{scope}:{nu}:cv_{seed}:{f}') if cf is None else (cf,-1)
      pos=np.where(fd==f)[0];pp[pos]=predict(nu,cc,te);fc[str(f)]={'coefficients':cc,'selected_fit_id':fid}
      for j in pos:cfs[j]=cc;fids[j]=fid
     assert np.all(np.isfinite(pp));v['CV'][str(seed)]={'metrics':metrics(pp,idx),'fold_coefficients':fc};save(scope,nu,model,f'cv_{seed}',pp,idx,fd,cfs,fids)
     if seed==SEEDS[0]:predcache[scope,nu,model]=pp
    print(scope,nu,model,'eta,Rc',v['coefficients'][0],10**v['coefficients'][1],'CV',v['CV'][str(SEEDS[0])]['metrics']['B_dex'],flush=True)
   if scope=='joint':
    for name,te in [('source_CLASH',I[source=='CLASH']),('source_XCOP',I[source=='XCOP']),('domain_disk',DI),('domain_cluster',CI)]:
     tr=I[~np.isin(system,system[te])];assert not set(system[tr])&set(system[te]);cc,fid=fit(nu,tr,f'joint:{nu}:{name}');pr=predict(nu,cc,te);out['transfers'][nu+':'+name]={'coefficients':cc,'selected_fit_id':fid,'train_rows':tr.tolist(),'test_rows':te.tolist(),'metrics':metrics(pr,te),'Rc_externally_fixed_400kpc':bool(name=='domain_cluster')};save('joint',nu,'fitted',name,pr,te,[name]*len(te),[cc]*len(te),[fid]*len(te))
  (R/'results.json').write_text(json.dumps(out,indent=2))
 # Frozen-prediction descriptive paired group bootstrap.
 bs={};rng=np.random.default_rng(20261029)
 def bootstrap_scores(pp,idx,draws):
  ee=(pp-target[idx])**2;dd=[]
  for d in sorted(set(domain[idx])):
   ss=[]
   for so in sorted(set(source[idx][domain[idx]==d])):
    gg=np.unique(system[idx][source[idx]==so]);means=np.array([np.mean(ee[system[idx]==g]) for g in gg]);ss.append(np.mean(means[draws[d,so]],axis=1))
   dd.append(np.mean(ss,axis=0))
  return np.sqrt(np.mean(dd,axis=0))
 for nu in ['simple','exp']:
  cases=[('disk_CLASH_frozen_minus_zero',DI,predcache['disk',nu,'CLASH_frozen'],predcache['disk',nu,'zero']),('disk_literature_minus_zero',DI,predcache['disk',nu,'literature'],predcache['disk',nu,'zero']),('disk_joint_minus_separate',DI,predcache['joint',nu,'fitted'][DI],predcache['disk',nu,'fitted']),('cluster_joint_minus_separate',CI,predcache['joint',nu,'fitted'][CI],predcache['cluster',nu,'fitted']),('joint_fitted_minus_zero',I,predcache['joint',nu,'fitted'],predcache['joint',nu,'zero'])]
  for name,idx,pa,pb in cases:
   draws={}
   for d in sorted(set(domain[idx])):
    for so in sorted(set(source[idx][domain[idx]==d])):
     ng=len(set(system[idx][source[idx]==so]));draws[d,so]=rng.integers(0,ng,size=(5000,ng))
   diff=bootstrap_scores(pa,idx,draws)-bootstrap_scores(pb,idx,draws)
   bs[nu+':'+name]={'difference_dex':metrics(pa,idx)['B_dex']-metrics(pb,idx)['B_dex'],'percentile_95_dex':np.percentile(diff,[2.5,97.5]).tolist()}
 out['descriptive_bootstrap']=bs
 # All-disk extrapolation diagnostics: preserve invalid source forces explicitly.
 fr=[];fm={}
 for nu in ['simple','exp']:
  for name,c in [('zero',[0.,np.log10(400.)]),('literature',[10.,np.log10(400.)]),('CLASH_frozen',FROZEN[nu])]+[(sc+'_all_fit',out['scopes'][sc][nu]['fitted']['coefficients']) for sc in ['disk','cluster','joint']]:
   valid=[];invalid=[]
   for z in full:
    gg=float(z['gas_accel_over_aref'])*1e-10;gs=float(z['stellar_accel_over_aref'])*1e-10;src=gs+(1+c[0])*gg;good=src>0;pred=float(10**logmond(np.array([src]),nu)[0]) if good else '';res=np.log10(pred/float(z['gobs_m_s2'])) if good else ''
    rr={'nu':nu,'model':name,'origin_row':int(z['row']),'system':z['galaxy'],'radius_kpc':float(z['radius_kpc']),'eta':c[0],'gbar_m_s2':float(z['gbar_m_s2']),'ggas_m_s2':gg,'gstar_m_s2':gs,'gsource_m_s2':src,'valid_positive_source':good,'gobs_m_s2':float(z['gobs_m_s2']),'gpred_m_s2':pred,'log10_residual':res};fr.append(rr);(valid if good else invalid).append(rr)
   fm[nu+':'+name]={'total_rows':len(full),'invalid_rows':len(invalid),'invalid_systems':sorted({z['system'] for z in invalid}),'invalid_origin_rows':[z['origin_row'] for z in invalid],'valid_conditional_row_RMS_dex':float(np.sqrt(np.mean([z['log10_residual']**2 for z in valid]))),'valid_conditional_median_acceleration_error_percent':float(np.median([abs(10**z['log10_residual']-1) for z in valid])*100),'not_a_whole_sample_score_when_invalid':bool(invalid)}
 out['full_disk_diagnostics']=fm
 writecsv(R/'predictions.csv',TABLE);writecsv(R/'full_disk_predictions.csv',fr)
 out['fit_sets']=sum(z['start_index']==0 for z in FITLOG);out['optimizer_attempts']=len(FITLOG);out['optimizer_failures']=[z for z in FITLOG if not z['success']];out['prediction_rows']=len(TABLE);out['prediction_tables']=len({(z['scope'],z['nu'],z['model'],z['split']) for z in TABLE});out['full_disk_diagnostic_rows']=len(fr)
 (R/'results.json').write_text(json.dumps(out,indent=2));print('DONE',out['fit_sets'],out['optimizer_attempts'],out['prediction_tables'],len(TABLE),'failures',len(out['optimizer_failures']),flush=True)

if __name__=='__main__':run()
