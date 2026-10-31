"""Frozen-protocol cluster gas-density repair. See PROTOCOL.txt and README.txt.
The derived CLASH parameter CSV is sufficient: no external source code is run.
"""
from pathlib import Path
import csv,json,hashlib,time,sys
import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss

R=Path(__file__).resolve().parent
sys.path.insert(0,str(R.parent/'dm_fit/deps'))
from astropy.io import fits
G=6.67430e-11;MSUN=1.988409870698051e30;KPC=3.085677581491367e19
A0=1.2e-10
SRC=R.parent/'dm_unified_regular/augmented_input.csv'
FITLOG=[]

def writecsv(path,rr):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)

def literal_parameters():
    keys=['ne0','r0','alpha','re0','beta0','ne1','re1','beta1']
    return {z['system']:[float(z[k]) for k in keys] for z in csv.DictReader((R/'clash_gas_parameters.csv').open())}

def rho_clash(r,c):
    n0,r0,al,re,beta,n1,re1,beta1=c
    z=n0*1e5*np.exp(-al*np.log(r/r0)-1.5*beta*np.log1p((r/re)**2))
    if n1>0:z+=n1*1e5*np.exp(-1.5*beta1*np.log1p((r/re1)**2))
    return z

def nodes_clash(r,c,order=96):
    q,w=leggauss(order);t=(q+1)/2;w=w/2;s=r*t*t
    return s,w*(2*r*t)*4*np.pi*s*s*rho_clash(s,c)

def nodes_xcop(r,rs,ms,order=12,inner_order=64):
    # Piecewise power-law enclosed mass; constant-density unresolved core.
    q,w=leggauss(inner_order);t=(q+1)/2;w=w/2;r0=rs[0];m0=ms[0]
    ss=[r0*t*t];ww=[w*6*m0*t**5]
    q,w=leggauss(order)
    for lo,hi,ml,mh in zip(rs[:-1],rs[1:],ms[:-1],ms[1:]):
        if lo>=r:break
        upper=min(r,hi);p=np.log(mh/ml)/np.log(hi/lo)
        a,b=np.log(lo),np.log(upper);z=(q+1)*(b-a)/2+a;s=np.exp(z)
        ss.append(s);ww.append(w*(b-a)/2*p*ml*(s/lo)**p)
        if hi>=r:break
    return np.concatenate(ss),np.concatenate(ww)

def prepare():
    params=literal_parameters();rr=list(csv.DictReader(SRC.open()));oldxc={}
    for x in csv.DictReader((R.parent/'dm_fixed_radius/fixed_radius_inputs.csv').open()):
        oldxc[(x['cluster'],float(x['radius_kpc']))]={k:(v if k=='cluster' else float(v)) for k,v in x.items()}
    rows=[];quadr=[];profiles=[];checks=[];skipped=[];sourcefiles={str(SRC.relative_to(R.parent)):SRC}
    sourcefiles['clash_gas_parameters.csv']=R/'clash_gas_parameters.csv'
    for i,z in enumerate(rr):
        if z['source'] not in ['Tian2020_CDS_fig2','XCOP']:continue
        name=z['system'].replace('CLASH_','').replace('XCOP_','');r=float(z['radius_kpc']);gb=float(z['gbar_m_s2']);go=float(z['gobs_m_s2'])
        if z['source']=='Tian2020_CDS_fig2':
            if name not in params:skipped.append({'system':z['system'],'row':i,'reason':'No gas-profile parameters in published16-host archive'});continue
            s,w=nodes_clash(r,params[name]);s2,w2=nodes_clash(r,params[name],192)
            checks.append({'system':z['system'],'r_kpc':r,'relative_mass_quadrature_change':float(abs(w2.sum()/w.sum()-1))})
            src='CLASH';rmin=0.;mmin=0.;old={}
        else:
            src='XCOP';path=R.parent/'dm_fit/data/xcop'/(name+'_fgas_profile.fits');sourcefiles[str(path.relative_to(R.parent))]=path
            with fits.open(path) as f:
                rs=np.asarray(f['FGAS'].data['RADIUS'],float)*f['FGAS'].header['R500'];ms=np.asarray(f['FGAS'].data['MGAS'],float)
            assert np.all(np.diff(rs)>0) and np.all(np.diff(ms)>0)
            assert rs[0]<r<rs[-1]
            s,w=nodes_xcop(r,rs,ms);s2,w2=nodes_xcop(r,rs,ms,24,128)
            checks.append({'system':z['system'],'r_kpc':r,'relative_mass_quadrature_change':float(abs(w2.sum()/w.sum()-1))})
            if not any(a['system']==z['system'] for a in profiles):
                profiles.extend({'system':z['system'],'radius_kpc':float(a),'mgas_msun':float(b)} for a,b in zip(rs,ms))
            rmin=float(rs[0]);mmin=float(ms[0]);old=oldxc[name,r]
        mb=gb*(r*KPC)**2/G/MSUN
        row={'row':len(rows),'broad_row':i,'system':z['system'],'source':src,'radius_kpc':r,'gbar_m_s2':gb,'gobs_m_s2':go,'mgas_newprofile_msun':float(w.sum()),'mgas_over_old_mb':float(w.sum()/mb),'inner_radius_kpc':rmin,'inner_mass_msun':mmin}
        for key in ['gbar_no_stars','gbar_full_projected_stars','gobs_forward','gobs_nfw','gobs_einasto']:row[key]=old.get(key,gb if key.startswith('gbar') else go)
        rows.append(row);quadr.append({'row':row['row'],'radius_nodes_kpc':s.tolist(),'mass_weights_msun':w.tolist(),'check_nodes_kpc':s2.tolist(),'check_weights_msun':w2.tolist()})
    writecsv(R/'inputs.csv',rows);writecsv(R/'gas_profiles.csv',profiles)
    writecsv(R/'clash_gas_parameters.csv',[dict(system=k,ne0=c[0],r0=c[1],alpha=c[2],re0=c[3],beta0=c[4],ne1=c[5],re1=c[6],beta1=c[7]) for k,c in params.items()])
    (R/'quadrature.json').write_text(json.dumps(quadr))
    prov={'source_urls':{'CLASH_gas':'https://zenodo.org/records/15299349/files/ClusterInfo.py?download=1','paper':'https://arxiv.org/html/2410.02612v2','published_equation':17},'source_hashes':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in sourcefiles.items()},'missing_gas_exclusions':skipped,'quadrature_checks':checks,'counts':{s:{'rows':sum(r['source']==s for r in rows),'systems':len({r['system'] for r in rows if r['source']==s})} for s in ['CLASH','XCOP']},'gas_old_baryon_consistency_CLASH':{'median_ratio':float(np.median([r['mgas_over_old_mb'] for r in rows if r['source']=='CLASH'])),'max_ratio':max(r['mgas_over_old_mb'] for r in rows if r['source']=='CLASH'),'rows_gas_exceeds_old_total_baryons':sum(r['mgas_over_old_mb']>1 for r in rows if r['source']=='CLASH')}}
    prov['CLASH_original_archive_provenance']=json.loads((R/'source_provenance.json').read_text())
    (R/'provenance.json').write_text(json.dumps(prov,indent=2))
    return rows,quadr,prov

def logmond(g,a,kind):
    if kind=='simple':p=.5*g+np.sqrt(.25*g*g+a*g)
    else:p=g/(-np.expm1(-np.sqrt(g/a)))
    return np.log10(p)

def main():
    global rows,quadr,sysnames,source,rad,bar,target,nodes,mweights,I
    rows,quadr,prov=prepare();n=len(rows);I=np.arange(n)
    sysnames=np.array([r['system'] for r in rows]);source=np.array([r['source'] for r in rows]);rad=np.array([r['radius_kpc'] for r in rows]);bar=np.array([r['gbar_m_s2'] for r in rows]);target=np.array([r['gobs_m_s2'] for r in rows]);maxn=max(len(v['mass_weights_msun']) for v in quadr);nodes=np.zeros((n,maxn));mweights=nodes.copy()
    for j,v in enumerate(quadr):nn=len(v['mass_weights_msun']);nodes[j,:nn]=v['radius_nodes_kpc'];mweights[j,:nn]=v['mass_weights_msun']
    allpred=[];result={'protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((R/'inputs.csv').read_bytes()).hexdigest(),'quadrature_sha256':hashlib.sha256((R/'quadrature.json').read_bytes()).hexdigest(),'provenance':prov,'fits':FITLOG,'cases':{}}
    models=[f'{nu}_{mo}' for nu in ['simple','exp'] for mo in ['fixed','frozen_gas','refit_a0','repair']]
    def weights(idx):
        out=np.zeros(len(idx));srcs=np.unique(source[idx])
        for so in srcs:
            gs=np.unique(sysnames[idx][source[idx]==so])
            for g in gs:
                z=sysnames[idx]==g;out[z]=1/(len(srcs)*len(gs)*sum(z))
        return out
    def pred(model,c,idx,bb):
        nu,mo=model.split('_',1)
        if mo in ['fixed','refit_a0']:return logmond(bb[idx],A0 if mo=='fixed' else 1e-10*10**c[0],nu)
        eta,rc=(10.,400.) if mo=='frozen_gas' else 10**np.array(c)
        extra=eta*np.sum(mweights[idx]*np.exp(-nodes[idx]/rc),axis=1)
        g=bb[idx]+G*MSUN*extra/(rad[idx]*KPC)**2
        return logmond(g,A0,nu)
    def fit(model,idx,bb,yy,tag):
        if model.endswith('fixed') or model.endswith('frozen_gas'):return []
        repair=model.endswith('repair');starts=[[1.,np.log10(400.)],[0.,2.],[1.,3.],[2.,1.5]] if repair else [[-1.],[0.],[1.]];bounds=([-3.,1.],[3.,4.]) if repair else ([-3.],[3.]);rr=[]
        for j,start in enumerate(starts):
            r=least_squares(lambda c:np.sqrt(weights(idx))*(pred(model,c,idx,bb)-yy[idx]),start,bounds=bounds,max_nfev=2000,ftol=1e-11,xtol=1e-11,gtol=1e-11)
            record={'tag':tag,'model':model,'start':j,'train_rows':idx.tolist(),'success':bool(r.success),'coefficients':r.x.tolist(),'cost':float(r.cost),'active_bounds':r.active_mask.tolist(),'nfev':r.nfev};FITLOG.append(record);rr.append(r)
        ok=[r for r in rr if r.success];assert ok
        return min(ok,key=lambda z:z.cost).x.tolist()
    def metrics(pp,yy,idx):
        e=pp-yy[idx]
        return {'B_dex':float(np.sqrt(np.sum(weights(idx)*e*e))),'median_acceleration_error_percent':float(np.median(abs(10**e-1))*100),'rows':len(idx),'systems':len(set(sysnames[idx])),'sources':{so:{'RMS_dex':float(np.sqrt(np.mean(e[source[idx]==so]**2))),'median_acceleration_error_percent':float(np.median(abs(10**e[source[idx]==so]-1))*100)} for so in np.unique(source[idx])}}
    def save(case,model,split,pp,bb,yy,idx,folds,coef):
        for j,i in enumerate(idx):allpred.append({'case':case,'model':model,'split':split,'row':int(i),'system':sysnames[i],'source':source[i],'fold':folds[j],'coefficients_json':json.dumps(coef[j]),'radius_kpc':rad[i],'gbar_m_s2':bb[i],'gobs_m_s2':10**yy[i],'gpred_m_s2':10**pp[j],'log10_residual':pp[j]-yy[i]})
    cases=[('combined_forward_mid',I,bar,np.log10(target)),('CLASH_forward_mid',I[source=='CLASH'],bar,np.log10(target))]
    for tar in ['forward','nfw','einasto']:
        for bkey in ['gbar_m_s2','gbar_no_stars','gbar_full_projected_stars']:
            cases.append((f'XCOP_{tar}_{bkey}',I[source=='XCOP'],np.array([r[bkey] for r in rows]),np.log10([r['gobs_'+tar] for r in rows])))
    for cname,idx,bb,yy in cases:
        out={};result['cases'][cname]=out
        for model in models:
            cc=fit(model,idx,bb,yy,cname+':all');pp=pred(model,cc,idx,bb);d={'coefficients':cc,'all':metrics(pp,yy,idx)};save(cname,model,'all',pp,bb,yy,idx,['all']*len(idx),[cc]*len(idx));out[model]=d
            pp=np.full(len(idx),np.nan);folds=[];coefs=[];cd={}
            for host in np.unique(sysnames[idx]):
                train=idx[sysnames[idx]!=host];test=idx[sysnames[idx]==host];assert host not in set(sysnames[train]);cf=fit(model,train,bb,yy,cname+':host:'+host);pp[sysnames[idx]==host]=pred(model,cf,test,bb);cd[host]=cf
            folds=[sysnames[i] for i in idx];coefs=[cd[v] for v in folds];d['host_holdout']=metrics(pp,yy,idx);d['host_coefficients']=cd;save(cname,model,'host_holdout',pp,bb,yy,idx,folds,coefs)
            if cname.startswith('combined'):
                d['source_holdout']={}
                for so in ['CLASH','XCOP']:
                    train=idx[source[idx]!=so];test=idx[source[idx]==so];cf=fit(model,train,bb,yy,cname+':source:'+so);pp=pred(model,cf,test,bb);d['source_holdout'][so]={'coefficients':cf,'metrics':metrics(pp,yy,test)};save(cname,model,'source_'+so,pp,bb,yy,test,[so]*len(test),[cf]*len(test))
        (R/'results.json').write_text(json.dumps(result,indent=2));print(cname,{m:round(v['host_holdout']['B_dex'],6) for m,v in out.items()},flush=True)
    writecsv(R/'predictions.csv',allpred);result['prediction_rows']=len(allpred);result['prediction_tables']=len(set((r['case'],r['model'],r['split']) for r in allpred));result['fit_calls']=sum(r['start']==0 for r in FITLOG);result['optimizer_attempts']=len(FITLOG);result['optimizer_failures']=sum(not r['success'] for r in FITLOG)
    (R/'results.json').write_text(json.dumps(result,indent=2));print('DONE',result['prediction_rows'],result['prediction_tables'],result['fit_calls'],result['optimizer_attempts'],flush=True)

if __name__=='__main__':main()
