from pathlib import Path
import csv, json, hashlib
import numpy as np

R=Path(__file__).resolve().parent;T=R.parent/'dm_unified_regular';P=R.parent/'dm_unified_repair'
J=json.loads((T/'results/results.json').read_text());prior=json.loads((P/'results/results.json').read_text())
manifest=json.loads((T/'pre_fit_manifest.json').read_text())
for fn,expected in manifest.items():assert hashlib.sha256((T/fn).read_bytes()).hexdigest()==expected
assert J['fit_script_sha256']==manifest['run.py'] and J['protocol_sha256']==manifest['PROTOCOL.txt']
spec=J['spec']['mixed_regular'];assert spec['lo']==[-3.,-5.,-5.,0.,0.,-2.,0.] and spec['hi']==[3.,3.,3.,.95,.95,0.,4.]
assert len(spec['starts'])==3 and len(J['fit_log'])==144 and J['fit_calls']==48
assert all(v['success'] for v in J['fit_log']) and not J['optimizer_failures']
def direct(c,gb,r):
    l=10**c[5];rc=10**c[6];m=(gb/1e-10)*r*r;s=(gb/1e-10)*r
    inverse=10**c[1]*m**(-c[3])*(r**4/(r**4+l**4))/(1+(gb/1e-8)**2)
    positive=10**c[2]*s**c[4]/(1+r/rc)
    ae=1e-10*(10**c[0]+inverse+positive)
    return np.log10(gb/(-np.expm1(-np.sqrt(gb/ae)))/1e-10)

checked=[];maxerr=0.;boot_checked=0;probe=[];diag={};selected_bound=[]
origbase=list(csv.DictReader((R.parent/'dm_broad_joint/combined_input.csv').open()))
xcop=list(csv.DictReader((R.parent/'dm_fixed_radius/fixed_radius_inputs.csv').open()))
for dataset in ['original','augmented']:
    f=T/(dataset+'_input.csv');assert f.read_bytes()==(P/f.name).read_bytes();assert hashlib.sha256(f.read_bytes()).hexdigest()==J['input_hashes'][f.name]
    assert J['source_hashes']['dm_unified_repair/'+f.name]==J['input_hashes'][f.name]
    rows=list(csv.DictReader(f.open()));n=len(rows);I=np.arange(n);assert n==(2910 if dataset=='original' else 2931)
    if dataset=='original':assert rows==origbase
    else:
        for old,new in zip(origbase,rows[:2910]):
            expect=dict(old)
            if expect['population']=='CLASH':expect['population']='clusters'
            assert new==expect
        for raw,new in zip(xcop,rows[2910:]):
            assert new['system']==new['group']=='XCOP_'+raw['cluster']
            assert new['gbar_m_s2']==raw['gbar'] and new['gobs_m_s2']==raw['gobs_forward'] and new['radius_kpc']==raw['radius_kpc'] and new['source']=='XCOP'
    d=J['models_by_dataset'][dataset];model=d['models']['mixed_regular']
    assert d['comparison_source_sha256']==hashlib.sha256((P/'results/results.json').read_bytes()).hexdigest()
    for k in J['comparators_reused']:assert d['comparator_metrics_reused'][k]==prior['models_by_dataset'][dataset]['models'][k]
    pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows]);source=np.array(['CLASH' if r['source']=='Tian2020_CDS_fig2' else r['source'] for r in rows])
    gb=np.array([float(r['gbar_m_s2']) for r in rows]);rad=np.array([float(r['radius_kpc']) for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows]);PS=sorted(set(pop))
    assert len(set(group))==(278 if dataset=='original' else 285)
    def B(pr,ind):return float(np.sqrt(np.mean([np.mean([np.mean((pr[(pop[ind]==p)&(group[ind]==g)]-y[ind][(pop[ind]==p)&(group[ind]==g)])**2) for g in set(group[ind][pop[ind]==p])]) for p in set(pop[ind]) ])))
    def match_metrics(pr,ind,met):
        assert abs(B(pr,ind)-met['B_dex'])<1e-11
        for key,vec in [('populations',pop),('sources',source)]:
            for name,m in met[key].items():
                mask=vec[ind]==name;idx=ind[mask];err=pr[mask]-y[idx]
                assert m['rows']==len(idx) and m['groups']==len(set(group[idx]))
                assert abs(B(pr[mask],idx)-m['group_RMS_dex'])<1e-11
                assert abs(np.sqrt(np.mean(err**2))-m['RMS_dex'])<1e-11
                assert abs(100*np.median(abs(10**err-1))-m['median_acceleration_error_percent'])<1e-9
                assert abs(100*np.median(abs(10**(err/2)-1))-m['median_speed_error_percent'])<1e-9
    for file in sorted((T/'results'/dataset).glob('mixed_regular_*.csv')):
        rr=list(csv.DictReader(file.open()));ind=np.array([int(r['row']) for r in rr]);fd=np.array([r['fold'] for r in rr]);pred=np.log10([float(r['gpred_m_s2'])/1e-10 for r in rr]);tail=file.stem[len('mixed_regular_'):];expected=np.empty(len(ind));assert len(ind)>0
        for rec,i in zip(rr,ind):assert rec['system']==rows[i]['system'] and rec['group']==group[i] and rec['source']==source[i]
        if tail=='all':
            expected=direct(model['coefficients'],gb[ind],rad[ind]);met=model['training'];selectors=[('all',I,model['coefficients'])]
        else:
            selectors=[]
            for fold in set(fd):
                mask=fd==fold;te=ind[mask]
                if tail.startswith('cv_'):
                    seed=tail[3:];c=model['CV'][seed]['coefficients'][int(fold)];met=model['CV'][seed]['metrics'];tr=np.setdiff1d(I,te);tag=f'cv_{seed}_{fold}'
                elif tail=='population_holdout':
                    assert set(pop[te])=={fold};c=model['population_holdout']['coefficients'][fold];met=model['population_holdout']['metrics'];tr=I[pop!=fold];tag='population_'+fold
                elif tail.startswith('host_'):
                    h=fold;names=['Milky_Way','G/L'] if h=='MW' else ['M31','A/L'];assert set(te)==set(I[(pop=='dSph')&np.isin(env,names)]) and len(te)==(26 if h=='MW' else 18)
                    c=model['host_holdout'][h]['coefficients'];met=model['host_holdout'][h]['metrics'];tr=np.setdiff1d(I,te);tag='host_'+h
                else:
                    assert tail=='source_'+fold;assert set(te)==set(I[source==fold]) and len(te)==(21 if fold=='XCOP' else 84)
                    c=model['source_holdout'][fold]['coefficients'];met=model['source_holdout'][fold]['metrics'];tr=np.setdiff1d(I,te);tag='source_'+fold
                assert not set(group[tr])&set(group[te]);expected[mask]=direct(c,gb[te],rad[te]);selectors.append((tag,tr,c))
        for tag,tr,c in selectors:
            fits=[v for v in J['fit_log'] if v['dataset']==dataset and v['split']==tag];assert len(fits)==3 and all(v['train_rows']==len(tr) for v in fits)
            selected=min(fits,key=lambda z:z['cost']);assert np.array_equal(np.array(c),np.array(selected['coefficients']))
            assert selected['success'];selected_bound.append({'dataset':dataset,'split':tag,'active_bounds':selected['active_bounds']})
        err=float(np.max(abs(expected-pred)));maxerr=max(maxerr,err);assert err<1e-11
        assert max(abs(pred-y[ind]-np.array([float(r['log10_residual']) for r in rr])))<1e-12
        match_metrics(pred,ind,met);checked.append(dataset+'/'+file.name)
    # Reproduce all paired frozen-prediction resampling intervals.
    for split in ['cv_20261027','population_holdout']:
        files={'mixed_regular':T/'results'/dataset/('mixed_regular_'+split+'.csv')}
        files.update({k:P/'results'/dataset/(k+'_'+split+'.csv') for k in J['comparators_reused']})
        baseline_rows=list(csv.DictReader(files['mixed_regular'].open()));groups=sorted(set(group));lookup={g:i for i,g in enumerate(groups)}
        draws=np.random.default_rng(20261028).poisson(1,size=(5000,len(groups)));scores={};point={}
        for k,file in files.items():
            rr=list(csv.DictReader(file.open()));assert all(a['row']==b['row'] and a['group']==b['group'] and a['fold']==b['fold'] for a,b in zip(baseline_rows,rr))
            pred=np.log10([float(r['gpred_m_s2'])/1e-10 for r in rr]);point[k]=B(pred,I);parts=[]
            for p in PS:
                gs=sorted(set(group[pop==p]));v=np.array([np.mean((pred[(pop==p)&(group==g)]-y[(pop==p)&(group==g)])**2) for g in gs]);w=draws[:,[lookup[g] for g in gs]];parts.append((w@v)/w.sum(axis=1))
            scores[k]=np.sqrt(np.mean(parts,axis=0))
        for k in J['comparators_reused']:
            a=d['bootstrap_'+split]['mixed_regular minus '+k];e=scores['mixed_regular']-scores[k]
            assert abs(point['mixed_regular']-point[k]-a['point_difference_dex'])<1e-11
            assert np.max(abs(np.percentile(e,[2.5,97.5])-a['percentile_95_dex']))<1e-11
            assert abs(np.mean(e<0)-a['fraction_negative_descriptive_only'])<1e-12
            boot_checked+=1
    # Finite stellar probes and fixed-radius response derivative on an explicit grid.
    GM=1.32712440018e20;AU=149597870700.;KPC=3.0856775814913673e19;au=np.array([.4,1,10,30,1000,10000.]);r=au*AU;g=GM/r**2;c=model['coefficients'];lp=direct(c,g,r/KPC);ratio=10**lp*1e-10/g
    for j,aa in enumerate(au):probe.append({'dataset':dataset,'radius_AU':aa,'gbar_m_s2':g[j],'predicted_over_Newtonian':ratio[j]})
    gx=np.linspace(-15,15,601);gu=np.linspace(-12,6,73);lp=direct(c,1e-10*10**gx[:,None],10**gu[None,:]);slope=np.diff(lp,axis=0)/(gx[1]-gx[0]);ij=np.unravel_index(np.argmin(slope),slope.shape)
    diag[dataset]={'minimum_fixed_radius_dlogg_dloggbar':float(slope[ij]),'grid_x':float(gx[ij[0]]),'grid_u':float(gu[ij[1]]),'negative_beyond_1e-8_steps':int(np.sum(slope < -1e-8))}

with (R/'regular_stellar_probes.csv').open('w') as f:
    wr=csv.DictWriter(f,fieldnames=list(probe[0]));wr.writeheader();wr.writerows(probe)
out={'status':'passed','prediction_tables_checked':len(checked),'bootstrap_comparisons_checked':boot_checked,'maximum_direct_equation_error_dex':maxerr,'source_results_sha256':hashlib.sha256((T/'results/results.json').read_bytes()).hexdigest(),'checks':['pre-fit script and protocol hashes','fixed lower length bound10pc','original inputs unchanged','augmented rows preserve original sources and append21XCOP entries','direct formula predictions','population/source/group/median metrics','three starts per selected fit and minimum-cost selection','group/source/host exclusions','baseline comparator inputs and folds','5000pairedresample24intervals'],'physical_audit_note':'See regularity_audit.txt: analytic regular-force and nonnegative isolated-point-source effective-density proofs, incomplete field-theory interpretation, and imposed stellar-screening lower bound.','monotonicity_grid':{'x_range':[-15,15],'x_points':601,'u_range':[-12,6],'u_points':73,'results':diag},'selected_bound_contacts':selected_bound,'files':checked}
(R/'regular_verification.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['selected_bound_contacts','files']},indent=2));print('Stellar probes:',probe)
