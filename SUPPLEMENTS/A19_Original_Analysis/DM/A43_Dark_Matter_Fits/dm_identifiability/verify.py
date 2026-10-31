from pathlib import Path
import csv, json, hashlib
import numpy as np

R=Path(__file__).resolve().parent;O=R/'results'
D=R.parent/'dm_broad_joint/combined_input.csv'
rows=list(csv.DictReader(D.open()));n=len(rows);I=np.arange(n)
J=json.loads((O/'results.json').read_text())
assert hashlib.sha256(D.read_bytes()).hexdigest()==J['input_sha256']
assert hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest()==J['protocol_sha256']
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows])
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows]);u=np.log10([float(r['radius_kpc']) for r in rows]);y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows])
def rarity(i,a):return np.log10((10**x[i])/(1-np.exp(-np.sqrt(10**(x[i]-a)))))
def weight(i):
    d={p:sorted(set(group[i][pop[i]==p])) for p in set(pop[i])}
    c={(p,g):sum((pop[i]==p)&(group[i]==g)) for p in d for g in d[p]}
    return np.array([1/(len(d)*len(d[pop[j]])*c[(pop[j],group[j])]) for j in i])
def score(pred,ind):
    return float(np.sqrt(np.mean([np.mean([np.mean((pred[(pop[ind]==p)&(group[ind]==g)]-y[ind][(pop[ind]==p)&(group[ind]==g)])**2) for g in set(group[ind][pop[ind]==p])]) for p in set(pop[ind]) ])))
checked=[]
for file in sorted(O.glob('*.csv')):
    rr=list(csv.DictReader(file.open()))
    if file.name.startswith('nearest_'):
        cross=file.name=='nearest_cross_population.csv'
        assert len(rr)==n
        for rec in rr:
            i=int(rec['row']);j=int(rec['neighbor_row']);assert group[i]!=group[j]
            if cross:assert pop[i]!=pop[j]
            keep=group!=group[i]
            if cross:keep&=pop!=pop[i]
            dist=np.hypot(x-x[i],u-u[i]);expected=np.min(dist[keep])
            assert abs(np.hypot(x[i]-x[j],u[i]-u[j])-expected)<1e-12
            assert abs(abs(y[i]-y[j])-float(rec['delta_y']))<1e-12
            assert abs(10**abs(y[i]-y[j])-float(rec['target_factor']))<1e-10
        checked.append(file.name);continue
    ind=np.array([int(r['row']) for r in rr]);fold=np.array([r['fold'] for r in rr]);pred=np.log10([float(r['gpred_m_s2'])/1e-10 for r in rr]);is_kernel=file.name.startswith('kernel_')
    assert len(set(ind))==len(ind)
    assert max(abs(pred-y[ind]-np.array([float(r['log10_residual']) for r in rr])))<1e-12
    for f in set(fold):
        test=ind[fold==f]
        if 'group_holdout' in file.name:
            tr=np.setdiff1d(I,test);sel=J['kernel']['group_holdout']['selection'][int(f)]
        elif 'population_holdout' in file.name:
            tr=I[~np.isin(group,list(set(group[test])))];sel=J['kernel']['population_holdout']['selection'][f]
        else:
            host=f;names=['Milky_Way','G/L'] if host=='MW' else ['M31','A/L'];excluded=I[(pop=='dSph')&np.isin(env,names)]
            assert set(excluded)==set(test) and len(test)==(26 if host=='MW' else 18)
            tr=I[~np.isin(group,list(set(group[test])))];sel=J['kernel']['host_holdout'][host]['selection']
        assert not set(group[tr])&set(group[test])
        h=sel['chosen_bandwidth'];a=sel['training_RAR_loga'];assert str(h)==min(sel['inner_RMS_dex'],key=sel['inner_RMS_dex'].get)
        expected=rarity(test,a)
        if is_kernel:
            delta=(x[test,None]-x[None,tr])**2+(u[test,None]-u[None,tr])**2
            ww=weight(tr)[None,:]*np.exp(-(delta-np.min(delta,axis=1,keepdims=True))/(2*h*h));ww/=ww.sum(axis=1,keepdims=True)
            expected+=np.sum(ww*(y[tr]-rarity(tr,a))[None,:],axis=1)
        assert np.max(abs(expected-pred[fold==f]))<1e-11
    if 'group_holdout' in file.name:metric=J['kernel']['group_holdout']
    elif 'population_holdout' in file.name:metric=J['kernel']['population_holdout']
    else:metric=J['kernel']['host_holdout'][fold[0]]
    assert abs(score(pred,ind)-metric['kernel' if is_kernel else 'RAR']['B_dex'])<1e-11
    checked.append(file.name)
demo=J['curl_demo'];assert abs(demo['analytic_curl_z']+5**(-1.25))<1e-14
assert abs(demo['finite_difference_curl_z']['1e-05']-demo['analytic_curl_z'])<1e-9
v={'status':'passed','tables_checked':len(checked),'files':checked,'checks':['input and protocol hashes','input-only nearest neighbors from distinct groups','kernel and baseline predictions evaluated independently','equal population/group metric','unsplit outer validation groups','complete 26/18 host exclusions','bandwidth minimizes recorded inner score','analytic curl and finite differences']}
(R/'verification.json').write_text(json.dumps(v,indent=2));print(json.dumps(v,indent=2))
