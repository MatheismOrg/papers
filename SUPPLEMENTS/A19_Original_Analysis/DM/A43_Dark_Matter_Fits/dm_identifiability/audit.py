from pathlib import Path
import csv, json, hashlib
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist

R=Path(__file__).resolve().parent
DATA=R.parent/'dm_broad_joint/combined_input.csv'
O=R/'results';O.mkdir(exist_ok=True)
rows=list(csv.DictReader(DATA.open()))
n=len(rows);I=np.arange(n)
x=np.log10([float(r['gbar_m_s2'])/1e-10 for r in rows])
u=np.log10([float(r['radius_kpc']) for r in rows])
y=np.log10([float(r['gobs_m_s2'])/1e-10 for r in rows])
X=np.column_stack([x,u])
pop=np.array([r['population'] for r in rows]);group=np.array([r['group'] for r in rows]);env=np.array([r['environment'] for r in rows])
POPS=sorted(set(pop));D2=cdist(X,X,'sqeuclidean')
HS=[.15,.30,.60,1.20]
WARN=[];COUNTS=0

def weights(idx):
    w=np.empty(len(idx));ps=np.unique(pop[idx])
    for p in ps:
        gs=np.unique(group[idx][pop[idx]==p])
        for g in gs:
            m=(pop[idx]==p)&(group[idx]==g)
            w[m]=1/(len(ps)*len(gs)*sum(m))
    return w

def B(pr,idx):return float(np.sqrt(np.sum(weights(idx)*(pr-y[idx])**2)))

def metrics(pr,idx):
    return {'B_dex':B(pr,idx),'population':{p:{'group_RMS_dex':B(pr[pop[idx]==p],idx[pop[idx]==p]),'median_acceleration_error_percent':float(100*np.median(abs(10**(pr[pop[idx]==p]-y[idx[pop[idx]==p]])-1))),'median_speed_error_percent':float(100*np.median(abs(10**((pr[pop[idx]==p]-y[idx[pop[idx]==p]])/2)-1)))} for p in np.unique(pop[idx])}}

def folds(idx,seed,nf):
    rng=np.random.default_rng(seed);mp={}
    for p in sorted(set(pop[idx])):
        gs=[g for g in np.unique(group[idx][pop[idx]==p]) if g not in mp]
        ct=np.array([sum(v==j for g,v in mp.items() if g in group[idx][pop[idx]==p]) for j in range(nf)])
        for g in rng.permutation(gs):
            j=int(np.argmin(ct));mp[g]=j;ct[j]+=1
    return np.array([mp[g] for g in group[idx]])

def rar(idx,a):return x[idx]-np.log10(-np.expm1(-10**(.5*(x[idx]-a))))

def fit(idx):
    global COUNTS
    COUNTS+=1
    res=least_squares(lambda c:np.sqrt(weights(idx))*(rar(idx,c[0])-y[idx]),[.1],bounds=([-3.],[3.]),ftol=1e-10,xtol=1e-10,gtol=1e-10)
    if not res.success:WARN.append(res.message)
    return float(res.x[0])

def kernel(tr,te,h,a):
    assert not set(group[tr])&set(group[te])
    logw=-D2[np.ix_(te,tr)]/(2*h*h)+np.log(weights(tr))[None,:]
    logw-=np.max(logw,axis=1,keepdims=True)
    w=np.exp(logw);w/=np.sum(w,axis=1,keepdims=True)
    return rar(te,a)+w@(y[tr]-rar(tr,a))

def choose(tr,seed):
    ff=folds(tr,seed,3);ip={h:np.empty(len(tr)) for h in HS};aa=[]
    for j in range(3):
        tri=tr[ff!=j];tei=tr[ff==j];a=fit(tri);aa.append(a)
        for h in HS:ip[h][ff==j]=kernel(tri,tei,h,a)
    scores={h:B(ip[h],tr) for h in HS}
    return min(scores,key=scores.get),{'inner_RMS_dex':scores,'inner_RAR_loga':aa}

def savepred(name,idx,pr,fd):
    with (O/(name+'.csv')).open('w') as f:
        wr=csv.writer(f);wr.writerow(['row','population','system','group','fold','gbar_m_s2','radius_kpc','gobs_m_s2','gpred_m_s2','log10_residual'])
        for j,i in enumerate(idx):wr.writerow([i,pop[i],rows[i]['system'],group[i],fd[j],float(10**x[i]*1e-10),float(10**u[i]),float(10**y[i]*1e-10),float(10**pr[j]*1e-10),float(pr[j]-y[i])])

J={'input_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((R/'PROTOCOL.txt').read_bytes()).hexdigest(),'bandwidths_dex':HS,'nearest':{},'kernel':{}}

# Input-only nearest matching, independent of target differences.
WW=weights(I)
for cross in [False,True]:
    dd=D2.copy();dd[group[:,None]==group[None,:]]=np.inf
    if cross:dd[pop[:,None]==pop[None,:]]=np.inf
    nn=np.argmin(dd,axis=1);dx=abs(x-x[nn]);du=abs(u-u[nn]);dy=abs(y-y[nn]);dist=np.sqrt(dd[I,nn])
    label='cross_population' if cross else 'distinct_group'
    tab=[]
    for i,j in enumerate(nn):
        tab.append({'row':i,'neighbor_row':int(j),'population':pop[i],'system':rows[i]['system'],'group':group[i],'neighbor_population':pop[j],'neighbor_system':rows[j]['system'],'neighbor_group':group[j],'delta_x':float(dx[i]),'delta_u':float(du[i]),'distance_dex':float(dist[i]),'delta_y':float(dy[i]),'target_factor':float(10**dy[i]),'gbar_m_s2':rows[i]['gbar_m_s2'],'neighbor_gbar_m_s2':rows[j]['gbar_m_s2'],'radius_kpc':rows[i]['radius_kpc'],'neighbor_radius_kpc':rows[j]['radius_kpc'],'gobs_m_s2':rows[i]['gobs_m_s2'],'neighbor_gobs_m_s2':rows[j]['gobs_m_s2'],'loggobs_err_dex':rows[i].get('loggobs_err_dex',''),'neighbor_loggobs_err_dex':rows[j].get('loggobs_err_dex','')})
    with (O/('nearest_'+label+'.csv')).open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(tab[0]));wr.writeheader();wr.writerows(tab)
    rec={}
    for tol in [.05,.1,.2]:
        m=(dx<=tol)&(du<=tol);wm=WW[m]/sum(WW[m]);unique_pairs={tuple(sorted([int(i),int(nn[i])])) for i in I[m]}
        rec[str(tol)]={'rows_with_match':int(sum(m)),'groups_with_match':len(set(group[m])),'unordered_pairs':len(unique_pairs),'balanced_coverage':float(sum(WW[m])),'conditional_balanced_RMS_pair_difference_dex':float(np.sqrt(sum(wm*dy[m]**2))),'conditional_balanced_fraction_factor2':float(sum(wm*(dy[m]>=np.log10(2)))),'conditional_balanced_fraction_factor3':float(sum(wm*(dy[m]>=np.log10(3)))),'population_coverage':{p:float(sum(WW[(pop==p)&m])/sum(WW[pop==p])) for p in POPS}}
    candidates=[t for t in tab if t['delta_x']<=.1 and t['delta_u']<=.1]
    seen=set();ex=[]
    for t in sorted(candidates,key=lambda z:-z['delta_y']):
        key=tuple(sorted([t['group'],t['neighbor_group']]))
        if key not in seen:seen.add(key);ex.append(t)
        if len(ex)>=12:break
    J['nearest'][label]={'tolerances':rec,'largest_matched_examples':ex}

# Nested group CV and population/host transfer.
ff=folds(I,20261027,5);pr=np.empty(n);base=np.empty(n);selection=[]
for j in range(5):
    tr=I[ff!=j];te=I[ff==j];h,rec=choose(tr,39173+j);a=fit(tr);pr[te]=kernel(tr,te,h,a);base[te]=rar(te,a)
    selection.append(dict(outer_fold=j,chosen_bandwidth=h,training_RAR_loga=a,**rec))
    print('CV',j,h,flush=True)
J['kernel']['group_holdout']={'kernel':metrics(pr,I),'RAR':metrics(base,I),'selection':selection}
savepred('kernel_group_holdout',I,pr,ff);savepred('RAR_group_holdout',I,base,ff)

pr=np.empty(n);base=np.empty(n);selection={};ffp=np.empty(n,dtype=object)
for j,p in enumerate(POPS):
    te=I[pop==p];tg=set(group[te]);tr=I[np.array([g not in tg for g in group])];h,rec=choose(tr,50000+j);a=fit(tr);pr[te]=kernel(tr,te,h,a);base[te]=rar(te,a);ffp[te]=p
    selection[p]=dict(chosen_bandwidth=h,training_RAR_loga=a,**rec)
    print('Population',p,h,flush=True)
J['kernel']['population_holdout']={'kernel':metrics(pr,I),'RAR':metrics(base,I),'selection':selection}
savepred('kernel_population_holdout',I,pr,ffp);savepred('RAR_population_holdout',I,base,ffp)

J['kernel']['host_holdout']={}
for j,hst in enumerate(['MW','M31']):
    names=['Milky_Way','G/L'] if hst=='MW' else ['M31','A/L']
    te=I[(pop=='dSph')&np.isin(env,names)];assert len(te)==(26 if hst=='MW' else 18)
    testgroups=set(group[te]);tr=I[np.array([g not in testgroups for g in group])]
    assert not set(group[tr])&testgroups
    h,rec=choose(tr,61000+j);a=fit(tr);pr=kernel(tr,te,h,a);base=rar(te,a)
    J['kernel']['host_holdout'][hst]={'kernel':metrics(pr,te),'RAR':metrics(base,te),'selection':dict(chosen_bandwidth=h,training_RAR_loga=a,**rec),'excluded_rows':len(te)}
    savepred('kernel_host_'+hst,te,pr,np.repeat(hst,len(te)));savepred('RAR_host_'+hst,te,base,np.repeat(hst,len(te)))
    print('Host',hst,h,flush=True)

# Analytic and central-difference conservative-field test.
def gn(a,b,xx,yy):return np.array([-a*xx,-b*yy])
def ga(a,b,xx,yy):
    v=gn(a,b,xx,yy);return np.linalg.norm(v)**(-.5)*v
a,b,xx,yy=1.,2.,1.,1.;s=np.sqrt(a*a*xx*xx+b*b*yy*yy)
analytic=.5*a*b*(a-b)*xx*yy/s**2.5
num={}
for eps in [1e-3,1e-4,1e-5]:
    curl=(ga(a,b,xx+eps,yy)[1]-ga(a,b,xx-eps,yy)[1])/(2*eps)-(ga(a,b,xx,yy+eps)[0]-ga(a,b,xx,yy-eps)[0])/(2*eps)
    num[str(eps)]=float(curl)
J['curl_demo']={'potential':'PhiN=(a*x^2+b*y^2)/2','nu':'s^(-1/2), s=|gN|, a0=1','a':a,'b':b,'x':xx,'y':yy,'analytic_curl_z':float(analytic),'finite_difference_curl_z':num,'formula':'a*b*(a-b)*x*y/(2*(a*a*x*x+b*b*y*y)^(5/4))','interpretation':'Nonzero curl prevents this algebraically rescaled field being a gradient of one scalar potential. It does not apply to properly solved AQUAL/QUMOND.'}
J['fit_calls']=COUNTS;J['optimizer_warnings']=WARN
(O/'results.json').write_text(json.dumps(J,indent=2))
print('DONE',COUNTS,'fits',WARN,flush=True)
print(json.dumps({k:{m:v[m]['B_dex'] for m in ['kernel','RAR']} for k,v in J['kernel'].items() if k!='host_holdout'},indent=2))
