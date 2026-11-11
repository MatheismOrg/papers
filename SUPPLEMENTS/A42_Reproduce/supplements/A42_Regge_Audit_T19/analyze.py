"""Independent diagnostics; does not import the original analysis driver."""
from pathlib import Path
import json, math
import numpy as np
D=Path(__file__).resolve().parent

def batches(x):
    ans=[];b=1
    while len(x)//b>=16:
        n=len(x)//b; m=x[:n*b].reshape(n,b).mean(1)
        ans.append(dict(batch=b,batches=n,se=float(m.std(ddof=1)/np.sqrt(n))))
        b*=2
    return ans

def ess(x):
    z=x-x.mean(); n=len(z)
    f=np.fft.rfft(z,n=2*n); ac=np.fft.irfft(f*f.conj(),n=2*n)[:n]; ac/=np.arange(n,0,-1)
    if ac[0]<=0:return float(n)
    ac/=ac[0];total=0.;previous=float('inf')
    for i in range(1,n-1,2):
        pair=float(ac[i]+ac[i+1])
        if pair<=0:break
        pair=min(pair,previous);previous=pair;total+=pair
    return float(min(n,n/(1+2*total)))

def rhat(chains):
    # Ordinary split R-hat. Useful rejection diagnostic; not rank-normalized.
    n=min(map(len,chains))//2
    z=np.array([part for x in chains for part in (x[:n],x[-n:])])
    W=z.var(1,ddof=1).mean(); B=n*z.mean(1).var(ddof=1)
    return float(np.sqrt(((n-1)*W/n+B/n)/W))

def load(tag):
    a=np.loadtxt(D/'runs'/f'{tag}.dat');return a[a[:,1]==1]

def main():
    plan=json.loads((D/'plan.json').read_text());out=[];arrays={}
    for j in plan['main']:
        a=load(j['tag']);assert len(a)==j['measure'],j['tag'];x=a[:,4];arrays[j['tag']]=x
        h=len(x)//2;bc=batches(x);b1=batches(x[:h]);b2=batches(x[h:])
        se=max(z['se'] for z in bc);s1=max(z['se'] for z in b1);s2=max(z['se'] for z in b2)
        scale=(2*a[:,2]-j['k']*a[:,3])/(15*j['L']**4)
        out.append(dict(**j,mean_R=float(x.mean()),batch_se=se,ess_R=ess(x),first_half_R=float(x[:h].mean()),last_half_R=float(x[h:].mean()),
          drift_z_descriptive=float((x[h:].mean()-x[:h].mean())/np.hypot(s1,s2)),batch_curve=bc,
          scaling_identity_ratio=float(scale.mean()),min_R=float(x.min()),max_R=float(x.max()),
          mean_deficit_std=float(a[:,21].mean()),mean_maxq=float(a[:,9].mean()),max_incremental_drift=float(a[:,15].max())))
    groups=[]
    for L in [3,4]:
      for k in [.05,.065]:
        rows=[r for r in out if r['L']==L and r['k']==k];rh=rhat([arrays[r['tag']] for r in rows]);cats={}
        for start in ['flat','rough']:
            rr=[r for r in rows if r['start']==start];means=np.array([r['mean_R'] for r in rr]);cats[start]=dict(seed_means=list(means),mean_R=float(means.mean()),between_seed_se=float(means.std(ddof=1)/np.sqrt(len(means))),split_rhat=rhat([arrays[r['tag']] for r in rr]))
        den=np.hypot(cats['flat']['between_seed_se'],cats['rough']['between_seed_se'])
        contrast=cats['rough']['mean_R']-cats['flat']['mean_R']
        groups.append(dict(L=L,k=k,split_rhat=rh,by_start=cats,start_difference=contrast,start_difference_se=float(den),start_difference_z_descriptive=float(contrast/den),max_absolute_drift_z=max(abs(r['drift_z_descriptive']) for r in rows),passes_rhat_screen=rh<1.05))
    result=dict(description='New execution of unchanged source; full per-sweep histories retained. R-hat/ESS do not prove equilibrium. Drift/start z values are descriptive, not calibrated Gaussian significance.',runs=out,groups=groups)
    (D/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    for g in groups: print('L',g['L'],'k',g['k'],'Rhat',round(g['split_rhat'],4),'means', {s:round(g['by_start'][s]['mean_R'],5) for s in ['flat','rough']},'max drift z',round(g['max_absolute_drift_z'],2))
if __name__=='__main__':main()
