"""T25 diagnostics for the a=0.005 rerun of the T19 audit design.
The statistics functions (batches, ess, rhat) are copied verbatim from the T19 analyze.py."""
from pathlib import Path
import json, math, sys
import numpy as np
D = Path(__file__).resolve().parent

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

def load(runs, tag):
    a=np.loadtxt(runs/f'{tag}.dat'); return a, a[a[:,1]==1]

# T19 (a = 0) group results, copied from A42_Regge_Audit_T19/AUDIT_REPORT.txt
T19 = {(3,0.05):(-1.0470,-0.6731,1.573),(3,0.065):(-0.6022,44.1443,13.032),
       (4,0.05):(-1.0655,31.2166,5.498),(4,0.065):(-0.6793,45.4855,17.716)}

def main(planfile='plan.json', runsdir='runs', outfile='results.json'):
    plan=json.loads((D/planfile).read_text()); runs=D/runsdir; out=[]; arrays={}; full={}
    for j in plan['main']:
        allrows,a=load(runs,j['tag']); assert len(a)==j['measure'],j['tag']; x=a[:,4]; arrays[j['tag']]=x; full[j['tag']]=allrows
        h=len(x)//2; bc=batches(x); b1=batches(x[:h]); b2=batches(x[h:])
        se=max(z['se'] for z in bc); s1=max(z['se'] for z in b1); s2=max(z['se'] for z in b2)
        NE=15*j['L']**4; scale=(2*a[:,2]-j['k']*a[:,3])/NE
        th=allrows[allrows[:,1]==0][:,4]
        # first thermal sweep (in 100-sweep blocks) whose mean is within the measured band of this chain
        band=x.mean()+3*x.std(); blocks=[th[i:i+100].mean() for i in range(0,len(th),100)]
        settle=next((100*i for i,m in enumerate(blocks) if m<=band),None)
        out.append(dict(**j,mean_R=float(x.mean()),batch_se=se,ess_R=ess(x),first_half_R=float(x[:h].mean()),last_half_R=float(x[h:].mean()),
          drift_z_descriptive=float((x[h:].mean()-x[:h].mean())/np.hypot(s1,s2)),
          scaling_identity_ratio=float(scale.mean()),min_R=float(x.min()),max_R=float(x.max()),
          mean_deficit_std=float(a[:,21].mean()),mean_maxq=float(a[:,9].mean()),mean_q=float(a[:,5].mean()),
          max_incremental_drift=float(allrows[:,15].max()),
          mean_sumR2_per_hinge=float(a[:,27].mean()/(50*j['L']**4)) if a.shape[1]>27 else None,
          # completing the square per triangle: I = sumV - k sumDA + a sumR2 >= (1 - k^2/(4a)) sumV (lambda = 1)
          bound_coefficient=(1-j['k']**2/(4*j['a'])) if j.get('a') else None,
          bound_min_margin=float(((allrows[:,2]-j['k']*allrows[:,3]+j['a']*allrows[:,27])/allrows[:,2]-(1-j['k']**2/(4*j['a']))).min()) if j.get('a') and allrows.shape[1]>27 else None,
          first_thermal_block_in_band=settle))
    groups=[]
    keys=sorted(set((r['L'],r['k']) for r in out))
    for L,k in keys:
        rows=[r for r in out if r['L']==L and r['k']==k]; rh=rhat([arrays[r['tag']] for r in rows]); cats={}
        for start in ['flat','rough']:
            rr=[r for r in rows if r['start']==start]; means=np.array([r['mean_R'] for r in rr])
            cats[start]=dict(seed_means=list(means),mean_R=float(means.mean()),between_seed_se=float(means.std(ddof=1)/np.sqrt(len(means))),
                             split_rhat=rhat([arrays[r['tag']] for r in rr]))
        den=np.hypot(cats['flat']['between_seed_se'],cats['rough']['between_seed_se']); contrast=cats['rough']['mean_R']-cats['flat']['mean_R']
        g=dict(L=L,k=k,split_rhat=rh,by_start=cats,start_difference=contrast,start_difference_se=float(den),
               start_difference_z_descriptive=float(contrast/den),max_absolute_drift_z=max(abs(r['drift_z_descriptive']) for r in rows),
               passes_rhat_screen=rh<1.05,max_rough_settle_sweeps=max((r['first_thermal_block_in_band'] if r['first_thermal_block_in_band'] is not None else 10**9) for r in rows if r['start']=='rough'))
        if (L,k) in T19: g['T19_a0']=dict(flat_mean_R=T19[(L,k)][0],rough_mean_R=T19[(L,k)][1],split_rhat=T19[(L,k)][2])
        groups.append(g)
    result=dict(description='T25 rerun of the T19 audit design with a=0.005. R-hat/ESS do not prove equilibrium; drift/start z are descriptive.',
                runs=out,groups=groups)
    (D/outfile).write_text(json.dumps(result,indent=2)+'\n')
    for g in groups:
        print(f"L={g['L']} k={g['k']}: split R-hat {g['split_rhat']:.4f} (flat {g['by_start']['flat']['split_rhat']:.4f}, rough {g['by_start']['rough']['split_rhat']:.4f})"
              f"  mean R flat {g['by_start']['flat']['mean_R']:.4f} rough {g['by_start']['rough']['mean_R']:.4f}  diff z {g['start_difference_z_descriptive']:.2f}"
              f"  max drift z {g['max_absolute_drift_z']:.2f}  rough settle <= {g['max_rough_settle_sweeps']} sweeps  pass={g['passes_rhat_screen']}")
    print('scale identity ratios', min(r['scaling_identity_ratio'] for r in out), max(r['scaling_identity_ratio'] for r in out))
    print('max incremental drift', max(r['max_incremental_drift'] for r in out))
    print('ESS range', min(r['ess_R'] for r in out), max(r['ess_R'] for r in out))
if __name__=='__main__':
    main(*sys.argv[1:])
