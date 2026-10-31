"""Independent derivatives/direct covariance checks before fitting."""
import json,copy
import numpy as np
import run as r
p=[11.4672833436,.993748658,-.635500615,0.]
base={'base_RAR':[.015701365],'base_NFW':p[:3],'weight':.532087775}
h=1e-5;errors={}
for model,pars in [('NFW_response_mc',p+[-1.2,.2]),('RAR_nuisance',[.015701365,-1.2]),('Blend_offset',[-1.2,.1])]:
 c=r.central(model,pars,base)
 for variable in ['D','Y']:
  ys=[]
  for z in [-h,h]:
   d=copy.deepcopy(r.D);factor=10**z
   if variable=='D':
    d['r']*=factor;d['re']*=factor;d['mp']*=factor**2;d['ms']*=factor**2
   else:
    d['gb']=d['gb']-d['star']+d['star']*factor;d['star']*=factor;d['mp']=d['mp']-d['ms']+d['ms']*factor;d['ms']*=factor
   val=r.central(model,pars,base,d)['log']+(z if variable=='D' else 0)
   ys.append(val)
  num=(ys[1]-ys[0])/(2*h);error=float(np.max(abs(num[r.MASK]-c['J'+variable][r.MASK])));errors[model+'_'+variable]=error
for variable,kw in [('M','delta_m'),('C','delta_c')]:
 num=(r.nfw(p,**{kw:h})['log']-r.nfw(p,**{kw:-h})['log'])/(2*h)
 errors['NFW_'+variable]=float(np.max(abs(num[r.MASK]-r.nfw(p)['J'+variable][r.MASK])))
assert max(errors.values())<1e-7,errors
largest=0
for model in r.MODELS:
 _,starts=r.spec(model);pars=starts[0];c,v,U=r.distribution(model,pars,base);scores=r.score_arrays(model,pars,base)
 for j in [0,20,70,100,148]:
  n=r.COUNT[j];u=U[j,:n];C=np.diag(v[j,:n])+u@u.T;e=c['log'][j,:n]-np.log10(r.D['go'][j,:n])
  sign,ld=np.linalg.slogdet(C);assert sign>0
  val=.5*(n*np.log(2*np.pi)+ld+e@np.linalg.solve(C,e))/n
  largest=max(largest,abs(val-scores['nll'][j]))
assert largest<1e-10
old=json.loads((r.OLD/'results/all_record.json').read_text());pr=10**r.flatten(r.nfw(old['models']['NFW']['fit']['p'])['log']);ref=np.array(old['models']['NFW']['pred'])
diff=float(np.max(abs(np.log10(pr/ref))));assert diff<1e-10
out=dict(derivative_errors=errors,max_direct_covariance_score_difference=largest,baseline_max_log_prediction_difference=diff,galaxies=r.NG,rows=len(r.DF),max_fractional_distance_error=float(max(r.D['sd'][:,0])*r.LN),passed=True)
(r.OUT/'audit_pre.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
