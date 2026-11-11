"""Known-behavior checks for the new diagnostics, independent of Regge outcomes."""
import numpy as np,json
from pathlib import Path
from analyze import rhat,ess,batches
rng=np.random.default_rng(57103)
x=[rng.normal(size=10000) for _ in range(6)]
a=rhat(x);b=rhat([z+(0 if i<3 else 10) for i,z in enumerate(x)])
e=[ess(z) for z in x]
assert .99<a<1.01 and b>4 and min(e)>7000
# Known independent-normal standard error: sigma/sqrt(N)=0.01.
s=batches(x[0])[0]['se'];assert .009<s<.011
out=dict(iid_split_rhat=a,separated_split_rhat=b,iid_ess=e,iid_unbatched_se=s,passed=True)
Path(__file__).with_name('diagnostic_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
