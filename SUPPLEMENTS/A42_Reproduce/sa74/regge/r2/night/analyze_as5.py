"""SA74 night stages AS4 and AS5: start independence near the crossing for small curvature-squared couplings, on 4^4 with three flat
starts and three rough starts (the a = 0, k = 0.1 preparations): a = 0.002 at k = 0.175 (AS4) and a = 0.003 at k = 0.20 (AS5).
Split R-hat of the curvature and of the volume over all six chains (screen 1.05, as in T19/T25), flat and rough means, and the
difference in units of its standard error from the seed-to-seed spread. Writes as5_results.json."""
import json, sys, numpy as np
from pathlib import Path
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import batches, ess, rhat
RUNS = Path('/home/claude/sa74/regge/r2/audit_a005/runs')
out = {}
for stage, a, k in (('AS4', '0.002', '0.175'), ('AS5', '0.003', '0.200')):
    grp = {'flat': [], 'rough': []}
    for fn in sorted(RUNS.glob(f'{stage}_a{a}_L4_k{k}_*_s*.dat')):
        if '# done' not in fn.read_text()[-20000:]: continue
        d = np.loadtxt(fn); m = d[d[:, 1] == 1]
        grp['rough' if '_rough_' in fn.name else 'flat'].append(m)
    if len(grp['flat']) < 3 or len(grp['rough']) < 3: print(stage, 'incomplete', {s: len(v) for s, v in grp.items()}); continue
    allm = grp['flat'] + grp['rough']
    fm = np.array([m[:, 4].mean() for m in grp['flat']]); rm = np.array([m[:, 4].mean() for m in grp['rough']])
    z = float((fm.mean() - rm.mean()) / np.sqrt(fm.var(ddof=1) / 3 + rm.var(ddof=1) / 3))
    out[stage] = dict(a=float(a), k=float(k), L=4, rhat_R=float(rhat([m[:, 4] for m in allm])), rhat_V=float(rhat([m[:, 2] for m in allm])),
                      flat_mean=float(fm.mean()), rough_mean=float(rm.mean()), z=z, ess=[float(ess(m[:, 4])) for m in allm],
                      flat=[float(x) for x in fm], rough=[float(x) for x in rm])
    o = out[stage]
    print(f"{stage} a={a} k={k}: split R-hat(R) {o['rhat_R']:.3f}, R-hat(V) {o['rhat_V']:.3f}; flat {o['flat_mean']:+.4f} rough {o['rough_mean']:+.4f} z={z:+.2f}; ESS {min(o['ess']):.0f}-{max(o['ess']):.0f}")
json.dump(dict(description=__doc__, **out), open('/home/claude/sa74/regge/r2/night/as5_results.json', 'w'), indent=1)
