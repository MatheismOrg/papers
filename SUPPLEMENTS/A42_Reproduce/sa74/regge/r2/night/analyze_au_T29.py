"""SA74 night stages AU1 and AU2: the curvature-squared ensemble at a = 0.02 (Hamber 1993: k_c = 0.4244(36), delta = 0.628(20)).
For each finished flat-start chain: average curvature R with a batch-means error, the curvature fluctuation chi = 4 var(sum dA)/<sum V>
with a five-block error, the fluctuation slope cov(R, sum dA) (= d<R>/dk in equilibrium) and tau = n/ESS; then the zero crossing of
R on each grid by linear interpolation and the fitted power of chi in L at each k measured on three or more grids. Writes au_results.json."""
import json, re, sys, glob
from pathlib import Path
import numpy as np
AUD = Path('/home/claude/sa74/regge/r2/audit_a005'); sys.path.insert(0, str(AUD)); from analyze_r2 import batches, ess
rows = []
for fn in sorted(glob.glob(str(AUD / 'runs' / 'AU[12]_a0.02_L*_k*_flat.dat'))):
    if '# done' not in open(fn).read()[-20000:]: continue
    a = np.loadtxt(fn); m = a[a[:, 1] == 1]; x = m[:, 4]
    L = int(re.search(r'_L(\d+)_', fn).group(1)); k = float(re.search(r'_k([0-9.]+)_', fn).group(1))
    bl = np.array_split(np.arange(len(x)), 5); cb = [4 * m[b, 3].var() / m[b, 2].mean() for b in bl]
    rows.append(dict(chain=Path(fn).stem, L=L, k=k, n=len(x), R=float(x.mean()), R_se=float(max(b['se'] for b in batches(x))),
                     chi=float(4 * m[:, 3].var() / m[:, 2].mean()), chi_se=float(np.std(cb, ddof=1) / np.sqrt(5)),
                     slope=float(np.cov(x, m[:, 3])[0, 1]), tau=float(len(x) / ess(x)), l2=float(m[:, 5].mean())))
cross = {}
for L in sorted(set(r['L'] for r in rows)):
    p = sorted((r['k'], r['R']) for r in rows if r['L'] == L)
    cross[L] = next((k1 + (k2 - k1) * (-r1) / (r2 - r1) for (k1, r1), (k2, r2) in zip(p, p[1:]) if r1 < 0 <= r2), None)
powers = {}
for k in sorted(set(r['k'] for r in rows)):
    q = sorted((r['L'], r['chi']) for r in rows if r['k'] == k)
    if len(q) >= 3: powers[k] = dict(L=[a for a, _ in q], chi=[b for _, b in q], power=float(np.polyfit(np.log([a for a, _ in q]), np.log([b for _, b in q]), 1)[0]))
json.dump(dict(description=__doc__, chains=rows, crossing=cross, chi_power=powers), open('/home/claude/sa74/regge/r2/night/au_results.json', 'w'), indent=1)
for r in rows: print(f"L={r['L']} k={r['k']:.3f} R={r['R']:+.4f}±{r['R_se']:.4f} chi={r['chi']:.3f}±{r['chi_se']:.3f} dR/dk={r['slope']:.3f} tau={r['tau']:.0f}")
print('crossings', cross); print('chi powers', {k: round(v['power'], 2) for k, v in powers.items()})
