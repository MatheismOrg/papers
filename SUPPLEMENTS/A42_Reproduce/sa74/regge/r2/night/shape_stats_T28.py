"""SA74 night: distribution and memory of the curvature near the zero crossing, a = 0.005, flat starts.
Per chain: integrated autocorrelation time of R (n / ESS), excess kurtosis of R and of the curvature sum (a bimodal two-state
distribution drives it toward -2; Gaussian fluctuations give 0), and a two-peak check (a histogram dip between two maxima).
Chains: cloud stages (FSS, FSS2, the 4^4 and 3^4 scans are skipped: different lengths) and Mac N/Q flat chains (rough starts excluded)."""
import re, sys, glob, json, numpy as np
from pathlib import Path
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import ess
RUNS = Path('/home/claude/sa74/regge/r2/audit_a005/runs'); MAC = Path('/home/claude/sa74/NIGHT/mac_pull')
def kurt(x): z = x - x.mean(); return float((z ** 4).mean() / (z ** 2).mean() ** 2 - 3)
def twopeak(x, bins=30):
    h, e = np.histogram(x, bins=bins); h = np.convolve(h, np.ones(3) / 3, mode='same')
    pk = [i for i in range(1, bins - 1) if h[i] >= h[i - 1] and h[i] >= h[i + 1] and h[i] > 0.25 * h.max()]
    if len(pk) < 2: return 0.0
    a, b = pk[0], pk[-1]; dip = h[a:b + 1].min(); return float(1 - dip / min(h[a], h[b]))
rows = []
for fn in sorted(glob.glob(str(RUNS / 'FSS*_L*_k*_flat.dat'))):
    if '# done' not in open(fn).read()[-20000:]: continue
    m = re.match(r'FSS2?_L(\d+)_k([0-9.]+)_flat', Path(fn).stem); d = np.loadtxt(fn); rows.append((int(m.group(1)), float(m.group(2)), Path(fn).stem, d[d[:, 1] == 1]))
for fn in sorted(MAC.glob('[NQ]_L*_k*.dat')):
    m = re.match(r'[NQ]_L(\d+)_k([0-9.]+)$', fn.stem)
    if not m: continue
    d = np.loadtxt(fn, comments='#'); mm = d[d[:, 1] == 1]
    if len(mm) >= 800: rows.append((int(m.group(1)), float(m.group(2)), fn.stem, mm))
out = []
for L, k, name, m in sorted(rows, key=lambda r: (r[1], r[0])):
    if not (0.239 <= k <= 0.261): continue
    x = m[:, 4]; s = m[:, 3]; n = len(x); e = ess(x)
    out.append(dict(L=L, k=k, chain=name, n=n, tau=float(n / e), kurt_R=kurt(x), kurt_sum=kurt(s), twopeak=twopeak(x)))
json.dump(out, open('/home/claude/sa74/regge/r2/night/shape_stats.json', 'w'), indent=1)
print(f"{'k':>6} {'L':>2} {'chain':<22} {'n':>5} {'tau':>6} {'kurt R':>7} {'kurt sum':>8} {'2-peak':>6}")
for r in out: print(f"{r['k']:.4f} {r['L']:>2} {r['chain']:<22} {r['n']:5d} {r['tau']:6.0f} {r['kurt_R']:+7.2f} {r['kurt_sum']:+8.2f} {r['twopeak']:6.2f}")
