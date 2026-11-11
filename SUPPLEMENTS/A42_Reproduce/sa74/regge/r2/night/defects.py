"""SA74 night: what keeps a rough start above the flat start on 8^4 (plans Q1, Q3, Q6). For each 8^4 chain near the crossing, the
longest squared edge (column 10, maxq) and R in twelve equal blocks of the whole chain (continuations appended); and, from each
chain's final configuration (mac_pull/cfg/*.cfg: int L, double, then the 15 L^4 squared edge lengths), the number of edges longer
than four times the mean squared length and the number of vertices they touch; and, in flat starts near the crossing on every grid,
the correlation of 100-sweep means of the longest edge and of R, with the shift in R per unit of longest edge. Writes defects.json."""
import json, struct, glob
from pathlib import Path
import numpy as np
MAC = Path('/home/claude/sa74/NIGHT/mac_pull')
CH = {'N_L8_k0.245': ['N_L8_k0.245'], 'Q_L8_k0.245_r2': ['Q_L8_k0.245_r2'], 'Q_L8_k0.245_r3': ['Q_L8_k0.245_r3'],
      'Q_L8_k0.245_rough': ['Q_L8_k0.245_rough', 'Q_L8_k0.245_rough_cont'], 'Q_L8_k0.245_rough2': ['Q_L8_k0.245_rough2'],
      'Q_L8_k0.2475': ['Q_L8_k0.2475'], 'Q_L8_k0.2475_rough': ['Q_L8_k0.2475_rough'],
      'N_L8_k0.250': ['N_L8_k0.250'], 'Q_L8_k0.250_rough': ['Q_L8_k0.250_rough', 'Q_L8_k0.250_rough_cont'], 'Q_L8_k0.250_rough2': ['Q_L8_k0.250_rough2'],
      'N_L8_k0.255': ['N_L8_k0.255'], 'Q_L8_k0.255_rough': ['Q_L8_k0.255_rough']}
out = dict(description=__doc__, series={}, configs={})
for name, parts in CH.items():
    arrs = [np.loadtxt(MAC / f'{p}.dat', comments='#') for p in parts if (MAC / f'{p}.dat').exists()]
    arrs = [a for a in arrs if a.ndim == 2 and len(a) > 24]
    if not arrs: continue
    a = np.vstack(arrs); b = np.array_split(np.arange(len(a)), 12)
    out['series'][name] = dict(rows=len(a), R=[float(a[i, 4].mean()) for i in b], maxq=[float(a[i, 9].mean()) for i in b],
                               R_last_third=float(a[2 * len(a) // 3:, 4].mean()), maxq_last_third=float(a[2 * len(a) // 3:, 9].mean()))
for fn in sorted((MAC / 'cfg').glob('*.cfg')):
    raw = fn.read_bytes(); L = struct.unpack('<i', raw[:4])[0]; q = np.frombuffer(raw[12:], dtype='<f8')
    if len(q) != 15 * L ** 4: continue
    long = np.nonzero(q > 4 * q.mean())[0]
    out['configs'][fn.stem] = dict(L=L, mean_q=float(q.mean()), max_q=float(q.max()), n_long=int(len(long)), vertices=int(len(set((long // 15).tolist()))))
# flat starts near the crossing (k = 0.2425 to 0.2575, a = 0.005): 100-sweep means of the longest squared edge and of R, within each chain
src = open('/home/claude/sa74/regge/r2/night/reweight.py').read().split('chains = []')[0]; exec(src)   # the list `runs` of flat-start chain files
byL = {}
for L, k, name, fn in runs:
    if not (0.2425 <= k <= 0.2575) or L < 4: continue
    a = np.loadtxt(fn, comments='#')
    if a.ndim != 2 or a.shape[1] < 28: continue
    m = a[a[:, 1] == 1]
    if len(m) < 800: continue
    n = len(m) // 100; byL.setdefault(L, []).append((name, k, m[:n * 100, 4].reshape(n, 100).mean(1), m[:n * 100, 9].reshape(n, 100).mean(1)))
out['flat_correlation'] = {}
for L, cs in sorted(byL.items()):
    x = np.concatenate([q - q.mean() for _, _, r, q in cs]); y = np.concatenate([r - r.mean() for _, _, r, q in cs])
    out['flat_correlation'][str(L)] = dict(chains=len(cs), blocks=int(len(x)), corr=float(np.corrcoef(x, y)[0, 1]), dR_per_maxq=float((x * y).sum() / (x * x).sum()))
    print(f"flat L={L:2d}: corr(longest edge, R) of 100-sweep means {out['flat_correlation'][str(L)]['corr']:+.2f}, dR per unit {out['flat_correlation'][str(L)]['dR_per_maxq']:+.4f} ({len(cs)} chains)")
json.dump(out, open('/home/claude/sa74/regge/r2/night/defects.json', 'w'), indent=1)
for n, s in out['series'].items(): print(f"{n:22s} rows {s['rows']:6d} last third: R {s['R_last_third']:+.4f} maxq {s['maxq_last_third']:6.1f} | maxq by twelfths", ' '.join(f'{v:.0f}' for v in s['maxq']))
for n, c in out['configs'].items(): print(f"cfg {n:26s} max q {c['max_q']:5.1f} edges > 4x mean: {c['n_long']:3d} on {c['vertices']:3d} vertices")
