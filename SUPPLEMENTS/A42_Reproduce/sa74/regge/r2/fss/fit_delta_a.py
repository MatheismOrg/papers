"""Night stage AS2: average curvature against k on 4^4 for a = 0.001 and 0.002 (flat starts), with the a = 0.005 4^4 points
for comparison; fit Hamber's form R = -A (k_c - k)^delta (Hamber's definition, R scaled by the mean squared edge length) to the
negative points nearest the crossing, and report the zero crossing by linear interpolation."""
import json, glob, re, numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
import sys
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import batches, ess
RUNS = Path('/home/claude/sa74/regge/r2/audit_a005/runs')
def pts(a):
    out = []
    for fn in sorted(RUNS.glob(f'AS2_a{a}_L4_k*_flat.dat')):
        if '# done' not in fn.read_text()[-20000:]: continue
        k = float(re.search(r'_k([0-9.]+)_flat', fn.name).group(1)); d = np.loadtxt(fn); m = d[d[:, 1] == 1]; x = m[:, 4]
        RH = m[:, 5].mean() * 2 * m[:, 3].mean() / m[:, 2].mean(); se = max(b['se'] for b in batches(x)) * m[:, 5].mean()
        out.append(dict(k=k, R=float(x.mean()), RH=float(RH), RH_se=float(se), ess=float(ess(x)), chi=float(4 * m[:, 3].var() / m[:, 2].mean())))
    return sorted(out, key=lambda r: r['k'])
def f(k, A, kc, d): return -A * np.clip(kc - k, 1e-9, None) ** d
res = {}
for a in ('0.001', '0.002'):
    P = pts(a); res[a] = dict(points=P)
    neg = [p for p in P if p['RH'] < 0]; pos = [p for p in P if p['RH'] >= 0]
    if neg and pos:
        p1, p2 = neg[-1], min(pos, key=lambda p: p['k'])
        res[a]['crossing'] = p1['k'] + (p2['k'] - p1['k']) * (-p1['RH']) / (p2['RH'] - p1['RH'])
    win = [p for p in neg if p['k'] >= neg[-1]['k'] - 0.061]
    if len(win) >= 4:
        k = np.array([p['k'] for p in win]); y = np.array([p['RH'] for p in win]); e = np.array([max(p['RH_se'], 1e-3) for p in win])
        lo = k.max() + 2e-4
        p0, _ = curve_fit(f, k, y, p0=[3, lo + 0.01, 0.6], sigma=e, absolute_sigma=True, bounds=([0, lo, 0.05], [100, 1.0, 3]), maxfev=20000)
        rng = np.random.default_rng(1); bs = []
        for _ in range(300):
            try: bs.append(curve_fit(f, k, y + rng.normal(0, e), p0=p0, sigma=e, absolute_sigma=True, bounds=([0, lo, 0.05], [100, 1.0, 3]), maxfev=20000)[0])
            except Exception: pass
        bs = np.array(bs)
        res[a]['fit'] = dict(ks=k.tolist(), A=float(p0[0]), kc=float(p0[1]), delta=float(p0[2]), kc_sd=float(bs[:, 1].std()), delta_sd=float(bs[:, 2].std()),
                             chi2=float((((y - f(k, *p0)) / e) ** 2).sum()), dof=len(k) - 3)
json.dump(res, open('/home/claude/sa74/regge/r2/fss/fit_delta_a.json', 'w'), indent=1)
for a, v in res.items():
    print('a =', a, 'crossing', round(v.get('crossing', float('nan')), 4), 'fit', {k: (round(x, 4) if isinstance(x, float) else x) for k, x in v.get('fit', {}).items() if k != 'ks'})
    for p in v['points']: print(f"   k={p['k']:.3f} R={p['R']:8.4f} RH={p['RH']:8.4f} ess={p['ess']:.0f} chi={p['chi']:.2f}")
