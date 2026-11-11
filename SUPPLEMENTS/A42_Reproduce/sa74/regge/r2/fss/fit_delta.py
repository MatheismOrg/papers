"""Fit Hamber's smooth-phase form R(k) = -A (k_c - k)^delta to the pooled curvature means (a = 0.005), all sizes together
and per size, using 0.15 <= k <= 0.24. Errors by a parametric bootstrap over the point errors. Compare: Hamber (1993),
a = 0.005, lattices 4^4 to 16^4: k_c = 0.2443(11), delta = 0.626(25)."""
import json, numpy as np
from scipy.optimize import curve_fit
S = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results.json'))['summary']
def f(k, A, kc, d): return -A * np.clip(kc - k, 1e-9, None) ** d
def fit(rows, nboot=400, seed=1, key='R'):
    k = np.array([r['k'] for r in rows]); y = np.array([r[key] for r in rows]); e = np.array([max(r[key + '_se'], 1e-3) for r in rows])
    p, _ = curve_fit(f, k, y, p0=[0.5, 0.25, 0.6], sigma=e, absolute_sigma=True, bounds=([0, 0.2402, 0.05], [50, 0.40, 3.0]), maxfev=20000)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(nboot):
        try:
            q, _ = curve_fit(f, k, y + rng.normal(0, e), p0=p, sigma=e, absolute_sigma=True, bounds=([0, 0.2402, 0.05], [50, 0.40, 3.0]), maxfev=20000)
            bs.append(q)
        except Exception: pass
    bs = np.array(bs); chi2 = float((((y - f(k, *p)) / e) ** 2).sum())
    return dict(n=len(k), A=float(p[0]), kc=float(p[1]), delta=float(p[2]), kc_sd=float(bs[:, 1].std()), delta_sd=float(bs[:, 2].std()),
                chi2=chi2, dof=len(k) - 3, ks=k.tolist())
sel = [r for r in S if 0.15 <= r['k'] <= 0.24 + 1e-9]
out = {'all sizes, R = <2 sum dA / sum V>': fit(sel), 'all sizes, Hamber R = <l^2><2 sum dA>/<sum V>': fit(sel, key='RH')}
for L in sorted(set(r['L'] for r in sel)):
    rows = [r for r in sel if r['L'] == L]
    if len(rows) >= 4: out[f'L={L}, Hamber R'] = fit(rows, key='RH')
json.dump(out, open('/home/claude/sa74/regge/r2/fss/fit_delta.json', 'w'), indent=1)
for k, v in out.items():
    print(f"{k:46s} n={v['n']} k_c={v['kc']:.4f}+-{v['kc_sd']:.4f} delta={v['delta']:.3f}+-{v['delta_sd']:.3f} chi2/dof={v['chi2']:.1f}/{v['dof']}")
