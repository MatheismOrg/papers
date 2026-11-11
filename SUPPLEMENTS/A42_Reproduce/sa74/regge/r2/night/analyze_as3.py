"""SA74 night stage AS3 analysis (a = 0.001, the smallest curvature-squared coupling tested that still mixes).
(1) Start independence at k = 0.14 on 4^4: three flat and three rough starts (the a = 0, k = 0.1 preparations); split R-hat of
    the curvature and of the volume, screen 1.05 as in T19/T25, and flat against rough means.
(2) Average curvature against k on 4^4 (AS2 + AS3) and 6^4 (AS3): the zero crossing by linear interpolation, and Hamber's form
    R = -A (k_c - k)^delta (Hamber's R, scaled by the mean squared edge length) fitted to the negative points nearest the crossing.
Statistics functions from analyze_r2 (copied from the T19 analysis)."""
import json, re, sys, numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import batches, ess, rhat
RUNS = Path('/home/claude/sa74/regge/r2/audit_a005/runs')
def chain(fn):
    d = np.loadtxt(fn); m = d[d[:, 1] == 1]; x = m[:, 4]
    RH = m[:, 5].mean() * 2 * m[:, 3].mean() / m[:, 2].mean()
    return dict(name=fn.stem, n=len(x), R=float(x.mean()), se=float(max(b['se'] for b in batches(x))), ess=float(ess(x)),
                RH=float(RH), RH_se=float(max(b['se'] for b in batches(x)) * m[:, 5].mean()), chi=float(4 * m[:, 3].var() / m[:, 2].mean()),
                R_first=float(x[:len(x) // 2].mean()), R_last=float(x[len(x) // 2:].mean())), m
out = {}
# (1) start independence at k = 0.14, 4^4
grp = {}
for fn in sorted(RUNS.glob('AS3_a0.001_L4_k0.140_*_s*.dat')):
    c, m = chain(fn); st = 'rough' if '_rough_' in fn.name else 'flat'; grp.setdefault(st, []).append((c, m))
allm = [m for st in grp for c, m in grp[st]]
out['start_independence'] = dict(k=0.14, L=4, a=0.001,
    rhat_R=rhat([m[:, 4] for m in allm]), rhat_V=rhat([m[:, 2] for m in allm]),
    flat=[c for c, m in grp['flat']], rough=[c for c, m in grp['rough']])
fm = np.array([c['R'] for c, m in grp['flat']]); rm = np.array([c['R'] for c, m in grp['rough']])
out['start_independence']['flat_mean'] = float(fm.mean()); out['start_independence']['rough_mean'] = float(rm.mean())
out['start_independence']['z'] = float((fm.mean() - rm.mean()) / np.sqrt(fm.var(ddof=1) / len(fm) + rm.var(ddof=1) / len(rm)))
# (2) k-scans
def f(k, A, kc, d): return -A * np.clip(kc - k, 1e-9, None) ** d
for L, pats in ((4, ['AS2_a0.001_L4_k*_flat.dat', 'AS3_a0.001_L4_k0.1[35]0_flat.dat', 'AS3_a0.001_L4_k0.140_flat_s*.dat']),
                (6, ['AS3_a0.001_L6_k*_flat.dat'])):
    pts = {}
    for pat in pats:
        for fn in sorted(RUNS.glob(pat)):
            k = round(float(re.search(r'_k([0-9.]+)_', fn.name).group(1)), 3); c, m = chain(fn); pts.setdefault(k, []).append(c)
    P = []
    for k, cs in sorted(pts.items()):
        if len(cs) > 2:
            RH = np.mean([c['RH'] for c in cs]); se = np.std([c['RH'] for c in cs], ddof=1) / np.sqrt(len(cs))
        else:
            RH = cs[0]['RH']; se = cs[0]['RH_se']
        P.append(dict(k=k, chains=len(cs), RH=float(RH), RH_se=float(se), R=float(np.mean([c['R'] for c in cs])),
                      ess_min=float(min(c['ess'] for c in cs)), chi=float(np.mean([c['chi'] for c in cs])),
                      drift=float(np.mean([c['R_last'] - c['R_first'] for c in cs]))))
    res = dict(points=P)
    neg = [p for p in P if p['RH'] < 0]; pos = [p for p in P if p['RH'] >= 0]
    if neg and pos:
        p1 = max(neg, key=lambda p: p['k']); p2 = min(pos, key=lambda p: p['k'])
        res['crossing'] = p1['k'] + (p2['k'] - p1['k']) * (-p1['RH']) / (p2['RH'] - p1['RH'])
    for width in (0.041, 0.061):
        win = [p for p in neg if p['k'] >= max(q['k'] for q in neg) - width]
        if len(win) < 4: continue
        k = np.array([p['k'] for p in win]); y = np.array([p['RH'] for p in win]); e = np.array([max(p['RH_se'], 1e-3) for p in win])
        lo = k.max() + 2e-4
        try:
            p0, _ = curve_fit(f, k, y, p0=[3, lo + 0.01, 0.6], sigma=e, absolute_sigma=True, bounds=([0, lo, 0.05], [100, 1.0, 3]), maxfev=20000)
        except Exception as ex:
            res[f'fit_w{width}'] = dict(error=str(ex)); continue
        rng = np.random.default_rng(1); bs = []
        for _ in range(300):
            try: bs.append(curve_fit(f, k, y + rng.normal(0, e), p0=p0, sigma=e, absolute_sigma=True, bounds=([0, lo, 0.05], [100, 1.0, 3]), maxfev=20000)[0])
            except Exception: pass
        bs = np.array(bs)
        res[f'fit_w{width}'] = dict(ks=k.tolist(), A=float(p0[0]), kc=float(p0[1]), delta=float(p0[2]), kc_sd=float(bs[:, 1].std()),
                                    delta_sd=float(bs[:, 2].std()), chi2=float((((y - f(k, *p0)) / e) ** 2).sum()), dof=len(k) - 3)
    out[f'L{L}'] = res
json.dump(out, open('/home/claude/sa74/regge/r2/night/as3_results.json', 'w'), indent=1)
si = out['start_independence']
print(f"k=0.14 4^4 a=0.001: R-hat(R) {si['rhat_R']:.3f}  R-hat(V) {si['rhat_V']:.3f}  flat {si['flat_mean']:+.4f}  rough {si['rough_mean']:+.4f}  z {si['z']:+.2f}")
for st in ('flat', 'rough'):
    for c in si[st]: print(f"   {c['name']:<34} R={c['R']:+.4f}±{c['se']:.4f} ess={c['ess']:.0f} first={c['R_first']:+.4f} last={c['R_last']:+.4f}")
for L in (4, 6):
    r = out[f'L{L}']; print(f"L={L}: crossing {r.get('crossing', float('nan')):.4f}")
    for w in ('fit_w0.041', 'fit_w0.061'):
        if w in r: print('   ', w, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r[w].items() if k != 'ks'})
    for p in r['points']: print(f"   k={p['k']:.3f} n={p['chains']} RH={p['RH']:+.4f}±{p['RH_se']:.4f} R={p['R']:+.4f} ess_min={p['ess_min']:.0f} chi={p['chi']:.2f} drift={p['drift']:+.4f}")
