"""SA74 night: single-histogram reweighting in k of every a = 0.005 flat-start chain (Ferrenberg and Swendsen, 1988).
The engine samples exp(-I) with I = lambda sum V - k sum dA + a sum (dA)^2/V, so a chain at k samples k' with weights
exp((k' - k) sum dA) (column 3 = sum dA). For each grid L the fluctuation chi(k') = 4 var(sum dA)/<sum V> is reweighted
from every chain whose reweighting efficiency (sum w)^2 / (n sum w^2) stays above 0.5, averaged over chains, and its
largest value over k' is compared across L. At a critical point that peak grows as L^(2/nu - 4) (about L^0.94 for
delta = 0.62); at a first-order transition as L^4. Measurement sweeps only; the second half of each chain if SKIPFRAC=0.5.
Also checks d<R>/dk from the fluctuations, cov(R, sum dA), against finite differences between neighboring k."""
import json, re, glob, os, sys
from pathlib import Path
import numpy as np
R2 = Path('/home/claude/sa74/regge/r2'); AUD = R2 / 'audit_a005'; MAC = Path('/home/claude/sa74/NIGHT/mac_pull')
SKIP = float(os.environ.get('SKIPFRAC', '0')); MINROWS = int(os.environ.get('MINROWS', '800')); EFF = float(os.environ.get('EFF', '0.5'))
runs = []
# T32: glob returns files in the file system's order, and the bootstrap in multihist.py draws chains in the order they are read, so its errors
# depend on that order. chain_order_T30.txt records the order of the paper's run; files it does not list follow, sorted by name.
_rank = {n: i for i, n in enumerate((AUD / 'chain_order_T30.txt').read_text().split())} if (AUD / 'chain_order_T30.txt').exists() else {}
for fn in sorted(glob.glob(str(AUD / 'runs' / '*.dat')), key=lambda f: (_rank.get(os.path.basename(f), len(_rank)), os.path.basename(f))):
    name = Path(fn).stem
    m = re.match(r'(?:FSS\d?_|HY1_)?L(\d+)_k([0-9.]+)_(flat|rough)', name)
    if not m or name.startswith('prepare') or m.group(3) != 'flat': continue
    if '# done' not in open(fn).read()[-20000:]: continue
    runs.append((int(m.group(1)), float(m.group(2)), name, fn))
for fn in sorted(MAC.glob('[NQ]_L*_k*.dat')):
    m = re.match(r'[NQ]_L(\d+)_k([0-9.]+)(_rough)?(_r\d+)?$', fn.stem)
    if m and not m.group(3): runs.append((int(m.group(1)), float(m.group(2)), 'mac_' + fn.stem, str(fn)))
for fn in sorted(glob.glob(str(R2 / 'scan' / 'L3_k*_flat.dat'))) + sorted(glob.glob(str(R2 / 'scan4' / 'L4_k*_flat.dat'))):   # T32: sorted
    m = re.match(r'L(\d)_k([0-9.]+)_flat', Path(fn).stem); runs.append((int(m.group(1)), float(m.group(2)), 'scan_' + Path(fn).stem, fn))
chains = []
for L, k, name, fn in runs:
    a = np.loadtxt(fn, comments='#')
    if a.ndim != 2 or a.shape[1] < 28: continue
    a = a[a[:, 1] == 1]
    if len(a) < MINROWS: continue
    a = a[int(SKIP * len(a)):]
    chains.append(dict(L=L, k=round(k, 4), name=name, V=a[:, 2], DA=a[:, 3], R=a[:, 4]))
grid = np.round(np.arange(0.2300, 0.2701, 0.0005), 4)
def rew(c, kp):
    x = (kp - c['k']) * (c['DA'] - c['DA'].mean()); w = np.exp(x - x.max()); w /= w.sum()
    eff = 1.0 / (len(w) * (w ** 2).sum())
    mDA = (w * c['DA']).sum(); vDA = (w * (c['DA'] - mDA) ** 2).sum(); mV = (w * c['V']).sum()
    return eff, 4 * vDA / mV, (w * c['R']).sum()
out = dict(description=__doc__, skipfrac=SKIP, eff_min=EFF, per_L={})
for L in sorted(set(c['L'] for c in chains)):
    cs = [c for c in chains if c['L'] == L]
    curve = []
    for kp in grid:
        vals = [rew(c, kp) for c in cs]; ok = [(ch, R) for (e, ch, R) in vals if e >= EFF]
        if ok: curve.append(dict(k=float(kp), chi=float(np.mean([v[0] for v in ok])), R=float(np.mean([v[1] for v in ok])), chains=len(ok)))
    if not curve: continue
    pk = max(curve, key=lambda r: r['chi'])
    # fluctuation slope against finite differences of the measured points
    pts = {}
    for c in cs: pts.setdefault(c['k'], []).append(c)
    ks = sorted(pts)
    slope_fd = [dict(k_mid=(k1 + k2) / 2, fd=float((np.mean([c['R'].mean() for c in pts[k2]]) - np.mean([c['R'].mean() for c in pts[k1]])) / (k2 - k1)),
                     fluct=float(np.mean([np.cov(c['R'], c['DA'])[0, 1] for c in pts[k1] + pts[k2]])))
                for k1, k2 in zip(ks, ks[1:]) if k2 - k1 <= 0.0051]
    out['per_L'][str(L)] = dict(chains=len(cs), k_points=ks, k_covered=[curve[0]['k'], curve[-1]['k']], n_grid=len(curve),
                                peak=pk, curve=curve, slope_check=slope_fd)
Ls = sorted(int(L) for L in out['per_L']); P = [out['per_L'][str(L)]['peak']['chi'] for L in Ls]
sel = [(L, p) for L, p in zip(Ls, P) if out['per_L'][str(L)]['k_covered'][0] <= 0.245 and out['per_L'][str(L)]['k_covered'][1] >= 0.2525]
if len(sel) >= 3:
    x = np.log([s[0] for s in sel]); y = np.log([s[1] for s in sel]); A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ coef; s2 = (r ** 2).sum() / max(1, len(x) - 2)
    out['peak_power'] = dict(L=[s[0] for s in sel], power=float(coef[0]), se=float(np.sqrt(s2 * np.linalg.inv(A.T @ A)[0, 0])))
fn = R2 / 'night' / ('reweight_skip%02d.json' % int(100 * SKIP)); fn.write_text(json.dumps(out, indent=1) + '\n')
for L in Ls:
    d = out['per_L'][str(L)]; pk = d['peak']
    print(f"L={L:2d} chains={d['chains']:2d} k covered {d['k_covered'][0]:.4f}-{d['k_covered'][1]:.4f} ({d['n_grid']} grid pts)  peak chi={pk['chi']:.2f} at k={pk['k']:.4f} (R there {pk['R']:+.4f}, {pk['chains']} chains)")
    for s in d['slope_check']: print(f"      dR/dk at {s['k_mid']:.5f}: finite difference {s['fd']:6.2f}  fluctuation cov(R, sum dA) {s['fluct']:6.2f}")
if 'peak_power' in out: print('peak chi ~ L^%.2f +- %.2f over L=%s' % (out['peak_power']['power'], out['peak_power']['se'], out['peak_power']['L']))
