"""SA73 T56 exploration (no protocol): cold gas weighted mu times in the baryon term, mu fitted on the disks with A, B, C.
M_ref = (M* + mu M_gas)/f_b + 10^A (M*/1e10)^B (R_e/3)^C (+ q M_E, floor, as the final formula). (a) scan of mu on all disks;
(b) the benchmark's 15 disk training folds: fit on the training disks, score the held-out disks and, never fitted, the dwarfs
and early types. Usage: python3 mu56.py <dmfit dir>"""
import sys, json, csv; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows())}
d = P['disk']; dw = P['dwarfs']
def gpred(a, mu):
    isd = a['kind'][0] == 'disk'; mp = a['ms'] + mu * (a['mp'] - a['ms']) if isd else a['mp']
    m = G43.mass(th, a['ms'], mp, a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.); rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def hbsel(a, e, sel=None):
    u, i = np.unique(a['name'], return_inverse=True); ms = np.bincount(i, e ** 2) / np.bincount(i)
    if sel is not None: ms = ms[np.isin(u, sel)]
    return float(np.sqrt(np.mean(ms))), ms
def dwslope(mu):
    e = np.log10(dw['go'] / gpred(dw, mu)); X = np.c_[np.ones(len(e)), np.log10(dw['ms']), np.log10(dw['re'])]; return float(np.linalg.lstsq(X, e, rcond=None)[0][2]), float(np.mean(e))
def fit(train, free_mu, mu_fixed=1.):
    k = np.isin(d['name'], train); un, inv, ct = np.unique(d['name'][k], return_inverse=True, return_counts=True); hw = 1 / np.sqrt(len(un) * ct[inv])
    sub = {n: (v[k] if isinstance(v, np.ndarray) and len(v) == len(k) else v) for n, v in d.items()}
    def f(p):
        B.FROZEN['abc'] = np.asarray(p[:3]); return np.log10(sub['go'] / gpred(sub, abs(p[3]) if free_mu else mu_fixed)) * hw
    best = None
    for st in ([1.], [4.]) if free_mu else ([],):
        r = least_squares(f, list(ABC0) + st)
        if best is None or r.cost < best.cost: best = r
    return best.x[:3], (abs(best.x[3]) if free_mu else mu_fixed)
names = np.unique(d['name'])
print('(a) scan on all disks')
for mu in (1., 2., 3., 4., 5., 6., 8., 12.):
    abc, _ = fit(names, False, mu); B.FROZEN['abc'] = abc; s, mn = dwslope(mu)
    print(f'  mu {mu:4.1f}: disks {hbsel(d, np.log10(d["go"]/gpred(d, mu)))[0]:.4f} dwarfs {hbsel(dw, np.log10(dw["go"]/gpred(dw, mu)))[0]:.4f} (mean {mn:+.3f}, slope {s:+.2f}) ETG {hbsel(P["ETG"], np.log10(P["ETG"]["go"]/gpred(P["ETG"], mu)))[0]:.4f}  A,B,C {np.round(abc, 3).tolist()}', flush=True)
print('(b) the 15 training folds')
folds = list(csv.DictReader(open(H.PKG / 'data/SPARC/A19_Disk_Folds.csv'))); cols = [c for c in folds[0] if c.lower().startswith('fold')]
key = [c for c in folds[0] if c.lower() in ('galaxy', 'name')][0]
acc = {m: dict(disk=[], dw=[], etg=[], mu=[], abc=[]) for m in ('base', 'mu')}
for c in cols:
    for fo in sorted({r[c] for r in folds}):
        test = [r[key] for r in folds if r[c] == fo and r[key] in names]; train = [n for n in names if n not in test]
        for m in ('base', 'mu'):
            abc, mu = fit(train, m == 'mu'); B.FROZEN['abc'] = abc
            acc[m]['disk'].append(hbsel(d, np.log10(d['go'] / gpred(d, mu)), test)[1]); acc[m]['dw'].append(hbsel(dw, np.log10(dw['go'] / gpred(dw, mu)))[1])
            acc[m]['etg'].append(hbsel(P['ETG'], np.log10(P['ETG']['go'] / gpred(P['ETG'], mu)))[1]); acc[m]['mu'].append(mu); acc[m]['abc'].append(abc)
for m in ('base', 'mu'):
    a = acc[m]; print(f"  {m:5s} held-out disks {np.sqrt(np.mean(np.concatenate(a['disk']))):.4f}  dwarfs {np.sqrt(np.mean(np.mean(a['dw'], axis=0))):.4f}  ETG {np.sqrt(np.mean(np.mean(a['etg'], axis=0))):.4f}",
                      ' mu: median %.2f range %.2f to %.2f' % (np.median(a['mu']), min(a['mu']), max(a['mu'])), ' B: %.2f to %.2f  C: %.2f to %.2f' % (min(x[1] for x in a['abc']), max(x[1] for x in a['abc']), min(x[2] for x in a['abc']), max(x[2] for x in a['abc'])))
json.dump({m: dict(mu=[float(x) for x in acc[m]['mu']], abc=[np.round(x, 4).tolist() for x in acc[m]['abc']]) for m in acc}, open('res56_mu_folds.json', 'w'))
B.FROZEN['abc'] = ABC0
