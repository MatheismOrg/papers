"""SA73 T61, after the runs (not in the protocol; a check in MOND's favor): the outside checks (d) and (e) against MOND's
relation with its constant refitted on all points of the 149 benchmark disks with T61's error weights.
Usage: PYTHONPATH=<dmfit>:<t59 folder> python3 outside61m.py <t59 folder>"""
import sys, json, numpy as np
from scipy.optimize import minimize_scalar
src = open('outside61.py').read().split("\nrng = np.random.default_rng(611)")[0]; exec(compile(src, 'outside61_defs', 'exec'))
bench = [g['name'] for g in H.disks()]; tr = M.pick(M.ALL, np.isin(M.ALL['name'], bench)); nw = M.weights(tr)
def g_mond(a, a0): return a['gb'] * C.nu(a['gb'] / a0)
o = minimize_scalar(lambda la: float(np.sum((nw * np.log10(tr['go'] / g_mond(tr, 10 ** la))) ** 2)), bounds=(-10.6, -9.4), method='bounded'); a0 = 10 ** o.x; print('a0 refitted on all points of the 149 disks: %.4g' % a0)
rng = np.random.default_rng(613); out = dict(a0=a0); B.FROZEN['abc'] = ABC; M.U['t'], M.U['k'] = TK
def score(a, lab):
    e = dict(T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), MOND_refit=np.log10(a['go'] / g_mond(a, a0)))
    un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi); f = {m: np.bincount(gi, v ** 2) / ct for m, v in e.items()}; idx = rng.integers(0, len(un), (10000, len(un)))
    dd = np.sqrt(f['T61'][idx].mean(1)) - np.sqrt(f['MOND_refit'][idx].mean(1))
    r = dict(galaxies=len(un), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()}, T61_minus_MOND_refit_boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f['T61'] < f['MOND_refit'])))
    print(lab, json.dumps(r)); return r
test = [n for n, g in gal.items() if g['sparc'] is None]
out['d'] = score(rows(test, 9.), '(d) Iorio outside, every point'); A = M.ALL; out['e'] = score(M.pick(A, ~np.isin(A['name'], bench)), '(e) SPARC outside the benchmark, every point')
json.dump(out, open('res61_outside_mond_refit.json', 'w'), indent=1)
