"""SA73 T61, after the held-out run (not in the protocol; a check in MOND's favor): MOND's relation with its constant a0
refitted on each partition's training disks, on all their points with the same error weights as T61's two numbers, then
scored on the held-out disks as heldout61.py scores T61. Usage: PYTHONPATH=<dmfit> python3 heldout61m.py <res61 primary json>"""
import sys, json, numpy as np
from scipy.optimize import minimize_scalar
import harness_T61 as M
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B
recs = [r for r in json.load(open(sys.argv[1]))['records'] if r['lens'] == 0]; A = M.ALL; bench = sorted({n for r in recs for n in r['losses']['disk']}); assert len(bench) == 149
def g_mond(a, a0): return a['gb'] * C.nu(a['gb'] / a0)
acc = {}; a0s = []
for r in recs:
    held = list(r['losses']['disk']); tr = M.pick(A, np.isin(A['name'], [n for n in bench if n not in held])); nw = M.weights(tr)
    o = minimize_scalar(lambda la: float(np.sum((nw * np.log10(tr['go'] / g_mond(tr, 10 ** la))) ** 2)), bounds=(-10.6, -9.4), method='bounded'); a0 = 10 ** o.x; a0s.append(a0)
    a = M.pick(A, np.isin(A['name'], held)); th = np.array(r['fit']['th']); B.FROZEN['abc'] = np.array(r['fit']['abc']); M.U['t'], M.U['k'] = r['fit']['t'], r['fit']['kappa']
    e = dict(T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), MOND_refit=np.log10(a['go'] / g_mond(a, a0)))
    for s, k in dict(kept=a['fe'] <= .1, dropped=a['fe'] > .1, all=np.ones(len(a['r']), bool)).items():
        for n in np.unique(a['name'][k]):
            kk = k & (a['name'] == n)
            for m, v in e.items(): acc.setdefault((s, m, n), []).append(float(np.mean(v[kk] ** 2)))
rng = np.random.default_rng(612); out = dict(a0_range=[float(min(a0s)), float(max(a0s))], a0_median=float(np.median(a0s)))
print('a0 refitted: %.3g to %.3g (median %.3g)' % (min(a0s), max(a0s), np.median(a0s)))
for s in ('kept', 'dropped', 'all'):
    hosts = sorted({n for (ss, m, n) in acc if ss == s and m == 'T61'}); f = {m: np.array([np.mean(acc[(s, m, n)]) for n in hosts]) for m in ('T61', 'MOND', 'MOND_refit')}
    idx = rng.integers(0, len(hosts), (10000, len(hosts))); dd = np.sqrt(f['T61'][idx].mean(1)) - np.sqrt(f['MOND_refit'][idx].mean(1))
    out[s] = dict(galaxies=len(hosts), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()}, T61_minus_MOND_refit_boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), p_not_lower=round(float(np.mean(dd >= 0)), 4))
    print(s, json.dumps(out[s]))
json.dump(out, open('res61_heldout_mond_refit.json', 'w'), indent=1)
