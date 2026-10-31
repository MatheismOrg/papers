"""A19 supplement: post-run predictive checks of protocol 4's 85 dwarfs, added in turn 13 after an outside review.
Not a protocol and not a blind test: the dwarfs, the models and protocol 4's likelihood are those already scored.
(1) Leave one dwarf out: each model's scatter is fitted on the other 84 dwarfs with protocol 4's likelihood
    (a measured value with its error, an upper limit as a ceiling) and the held-out dwarf is scored; the 85 scores add.
(2) The same with an offset fitted as well as the scatter, which lets every model shift its predictions.
(3) A bootstrap over dwarfs of the summed held-out differences (it resamples the saved scores, so it leaves out
    the overlap of the training sets and any error shared by dwarfs of one host).
(4) Root-mean-square errors on the 17 dwarfs of protocol 2 scored blind under protocol 1.
The outside review's own script (a separate likelihood code and optimizer) gave 35.71, 8.81, 21.8 to 47.8 and
0.4047, 0.7762, 0.3974 for these; this script checks them with protocol 4's own functions.
Usage: python3 A19_Predictive_Audit.py dwarf_all.csv blind_test_dir out.json"""
import sys, os, json, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_Combined_v4 as V4

def fit_offset(res, e, ul):
    f = lambda p: -np.sum(V4.lnl_terms(res, e, ul, min(max(p[1], 0.01), 3.0), p[0]))
    best = min((minimize(f, x0, method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-10, 'maxiter': 6000})
                for x0 in ([0.0, 0.3], [0.4, 0.6], [-0.4, 0.6])), key=lambda o: o.fun)
    return float(best.x[0]), float(min(max(best.x[1], 0.01), 3.0))

def main(lvdb, bt, out):
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    rows, src = V4.load_rows(lvdb, bt); A = V4.arrays(rows, rules); n = len(rows)
    R = {m: A['lo'] - A['logg'][m] for m in V4.MODELS}; e, ul = A['e'], A['ul']
    res = dict(n=n, n_measured=int((~ul).sum()), n_limits=int(ul.sum()), status='post-run check, not a protocol or a blind test',
               in_sample={}, loo_scatter={}, loo_offset_and_scatter={})
    per = {}
    for m in V4.MODELS:
        s, L = V4.fit_s(R[m], e, ul); res['in_sample'][m] = dict(s_dex=round(s, 4), lnL=round(L, 2))
        a, b, ss, oo = [], [], [], []
        for j in range(n):
            keep = np.arange(n) != j
            s1, _ = V4.fit_s(R[m], e, ul, keep); ss.append(s1)
            a.append(float(V4.lnl_terms(R[m][j:j + 1], e[j:j + 1], ul[j:j + 1], s1)[0]))
            o1, s2 = fit_offset(R[m][keep], e[keep], ul[keep]); oo.append(o1)
            b.append(float(V4.lnl_terms(R[m][j:j + 1], e[j:j + 1], ul[j:j + 1], s2, o1)[0]))
        per[m] = dict(scatter=np.array(a), offset=np.array(b))
        res['loo_scatter'][m] = dict(score=round(float(np.sum(a)), 2), scatter_range=[round(min(ss), 4), round(max(ss), 4)])
        res['loo_offset_and_scatter'][m] = dict(score=round(float(np.sum(b)), 2), offset_range=[round(min(oo), 3), round(max(oo), 3)])
    rng = np.random.default_rng(20261031); ii = rng.integers(0, n, (20000, n))
    for kind in ('scatter', 'offset'):
        key = 'loo_scatter' if kind == 'scatter' else 'loo_offset_and_scatter'
        res[key + '_stellar_minus'] = {}
        for m in V4.MODELS:
            if m == 'stellar_rule': continue
            d = per['stellar_rule'][kind] - per[m][kind]
            res[key + '_stellar_minus'][m] = dict(difference=round(float(d.sum()), 2),
                                                 bootstrap_95=[round(float(x), 2) for x in np.percentile(d[ii].sum(axis=1), [2.5, 97.5])],
                                                 dwarfs_favoring_stellar=int(np.sum(d > 0)))
    j2 = json.load(open(os.path.join(bt, 'A19_OutOfSample_v2_Results.json'))); k17 = set(j2['blind_v1_subset']['keys'])
    q = np.array([x['key'] in k17 for x in rows]); assert q.sum() == 17
    res['blind_17_rmse_dex'] = {m: round(float(np.sqrt(np.mean(R[m][q] ** 2))), 4) for m in V4.MODELS}
    res['per_dwarf_loo_scatter'] = {x['key']: {m: round(float(per[m]['scatter'][i]), 3) for m in V4.MODELS} for i, x in enumerate(rows)}
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'per_dwarf_loo_scatter'}, indent=1))

if __name__ == '__main__':
    main(*sys.argv[1:4])
