"""T56 against T48, held out (PROTOCOL_T56 pass marks). Usage: python3 analyze56.py (from the dmfit folder)"""
import json, numpy as np
TOL = {'disk': .01, 'dwarfs': .02, 'ETG': .05, 'MaNGA': .05, 'CLASH_centre': .05, 'CLASH_outer': .05}
rng = np.random.default_rng(56); out = {}
for mask in ('primary', 'strict'):
    a = json.load(open(f'res48_MOGF0IESB_{mask}.json')); b = json.load(open(f'res56_MOGF0IESBU_{mask}.json'))
    def pooled(recs, p):
        d = {}
        for r in recs:
            for h, v in r['losses'][p].items(): d.setdefault(h, []).append(v)
        return {h: float(np.mean(v)) for h, v in d.items()}
    row = {}; ok = True
    for p, tol in TOL.items():
        x, y = pooled(a['records'], p), pooled(b['records'], p); hs = sorted(x); xv = np.array([x[h] for h in hs]); yv = np.array([y[h] for h in hs])
        s0, s1 = np.sqrt(xv.mean()), np.sqrt(yv.mean()); idx = rng.integers(0, len(hs), (4000, len(hs))); dd = np.sqrt(yv[idx].mean(1)) - np.sqrt(xv[idx].mean(1))
        row[p] = dict(T48=round(float(s0), 4), T56=round(float(s1), 4), change_pct=round(float(100 * (s1 / s0 - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), within_tol=bool(s1 <= s0 * (1 + tol)))
        ok &= row[p]['within_tol']
    l0 = sum(r['test_lensing_chi2'] for r in a['records']) / 15; l1 = sum(r['test_lensing_chi2'] for r in b['records']) / 15
    row['lensing'] = dict(T48=round(l0, 2), T56=round(l1, 2), change_pct=round(100 * (l1 / l0 - 1), 2), within_tol=bool(l1 <= l0 * 1.05)); ok &= row['lensing']['within_tol']
    mus = [r['fit']['mu'] for r in b['records']]; abc = np.array([r['fit']['abc'] for r in b['records']]); bh = sum(bool(r['fit']['bound_hits']) for r in b['records'])
    row['fit'] = dict(mu_median=round(float(np.median(mus)), 2), mu_range=[round(min(mus), 2), round(max(mus), 2)], B_range=[round(float(abc[:, 1].min()), 2), round(float(abc[:, 1].max()), 2)],
                      C_range=[round(float(abc[:, 2].min()), 2), round(float(abc[:, 2].max()), 2)], logMt_median=round(float(np.median([r['fit']['th'][0] for r in b['records']])), 3), records_with_bound_hits=bh)
    row['all_within_tolerance'] = bool(ok); row['dwarfs_better_boot'] = bool(row['dwarfs']['boot95'][1] < 0); row['PASS'] = bool(ok and row['dwarfs_better_boot'])
    out[mask] = row; print(mask); [print('  ', k, v) for k, v in row.items()]
out['PASS_both_masks'] = bool(out['primary']['PASS'] and out['strict']['PASS']); print('PASS on both masks:', out['PASS_both_masks'])
json.dump(out, open('res56_analysis.json', 'w'), indent=1)
