"""T57 against T48, held out (PROTOCOL_T57 pass marks). Usage: python3 analyze57.py (from the dmfit folder)"""
import json, numpy as np
TOL = {'disk': .01, 'dwarfs': .02, 'ETG': .05, 'MaNGA': .05, 'CLASH_centre': .05, 'CLASH_outer': .05}
out = {}
for cfg in ('MOGF0IESBD', 'MOGF0IESBD5'):
    for mask in ('primary', 'strict'):
        rng = np.random.default_rng(57)
        a = json.load(open(f'res48_MOGF0IESB_{mask}.json')); b = json.load(open(f'res57_{cfg}_{mask}.json'))
        def pooled(recs, p):
            d = {}
            for r in recs:
                for h, v in r['losses'][p].items(): d.setdefault(h, []).append(v)
            return {h: float(np.mean(v)) for h, v in d.items()}
        row = {}; ok = True
        for p, tol in TOL.items():
            x, y = pooled(a['records'], p), pooled(b['records'], p); hs = sorted(x); xv = np.array([x[h] for h in hs]); yv = np.array([y[h] for h in hs])
            s0, s1 = np.sqrt(xv.mean()), np.sqrt(yv.mean()); idx = rng.integers(0, len(hs), (10000, len(hs))); dd = np.sqrt(yv[idx].mean(1)) - np.sqrt(xv[idx].mean(1))
            row[p] = dict(T48=round(float(s0), 4), T57=round(float(s1), 4), pct=round(float(100 * (s1 / s0 - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), p_not_better=round(float(np.mean(dd >= 0)), 4), within_tol=bool(s1 <= s0 * (1 + tol)))
            ok &= row[p]['within_tol']
            if p == 'disk': row['disk_n_better'] = [int(np.sum(yv < xv)), len(hs)]
        l0 = sum(r['test_lensing_chi2'] for r in a['records']) / 15; l1 = sum(r['test_lensing_chi2'] for r in b['records']) / 15
        row['lensing'] = dict(T48=round(l0, 2), T57=round(l1, 2), pct=round(100 * (l1 / l0 - 1), 2), within_tol=bool(l1 <= l0 * 1.05)); ok &= row['lensing']['within_tol']
        ts = np.array([[r['fit']['t'], r['fit']['s']] for r in b['records']]); abc = np.array([r['fit']['abc'] for r in b['records']])
        row['fit'] = dict(t_median=round(float(np.median(ts[:, 0])), 3), t_range=np.round([ts[:, 0].min(), ts[:, 0].max()], 3).tolist(), s_median=round(float(np.median(ts[:, 1])), 3), s_range=np.round([ts[:, 1].min(), ts[:, 1].max()], 3).tolist(),
                          C_range=np.round([abc[:, 2].min(), abc[:, 2].max()], 2).tolist(), bound_hits=sum(bool(r['fit']['bound_hits']) for r in b['records']))
        row['all_within_tolerance'] = bool(ok); d = row['disk']; row['disk_gain_beyond_1pct'] = bool(d['pct'] < -1.); row['disk_boot_below_zero'] = bool(d['boot95'][1] < 0)
        row['PASS'] = bool(ok and row['disk_gain_beyond_1pct'] and row['disk_boot_below_zero']); out[f'{cfg}|{mask}'] = row
        print(cfg, mask); [print('  ', k, v) for k, v in row.items()]
json.dump(out, open('res57_analysis.json', 'w'), indent=1)
