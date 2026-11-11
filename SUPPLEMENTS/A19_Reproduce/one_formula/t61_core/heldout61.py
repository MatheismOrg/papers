"""PROTOCOL_T61, marks (b) and (c): with each partition's A, B, C, t, kappa, the held-out disks' dropped points (velocity
error above 10 percent) and all points, against T48 (the same A, B, C), MOND's relation and standard halos. Each disk is
held out once per fold assignment; its three mean squared errors are averaged. Host-balanced rms, paired host bootstrap.
Usage: PYTHONPATH=<dmfit> python3 heldout61.py <res61_MOGF0IESBC_primary.json>"""
import sys, json, numpy as np
import harness_T61 as M
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B
recs = [r for r in json.load(open(sys.argv[1]))['records'] if r['lens'] == 0]; assert len(recs) == 15
A = M.ALL; acc = {}
for r in recs:
    names = list(r['losses']['disk']); a = M.pick(A, np.isin(A['name'], names)); th = np.array(r['fit']['th']); B.FROZEN['abc'] = np.array(r['fit']['abc'])
    M.U['t'], M.U['k'] = r['fit']['t'], r['fit']['kappa']
    e = dict(T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)), T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
    V = np.sqrt(a['go'] * a['r'] / H.ACC); sets = dict(kept=a['fe'] <= .1, dropped=a['fe'] > .1, all=np.ones(len(V), bool), slow_dropped=(a['fe'] > .1) & (V < 30))
    for s, k in sets.items():
        for n in np.unique(a['name'][k]):
            kk = k & (a['name'] == n)
            for m, v in e.items(): acc.setdefault((s, m, n), []).append(float(np.mean(v[kk] ** 2)))
            acc.setdefault((s, 'bias_T61', n), []).append(float(np.mean(e['T61'][kk]))); acc.setdefault((s, 'bias_MOND', n), []).append(float(np.mean(e['MOND'][kk])))
rng = np.random.default_rng(61); out = dict(fits=dict(t=[round(r['fit']['t'], 4) for r in recs], kappa=[round(r['fit']['kappa'], 4) for r in recs], bound_hits=sum(len(r['fit']['bound_hits']) for r in recs)))
for s in ('kept', 'dropped', 'all', 'slow_dropped'):
    hosts = sorted({n for (ss, m, n) in acc if ss == s and m == 'T61'}); assert all(len(acc[(s, 'T61', n)]) == 3 for n in hosts)
    f = {m: np.array([np.mean(acc[(s, m, n)]) for n in hosts]) for m in ('T61', 'T48', 'MOND', 'STD')}; idx = rng.integers(0, len(hosts), (10000, len(hosts)))
    o = dict(galaxies=len(hosts), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()}, mean_T61=round(float(np.mean([np.mean(acc[(s, 'bias_T61', n)]) for n in hosts])), 4), mean_MOND=round(float(np.mean([np.mean(acc[(s, 'bias_MOND', n)]) for n in hosts])), 4))
    for ref in ('T48', 'MOND', 'STD'):
        dd = np.sqrt(f['T61'][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
        o['T61_minus_' + ref] = dict(diff=round(float(np.sqrt(f['T61'].mean()) - np.sqrt(f[ref].mean())), 4), pct=round(float(100 * (np.sqrt(f['T61'].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(),
                                     p_not_lower=round(float(np.mean(dd >= 0)), 4), n_lower=int(np.sum(f['T61'] < f[ref])), n_higher=int(np.sum(f['T61'] > f[ref])))
    dd = np.sqrt(f['T48'][idx].mean(1)) - np.sqrt(f['MOND'][idx].mean(1)); o['T48_minus_MOND_boot95'] = np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist()
    out[s] = o; print(s, json.dumps(o), flush=True)
json.dump(out, open('res61_heldout.json', 'w'), indent=1)
