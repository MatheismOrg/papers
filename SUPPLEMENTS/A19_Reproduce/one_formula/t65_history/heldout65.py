"""PROTOCOL_T65, Test 3, marks (b) and (c): the core under the gate "any cold gas" (H, harness_T65.py) on the held-out disks'
dropped points and all points, with each partition's A, B, C, t, kappa, against T48, the final formula's gas gate (T61,
that partition's t and kappa) and MOND's relation. Each disk is held out once per fold assignment; its three mean squared
errors are averaged. Host-balanced rms, paired host bootstrap (10,000 draws).
Usage: PYTHONPATH=../dmfit python3 heldout65.py <res65 primary> <res61 primary>"""
import sys, json, numpy as np
import harness_T65 as M5
import harness_T48 as M8
import compare46 as C
H = M5.H; B = M5.B; M61 = M5.M61
def recs(p): return {(r['assignment'], r['fold']): r for r in json.load(open(p))['records'] if r['lens'] == 0}
R5, R61 = recs(sys.argv[1]), recs(sys.argv[2]); assert len(R5) == len(R61) == 15
A = M5.ALL; acc = {}; MODELS = ('H', 'T61', 'T48', 'MOND')
for key in sorted(R5):
    r5, r6 = R5[key], R61[key]; names = list(r5['losses']['disk']); assert set(names) == set(r6['losses']['disk']) and np.allclose(r5['fit']['abc'], r6['fit']['abc'])
    a = M5.pick(A, np.isin(A['name'], names)); B.FROZEN['abc'] = np.array(r5['fit']['abc']); M5.U['t'], M5.U['k'] = r5['fit']['t'], r5['fit']['kappa']; M61.U['t'], M61.U['k'] = r6['fit']['t'], r6['fit']['kappa']
    e = dict(H=np.log10(a['go'] / H.gpred(M5.MODEL, np.array(r5['fit']['th']), a)), T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', np.array(r6['fit']['th']), a)), T48=np.log10(a['go'] / H.gpred(M8.MODEL, np.array(r5['fit']['th']), a)), MOND=np.log10(a['go'] / C.g_mond(a)))
    for s, k in dict(kept=a['fe'] <= .1, dropped=a['fe'] > .1, all=np.ones(len(a['r']), bool)).items():
        for n in np.unique(a['name'][k]):
            kk = k & (a['name'] == n)
            for m, v in e.items(): acc.setdefault((s, m, n), []).append(float(np.mean(v[kk] ** 2)))
rng = np.random.default_rng(653); out = dict(fits=dict(t=[round(R5[k]['fit']['t'], 4) for k in sorted(R5)], kappa=[round(R5[k]['fit']['kappa'], 4) for k in sorted(R5)], bound_hits=sum(len(R5[k]['fit']['bound_hits']) for k in R5),
                                                   loss_over_gas_gate=[round(R5[k]['fit']['disk_loss'] / R61[k]['fit']['disk_loss'], 4) for k in sorted(R5)]))
for s in ('kept', 'dropped', 'all'):
    hosts = sorted({n for (ss, m, n) in acc if ss == s and m == 'H'}); assert all(len(acc[(s, 'H', n)]) == 3 for n in hosts)
    f = {m: np.array([np.mean(acc[(s, m, n)]) for n in hosts]) for m in MODELS}; idx = rng.integers(0, len(hosts), (10000, len(hosts))); o = dict(galaxies=len(hosts), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()})
    for ref in ('T48', 'MOND', 'T61'):
        dd = np.sqrt(f['H'][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
        o['H_minus_' + ref] = dict(diff=round(float(np.sqrt(f['H'].mean()) - np.sqrt(f[ref].mean())), 4), pct=round(float(100 * (np.sqrt(f['H'].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f['H'] < f[ref])), n_higher=int(np.sum(f['H'] > f[ref])))
    out[s] = o; print(s, json.dumps(o), flush=True)
print(json.dumps(out['fits'])); json.dump(out, open('res65_heldout.json', 'w'), indent=1)
