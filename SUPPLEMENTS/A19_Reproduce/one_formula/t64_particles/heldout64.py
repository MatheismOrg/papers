"""PROTOCOL_T64, marks (b) and (c) for the fuzzy and warm models: with each partition's A, B, C and particle number, the
held-out disks' dropped points (velocity error above 10 percent) and all points, against T48 (the same A, B, C), the
formula with its gas-gated core (T61, that partition's t and kappa), MOND's relation and standard halos. Each disk is held
out once per fold assignment; its three mean squared errors are averaged. Host-balanced rms, paired host bootstrap.
Usage: PYTHONPATH=../dmfit python3 heldout64.py <res64 F primary> <res64 W primary> <res61 primary>"""
import sys, json, numpy as np
import harness_T64 as M
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B; M61 = M.M61
def recs(p): return {(r['assignment'], r['fold']): r for r in json.load(open(p))['records'] if r['lens'] == 0}
RF, RW, R61 = recs(sys.argv[1]), recs(sys.argv[2]), recs(sys.argv[3]); assert len(RF) == len(RW) == len(R61) == 15
A = M.ALL; acc = {}; MODELS = ('F', 'W', 'T61', 'T48', 'MOND', 'STD')
for key in sorted(RF):
    rf, rw, r6 = RF[key], RW[key], R61[key]; names = list(rf['losses']['disk']); assert set(names) == set(rw['losses']['disk']) == set(r6['losses']['disk'])
    assert np.allclose(rf['fit']['abc'], rw['fit']['abc']) and np.allclose(rf['fit']['abc'], r6['fit']['abc'])
    a = M.pick(A, np.isin(A['name'], names)); B.FROZEN['abc'] = np.array(rf['fit']['abc']); th = np.array(rf['fit']['th'])
    M.V['lm22'] = rf['fit']['par']; M.V['lhm'] = rw['fit']['par']; M61.U['t'], M61.U['k'] = r6['fit']['t'], r6['fit']['kappa']
    e = dict(F=np.log10(a['go'] / H.gpred('MOGF0IESBF', th, a)), W=np.log10(a['go'] / H.gpred('MOGF0IESBW', np.array(rw['fit']['th']), a)), T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', np.array(r6['fit']['th']), a)),
             T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
    Vc = np.sqrt(a['go'] * a['r'] / H.ACC); sets = dict(kept=a['fe'] <= .1, dropped=a['fe'] > .1, all=np.ones(len(Vc), bool), slow_dropped=(a['fe'] > .1) & (Vc < 30))
    for s, k in sets.items():
        for n in np.unique(a['name'][k]):
            kk = k & (a['name'] == n)
            for m, v in e.items(): acc.setdefault((s, m, n), []).append(float(np.mean(v[kk] ** 2)))
rng = np.random.default_rng(64)
out = dict(fits=dict(lm22=[round(RF[k]['fit']['par'], 4) for k in sorted(RF)], lhm=[round(RW[k]['fit']['par'], 4) for k in sorted(RW)],
                     fuzzy_bound_hits=sum('lm22' in RF[k]['fit']['bound_hits'] for k in RF), warm_bound_hits=sum('lhm' in RW[k]['fit']['bound_hits'] for k in RW),
                     fuzzy_loss_over_cold=[round(RF[k]['fit']['disk_loss'] / RF[k]['fit']['disk_loss_cold'], 5) for k in sorted(RF)], warm_loss_over_cold=[round(RW[k]['fit']['disk_loss'] / RW[k]['fit']['disk_loss_cold'], 5) for k in sorted(RW)]))
for s in ('kept', 'dropped', 'all', 'slow_dropped'):
    hosts = sorted({n for (ss, m, n) in acc if ss == s and m == 'F'}); assert all(len(acc[(s, 'F', n)]) == 3 for n in hosts)
    f = {m: np.array([np.mean(acc[(s, m, n)]) for n in hosts]) for m in MODELS}; idx = rng.integers(0, len(hosts), (10000, len(hosts)))
    o = dict(galaxies=len(hosts), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()})
    for mod in ('F', 'W'):
        for ref in ('T48', 'MOND', 'T61'):
            dd = np.sqrt(f[mod][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
            o[mod + '_minus_' + ref] = dict(diff=round(float(np.sqrt(f[mod].mean()) - np.sqrt(f[ref].mean())), 4), pct=round(float(100 * (np.sqrt(f[mod].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(),
                                            n_lower=int(np.sum(f[mod] < f[ref])), n_higher=int(np.sum(f[mod] > f[ref])))
    out[s] = o; print(s, json.dumps(o), flush=True)
print(json.dumps(out['fits'])); json.dump(out, open('res64_heldout.json', 'w'), indent=1)
