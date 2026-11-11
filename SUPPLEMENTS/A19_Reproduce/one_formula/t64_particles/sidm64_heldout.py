"""PROTOCOL_T64, the self-interacting model: in each of the 15 disk partitions (3 fold assignments x 5 folds) the
cross-section s0 is chosen from sidm64.py's grid on ALL points of the training disks with T61's weights (the same loss
that sets the fuzzy and warm numbers), and every population's held-out systems are scored at that s0 and at s0 = 0 (the
same halos with no self-interaction). Disks: benchmark points, dropped points, all points. Host-balanced rms, each host's
three held-out mean squared errors averaged; paired host bootstrap against s0 = 0 and MOND's relation.
Usage: PYTHONPATH=../dmfit python3 sidm64_heldout.py <res64_sidm_points.json> <res48 or res61 primary, for the partitions>"""
import sys, json, numpy as np
import harness_T64 as M
import compare46 as C
H = M.H; S0 = M.M61.S0
res = json.load(open(sys.argv[1])); recs = {(r['assignment'], r['fold']): r for r in json.load(open(sys.argv[2]))['records'] if r['lens'] == 0}; assert len(recs) == 15
grid = sorted(res[0][2]['by_s0'], key=float); bench = {g['name'] for g in H.disks()}
D = {n: o for p, n, o in res if p == 'disk'}; EXT = {p: {n: o for pp, n, o in res if pp == p} for p in ('dwarfs', 'ETG', 'MaNGA', 'CLASH')}
A = M.ALL; emond = np.log10(A['go'] / C.g_mond(A)); MOND = {n: emond[A['name'] == n] for n in D}
for n in D: assert len(MOND[n]) == len(D[n]['fe']) and np.allclose(np.sort(A['r'][A['name'] == n]), np.sort(D[n]['r']))
def wmse(o, s):
    fe = np.array(o['fe']); pw = 1 / ((2 / np.log(10) * fe) ** 2 + S0 ** 2); er = np.array(o['by_s0'][s]['er']); return float(np.sum(pw * er ** 2) / np.sum(pw))
def mse(er, k=None): er = np.asarray(er); return float(np.mean((er if k is None else er[k]) ** 2))
acc = {}; chosen = []; losses = []
for key in sorted(recs):
    r = recs[key]; test = set(r['losses']['disk']); train = sorted(bench - test); assert test <= bench and len(train) + len(test) == len(bench)
    L = {s: float(np.mean([wmse(D[n], s) for n in train])) for s in grid}; sb = min(grid, key=lambda s: L[s]); chosen.append(float(sb)); losses.append({s: round(L[s], 6) for s in grid})
    for n in test:
        fe = np.array(D[n]['fe']); sets = dict(kept=fe <= .1, dropped=fe > .1, all=np.ones(len(fe), bool))
        for st, k in sets.items():
            if not k.any(): continue
            acc.setdefault(('disk_' + st, 'S', n), []).append(mse(D[n]['by_s0'][sb]['er'], k)); acc.setdefault(('disk_' + st, 'S0', n), []).append(mse(D[n]['by_s0'][grid[0]]['er'], k)); acc.setdefault(('disk_' + st, 'MOND', n), []).append(mse(MOND[n], k))
    for p in ('dwarfs', 'ETG', 'MaNGA', 'CLASH'):
        for n in r['losses'][p]:
            o = EXT[p][n]; kinds = np.array(o['kind'])
            for lab, k in ([(p, None)] if p in ('dwarfs', 'MaNGA') else [(p, None)] + [(kk, kinds == kk) for kk in sorted(set(kinds))]):
                acc.setdefault((lab, 'S', n), []).append(mse(o['by_s0'][sb]['er'], k)); acc.setdefault((lab, 'S0', n), []).append(mse(o['by_s0'][grid[0]]['er'], k))
rng = np.random.default_rng(642); out = dict(grid=[float(s) for s in grid], s0_chosen=chosen, s0_counts={s: chosen.count(float(s)) for s in grid}, training_losses=losses,
                                             flags={f: sum(o['by_s0'][s]['flag'] == f for p, n, o in res for s in grid) for f in ('ok', 'nfw', 'fail')})
for lab in sorted({k[0] for k in acc}):
    hosts = sorted({n for (l, m, n) in acc if l == lab and m == 'S'}); assert all(len(acc[(lab, 'S', n)]) == 3 for n in hosts), lab
    f = {m: np.array([np.mean(acc[(lab, m, n)]) for n in hosts]) for m in ('S', 'S0', 'MOND') if (lab, m, hosts[0]) in acc}; idx = rng.integers(0, len(hosts), (10000, len(hosts)))
    o = dict(hosts=len(hosts), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()})
    for ref in [m for m in ('S0', 'MOND') if m in f]:
        dd = np.sqrt(f['S'][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
        o['S_minus_' + ref] = dict(pct=round(float(100 * (np.sqrt(f['S'].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f['S'] < f[ref])), n_higher=int(np.sum(f['S'] > f[ref])))
    out[lab] = o; print(lab, json.dumps(o), flush=True)
# all data, every grid value (descriptive)
tab = {}
for s in grid:
    row = dict(loss=float(np.mean([wmse(D[n], s) for n in sorted(bench)])))
    for st in ('kept', 'dropped', 'all'):
        v = []
        for n in sorted(bench):
            fe = np.array(D[n]['fe']); k = dict(kept=fe <= .1, dropped=fe > .1, all=np.ones(len(fe), bool))[st]
            if k.any(): v.append(mse(D[n]['by_s0'][s]['er'], k))
        row['disk_' + st] = float(np.sqrt(np.mean(v)))
    for p in ('dwarfs', 'ETG', 'MaNGA'): row[p] = float(np.sqrt(np.mean([mse(o['by_s0'][s]['er']) for o in EXT[p].values()])))
    for kk in ('CLASH_centre', 'CLASH_outer'): row[kk] = float(np.sqrt(np.mean([mse(o['by_s0'][s]['er'], np.array(o['kind']) == kk) for o in EXT['CLASH'].values()])))
    out_ = [n for n in D if n not in bench]; row['sparc_outside_all'] = float(np.sqrt(np.mean([mse(D[n]['by_s0'][s]['er']) for n in out_])))
    if any(p == 'iorio' for p, n, o in res): row['iorio_all'] = float(np.sqrt(np.mean([mse(o['by_s0'][s]['er']) for p, n, o in res if p == 'iorio'])))
    row['fails'] = sum(o['by_s0'][s]['flag'] == 'fail' for p, n, o in res); tab[s] = {k: round(v, 4) for k, v in row.items()}
    print('all data s0 =', s, json.dumps(tab[s]), flush=True)
out['all_data'] = tab; out['s0_all_data'] = float(min(grid, key=lambda s: tab[s]['loss']))
print('chosen per partition:', chosen, '| all-data choice:', out['s0_all_data']); json.dump(out, open('res64_sidm_heldout.json', 'w'), indent=1)
