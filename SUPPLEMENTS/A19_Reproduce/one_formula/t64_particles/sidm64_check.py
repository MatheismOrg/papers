"""T64, after the protocol's run (a check, not a mark): where the self-interaction solver failed (a failed solve falls back
to the NFW halo, which flatters the model there), and the held-out scores with every host that failed at the chosen
cross-section left out of both sides. Also the core radii r1 at s0 = 1.
Usage: PYTHONPATH=../dmfit python3 sidm64_check.py res64_sidm_points.json ../dmfit/res48_MOGF0IESB_primary.json"""
import sys, json, numpy as np
import harness_T64 as M
H = M.H; S0 = M.M61.S0
res = json.load(open(sys.argv[1])); recs = {(r['assignment'], r['fold']): r for r in json.load(open(sys.argv[2]))['records'] if r['lens'] == 0}
grid = sorted(res[0][2]['by_s0'], key=float); bench = {g['name'] for g in H.disks()}; out = {}
fl = {}
for p, n, o in res:
    for s in grid:
        if o['by_s0'][s]['flag'] == 'fail': fl.setdefault(s, {}).setdefault(p, []).append(n)
out['failures'] = {s: {p: len(v) for p, v in d.items()} for s, d in fl.items()}; print('failures by cross-section and population:', json.dumps(out['failures']))
print('failed at s0 = 1:', json.dumps(fl.get('1.0', {})))
D = {n: o for p, n, o in res if p == 'disk'}; EXT = {p: {n: o for pp, n, o in res if pp == p} for p in ('dwarfs', 'ETG', 'MaNGA', 'CLASH')}
def wmse(o, s):
    fe = np.array(o['fe']); pw = 1 / ((2 / np.log(10) * fe) ** 2 + S0 ** 2); er = np.array(o['by_s0'][s]['er']); return float(np.sum(pw * er ** 2) / np.sum(pw))
def mse(er, k=None): er = np.asarray(er); return float(np.mean((er if k is None else er[k]) ** 2))
acc = {}
for key in sorted(recs):
    r = recs[key]; test = set(r['losses']['disk']); train = sorted(bench - test); sb = min(grid, key=lambda s: float(np.mean([wmse(D[n], s) for n in train])))
    for n in test:
        if D[n]['by_s0'][sb]['flag'] == 'fail': continue
        fe = np.array(D[n]['fe'])
        for st, k in dict(kept=fe <= .1, dropped=fe > .1, all=np.ones(len(fe), bool)).items():
            if k.any(): acc.setdefault(('disk_' + st, 'S', n), []).append(mse(D[n]['by_s0'][sb]['er'], k)); acc.setdefault(('disk_' + st, 'S0', n), []).append(mse(D[n]['by_s0'][grid[0]]['er'], k))
    for p in ('dwarfs', 'ETG', 'MaNGA', 'CLASH'):
        for n in r['losses'][p]:
            o = EXT[p][n]
            if o['by_s0'][sb]['flag'] == 'fail': continue
            kinds = np.array(o['kind'])
            for lab, k in ([(p, None)] if p != 'CLASH' else [('CLASH_centre', kinds == 'CLASH_centre')]):
                acc.setdefault((lab, 'S', n), []).append(mse(o['by_s0'][sb]['er'], k)); acc.setdefault((lab, 'S0', n), []).append(mse(o['by_s0'][grid[0]]['er'], k))
for lab in sorted({k[0] for k in acc}):
    hosts = sorted({n for (l, m, n) in acc if l == lab and m == 'S'}); f = {m: np.array([np.mean(acc[(lab, m, n)]) for n in hosts]) for m in ('S', 'S0')}
    out[lab] = dict(hosts=len(hosts), S=round(float(np.sqrt(f['S'].mean())), 4), S0=round(float(np.sqrt(f['S0'].mean())), 4), pct=round(float(100 * (np.sqrt(f['S'].mean() / f['S0'].mean()) - 1)), 2)); print('without failed hosts:', lab, json.dumps(out[lab]))
r1 = {p: [o['by_s0']['1.0']['r1'] for pp, n, o in res if pp == p and o['by_s0']['1.0']['flag'] == 'ok'] for p in ('disk', 'dwarfs', 'ETG', 'MaNGA', 'CLASH')}
out['r1_kpc_at_s0_1'] = {p: [round(float(np.percentile(v, q)), 2) for q in (16, 50, 84)] for p, v in r1.items()}; print('core radius r1 (kpc; 16th, 50th, 84th percentile) at s0 = 1:', json.dumps(out['r1_kpc_at_s0_1']))
json.dump(out, open('res64_sidm_check.json', 'w'), indent=1)
