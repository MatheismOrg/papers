"""T34 (PROTOCOL_T34.txt): sensitivity of T32 and T33 to the selection threshold and the scatter."""
import sys, csv, json
import numpy as np
import harness_T26 as M
F = M.F; H = M.H; K = H.K; B = __import__('harness_T15b')
sys.path.insert(0, str(H.PKG / 'code')); import A19_Clusters_v5 as V5
from colossus.lss import mass_function
from colossus.cosmology import cosmology
cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624}); cosmology.setCosmology('A19planck')
B.FROZEN['abc'] = abc = np.array([11.7481, .8178, -.7119]); th = np.array([11.1868]); HH = 0.671; GRID = np.arange(11.0, 16.5 + 1e-9, 0.005)
cl = H.cluster_rows()
def gdm(m, r, z):
    c = K.c_dm14(m, z); rv = (H.G * m / (100 * H.hz(z) ** 2)) ** (1 / 3)
    return H.G * (1 - H.FB) * m * H.ff(r * c / rv) / H.ff(c) / r ** 2 * H.ACC
def predict(a, msel, sig):
    pred = np.zeros(len(a['r']))
    for host in np.unique(a['name']):
        k = a['name'] == host; ms, mp, re, z = a['ms'][k][0], a['mp'][k][0], a['re'][k][0], a['z'][k][0]
        base = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
        q = float(F.Z.switch('MOG', th, np.array([ms]), np.array([mp]), np.array([re]))[0])
        prior = mass_function.massFunction(10 ** GRID * HH, z, mdef='200c', model='tinker08', q_out='dndlnM')
        like = np.exp(-.5 * ((np.log10(ms) - np.log10(V5.moster_mstar(10 ** GRID, z))) / sig) ** 2)
        w = prior * like * (10 ** GRID > msel); w /= w.sum(); j = w > 1e-8 * w.max()
        g = np.array([np.log10(a['gb'][k] + gdm(base + q * m, a['r'][k], z)) for m in 10 ** GRID[j]])
        pred[k] = w[j] @ g / w[j].sum()
    return pred
def score(a, p, kk):
    er = np.log10(a['go'][kk]) - p[kk]; hosts = a['name'][kk]
    return round(float(np.sqrt(np.mean([np.mean(er[hosts == h] ** 2) for h in np.unique(hosts)]))), 4)
res = {'CLASH': [], 'MaNGA': []}
c = H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG']); frozen_c = np.log10(H.gpred('MOGF0IE', th, c))
route = {r['host']: r['route'] for r in csv.DictReader(open(H.PKG / 'data/Clusters/umetsu2016_table2.csv'))}
assert set(np.unique(c['name'])) <= set(route), set(np.unique(c['name'])) - set(route)
for msel in (3e14, 4e14, 5e14, 7e14):
    for sig in (.15, np.hypot(.15, .10), .20):
        p = predict(c, msel, sig)
        d = dict(threshold=msel, sigma=round(float(sig), 3), centre=score(c, p, c['kind'] == 'CLASH_centre'), outer=score(c, p, c['kind'] == 'CLASH_outer'))
        if msel == 5e14 and abs(sig - .18) < .01:
            for rt in sorted(set(route.values())):
                kr = np.array([route[n] == rt for n in c['name']])
                d['by_route_' + rt] = dict(n=int(len(np.unique(c['name'][kr]))), centre=score(c, p, kr & (c['kind'] == 'CLASH_centre')), outer=score(c, p, kr & (c['kind'] == 'CLASH_outer')),
                                           frozen_centre=score(c, frozen_c, kr & (c['kind'] == 'CLASH_centre')), frozen_outer=score(c, frozen_c, kr & (c['kind'] == 'CLASH_outer')))
        res['CLASH'].append(d); print('CLASH', d, flush=True)
m = H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG'])
for lt in (13.8, 14.0, 14.2):
    for sig in (.15, np.hypot(.15, .10), .20):
        p = predict(m, 10 ** lt, sig); d = dict(log_threshold=lt, sigma=round(float(sig), 3), MaNGA=score(m, p, np.ones(len(m['r']), bool)))
        res['MaNGA'].append(d); print('MaNGA', d, flush=True)
res['robust'] = all(d['centre'] < .1545 and d['outer'] < .2073 for d in res['CLASH']) and all(d['MaNGA'] <= .1141 * 1.05 for d in res['MaNGA'])
print('robust:', res['robust']); json.dump(res, open('res34_selection_sensitivity.json', 'w'), indent=1)
