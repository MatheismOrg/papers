"""T33 (PROTOCOL_T33.txt): the MaNGA central galaxies predicted with the scatter-aware Moster inversion, with and without
the MaNGA BCG program's stated halo selection (M200c > 1e14 suns), against the frozen formula. All data, nothing fitted."""
import sys, json
import numpy as np
import harness_T26 as M
F = M.F; H = M.H; K = H.K; B = __import__('harness_T15b')
sys.path.insert(0, str(H.PKG / 'code')); import A19_Clusters_v5 as V5
from colossus.lss import mass_function
from colossus.cosmology import cosmology
cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624}); cosmology.setCosmology('A19planck')
B.FROZEN['abc'] = abc = np.array([11.7481, .8178, -.7119]); th = np.array([11.1868])
HH = 0.671; GRID = np.arange(11.0, 16.5 + 1e-9, 0.005); SIG = np.hypot(.15, .10); MSEL = 1e14
cl = H.cluster_rows(); a = H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG'])
def gdm(m, r, z):
    c = K.c_dm14(m, z); rv = (H.G * m / (100 * H.hz(z) ** 2)) ** (1 / 3)
    return H.G * (1 - H.FB) * m * H.ff(r * c / rv) / H.ff(c) / r ** 2 * H.ACC
frozen = np.log10(H.gpred('MOGF0IE', th, a))
pred = {'full_prior': np.zeros(len(a['r'])), 'selected': np.zeros(len(a['r']))}; post = {}
for host in np.unique(a['name']):
    k = a['name'] == host; ms, mp, re, z = a['ms'][k][0], a['mp'][k][0], a['re'][k][0], a['z'][k][0]
    base = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
    q = float(F.Z.switch('MOG', th, np.array([ms]), np.array([mp]), np.array([re]))[0])
    prior = mass_function.massFunction(10 ** GRID * HH, z, mdef='200c', model='tinker08', q_out='dndlnM')
    like = np.exp(-.5 * ((np.log10(ms) - np.log10(V5.moster_mstar(10 ** GRID, z))) / SIG) ** 2)
    for lab, sel in (('full_prior', np.ones_like(GRID)), ('selected', (10 ** GRID > MSEL).astype(float))):
        w = prior * like * sel; w /= w.sum(); j = w > 1e-8 * w.max()
        mh = base + q * 10 ** GRID[j]
        g = np.array([np.log10(a['gb'][k] + gdm(m, a['r'][k], z)) for m in mh])     # (grid, points)
        pred[lab][k] = w[j] @ g / w[j].sum()
        post.setdefault(host, {})[lab] = round(float(np.log10(np.sum(w[j] * 10 ** GRID[j]) / w[j].sum())), 2)
def score(p, kind):
    kk = a['kind'] == kind; er = np.log10(a['go'][kk]) - p[kk]; hosts = a['name'][kk]
    return round(float(np.sqrt(np.mean([np.mean(er[hosts == h] ** 2) for h in np.unique(hosts)]))), 4)
res = {lab: {kind: score(p, kind) for kind in ('MaNGA_BCG',)} for lab, p in [('frozen_mean_inversion', frozen), ('scatter_full_prior', pred['full_prior']), ('scatter_selected', pred['selected'])]}
um = {}
import csv
for r in csv.DictReader(open(H.PKG / 'data/Clusters/umetsu2016_table2.csv')): um = um or list(r.keys())
res['posterior_mean_log_M200c_selected'] = {h: v['selected'] for h, v in post.items()}
res['posterior_mean_log_M200c_full'] = {h: v['full_prior'] for h, v in post.items()}
fz, ss = res['frozen_mean_inversion'], res['scatter_selected']
res['selection_reading_holds'] = all(ss[k] <= fz[k] * 1.05 for k in fz)
print(json.dumps({k: v for k, v in res.items() if not k.startswith('posterior')}, indent=1)); print('posterior means (selected):', sorted(res['posterior_mean_log_M200c_selected'].values()))
json.dump(res, open('res33_manga_selection.json', 'w'), indent=1)
