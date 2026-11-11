"""FORECAST_T58: T57 against T48 on the SPARC galaxies excluded from the benchmark (Q = 3 or inclination < 30 deg), all
coefficients frozen. Usage: PYTHONPATH=<dmfit> python3 excluded58.py"""
import os, json, numpy as np
import harness_T57d as M7
import harness_T48 as M8
H = M7.H; B = M7.B; D = H.D
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); M7.U['t'] = 0.361; M7.U['s'] = 0.146
t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); rows = []
for g, r in sorted(t1.items()):
    if r['Q'] < 3 and r['Inc'] >= 30: continue
    d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= 0.1) & (vb2 > 0)
    Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']
    rows += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
a = H.arr(rows); un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi)
e8 = np.log10(a['go'] / H.gpred(M8.MODEL, th, a)); e7 = np.log10(a['go'] / H.gpred('MOGF0IESBD', th, a))
m8 = np.bincount(gi, e8 ** 2) / ct; m7 = np.bincount(gi, e7 ** 2) / ct; w = np.array([M7.gl(a['ms'][gi == i][:1], a['mp'][gi == i][:1], a['re'][gi == i][:1], 0.)[0] for i in range(len(un))])
rng = np.random.default_rng(58)
def cmp(sel, lab):
    x, y = m8[sel], m7[sel]; idx = rng.integers(0, len(x), (10000, len(x))); dd = np.sqrt(y[idx].mean(1)) - np.sqrt(x[idx].mean(1))
    r = dict(n=int(sel.sum()), T48=round(float(np.sqrt(x.mean())), 4), T57=round(float(np.sqrt(y.mean())), 4), pct=round(float(100 * (np.sqrt(y.mean() / x.mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(),
             p_not_better=round(float(np.mean(dd >= 0)), 4), n_better=int(np.sum(y < x)), n_worse=int(np.sum(y > x)))
    print(lab, r); return r
out = dict(galaxies=len(un), points=len(a['r']), all=cmp(np.ones(len(un), bool), 'all excluded galaxies'), with_term=cmp(w > 0, 'with a nonzero term (g L > 0)'), strong_term=cmp(w > np.median(w[w > 0]), 'term above its median'))
json.dump(out, open('res58_excluded.json', 'w'), indent=1)
