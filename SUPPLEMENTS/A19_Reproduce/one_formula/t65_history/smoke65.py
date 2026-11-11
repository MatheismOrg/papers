"""T65 smoke test: every function runs for every galaxy; prints only whether the predictions are finite and positive
(no comparison with any measurement)."""
import json, numpy as np
import lib65 as L
sf = json.load(open('sfh65_inputs.json')); ok = True
for lab, v in sf['read2019'].items():
    a, dsph = L.galaxy('../t59', lab); tsf = L.TU - v['t_trunc']; r = np.array([0.15])
    p = [float(L.reading(a, r, k, tsf=tsf, tk=(0.1, 0.5))[1][0]) for k in ('N', 'G', 'H1', 'H2')]; mm, rho, flag, r1 = L.sidm(a, r, dsph); p.append(float(rho[0]))
    good = all(np.isfinite(x) and x > 0 for x in p); ok &= good; print('%-9s rows %2d predictions finite and positive: %s (self-interaction solver: %s)' % (lab, len(a['r']), good, flag))
d = L.H.arr(L.H.dwarf_rows()); n = 0
for wn, nm in L.WEISZ.items():
    k = d['name'] == nm; a = {key: np.asarray(x)[k] for key, x in d.items()}; tsf = max(L.TU - 10 ** sf['weisz2015'][wn]['tau90'] / 1e9, 0.05)
    g = [float((a['gb'] + L.H.G * L.reading(a, a['r'], kk, tsf=tsf, tk=(0.1, 0.5))[0] / a['r'] ** 2 * L.H.ACC)[0]) for kk in ('G', 'H1', 'H2')]; n += all(np.isfinite(x) and x > 0 for x in g)
print('Weisz dwarfs with finite predictions:', n, 'of', len(L.WEISZ), '| all sixteen fine:', ok)
