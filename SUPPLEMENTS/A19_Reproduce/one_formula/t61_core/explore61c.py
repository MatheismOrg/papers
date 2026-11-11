"""SA73 T61 all-data looks, third pass: the core's gate. n = tanh(kappa g L) (gas share times diffuseness, as the mass
term) against n = tanh(kappa g) (gas share alone). Weighted fits on all points of the benchmark galaxies, as explore61b.py.
Usage: PYTHONPATH=<dmfit> python3 explore61c.py"""
import json, numpy as np
src = open('explore61b.py').read().split("\nout = {}")[0]; exec(compile(src, 'explore61b_defs', 'exec'))
gfr = np.clip(1 - a['ms'] / a['mp'], 0, 1)
def gp2(p, gate):
    t, kap = p; m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], 0.) * 10 ** (t * w); c = K.c_dm14(m, a['z'])
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c) * np.tanh(a['r'] / (1.75 * a['re'])) ** np.tanh(kap * (w if gate == 'gL' else gfr))
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
ka = SETS['all149']; un, gi = np.unique(a['name'][ka], return_inverse=True); nw = np.zeros(len(fe)); nw[ka] = np.sqrt(pw[ka] / np.bincount(gi, pw[ka])[gi] / len(un)); out = {}
for gate in ('gL', 'g'):
    best = None
    for s0 in ((.1, .6), (.2, 1.5), (0., 3.)):
        o = least_squares(lambda q: (nw * np.log10(a['go'] / gp2(q, gate)))[ka], np.array(s0), bounds=([-2, 0], [2, 30]), ftol=1e-10, xtol=1e-10)
        if best is None or o.cost < best.cost: best = o
    out[gate] = dict(p=best.x.tolist(), **rep('gate %s -> t %.3f kappa %.3f' % (gate, *best.x), gp2(best.x, gate)))
print('median g, L, gL over benchmark galaxies with gas:', np.round([np.median(gfr[SETS['kept149']]), np.median(w[SETS['kept149']] / np.maximum(gfr[SETS['kept149']], 1e-9)), np.median(w[SETS['kept149']])], 3))
json.dump(out, open('res61_explore_c.json', 'w'), indent=1)
