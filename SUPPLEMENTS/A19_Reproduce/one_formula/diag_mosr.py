"""All-data diagnostic (not held out): giant-galaxy term = published Moster (2013) halo mass x size term, or a gated Moster term."""
import numpy as np, harness_T17 as X
import diag_mosg as G
H = X.H; B = X.B; K = H.K
abc = G.abc; B.FROZEN['abc'] = abc
H.MODELS['MOSR'] = (['dA', 'C2'], [-4, -2], [3, 4], [[0., 2.5], [-1., 2.], [.5, 1.]])
H.MODELS['MOSG'] = (['lt'], [10., ], [12.5], [[11.0], [11.4], [11.8]])
_g, _e = H.gpred, H.esd
def mass(model, th, ms, mp, re, z):
    base = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
    if model == 'MOSR': return base + 10 ** th[0] * G.mos(ms, z) * (re / 10) ** th[1]
    return base + G.mos(ms, z) / (1 + (10 ** th[0] / ms) ** 4)
def gpred(model, th, a):
    if model not in ('MOSR', 'MOSG'): return _g(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in ('MOSR', 'MOSG'): return _e(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c; rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    d, ext, lens = G.d, G.ext, G.lens
    for model in ['SUBH', 'MOSR', 'MOSG']:
        if model == 'SUBH': th = np.array([12.0629, 2.63, 2.5098]); f = dict(bound_hits=[])
        else: f = H.fit(model, 0.1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
        ch = [float(np.sum(lens.white(model, th, [j]) ** 2)) for j in range(4)]
        print('%-5s all-data th' % model, np.round(th, 3), f['bound_hits'], '| disk %.4f' % G.sc(d, H.gpred(model, th, d)),
              ' '.join('%s %.4f' % (n, G.sc(a, H.gpred(model, th, a))) for n, a in ext.items()), '| lens %.2f' % sum(ch), np.round(ch, 1))
