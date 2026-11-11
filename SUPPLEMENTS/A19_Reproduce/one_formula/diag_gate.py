"""All-data diagnostic (not held out): rule + published Moster mass switched on by a stellar-mass and size gate."""
import numpy as np, sys
import diag_mosg as G
H = G.H; K = H.K; abc = G.abc
G.B.FROZEN['abc'] = abc
GATES = {
    'MOSH': (['h', 'gam'], [9., -3.], [30., 12.], [[11.5, 0.], [11.8, .6], [14., 4.]]),
    'MOSS8': (['logRt'], [0.], [2.5], [[1.], [1.3], [1.6]]),
    'MOSS4': (['logRt'], [0.], [2.5], [[1.], [1.3], [1.6]]),
}
for k, v in GATES.items(): H.MODELS[k] = v
def gate(model, th, ms, re):
    if model == 'MOSH': return 1 / (1 + 10 ** (-2 * (np.log10(ms) + th[1] * np.log10(re / 3) - th[0])))
    n = 8 if model == 'MOSS8' else 4
    return 1 / (1 + (10 ** th[0] / re) ** n)
def mass(model, th, ms, mp, re, z):
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + gate(model, th, ms, re) * G.mos(ms, z)
_g, _e = H.gpred, H.esd
def gpred(model, th, a):
    if model not in GATES: return _g(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in GATES: return _e(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c; rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    d, ext, lens = G.d, G.ext, G.lens
    for model in sys.argv[1:]:
        f = H.fit(model, 0.1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
        ch = [float(np.sum(lens.white(model, th, [j]) ** 2)) for j in range(4)]
        on = {n: float(np.mean(gate(model, th, a['ms'], a['re']) > .5)) for n, a in [('disk', d)] + list(ext.items())}
        print('%-6s th' % model, np.round(th, 3), f['bound_hits'], '| disk %.4f' % G.sc(d, H.gpred(model, th, d)),
              ' '.join('%s %.4f' % (n, G.sc(a, H.gpred(model, th, a))) for n, a in ext.items()), '| lens %.2f' % sum(ch), '| share of points with gate on', {k: round(v, 2) for k, v in on.items()})
