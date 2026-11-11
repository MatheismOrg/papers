"""All-data diagnostic: switch formulas (size: MOSS; cold gas: MOG) with contraction-only halo response (Cautun 2020 eq. 11,
F >= 1) at strength s fitted jointly with the switch coefficient; the paper's rule fixed as fitted on disks without contraction."""
import numpy as np, sys
import harness_T18 as T, harness_T17b as Y
H = T.H; B = T.B; K = H.K
F_ORIG, DF_ORIG = T.F, T.dF
Fc = lambda e: np.maximum(F_ORIG(e), 1.)
def dFc(e): return np.where(F_ORIG(e) > 1., DF_ORIG(e), 0.)
MODELS = {'MOSC': (['logRt', 's'], [0., 0.], [2.5, 1.5], [[1.15, .3], [1.15, .0], [1.3, .6]]),
          'MOGC': (['logMt', 's'], [10., 0.], [12.5, 1.5], [[11.2, .3], [11.0, .0], [11.4, .6]]),
          'MOG': (['logMt'], [10.], [12.5], [[11.2], [11.0], [11.4]])}
H.MODELS.update(MODELS)
def mass(model, th, ms, mp, re, z):
    abc = B.FROZEN['abc']; base = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
    if model == 'MOSC': s = 1 / (1 + (10 ** th[0] / re) ** 8)
    else: s = (ms / mp) ** 50 / (1 + (10 ** th[0] / ms) ** 4)
    return base + s * Y.moster(ms, z)
def strength(model, th): return 0. if model == 'MOG' else th[1]
_g, _e = H.gpred, H.esd
def gpred(model, th, a):
    if model not in MODELS: return _g(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Fc(eta) ** strength(model, th) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in MODELS: return _e(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    F0, d0 = T.F, T.dF; T.F, T.dF = Fc, dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, strength(model, th))
    finally: T.F, T.dF = F0, d0
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    import diag_mosg as G
    d, ext, lens = G.d, G.ext, G.lens
    B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
    for model in sys.argv[1:]:
        f = H.fit(model, 0.1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
        z = [float(np.sum(lens.white(model, th, [j]) ** 2)) for j in range(4)]
        print('%-5s all-data th %s %s | disk %.4f' % (model, np.round(th, 3), f['bound_hits'], G.sc(d, H.gpred(model, th, d))), ' '.join('%s %.4f' % (n, G.sc(a, H.gpred(model, th, a))) for n, a in ext.items()),
              '| lens %.1f' % sum(z), np.round(z, 1), flush=True)
