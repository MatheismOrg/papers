"""T23: T21's formulas (gas switch + halo floor 1e8 Msun, with or without contraction) with the steepness n of the mass
gate fitted instead of fixed at 4: q = (M*/Mp)^50 / (1 + (M_t/M*)^n), n in [1, 16].
MOGFn:  MOG + floor, coefficients log10 M_t, n                 (5 in all)
MOGCFn: MOGC + floor, coefficients log10 M_t, s, n             (6 in all)
Usage: python3 harness_T23.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T21 as F
Z = F.Z; H = Z.H; B = Z.B; K = H.K; T = Z.T; Y = Z.Y
NEW = {'MOGFn': False, 'MOGCFn': True}
H.MODELS['MOGFn'] = (['logMt', 'n'], [10., 1.], [12.5, 16.], [[11.2, 4.], [11.0, 2.], [11.3, 8.]])
H.MODELS['MOGCFn'] = (['logMt', 's', 'n'], [10., 0., 1.], [12.5, 1.5, 16.], [[11.2, .17, 4.], [11.0, 0., 2.], [11.3, .3, 8.]])
def mass(model, th, ms, mp, re, z):
    abc = B.FROZEN['abc']; n = th[-1]
    m = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + (ms / mp) ** 50 / (1 + (10 ** th[0] / ms) ** n) * Y.moster(ms, z)
    return (m ** 8 + 1e8 ** 8) ** (1 / 8)
def strength(model, th): return th[1] if NEW[model] else 0.
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in NEW: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Z.Fc(eta) ** strength(model, th) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in NEW: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, strength(model, th))
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
