"""T57 candidate (no protocol yet): T48's formula with cold gas making a light-selected galaxy's halo bigger and looser,
  M_ref = (M* + mu M_gas)/f_b + 10^A (M*/1e10)^B (R_e/3)^C + q M_E,   c = c_DM14(M_ref) x min(Sigma*/Sigma_p, 1)^(s g),
g = M_gas/M_p the gas share of the baryons, Sigma* = M*/(2 pi R_e^2), Sigma_p = 1e10/(2 pi 3^2) Msun/kpc^2 (the rule's pivot).
A, B, C, mu and s fitted on each fold's training disks. Gas-free systems (dwarf spheroidals, early types) and the
halo-selected clusters are untouched except through A, B, C. MOGF0IESBUL: 5 coefficients fitted to the benchmark."""
import sys
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as P8
P = P8.P; H = P.H; S = P.S; K = H.K; B = S.B; P44 = P.P; T = P44.T; Z = P44.Z; Y = P.Y
MODEL = 'MOGF0IESBUL'; H.MODELS[MODEL] = H.MODELS['MOGF0']
U = {'mu': 1.0, 's': 0.0, 'free': (True, True)}; _CACHE = {}; SP = 1e10 / (2 * np.pi * 9.)
LO, HI = [9, -.5, -2.5, 0., -1.], [14.5, 2, 1.5, 20., 1.]
def mass(th, ms, mp, re, z, thr):
    abc = B.FROZEN['abc']; thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    q = np.where(thr > 0, 1., Z.switch('MOG', th, ms, mp, re)); mpe = np.where(thr > 0, mp, ms + U['mu'] * (mp - ms))
    m = mpe / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + q * S.m_exp(ms, z, thr)
    return (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
def loose(ms, mp, re, thr):
    thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms)); g = np.where(thr > 0, 0., np.clip(1 - ms / mp, 0, 1))
    return np.minimum(ms / (2 * np.pi * re ** 2) / SP, 1.) ** (U['s'] * g)
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    thr = S.thresholds(a); m = mass(th, a['ms'], a['mp'], a['re'], a['z'], thr); c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * loose(a['ms'], a['mp'], a['re'], thr)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    m = mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m) * loose(st.ms, st.mp, st.re, 0.)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def fit_disks(dtr):
    key = (len(dtr['r']), float(np.sum(dtr['go'])), float(np.sum(dtr['r'])), U['free'])
    if key not in _CACHE:
        th0 = np.array([11.3]); fm, fs = U['free']
        def unpack(p): return np.asarray(p[:3]), (float(p[3]) if fm else 1.), (float(p[3 + fm]) if fs else 0.)
        def fun(p):
            B.FROZEN['abc'], U['mu'], U['s'] = unpack(p); return dtr['hw'] * np.log10(dtr['go'] / gpred(MODEL, th0, dtr)) / .15
        lo = LO[:3] + ([LO[3]] if fm else []) + ([LO[4]] if fs else []); hi = HI[:3] + ([HI[3]] if fm else []) + ([HI[4]] if fs else []); att = []
        for abc0, m0, s0 in (([11.7481, .8178, -.7119], 1., 0.), ([11.7481, .8178, -.7119], 4., .1), ([11.6, .9, -1.2], 6., .1)):
            st = abc0 + ([m0] if fm else []) + ([s0] if fs else [])
            o = least_squares(fun, st, bounds=(lo, hi), max_nfev=600, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        good = [x for x in att if x['success']] or att; best = min(good, key=lambda x: x['loss'])
        bh = [n for n, t, l, u in zip(['A', 'B', 'C'] + (['mu'] if fm else []) + (['s'] if fs else []), best['th'], lo, hi) if min(abs(t - l), abs(t - u)) < 1e-5]
        _CACHE[key] = (unpack(best['th']), best['loss'], att, bh)
    return _CACHE[key]
def fit(model, weight, dtr, etr, lens, groups):
    if model != MODEL: return _fit(model, weight, dtr, etr, lens, groups)
    (abc, mu, s), loss, att, bh = fit_disks(dtr); B.FROZEN['abc'] = np.array(abc); U['mu'] = mu; U['s'] = s
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in abc], mu=float(mu), s=float(s), disk_loss=loss)
H.fit, H.gpred, H.esd = fit, gpred, esd
