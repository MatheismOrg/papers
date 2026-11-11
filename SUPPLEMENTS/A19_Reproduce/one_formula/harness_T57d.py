"""T57 (PROTOCOL_T57.txt): T48's formula with bigger and looser halos for diffuse, gas-rich galaxies only. With
g = M_gas/M_p (gas share of the baryons), L = -log10 min(Sigma*/Sigma_p, 1), Sigma* = M*/(2 pi R_e^2),
Sigma_p = 1e10/(2 pi 3^2) Msun/kpc^2 (the rule's own pivot), for light-selected systems:
  M_ref -> M_ref x 10^(t g L),    c -> c x 10^(-s g L).
The term is zero for gas-free systems (dwarf spheroidals, early types), for dense ones (L = 0) and for the halo-selected
clusters. MOGF0IESBD: A, B, C as T48 (the paper's rule on the training disks), then t and s fitted on the same training
disks with A, B, C fixed. MOGF0IESBD5: A, B, C, t, s fitted together on the training disks. M_t as T45 in both."""
import sys
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as P8
P = P8.P; H = P.H; S = P.S; K = H.K; B = S.B; P44 = P.P; G43 = P44.G; T = P44.T; Z = P44.Z; Y = P.Y
MODELS = ('MOGF0IESBD', 'MOGF0IESBD5')
for m_ in MODELS: H.MODELS[m_] = H.MODELS['MOGF0']
U = {'t': 0.0, 's': 0.0}; _CACHE = {}; SP = 1e10 / (2 * np.pi * 9.)
def gl(ms, mp, re, thr):
    thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    return np.where(thr > 0, 0., np.clip(1 - ms / mp, 0, 1) * -np.log10(np.minimum(ms / (2 * np.pi * re ** 2) / SP, 1.)))
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model not in MODELS: return _gpred(model, th, a)
    thr = S.thresholds(a); w = gl(a['ms'], a['mp'], a['re'], thr)
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], thr) * 10 ** (U['t'] * w); c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * 10 ** (-U['s'] * w)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in MODELS: return _esd(model, th, st)
    w = gl(st.ms, st.mp, st.re, 0.); m = G43.mass(th, st.ms, st.mp, st.re, st.z, 0.) * 10 ** (U['t'] * w)
    c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m) * 10 ** (-U['s'] * w); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def fit_disks(model, dtr, abc0):
    key = (model, len(dtr['r']), float(np.sum(dtr['go'])), float(np.sum(dtr['r'])))
    if key not in _CACHE:
        th0 = np.array([11.3]); joint = model == 'MOGF0IESBD5'
        def unpack(p): return (np.asarray(p[:3]), float(p[3]), float(p[4])) if joint else (np.asarray(abc0), float(p[0]), float(p[1]))
        def fun(p):
            B.FROZEN['abc'], U['t'], U['s'] = unpack(p); return dtr['hw'] * np.log10(dtr['go'] / gpred(model, th0, dtr)) / .15
        lo = ([9, -.5, -2.5] if joint else []) + [-2., -2.]; hi = ([14.5, 2, 1.5] if joint else []) + [2., 2.]; att = []
        for t0, s0 in ((0., 0.), (.3, .15), (.6, .3)):
            st = (list(abc0) if joint else []) + [t0, s0]
            o = least_squares(fun, st, bounds=(lo, hi), max_nfev=600, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        good = [x for x in att if x['success']] or att; best = min(good, key=lambda x: x['loss'])
        names = (['A', 'B', 'C'] if joint else []) + ['t', 's']; bh = [n for n, v, l, u in zip(names, best['th'], lo, hi) if min(abs(v - l), abs(v - u)) < 1e-5]
        _CACHE[key] = (unpack(best['th']), best['loss'], att, bh)
    return _CACHE[key]
_RULE = {}
def fit(model, weight, dtr, etr, lens, groups):
    if model not in MODELS: return _fit(model, weight, dtr, etr, lens, groups)
    k = (len(dtr['r']), float(np.sum(dtr['go'])))
    if k not in _RULE: _RULE[k] = np.array(_fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])          # the paper's rule on these training disks, as chunk_split sets it
    (abc, t, s), loss, att, bh = fit_disks(model, dtr, _RULE[k]); B.FROZEN['abc'] = np.array(abc); U['t'] = t; U['s'] = s
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in abc], t=float(t), s=float(s), disk_loss=loss)
H.fit, H.gpred, H.esd = fit, gpred, esd
