"""T65 (PROTOCOL_T65.txt): the final formula's core (T61) with its gate changed from "how much cold gas" to "any cold gas".
T61: w = g L, g = M_gas/M_p, L = -log10 min(Sigma*/Sigma_p, 1). Here, MOGF0IESBH: w = s L with s = 1 for a light-selected
system that holds cold gas (M_p above M*) and s = 0 for one that holds none. Everything else is T61's: M_ref x 10^(t w),
M_DM(r) = M_NFW(r) x tanh(r / (1.75 R_e))^n, n = tanh(kappa w); lensing takes the mass term only; A, B, C as T48 on the
training disks' benchmark points; t and kappa on ALL points of the same training disks with T61's weights; M_t as T45.
The disks of the benchmark all hold cold gas and are all forming stars, so on them s = 1 stands for "still forming stars":
this is the gate a star-formation history would give the disks. Outside the harness (t65/lib65.py) s for a system without
gas is the share of cosmic time over which it formed stars, from published star-formation histories."""
import numpy as np
from scipy.optimize import least_squares
import harness_T61 as M61
H = M61.H; S = M61.S; K = M61.K; B = M61.B; P = M61.P; P44 = M61.P44; G43 = M61.G43; T = M61.T; Z = M61.Z
ALL = M61.ALL; pick = M61.pick; weights = M61.weights; ETA = M61.ETA; SP = M61.SP
MODEL = 'MOGF0IESBH'; H.MODELS[MODEL] = H.MODELS['MOGF0']
U = {'t': 0.0, 'k': 0.0}; _CACHE = {}
def diffuse(ms, re): return -np.log10(np.minimum(ms / (2 * np.pi * re ** 2) / SP, 1.))
def hl(ms, mp, re, thr):
    thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms)); return np.where((thr > 0) | (mp <= ms * (1 + 1e-9)), 0., diffuse(ms, re))
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    thr = S.thresholds(a); w = hl(a['ms'], a['mp'], a['re'], thr)
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], thr) * 10 ** (U['t'] * w); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    with np.errstate(all='ignore'): core = np.where(w > 0, np.tanh(a['r'] / (ETA * a['re'])) ** np.tanh(U['k'] * w), 1.)
    return a['gb'] + H.G * (1 - H.FB) * mt * core / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    w = hl(st.ms, st.mp, st.re, 0.); m = G43.mass(th, st.ms, st.mp, st.re, st.z, 0.) * 10 ** (U['t'] * w)
    c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def fit_disks(names, abc0):
    key = (tuple(sorted(names)), tuple(np.round(abc0, 9)))
    if key not in _CACHE:
        a = pick(ALL, np.isin(ALL['name'], list(names))); assert len(np.unique(a['name'])) == len(set(names))
        nw = weights(a); th0 = np.array([11.3]); B.FROZEN['abc'] = np.asarray(abc0)
        def fun(p):
            U['t'], U['k'] = float(p[0]), float(p[1]); return nw * np.log10(a['go'] / gpred(MODEL, th0, a))
        lo, hi = [-2., 0.], [2., 30.]; att = []
        for s0 in ((.1, .6), (.2, 1.5), (0., .2), (.05, .3)):
            o = least_squares(fun, list(s0), bounds=(lo, hi), max_nfev=600, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        good = [x for x in att if x['success']] or att; best = min(good, key=lambda x: x['loss'])
        bh = [n for n, v, l, u in zip(['t', 'kappa'], best['th'], lo, hi) if min(abs(v - l), abs(v - u)) < 1e-5]
        _CACHE[key] = ((float(best['th'][0]), float(best['th'][1])), best['loss'], att, bh, int(len(a['r'])))
    return _CACHE[key]
_RULE = {}
def fit(model, weight, dtr, etr, lens, groups):
    if model != MODEL: return _fit(model, weight, dtr, etr, lens, groups)
    k = (len(dtr['r']), float(np.sum(dtr['go'])))
    if k not in _RULE: _RULE[k] = np.array(_fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])
    abc = _RULE[k]; (t, kap), loss, att, bh, npts = fit_disks(np.unique(dtr['name']).tolist(), abc); B.FROZEN['abc'] = np.array(abc); U['t'] = t; U['k'] = kap
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in abc], t=float(t), kappa=float(kap), disk_loss=loss, disk_points=npts)
H.fit, H.gpred, H.esd = fit, gpred, esd
