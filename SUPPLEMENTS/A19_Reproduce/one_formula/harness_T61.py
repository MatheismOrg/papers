"""T61 (PROTOCOL_T61.txt): T48's formula with cores in diffuse, gas-rich galaxies. With g = M_gas/M_p (gas share of the
baryons), L = -log10 min(Sigma*/Sigma_p, 1), Sigma* = M*/(2 pi R_e^2), Sigma_p = 1e10/(2 pi 3^2) Msun/kpc^2 and w = g L
(the variable of T57), for light-selected systems:
  M_ref -> M_ref x 10^(t w),    M_DM(<r) = M_NFW(<r) x tanh(r / (1.75 R_e))^n,  n = tanh(kappa w)
(the coreNFW profile of Read, Agertz and Collins 2016, core radius 1.75 half-light radii, with the core's strength set by
w in place of their star-formation time). The term is zero for gas-free systems (dwarf spheroidals, early types), for
dense ones (L = 0) and for the halo-selected clusters. The lensing radii lie far outside the cores, so the lensing takes
only the mass term.
MOGF0IESBC: A, B, C as T48 (the paper's rule on the training disks' benchmark points); then t and kappa on ALL points of
the same training disks (positive radius, speed and baryonic V^2; no cut on the velocity error), each point weighted by
1 / (sigma_V^2 + sigma_0^2), sigma_V = (2 / ln 10) eV / V, sigma_0 = 0.131 dex, each galaxy's weights summing to the same.
M_t as T45."""
import os
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as P8
P = P8.P; H = P.H; S = P.S; K = H.K; B = S.B; P44 = P.P; G43 = P44.G; T = P44.T; Z = P44.Z; D = H.D
MODELS = ('MOGF0IESBC',)
for m_ in MODELS: H.MODELS[m_] = H.MODELS['MOGF0']
U = {'t': 0.0, 'k': 0.0}; ETA = 1.75; S0 = 0.131; SP = 1e10 / (2 * np.pi * 9.); _CACHE = {}
def gl(ms, mp, re, thr):
    thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    return np.where(thr > 0, 0., np.clip(1 - ms / mp, 0, 1) * -np.log10(np.minimum(ms / (2 * np.pi * re ** 2) / SP, 1.)))
def all_points(keep=lambda g, r: True):
    """every SPARC point with positive radius, speed and baryonic V^2, with its fractional velocity error ('fe')"""
    t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); rows = []; fe = []
    for g, r in sorted(t1.items()):
        if not keep(g, r): continue
        d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
        R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
        vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (vb2 > 0)
        Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']
        rows += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
        fe += list((eV / V)[k])
    a = H.arr(rows); a['fe'] = np.array(fe); return a
ALL = all_points()
def pick(a, k): return {n: v[k] for n, v in a.items()}
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model not in MODELS: return _gpred(model, th, a)
    thr = S.thresholds(a); w = gl(a['ms'], a['mp'], a['re'], thr)
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], thr) * 10 ** (U['t'] * w); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    with np.errstate(all='ignore'): core = np.where(w > 0, np.tanh(a['r'] / (ETA * a['re'])) ** np.tanh(U['k'] * w), 1.)
    return a['gb'] + H.G * (1 - H.FB) * mt * core / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in MODELS: return _esd(model, th, st)
    w = gl(st.ms, st.mp, st.re, 0.); m = G43.mass(th, st.ms, st.mp, st.re, st.z, 0.) * 10 ** (U['t'] * w)
    c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def weights(a):
    """per-point square-root weights: inverse variance, each galaxy's weights summing to 1 / (number of galaxies)"""
    un, gi = np.unique(a['name'], return_inverse=True); pw = 1 / ((2 / np.log(10) * a['fe']) ** 2 + S0 ** 2)
    return np.sqrt(pw / np.bincount(gi, pw)[gi] / len(un))
def fit_disks(model, names, abc0):
    key = (model, tuple(sorted(names)), tuple(np.round(abc0, 9)))
    if key not in _CACHE:
        a = pick(ALL, np.isin(ALL['name'], list(names))); assert len(np.unique(a['name'])) == len(set(names))
        nw = weights(a); th0 = np.array([11.3]); B.FROZEN['abc'] = np.asarray(abc0)
        def fun(p):
            U['t'], U['k'] = float(p[0]), float(p[1]); return nw * np.log10(a['go'] / gpred(model, th0, a))
        lo, hi = [-2., 0.], [2., 30.]; att = []
        for s0 in ((.1, .6), (.2, 1.5), (0., .2)):
            o = least_squares(fun, list(s0), bounds=(lo, hi), max_nfev=600, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        good = [x for x in att if x['success']] or att; best = min(good, key=lambda x: x['loss'])
        bh = [n for n, v, l, u in zip(['t', 'kappa'], best['th'], lo, hi) if min(abs(v - l), abs(v - u)) < 1e-5]
        _CACHE[key] = ((float(best['th'][0]), float(best['th'][1])), best['loss'], att, bh, int(len(a['r'])))
    return _CACHE[key]
_RULE = {}
def fit(model, weight, dtr, etr, lens, groups):
    if model not in MODELS: return _fit(model, weight, dtr, etr, lens, groups)
    k = (len(dtr['r']), float(np.sum(dtr['go'])))
    if k not in _RULE: _RULE[k] = np.array(_fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])          # the paper's rule on these training disks, as chunk_split sets it
    abc = _RULE[k]; (t, kap), loss, att, bh, npts = fit_disks(model, np.unique(dtr['name']).tolist(), abc); B.FROZEN['abc'] = np.array(abc); U['t'] = t; U['k'] = kap
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in abc], t=float(t), kappa=float(kap), disk_loss=loss, disk_points=npts)
H.fit, H.gpred, H.esd = fit, gpred, esd
