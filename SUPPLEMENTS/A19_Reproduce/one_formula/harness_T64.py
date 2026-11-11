"""T64 (PROTOCOL_T64.txt): the formula without its core (T48, MOGF0IESB) with, in EVERY halo, the inner structure that a
kind of dark-matter particle would give, one number each, in place of the gas-gated core of T61.
  MOGF0IESBF, fuzzy dark matter of particle mass m22 (in 1e-22 eV): the soliton of Schive et al. (2014, PRL 113, 261302),
     rho(r) = 1.9e7 m22^-2 (r_c/kpc)^-4 [1 + 0.091 (r/r_c)^2]^-8 Msun/kpc^3 (their eq. 3, physical units),
     r_c = 1.6 kpc m22^-1 a^(1/2) [zeta(z)/zeta(0)]^(-1/6) (M_h / 1e9 Msun)^(-1/3) (their eq. 7), M_h the formula's halo mass,
     inside r_t = RT r_c (RT = 3; 95 percent of the soliton's mass); outside, the formula's NFW halo, mass continuous:
     M_dm(r) = M_sol(r) for r <= r_t, and M_sol(r_t) + M_NFW(r) - M_NFW(r_t) beyond. Lensing unchanged (radii far outside).
  MOGF0IESBW, warm dark matter of half-mode mass M_hm (Msun): concentration c -> c (1 + 15 M_hm / M)^(-0.3)
     (Schneider et al. 2012, MNRAS 424, 684, Section 4.5), on the galaxies and in the lensing alike.
     Thermal relic mass: M_hm = 1.3e10 / h Msun (m / keV)^(-3.33) (their Table 1 and eq. 5).
A, B, C as T48 (the paper's rule on the training disks' benchmark points); then the particle number on ALL points of the
same training disks with the weights of T61 (harness_T61.weights). M_t as T45. Nothing else changes."""
import numpy as np
from scipy.optimize import least_squares
import harness_T61 as M61
H = M61.H; S = M61.S; K = M61.K; B = M61.B; P = M61.P; P44 = M61.P44; G43 = M61.G43; T = M61.T; Z = M61.Z
ALL = M61.ALL; pick = M61.pick; weights = M61.weights
MODELS = ('MOGF0IESBF', 'MOGF0IESBW')
for m_ in MODELS: H.MODELS[m_] = H.MODELS['MOGF0']
V = {'lm22': 6.0, 'lhm': 3.0}; RT = 3.0; BND = {'MOGF0IESBF': (-1., 6.), 'MOGF0IESBW': (3., 12.5)}; KEY = {'MOGF0IESBF': 'lm22', 'MOGF0IESBW': 'lhm'}
START = {'MOGF0IESBF': (0., 0.5, 1., 1.5, 2.5, 4.), 'MOGF0IESBW': (7., 9., 10., 11., 12.)}; _CACHE = {}
OM = .3175
def zeta(z):
    z = np.asarray(z, float); om = OM * (1 + z) ** 3 / (OM * (1 + z) ** 3 + 1 - OM); return (18 * np.pi ** 2 + 82 * (om - 1) - 39 * (om - 1) ** 2) / om
_X = np.linspace(0., 4., 40001); _I = 4 * np.pi * _X ** 2 * (1 + 0.091 * _X ** 2) ** -8; _F = np.concatenate([[0.], np.cumsum((_I[1:] + _I[:-1]) / 2 * np.diff(_X))])
def fsol(x):
    """soliton mass inside x = r / r_c, in units of rho_c r_c^3"""
    x = np.asarray(x, float); return np.where(x < 0.05, 4 * np.pi * (x ** 3 / 3 - 8 * 0.091 * x ** 5 / 5), np.interp(x, _X, _F))
def halo(th, a):
    thr = S.thresholds(a); m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], thr); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); return m, c, rv
def soliton(m, z, m22):
    rc = 1.6 / m22 * (1 + np.asarray(z, float)) ** -0.5 * (zeta(z) / zeta(0.)) ** (-1 / 6) * (m / 1e9) ** (-1 / 3); return rc, 1.9e7 * m22 ** -2 * rc ** -4
def mdm_fuzzy(m, c, rv, z, r, m22, rt=None):
    rt = RT if rt is None else rt; nfw = lambda x: (1 - H.FB) * m * H.ff(x * c / rv) / H.ff(c); rc, rhoc = soliton(m, z, m22); k = rhoc * rc ** 3
    return np.where(r <= rt * rc, k * fsol(r / rc), k * fsol(rt) + nfw(r) - nfw(rt * rc))
def wdm_fac(m, mhm): return (1 + 15 * mhm / m) ** -0.3
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model not in MODELS: return _gpred(model, th, a)
    m, c, rv = halo(th, a)
    if model == 'MOGF0IESBF': md = mdm_fuzzy(m, c, rv, a['z'], a['r'], 10 ** V['lm22'])
    else: c = c * wdm_fac(m, 10 ** V['lhm']); md = (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * md / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in MODELS: return _esd(model, th, st)
    if model == 'MOGF0IESBF': return _esd('MOGF0IESB', th, st)
    m = G43.mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m) * wdm_fac(m, 10 ** V['lhm'])
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3); T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def disk_loss(model, a, nw, x):
    V[KEY[model]] = float(x); r = nw * np.log10(a['go'] / gpred(model, np.array([11.3]), a)); return float(r @ r)
def fit_disks(model, names, abc0):
    key = (model, tuple(sorted(names)), tuple(np.round(abc0, 9)))
    if key not in _CACHE:
        a = pick(ALL, np.isin(ALL['name'], list(names))); assert len(np.unique(a['name'])) == len(set(names))
        nw = weights(a); th0 = np.array([11.3]); B.FROZEN['abc'] = np.asarray(abc0); lo, hi = BND[model]; att = []
        def fun(p):
            V[KEY[model]] = float(p[0]); return nw * np.log10(a['go'] / gpred(model, th0, a))
        for s0 in START[model]:
            o = least_squares(fun, [s0], bounds=([lo], [hi]), max_nfev=400, ftol=1e-10, xtol=1e-10, gtol=1e-10, diff_step=1e-3)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        grid = np.linspace(lo, hi, 81); gl_ = [disk_loss(model, a, nw, x) for x in grid]; j = int(np.argmin(gl_))      # a grid as well: the loss is flat at the cold end and can have two minima
        att.append(dict(th=[float(grid[j])], loss=float(gl_[j]), success=True, nfev=81))
        best = min(att, key=lambda x: x['loss']); cold = disk_loss(model, a, nw, hi if model == 'MOGF0IESBF' else lo)
        bh = [KEY[model]] if min(abs(best['th'][0] - lo), abs(best['th'][0] - hi)) < 1e-5 else []
        _CACHE[key] = (float(best['th'][0]), best['loss'], att, bh, int(len(a['r'])), cold)
    return _CACHE[key]
_RULE = {}
def fit(model, weight, dtr, etr, lens, groups):
    if model not in MODELS: return _fit(model, weight, dtr, etr, lens, groups)
    k = (len(dtr['r']), float(np.sum(dtr['go'])))
    if k not in _RULE: _RULE[k] = np.array(_fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])          # the paper's rule on these training disks, as chunk_split sets it
    abc = _RULE[k]; x, loss, att, bh, npts, cold = fit_disks(model, np.unique(dtr['name']).tolist(), abc); B.FROZEN['abc'] = np.array(abc); V[KEY[model]] = x
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(v) for v in abc], par=float(x), disk_loss=loss, disk_loss_cold=cold, disk_points=npts)
H.fit, H.gpred, H.esd = fit, gpred, esd
