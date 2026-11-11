"""T56 (PROTOCOL_T56.txt): T48's formula with the cold gas of light-selected galaxies weighted mu times in the baryon term,
  M_ref = (M* + mu M_gas)/f_b + 10^A (M*/1e10)^B (R_e/3)^C + q M_E(M*, z | S),   floor (M^8 + (1e8)^8)^(1/8),
mu fitted together with A, B, C on each fold's training disks (bounds of the paper's rule for A, B, C; 0 <= mu <= 20).
The hot gas of the halo-selected clusters keeps weight 1 and the switch q keeps the true M_p. Everything else as T48
(M_t on the red centrals of Mandelbaum et al. 2016, IMF trend, all-halo concentrations for the stacks, Moline factor).
MOGF0IESBU: 4 coefficients fitted to the benchmark (A, B, C, mu on the training disks)."""
import sys
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as P8
P = P8.P; H = P.H; S = P.S; K = H.K; B = S.B; P44 = P.P; T = P44.T; Z = P44.Z; Y = P.Y
MODEL = 'MOGF0IESBU'; H.MODELS[MODEL] = H.MODELS['MOGF0']
U = {'mu': 1.0}; _CACHE = {}
LO, HI = [9, -.5, -2.5, 0.], [14.5, 2, 1.5, 20.]
def mass(th, ms, mp, re, z, thr):
    abc = B.FROZEN['abc']; thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    q = np.where(thr > 0, 1., Z.switch('MOG', th, ms, mp, re)); mpe = np.where(thr > 0, mp, ms + U['mu'] * (mp - ms))
    m = mpe / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + q * S.m_exp(ms, z, thr)
    return (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    m = mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    m = mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def fit_disks(dtr):
    key = (len(dtr['r']), float(np.sum(dtr['go'])), float(np.sum(dtr['r'])))
    if key not in _CACHE:
        th0 = np.array([11.3])
        def fun(p):
            B.FROZEN['abc'] = np.asarray(p[:3]); U['mu'] = float(p[3]); return dtr['hw'] * np.log10(dtr['go'] / gpred(MODEL, th0, dtr)) / .15
        att = []
        for st in ([11.7481, .8178, -.7119, 1.], [11.7481, .8178, -.7119, 4.], [11.6, 1., -1.3, 4.]):
            o = least_squares(fun, st, bounds=(LO, HI), max_nfev=400, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
        good = [x for x in att if x['success']] or att; best = min(good, key=lambda x: x['loss'])
        bh = [n for n, t, l, u in zip(['A', 'B', 'C', 'mu'], best['th'], LO, HI) if min(abs(t - l), abs(t - u)) < 1e-5]
        _CACHE[key] = (best['th'], best['loss'], att, bh)
    return _CACHE[key]
def fit(model, weight, dtr, etr, lens, groups):
    if model != MODEL: return _fit(model, weight, dtr, etr, lens, groups)
    p, loss, att, bh = fit_disks(dtr); B.FROZEN['abc'] = np.array(p[:3]); U['mu'] = float(p[3])
    lmt, c2, edge = P.calibrate()
    return dict(th=[lmt], loss=c2, attempts=att, bound_hits=bh + (['logMt'] if edge else []), calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in p[:3]], mu=float(p[3]), disk_loss=loss)
H.fit, H.gpred, H.esd = fit, gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
