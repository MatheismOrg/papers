"""T22: T19's MOGC plus a hot gas halo holding a fraction f of each system's missing baryons, fb M_ref - Mp, with
density proportional to r^-1.5 inside r200 (the Milky Way's hot halo: beta = 0.50 +/- 0.03, holding at most about half of
its missing baryons; Miller and Bregman 2015). Dynamics: g_hot = G M_hot(<r)/r^2, M_hot(<r) = f (fb M_ref - Mp) (r/r200)^1.5.
Lensing: projected numerically (r = R cosh t quadrature, as T18). CLASH points get no hot term: their g_bar already
includes the intracluster gas. The contraction factor keeps T19's eta (galaxy baryons). The rule stays frozen on disks.
MOGH5:  f = 0.5 (Miller and Bregman's upper bound)              (5 coefficients, as MOGC)
MOGHf:  f fitted in [0, 1]                                        (6)
MOGFH5: MOGH5 + T21's floor M_min = 1e8 Msun                      (5)
Usage: python3 harness_T22.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T21 as F
Z = F.Z; H = Z.H; B = Z.B; K = H.K; T = Z.T; Y = Z.Y
HOT = {'MOGH5': ('MOGC', .5, False), 'MOGHf': ('MOGC', None, False), 'MOGFH5': ('MOGF', .5, True)}
H.MODELS['MOGH5'] = Z.MODELS['MOGC']; H.MODELS['MOGFH5'] = Z.MODELS['MOGC']
H.MODELS['MOGHf'] = (['logMt', 's', 'f'], [10., 0., 0.], [12.5, 1.5, 1.], [[11.2, .17, .5], [11.0, 0., 0.], [11.4, .4, 1.]])
def frac(model, th): return HOT[model][1] if HOT[model][1] is not None else th[2]
def base_mass(model, th, ms, mp, re, z):
    m = Z.mass('MOGC', th, ms, mp, re, z)
    return (m ** 8 + 1e8 ** 8) ** (1 / 8) if HOT[model][2] else m
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in HOT: return _gpred(model, th, a)
    m = base_mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    g = a['gb'] + H.G * (1 - H.FB) * mt * Z.Fc(eta) ** th[1] / a['r'] ** 2 * H.ACC
    mh = frac(model, th) * np.maximum(H.FB * m - a['mp'], 0.) * np.minimum(a['r'] / rv, 1.) ** 1.5
    clash = np.char.startswith(np.asarray(a['kind']).astype(str), 'CLASH')
    return g + np.where(clash, 0., H.G * mh / a['r'] ** 2 * H.ACC)
TN, TW = np.polynomial.legendre.leggauss(48); TMAX = 9.; TT = .5 * TMAX * (TN + 1); WT = .5 * TMAX * TW
CH, SH, TH = np.cosh(TT), np.sinh(TT), np.tanh(TT)
def esd(model, th, st):
    if model not in HOT: return _esd(model, th, st)
    m = base_mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, th[1])
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    e = e + hot_esd(frac(model, th) * np.maximum(H.FB * m - st.mp, 0.), rv, st.R)
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def hot_esd(tot, rv, R):
    """Excess surface density (Msun/pc^2) of a sphere of mass tot with density ~ r^-1.5 inside rv, zero outside;
    t-quadrature on [0, arccosh(rv/R)] so the edge is resolved."""
    R4 = R[..., None]; rv4 = rv[:, None, None, None]; tot4 = tot[:, None, None, None]
    tc = np.arccosh(np.maximum(rv4 / R4, 1.)); t = .5 * tc * (TN + 1); w = .5 * tc * TW
    r_ = R4 * np.cosh(t); dm = tot4 * 1.5 * r_ ** .5 / rv4 ** 1.5
    sigma = np.sum(dm / np.cosh(t) * w, axis=-1) / (2 * np.pi * R)
    mcyl = (tot4 * np.minimum(R4 / rv4, 1.) ** 1.5)[..., 0] + np.sum(dm * (1 - np.tanh(t)) * R4 * np.sinh(t) * w, axis=-1)
    return (mcyl / (np.pi * R ** 2) - sigma) / 1e12
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
