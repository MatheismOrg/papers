"""T19: one formula, two published halo relations joined by a physical property, optionally with halo contraction.
  M_ref = Mp/fb + 10^A (M*/1e10)^B (Re/3 kpc)^C + S x M_Moster(M*, z)
MOG, MOGC: S = (M*/Mp)^50 / (1 + (M_t/M*)^4): gas-poor (quenched) galaxies above M_t follow abundance matching
           (passive centrals have halos at least twice as massive as star-forming ones of the same stellar mass:
           Mandelbaum et al. 2016). One coefficient, log10 M_t in [10, 12.5]; 50 and 4 fixed.
MOSC:      S = 1 / (1 + (R_t/Re)^8) as T17b's MOSS.
MOGC, MOSC add contraction where baryons dominate (Cautun et al. 2020 eq. 11, contraction branch only, F >= 1):
  M_DM(<r) = (1 - fb) M_NFW(<r) max(F(eta), 1)^s, F = 0.45 + 0.38 (eta + 1.16)^0.53, eta = M_bar(<r)/(fb M_NFW(<r)),
  s in [0, 1.5] fitted with the switch coefficient. The expansion branch (eta < 1) is left out because g_bar omits hot gas.
Dynamics: M_bar(<r) = g_bar r^2/G. Lensing: lens M_p as a point mass, projected as in T18 (harness_T18.esd_dm).
The paper's rule A, B, C is frozen on each fold's training disks without contraction (as T17b), via harness_T15c.run.
Usage: python3 harness_T19.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T18 as T, harness_T17b as Y
H = T.H; B = T.B; K = H.K
F_ORIG, DF_ORIG = T.F, T.dF
def Fc(e): return np.maximum(F_ORIG(e), 1.)
def dFc(e): return np.where(F_ORIG(e) > 1., DF_ORIG(e), 0.)
MODELS = {'MOG': (['logMt'], [10.], [12.5], [[11.2], [11.0], [11.4]]),
          'MOGC': (['logMt', 's'], [10., 0.], [12.5, 1.5], [[11.2, .3], [11.0, 0.], [11.4, .6]]),
          'MOSC': (['logRt', 's'], [0., 0.], [2.5, 1.5], [[1.15, .3], [1.15, 0.], [1.3, .6]])}
H.MODELS.update(MODELS)
def switch(model, th, ms, mp, re):
    if model == 'MOSC': return 1 / (1 + (10 ** th[0] / re) ** 8)
    return (ms / mp) ** 50 / (1 + (10 ** th[0] / ms) ** 4)
def mass(model, th, ms, mp, re, z):
    abc = B.FROZEN['abc']
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + switch(model, th, ms, mp, re) * Y.moster(ms, z)
def strength(model, th): return 0. if model == 'MOG' else th[1]
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in MODELS: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Fc(eta) ** strength(model, th) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in MODELS: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Fc, dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, strength(model, th))
    finally: T.F, T.dF = F_ORIG, DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
