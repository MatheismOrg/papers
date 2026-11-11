"""T20: robustness of T19's gas switch to its fixed choices, and the two-way contraction variant.
S = (M*/Mp)^k / (1 + (M_t/M*)^n) x M_Moster as in harness_T19 (T19 used k = 50, n = 4).
MOGa: k = 20, n = 4.  MOGb: k = 100, n = 4.  MOGd: k = 50, n = 2.  MOGe: k = 50, n = 8.   (one coefficient, log10 M_t)
MOGC2: k = 50, n = 4, with the two-way Cautun et al. (2020) response F(eta)^s (expansion where eta < 1 kept), s fitted.
Everything else as harness_T19 (rule frozen on training disks, NFW, Dutton-Maccio c at each z, Moline factor).
Usage: python3 harness_T20.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T19 as Z
H = Z.H; B = Z.B; K = H.K; Y = Z.Y; T = Z.T
VAR = {'MOGa': (20, 4, None), 'MOGb': (100, 4, None), 'MOGd': (50, 2, None), 'MOGe': (50, 8, None), 'MOGC2': (50, 4, 'two-way')}
for k_, v in VAR.items():
    H.MODELS[k_] = Z.MODELS['MOGC'] if v[2] else Z.MODELS['MOG']
def mass(model, th, ms, mp, re, z):
    k, n, _ = VAR[model]; abc = B.FROZEN['abc']
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + (ms / mp) ** k / (1 + (10 ** th[0] / ms) ** n) * Y.moster(ms, z)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in VAR: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    s = th[1] if VAR[model][2] else 0.
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Z.F_ORIG(eta) ** s / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in VAR: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    s = th[1] if VAR[model][2] else 0.
    e = T.esd_dm(m, c, rv, st.mp, st.R, s)                 # T.F is the two-way Cautun function outside harness_T19.esd
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
