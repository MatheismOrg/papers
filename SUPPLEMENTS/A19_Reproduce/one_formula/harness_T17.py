"""T17: SUBH + published Moster redshift evolution of the high-mass branch, and/or gas-fraction concentration.
Imports harness_T16 (and through it the earlier harnesses) unchanged. Usage: python3 harness_T17.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T16 as T
H = T.H; B = T.B; C5 = H.C5
CACHE = {}
def dz(ms, z):
    """log10 M_Moster(M*, z) - log10 M_Moster(M*, 0) at fixed stellar mass (0 at z = 0)."""
    out = np.zeros(np.shape(ms))
    for i, (m, zz) in enumerate(zip(np.ravel(ms), np.ravel(z))):
        if zz <= 0: continue
        key = (float(m), float(zz))
        if key not in CACHE:
            CACHE[key] = float(np.log10(C5.moster_halo(m, zz)[0]) - np.log10(C5.moster_halo(m, 0.)[0]))
        out.flat[i] = CACHE[key]
    return out
NEW = {'MZ': (True, False), 'GASC': (False, True), 'MZGASC': (True, True)}     # (Moster evolution, gas concentration)
H.MODELS['MZ'] = H.MODELS['KNEESF']
H.MODELS['GASC'] = (['A2', 'B2', 'C2', 'f'], [9, 0, -2, -1], [16, 5, 3, 1], [[12.06, 2.63, 2.5, 0.], [12.06, 2.63, 2.5, -.3], [12.06, 2.63, 2.5, .3]])
H.MODELS['MZGASC'] = H.MODELS['GASC']
_gpred, _esd = H.gpred, H.esd
def parts(model, th, ms, mp, re, z):
    abc = B.FROZEN['abc']; mz, gas = NEW[model]
    branch = 10 ** th[0] * (ms / 1e11) ** th[1] * (re / 10) ** th[2]
    if mz: branch = branch * 10 ** dz(ms, z)
    m = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + branch
    lc = (th[3] * (mp - ms) / mp) if gas else 0.
    return m, lc
def gpred(model, th, a):
    if model not in NEW: return _gpred(model, th, a)
    m, lc = parts(model, th, a['ms'], a['mp'], a['re'], a['z']); c = H.K.c_dm14(m, a['z']) * 10 ** lc * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in NEW: return _esd(model, th, st)
    m, lc = parts(model, th, st.ms, st.mp, st.re, st.z); c = H.K.c_dm14(m, H.K.Z0) * 10 ** lc
    rc = H.K.rho_crit(H.K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c
    rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * H.K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
