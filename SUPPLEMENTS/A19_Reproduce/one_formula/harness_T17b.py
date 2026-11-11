"""T17b: the paper's rule (frozen on training disks) + the published Moster (2013) halo mass, switched on by size (MOSS)
or by stellar mass and size (MOSH); Dutton-Maccio concentration at each redshift; Moline factor for satellites.
Imports harness_T17 (and through it T16, T15c, T15b and the hashed T15 harness) unchanged.
Usage: python3 harness_T17b.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T17 as X
H = X.H; B = X.B; C5 = X.C5; K = H.K
MC = {}
def moster(ms, z):
    out = np.empty(np.size(ms))
    for i, (m, zz) in enumerate(zip(np.ravel(ms), np.ravel(z))):
        k = (float(m), float(zz))
        if k not in MC: MC[k] = float(C5.moster_halo(m, zz)[0])
        out[i] = MC[k]
    return out.reshape(np.shape(ms))
NEW = ('MOSS', 'MOSH')
H.MODELS['MOSS'] = (['logRt'], [0.], [2.5], [[1.], [1.3], [1.6]])
H.MODELS['MOSH'] = (['h', 'g'], [9., -3.], [30., 12.], [[11.5, 0.], [11.8, .6], [14., 4.]])
def switch(model, th, ms, re):
    if model == 'MOSS': return 1 / (1 + (10 ** th[0] / re) ** 8)
    return 1 / (1 + 10 ** (-2 * (np.log10(ms) + th[1] * np.log10(re / 3) - th[0])))
def mass(model, th, ms, mp, re, z):
    abc = B.FROZEN['abc']
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + switch(model, th, ms, re) * moster(ms, z)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in NEW: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in NEW: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)          # isolated lenses: no subhalo factor
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c
    rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
