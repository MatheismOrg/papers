"""T21: T19's formulas with a floor on the halo mass at the scale where the first galaxies could form.
M_ref -> (M_ref^8 + M_min^8)^(1/8): no galaxy sits in a halo lighter than M_min. A virial temperature of 1e4 K (atomic
cooling) corresponds to about 1e7 to 1e8 Msun at z of 20 to 10 (Johnson 2013, eq. 1); the smallest dwarfs formed before
reionization in such halos and have kept or grown their mass since, apart from tides.
MOGF:   MOGC (T19) + floor M_min = 1e8 Msun          (5 coefficients, as MOGC)
MOGF85: MOGC (T19) + floor M_min = 10^8.5 Msun       (5)
MOGF0:  MOG  (T19) + floor M_min = 1e8 Msun          (4)
The floor touches only the lightest dwarfs (rule masses below about 10^8.5); every other system is far heavier.
Usage: python3 harness_T21.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T19 as Z
H = Z.H; B = Z.B; K = H.K; Y = Z.Y; T = Z.T
FLOOR = {'MOGF': ('MOGC', 8.0), 'MOGF85': ('MOGC', 8.5), 'MOGF0': ('MOG', 8.0)}
for k_, (base, _) in FLOOR.items(): H.MODELS[k_] = Z.MODELS[base]
def mass(model, th, ms, mp, re, z):
    base, lf = FLOOR[model]
    m = Z.mass(base, th, ms, mp, re, z)
    return (m ** 8 + (10. ** lf) ** 8) ** (1 / 8)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in FLOOR: return _gpred(model, th, a)
    base = FLOOR[model][0]
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Z.Fc(eta) ** Z.strength(base, th) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in FLOOR: return _esd(model, th, st)
    base = FLOOR[model][0]
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, Z.strength(base, th))
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
