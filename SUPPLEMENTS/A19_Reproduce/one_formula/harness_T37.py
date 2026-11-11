"""T37 (PROTOCOL_T37.txt): the frozen T26 formula (MOGF0IE) with its high-mass branch replaced by the expected halo of a
galaxy given its starlight and its sample's selection (T31 to T33 built into the formula; no coefficient added):
  M_ref = M_p/f_b + 10^A (M*/1e10)^B (R_e/3 kpc)^C + q M_E(M*, z | S),  q as T19 (gas-poor and above M_t), floor as T21,
  M_E = sum_M M w(M) / sum_M w(M) on a grid in log10 M200c (9.0 to 16.5, step 0.005),
  w(M) = dn/dlnM (Tinker et al. 2008; colossus; A19 Planck cosmology as T32) x Gaussian in log10 M* about the Moster et al.
  (2013) M*(M, z) with sigma = sqrt(0.15^2 + 0.10^2) dex x S(M);
  S(M) = [M200c > 5e14] for CLASH (Postman et al. 2012), [M200c > 1e14] for the MaNGA BCGs (program aim, as T33), 1 otherwise.
  Below M* = 10^9.5 suns M_E is T26's mean inversion (the switch is below 1e-6 there).
IMF factors as T26 (imported first). Fits as T26: A, B, C on each fold's training disks (q ~ 0 for disks), log10 M_t
jointly at population weight 0.1. MOGF0IES: 4 coefficients.
Usage: python3 harness_T37.py <mask> <out> MOGF0IES@0.1   (or chunk_split.py harness_T37 ...);  describe24.py harness_T37 MOGF0IES"""
import sys
import numpy as np
import harness_T26 as M
F = M.F; H = M.H; Y = M.Y; Z = F.Z; K = H.K; T = Z.T; B = Z.B
sys.path.insert(0, str(H.PKG / 'code')); import A19_Clusters_v5 as V5
from colossus.lss import mass_function
from colossus.cosmology import cosmology
cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624})
cosmology.setCosmology('A19planck')
MODEL = 'MOGF0IES'
H.MODELS[MODEL] = H.MODELS['MOGF0']
HH = 0.671; GRID = np.arange(9.0, 16.5 + 1e-9, 0.005); MG = 10 ** GRID; SIG = np.hypot(.15, .10); LMSMIN = 9.5; LFLOOR = 8.0
SEL = {'CLASH_centre': 5e14, 'CLASH_outer': 5e14, 'MaNGA_BCG': 1e14}
_HMF, _ME = {}, {}
def hmf(z):
    k = round(float(z), 6)
    if k not in _HMF: _HMF[k] = mass_function.massFunction(MG * HH, k, mdef='200c', model='tinker08', q_out='dndlnM')
    return _HMF[k]
def m_exp1(ms, z, thr):
    k = (float(ms), float(z), float(thr))
    if k not in _ME:
        if ms < 10 ** LMSMIN: _ME[k] = float(Y.moster(np.array([ms]), np.array([z]))[0])
        else:
            ll = -.5 * ((np.log10(ms) - np.log10(V5.moster_mstar(MG, z))) / SIG) ** 2; sel = MG > thr
            w = np.where(sel, hmf(z) * np.exp(ll - ll[sel].max()), 0.)
            _ME[k] = float(np.sum(w * MG) / np.sum(w))
    return _ME[k]
def m_exp(ms, z, thr):
    ms, z, thr = np.broadcast_arrays(np.asarray(ms, float), np.asarray(z, float), np.asarray(thr, float))
    return np.array([m_exp1(a, b, c) for a, b, c in zip(ms.ravel(), z.ravel(), thr.ravel())]).reshape(ms.shape)
def thresholds(a): return np.array([SEL.get(k, 0.) for k in a['kind']]) if 'kind' in a else 0.
def mass(th, ms, mp, re, z, thr):
    abc = B.FROZEN['abc']
    m = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + Z.switch('MOG', th, ms, mp, re) * m_exp(ms, z, thr)
    return (m ** 8 + (10. ** LFLOOR) ** 8) ** (1 / 8)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    m = mass(th, a['ms'], a['mp'], a['re'], a['z'], thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    m = mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
