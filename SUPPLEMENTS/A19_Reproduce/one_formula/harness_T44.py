"""T44 (PROTOCOL_T44.txt): T43's formula with the KiDS lens stacks' halo concentrations scaled by Neto et al. (2007)'s
all-to-relaxed ratio of median concentrations, (4.67/5.26) (M200 h / 1e14)^(-0.01). Everything else as T43. MOGF0IESR: 4 coefficients."""
import sys
import numpy as np
import harness_T43 as G
S = G.S; H = S.H; K = H.K; T = S.T; Z = S.Z; Y = S.Y
MODEL = 'MOGF0IESR'
H.MODELS[MODEL] = H.MODELS['MOGF0']
def ratio_all_relaxed(m): return (4.67 / 5.26) * (m * S.HH / 1e14) ** (-0.01)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    return _gpred('MOGF0IESG' if model == MODEL else model, th, a)
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    m = S.mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * ratio_all_relaxed(m)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
