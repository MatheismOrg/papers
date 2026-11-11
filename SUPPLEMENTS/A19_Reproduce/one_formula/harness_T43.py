"""T43 (PROTOCOL_T43.txt): T37's formula with the switch q = 1 for samples selected by their halos (CLASH, MaNGA BCGs);
the gas and stellar-mass gate applies only to samples chosen by their light. MOGF0IESG: 4 coefficients as T37."""
import sys
import numpy as np
import harness_T37 as S
H = S.H; K = H.K; T = S.T; Z = S.Z; B = S.B; Y = S.Y
MODEL = 'MOGF0IESG'
H.MODELS[MODEL] = H.MODELS['MOGF0']
def mass(th, ms, mp, re, z, thr):
    abc = B.FROZEN['abc']; thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    q = np.where(thr > 0, 1., Z.switch('MOG', th, ms, mp, re))
    m = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + q * S.m_exp(ms, z, thr)
    return (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    m = mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    return _esd('MOGF0IES' if model == MODEL else model, th, st)       # lenses are chosen by their light: identical to T37
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
