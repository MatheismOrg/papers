"""T17 all-data diagnostic (not a protocol run): SUBH with the published Einasto halo (Dutton and Maccio 2014 shape and concentration) in place of NFW.
Imports harness_T17 (and through it T16, T15c, T15b and the hashed T15 harness) unchanged.
Usage: python3 harness_T17_einasto.py <mask> <out> cfg..."""
import sys
import numpy as np
import harness_T17 as X
import einasto as E
H = X.H; B = X.B
NEW = {'EIN': False, 'EINJ': True}         # name: refit the paper's A, B, C jointly?
H.MODELS['EIN'] = H.MODELS['KNEESF']
H.MODELS['EINJ'] = (['A', 'B', 'C', 'A2', 'B2', 'C2'], [9, -.5, -2.5, 9, 0, -2], [14.5, 2, 1.5, 16, 5, 3],
                    [[11.7481, .8178, -.7119, 12.06, 2.63, 2.5], [11.7481, .8178, -.7119, 12.5, 2.0, 1.], [11.6, .8, -.8, 12.0, 3.0, 2.0]])
_gpred, _esd = H.gpred, H.esd
def mass(model, th, ms, mp, re):
    abc, br = (th[:3], th[3:]) if NEW[model] else (B.FROZEN['abc'], th)
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + 10 ** br[0] * (ms / 1e11) ** br[1] * (re / 10) ** br[2]
def gpred(model, th, a):
    if model not in NEW: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re']); c = E.c_ein(m, a['z']) * a.get('fac', 1.); al = E.alpha(m, a['z'])
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * E.enc_ratio(x, c, al) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in NEW: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re); c = E.c_ein(m, H.K.Z0)              # isolated lenses: no subhalo factor
    al = np.clip(E.alpha(m, H.K.Z0), E.A_GRID[0], E.A_GRID[-1])                  # clip only beyond ~1e16.5 Msun
    rc = H.K.rho_crit(H.K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); r2 = rv / c
    norm = (1 - H.FB) * m / (4 * np.pi * r2 ** 2 * np.exp(E.log_I(c, al)))        # rho_-2 r_-2 for dark mass (1-fb) M inside r200
    xx = np.clip(st.R / r2[:, None, None], 10 ** E.LX_GRID[0], 10 ** E.LX_GRID[-1])
    e = norm[:, None, None] * E.D_ein(xx, np.broadcast_to(al[:, None, None], xx.shape)) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
