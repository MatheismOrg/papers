"""All-data diagnostic (not held out): the T17b switch formula (MOSS) with the Di Cintio et al. (2014) profile in place of NFW;
the paper's rule refitted on disks (no switch) with the same profile. Optional taper to NFW above a virial mass."""
import numpy as np, sys, time
import harness_T17b as Y, dc14 as DC
from scipy.optimize import least_squares
H = Y.H; B = Y.B; K = H.K
SET = {'taper': None, 'switch': True}
def halo(m, ms, z, fac=1.):
    c0 = K.c_dm14(m, z) * fac
    mv = DC.m200_to_mvir(m, c0, z); X = np.log10(ms / mv)
    p = DC.params(X)
    if SET['taper']: p = DC.blend(p, 1 / (1 + (mv / 10 ** SET['taper']) ** 4))
    a, b, g, cf = p
    return a, b, g, c0 * cf
def mass(th, a_):
    if SET['switch']: return Y.mass('MOSS', th, a_['ms'], a_['mp'], a_['re'], a_['z'])
    abc = B.FROZEN['abc']; return a_['mp'] / H.FB + 10 ** abc[0] * (a_['ms'] / 1e10) ** abc[1] * (a_['re'] / 3) ** abc[2]
def gp(th, d):
    m = mass(th, d); a, b, g, c = halo(m, d['ms'], d['z'], d.get('fac', 1.))
    rv = (H.G * m / (100 * H.hz(d['z']) ** 2)) ** (1 / 3); rs = rv / c * DC.rs_over_rm2(a, b, g)
    return d['gb'] + H.G * (1 - H.FB) * m * DC.I(d['r'] / rs, a, b, g) / DC.I(rv / rs, a, b, g) / d['r'] ** 2 * H.ACC
class St: pass
def es(th, st):
    d = dict(ms=st.ms, mp=st.mp, re=st.re, z=st.z); m = mass(th, d); a, b, g, c = halo(m, st.ms, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3); rs = rv / c * DC.rs_over_rm2(a, b, g)
    norm = (1 - H.FB) * m / (4 * np.pi * rs ** 2 * DC.I(rv / rs, a, b, g))
    e = norm[:, None, None] * DC.D(st.R / rs[:, None, None], a[:, None, None], b[:, None, None], g[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
from scipy.linalg import solve_triangular
def lens_chi2(th, lens):
    return [float(np.sum(solve_triangular(lens.L[j], lens.y[j] - es(th, lens.stacks[j])[lens.masks[j]], lower=True) ** 2)) for j in range(4)]
if __name__ == '__main__':
    import diag_mosg as G
    d, ext, lens = G.d, G.ext, G.lens
    th = np.array([1.14548])
    for taper in [None, 12.0, 12.5]:
        SET['taper'] = taper; SET['switch'] = False
        def fun(p):
            B.FROZEN['abc'] = p; return d['hw'] * np.log10(d['go'] / gp(th, d)) / .15
        t = time.time(); o = least_squares(fun, [11.7481, .8178, -.7119], bounds=([9, -.5, -2.5], [14.5, 2, 1.5])); B.FROZEN['abc'] = o.x
        SET['switch'] = True
        z = lens_chi2(th, lens)
        print('DC14 taper %s rule %s (%.0fs) | disk %.4f' % (taper, np.round(o.x, 3), time.time() - t, G.sc(d, gp(th, d))), ' '.join('%s %.4f' % (n, G.sc(a, gp(th, a))) for n, a in ext.items()),
              '| lens', np.round(z, 1), 'sum %.1f' % sum(z), flush=True)
