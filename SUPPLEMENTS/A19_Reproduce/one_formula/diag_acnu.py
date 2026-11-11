"""All-data diagnostic: contraction of strength nu, r_f/r_i = [M_i(r_i)/M_f(r_f)]^nu (nu = 1 is Blumenthal 1986, 0 is none),
on the T17b switch formula; the disk rule refitted for each nu."""
import numpy as np, harness_T17b as Y, diag_mosg as G
from scipy.optimize import least_squares
H = Y.H; B = Y.B; K = H.K
def g_dm(r, gb, mp, m, c, z, nu):
    rv = (H.G * m / (100 * H.hz(z) ** 2)) ** (1 / 3); md = (1 - H.FB) * m
    Mdm = lambda x: md * H.ff(x * c / rv) / H.ff(c)
    if nu == 0: return H.G * Mdm(r) / r ** 2 * H.ACC
    mb = np.minimum(gb * r ** 2 / (H.G * H.ACC), mp); eps = mp / md
    lo, hi = np.log(r) - 4, np.log(r) + 4
    for _ in range(80):
        mid = (lo + hi) / 2; ri = np.exp(mid)
        f = np.log(r / ri) - nu * np.log((1 + eps) * Mdm(ri) / (mb + Mdm(ri)))     # increasing in ri
        lo = np.where(f > 0, mid, lo); hi = np.where(f > 0, hi, mid)
    ri = np.exp((lo + hi) / 2)
    return H.G * Mdm(ri) / r ** 2 * H.ACC
def gp(th, a, nu):
    m = Y.mass('MOSS', th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    return a['gb'] + g_dm(a['r'], a['gb'], a['mp'], m, c, a['z'], nu)
d, ext = G.d, G.ext; th = np.array([1.14548])
for nu in [0., .1, .2, .3, .5, 1.]:
    def fun(p):
        B.FROZEN['abc'] = p; return d['hw'] * np.log10(d['go'] / gp(th, d, nu)) / .15
    o = least_squares(fun, [11.7481, .8178, -.7119], bounds=([9, -.5, -2.5], [14.5, 2, 1.5])); B.FROZEN['abc'] = o.x
    print('nu %.1f rule %s | disk %.4f' % (nu, np.round(o.x, 3), G.sc(d, gp(th, d, nu))), ' '.join('%s %.4f' % (n, G.sc(a, gp(th, a, nu))) for n, a in ext.items()), flush=True)
