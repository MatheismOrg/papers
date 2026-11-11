"""All-data diagnostic: Cautun et al. (2020, eq. 11) contraction, M_DM(<r) = M_DM,NFW(<r) [0.45 + 0.38 (eta + 1.16)^0.53],
eta = M_bar(<r) / (f_b M_tot,NFW(<r)), M_bar from g_bar; dynamics only; on the T17b switch formula, disk rule refitted."""
import numpy as np, harness_T17b as Y, diag_mosg as G
from scipy.optimize import least_squares
H = Y.H; B = Y.B; K = H.K
def gp(th, a, strength=1.):
    m = Y.mass('MOSS', th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    mb = a['gb'] * a['r'] ** 2 / (H.G * H.ACC); eta = mb / (H.FB * mt)
    fac = 0.45 + 0.38 * (eta + 1.16) ** 0.53
    fac = fac ** strength
    return a['gb'] + H.G * (1 - H.FB) * mt * fac / a['r'] ** 2 * H.ACC, eta, fac
d, ext = G.d, G.ext; th = np.array([1.14548])
for strength in [0., 1.]:
    def fun(p):
        B.FROZEN['abc'] = p; return d['hw'] * np.log10(d['go'] / gp(th, d, strength)[0]) / .15
    o = least_squares(fun, [11.7481, .8178, -.7119], bounds=([9, -.5, -2.5], [14.5, 2, 1.5])); B.FROZEN['abc'] = o.x
    print('Cautun x%.0f rule %s | disk %.4f' % (strength, np.round(o.x, 3), G.sc(d, gp(th, d, strength)[0])), ' '.join('%s %.4f' % (n, G.sc(a, gp(th, a, strength)[0])) for n, a in ext.items()), flush=True)
    if strength:
        for n, a in [('disk', d)] + list(ext.items()):
            g, eta, fac = gp(th, a, strength); print('   %-6s eta median %.2f (10-90%%: %.2f-%.2f)  DM factor median %.2f' % (n, np.median(eta), *np.percentile(eta, [10, 90]), np.median(fac)))
