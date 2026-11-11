"""All-data diagnostic (not held out): Blumenthal et al. (1986) adiabatic contraction of the NFW halo by the galaxy's
condensed baryons, for dynamics, on top of the T17b switch formula (MOSS). Baryons inside r from g_bar, capped at Mp."""
import numpy as np, harness_T17b as Y
H = Y.H; B = Y.B; K = H.K
def g_dm_ac(r, gb, mp, m, c, z):
    """DM acceleration at r after Blumenthal contraction; m = M_ref, dark mass (1-fb) m inside r200."""
    rv = (H.G * m / (100 * H.hz(z) ** 2)) ** (1 / 3); md = (1 - H.FB) * m
    Mdm = lambda x: md * H.ff(x * c / rv) / H.ff(c)
    mb = np.minimum(gb * r ** 2 / (H.G * H.ACC), mp); eps = mp / md
    lo, hi = np.log(r) - 4, np.log(r) + 4
    for _ in range(80):
        mid = (lo + hi) / 2; ri = np.exp(mid)
        f = ri * (1 + eps) * Mdm(ri) - r * (mb + Mdm(ri))
        lo = np.where(f < 0, mid, lo); hi = np.where(f < 0, hi, mid)
    ri = np.exp((lo + hi) / 2)
    return H.G * Mdm(ri) / r ** 2 * H.ACC, ri / r
def gpred_ac(th, a):
    m = Y.mass('MOSS', th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    g, ratio = g_dm_ac(a['r'], a['gb'], a['mp'], m, c, a['z'])
    return a['gb'] + g
if __name__ == '__main__':
    import diag_mosg as G
    from scipy.optimize import least_squares
    d, ext = G.d, G.ext
    def sc(a, g): return G.sc(a, g)
    B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([1.14548])
    print('MOSS no AC   | disk %.4f' % sc(d, H.gpred('MOSS', th, d)), ' '.join('%s %.4f' % (n, sc(a, H.gpred('MOSS', th, a))) for n, a in ext.items()))
    print('MOSS with AC | disk %.4f' % sc(d, gpred_ac(th, d)), ' '.join('%s %.4f' % (n, sc(a, gpred_ac(th, a))) for n, a in ext.items()))
    # refit the rule on disks with AC (no branch: switch off for all disks anyway)
    def fun(p):
        B.FROZEN['abc'] = p; return d['hw'] * np.log10(d['go'] / gpred_ac(th, d)) / .15
    o = least_squares(fun, [11.7481, .8178, -.7119], bounds=([9, -.5, -2.5], [14.5, 2, 1.5])); B.FROZEN['abc'] = o.x
    print('AC disk-only rule', np.round(o.x, 4), '| disk %.4f' % sc(d, gpred_ac(th, d)), ' '.join('%s %.4f' % (n, sc(a, gpred_ac(th, a))) for n, a in ext.items()))
    for n, a in [('disk', d)] + list(ext.items()):
        er = np.log10(a['go'] / gpred_ac(th, a)); print('   mean residual with AC and refit rule: %-6s %+.3f' % (n, er.mean()))
