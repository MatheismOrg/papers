"""All-data diagnostic (not held out): a type rule from published relations only. Gas-poor galaxies above 10^10.3 Msun get
at least twice the rule's halo mass (Mandelbaum et al. 2016), and the Moster mass where that is larger:
M_ref = Mp/fb + (1 - q) R + q max(10^delta R, M_Moster), R = the paper's rule term, q = (M*/Mp)^50 / (1 + (10^10.3/M*)^4)."""
import numpy as np, sys
import harness_T19 as Z, kids_split as KS
H = Z.H; B = Z.B; K = H.K; Y = Z.Y
DELTA = {'v': np.log10(2.)}
H.MODELS['TYPE'] = (['dummy'], [-1.], [1.], [[0.]])
def mass_type(ms, mp, re, z):
    abc = B.FROZEN['abc']; R = 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
    q = (ms / mp) ** 50 / (1 + (10 ** 10.3 / ms) ** 4)
    return mp / H.FB + (1 - q) * R + q * np.maximum(10 ** DELTA['v'] * R, Y.moster(ms, z))
_g, _e = H.gpred, H.esd
def gpred(model, th, a):
    if model != 'TYPE': return _g(model, th, a)
    m = mass_type(a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != 'TYPE': return _e(model, th, st)
    m = mass_type(st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c; rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    import diag_mosg as G
    d, ext, lens = G.d, G.ext, G.lens; S = KS.setup(H.disks())
    B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([0.])
    for dv in [0., np.log10(1.5), np.log10(2.), np.log10(3.)]:
        DELTA['v'] = dv
        z = [float(np.sum(lens.white('TYPE', th, [j]) ** 2)) for j in range(4)]; r = KS.check('TYPE', th, S)
        print('TYPE x%.1f | disk %.4f' % (10 ** dv, G.sc(d, H.gpred('TYPE', th, d))), ' '.join('%s %.4f' % (n, G.sc(a, H.gpred('TYPE', th, a))) for n, a in ext.items()),
              '| lens %.1f' % sum(z), np.round(z, 1), '| split %.3f, colour diff chi2 %.1f red %.1f blue %.1f' % (r['S_model'], r['colour']['chi2_difference'], r['colour']['chi2_early'], r['colour']['chi2_late']), flush=True)
