"""All-data diagnostic (not held out): the paper's rule plus the published Moster (2013) halo mass, switched on at high stellar mass."""
import numpy as np, harness_T17 as X
H = X.H; B = X.B; C5 = X.C5; K = H.K
abc = np.array([11.7481, .8178, -.7119])
MC = {}
def mos(ms, z):
    out = np.empty(np.size(ms))
    for i, (m, zz) in enumerate(zip(np.ravel(ms), np.ravel(z))):
        k = (float(m), float(zz))
        if k not in MC: MC[k] = C5.moster_halo(m, zz)[0]
        out[i] = MC[k]
    return out.reshape(np.shape(ms))
def mass(ms, mp, re, z, lt, k):
    g = 1 / (1 + (10 ** lt / ms) ** k)
    return mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + g * mos(ms, z)
def gp(a, lt, k):
    m = mass(a['ms'], a['mp'], a['re'], a['z'], lt, k); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def es(st, lt, k):
    m = mass(st.ms, st.mp, st.re, st.z, lt, k); c = K.c_dm14(m, K.Z0)
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c; rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
from scipy.linalg import solve_triangular
def sc(a, g):
    er = np.log10(a['go'] / g); return np.sqrt(np.mean([np.mean(er[a['name'] == n] ** 2) for n in np.unique(a['name'])]))
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
lens = H.Lens(gs, False)
B.FROZEN['abc'] = abc; th = np.array([12.0629, 2.63, 2.5098])
print('SUBH all-data      disk %.4f' % sc(d, H.gpred('SUBH', th, d)), ' '.join('%s %.4f' % (n, sc(a, H.gpred('SUBH', th, a))) for n, a in ext.items()),
      'lens %.2f' % np.sum(lens.white('SUBH', th, [0, 1, 2, 3]) ** 2))
for k in [2, 4, 8]:
    for lt in [10.6, 10.9, 11.2, 11.5]:
        ch = sum(float(np.sum(solve_triangular(lens.L[j], lens.y[j] - es(lens.stacks[j], lt, k)[lens.masks[j]], lower=True) ** 2)) for j in range(4))
        print('gate k=%d logMt=%.1f disk %.4f' % (k, lt, sc(d, gp(d, lt, k))), ' '.join('%s %.4f' % (n, sc(a, gp(a, lt, k))) for n, a in ext.items()), 'lens %.2f' % ch)
