"""SA73 T53 exploration (not pre-registered): how the MaNGA BCG residual trend with size responds to the split between
stars and halo. g_pred = f_s g_bar + f_h g_dm, the formula's halo term unchanged otherwise. Usage: python3 -I f5_split.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
a = H.arr([r for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG'])
gp = H.gpred(M.MODEL, th, a); gdm = gp - a['gb']; go = a['go']; lre = np.log10(a['re'])
fdm = gdm / gp
print('DM share of predicted g at r: median', round(float(np.median(fdm)), 3), 'range', np.round(np.percentile(fdm, [5, 95]), 3).tolist(),
      'corr with log R_e', round(float(np.corrcoef(lre, fdm)[0, 1]), 3))
e0 = np.log10(go / gp); print('benchmark rms', round(float(np.sqrt(np.mean(e0 ** 2))), 4), 'mean', round(float(np.mean(e0)), 3), 'corr', round(float(np.corrcoef(lre, e0)[0, 1]), 3))
res = {}
for ls in np.arange(-.40, .101, .05):
    for lh in np.arange(-.2, .61, .1):
        e = np.log10(go / (10 ** ls * a['gb'] + 10 ** lh * gdm))
        res[(round(ls, 2), round(lh, 2))] = (float(np.sqrt(np.mean(e ** 2))), float(np.mean(e)), float(np.corrcoef(lre, e)[0, 1]))
best = min(res, key=lambda k: res[k][0])
print('best (log f_s, log f_h)', best, 'rms %.4f mean %.3f corr %.3f' % res[best])
for ls in (0., -.1, -.2, -.3):
    row = [(lh, round(res[(round(ls, 2), round(lh, 2))][0], 3), round(res[(round(ls, 2), round(lh, 2))][2], 2)) for lh in np.round(np.arange(-.2, .61, .1), 2)]
    print('log f_s', ls, ' (log f_h, rms, corr):', row)
# stellar weight trend: f_s depending on the stellar surface density proxy log(M*/R_e^2)
lsd = np.log10(a['ms'] / a['re'] ** 2); print('corr(log M*/R_e^2, e0)', round(float(np.corrcoef(lsd, e0)[0, 1]), 3), ' corr(log g_bar, e0)', round(float(np.corrcoef(np.log10(a['gb']), e0)[0, 1]), 3))
