"""SA73 T56 exploration (all data, no protocol): size acting on halo concentration.
c = c_DM14(M_ref) x fac x 10^(-k delta), delta = log10 R_e - (p0 + p1 log10 M*), the size residual from the disks' own
size-mass relation (H.pop_coeff). A, B, C refitted on all disks at each k. Usage: python3 szc56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); cl = H.cluster_rows(); (p1, p0), _ = H.pop_coeff(gs)
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
     'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
d = P['disk']; dw = P['dwarfs']
print('disk size-mass relation: log Re = %.3f + %.3f log M*' % (p0, p1), '| dwarfs about it: mean delta %.3f sd %.3f' % (np.mean(np.log10(dw['re']) - p0 - p1 * np.log10(dw['ms'])), np.std(np.log10(dw['re']) - p0 - p1 * np.log10(dw['ms']))),
      '| disks sd %.3f' % np.std(np.log10(d['re']) - p0 - p1 * np.log10(d['ms'])))
def delta(a): return np.log10(a['re']) - p0 - p1 * np.log10(a['ms'])
def gpred(a, k, pops_k=True):
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * 10 ** (-k * delta(a))
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
def report(k):
    r = {p: round(hb(a, np.log10(a['go'] / gpred(a, k))), 4) for p, a in P.items()}
    e = np.log10(dw['go'] / gpred(dw, k)); X = np.c_[np.ones(len(e)), np.log10(dw['ms']), np.log10(dw['re'])]
    r['dw_mean'] = round(float(np.mean(e)), 3); r['dw_slope'] = round(float(np.linalg.lstsq(X, e, rcond=None)[0][2]), 3); return r
res = {}
for k in (-0.2, 0., .1, .2, .3, .4, .5, .7, 1.0):
    def f(p):
        B.FROZEN['abc'] = np.asarray(p); return np.log10(d['go'] / gpred(d, k)) * d['hw']
    fit = least_squares(f, ABC0, x_scale=[1, .3, .3]); B.FROZEN['abc'] = fit.x; r = report(k); r['abc'] = np.round(fit.x, 3).tolist(); res[k] = r
    print('k=%+.1f' % k, json.dumps(r), flush=True)
# free k on disks
def f4(p):
    B.FROZEN['abc'] = np.asarray(p[:3]); return np.log10(d['go'] / gpred(d, p[3])) * d['hw']
fit = least_squares(f4, list(ABC0) + [0.2], x_scale=[1, .3, .3, .2]); B.FROZEN['abc'] = fit.x[:3]; r = report(fit.x[3]); r['abc'] = np.round(fit.x[:3], 3).tolist()
print('disks choose k = %.3f:' % fit.x[3], json.dumps(r))
# what k the dwarfs alone prefer, with the final A, B, C
B.FROZEN['abc'] = ABC0
ks = np.arange(-.2, 1.21, .05); rm = [hb(dw, np.log10(dw['go'] / gpred(dw, k))) for k in ks]; print('dwarfs alone (final A, B, C) prefer k = %.2f, rms %.4f (k = 0: %.4f)' % (ks[int(np.argmin(rm))], min(rm), rm[list(np.round(ks, 2)).index(0.)]))
json.dump({str(k): v for k, v in res.items()}, open('res56_szc.json', 'w'), indent=1)
