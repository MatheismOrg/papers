"""SA73 T54 exploration (all data, no protocol): can size act on the halo's concentration and remove the dwarfs' trend with
size, out of sample? c = c_DM14(M) x fac x (Sigma*/1e8 Msun kpc^-2)^s, Sigma* = M*/(2 pi R_e^2); A, B, C (and s) fitted on
the disks only; dwarfs, early types and clusters are then predictions. Usage: python3 explore54_sizec.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
def gpred(a, s, only=None):
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); sig = a['ms'] / (2 * np.pi * a['re'] ** 2) / 1e8
    c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * sig ** s
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
     'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
d = P['disk']; dw = P['dwarfs']
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
def report(lab, s):
    r = {p: round(hb(a, np.log10(a['go'] / gpred(a, s))), 4) for p, a in P.items()}
    e = np.log10(dw['go'] / gpred(dw, s)); X = np.c_[np.ones(len(e)), np.log10(dw['ms']), np.log10(dw['re'])]
    r['dwarf_mean'] = round(float(np.mean(e)), 3); r['dwarf_slope_logRe'] = round(float(np.linalg.lstsq(X, e, rcond=None)[0][2]), 3)
    r['abc'] = np.round(B.FROZEN['abc'], 4).tolist(); r['s'] = round(float(s), 4); print(lab, json.dumps(r), flush=True); return r
res = {}
B.FROZEN['abc'] = ABC0; res['final formula'] = report('final formula (s = 0)', 0.)
def f4(p):
    B.FROZEN['abc'] = np.asarray(p[:3]); return np.log10(d['go'] / gpred(d, p[3])) * d['hw']
fit = least_squares(f4, list(ABC0) + [0.], x_scale=[1, .3, .3, .1]); B.FROZEN['abc'] = fit.x[:3]; res['A,B,C,s on disks'] = report('A, B, C, s fitted on disks', fit.x[3])
def f3(p):
    B.FROZEN['abc'] = np.array([p[0], p[1], 0.]); return np.log10(d['go'] / gpred(d, p[2])) * d['hw']
fit = least_squares(f3, [11.7, .8, .1], x_scale=[1, .3, .1]); B.FROZEN['abc'] = np.array([fit.x[0], fit.x[1], 0.]); res['A,B,s on disks (C = 0)'] = report('A, B, s fitted on disks, C = 0', fit.x[2])
# what the dwarfs would want (dwarfs only, disks' A, B, C): scan s
B.FROZEN['abc'] = ABC0
for s in (0.05, 0.1, 0.15, 0.2, 0.3): res[f'scan s={s}'] = report(f'final A, B, C with s = {s}', s)
json.dump(res, open('res54_sizec.json', 'w'), indent=1)
