"""SA73 T56 exploration (all data, no protocol): looser halos for diffuse disks, with and without a cold-gas gate, and what
each does to the other scales through A, B, C. c = c_DM14(M_ref) x fac x min(x, 1)^(s w), x = Sigma*/Sigma_p,
w = 1 (no gate), g (gas share of baryons) or [g > 0]. Fits: disks only, or disks plus dwarfs at weight 0.1.
Usage: python3 gate56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); cl = H.cluster_rows(); SP = 1e10 / (2 * np.pi * 9.)
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
     'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
d = P['disk']; dw = P['dwarfs']
def gate(a, mode):
    g = np.clip(1 - a['ms'] / a['mp'], 0, 1) if a['kind'][0] == 'disk' else np.zeros(len(a['ms']))   # only the disks carry cold gas in the benchmark rows
    return {'none': np.ones(len(g)), 'gas share': g, 'has gas': (g > 0).astype(float)}[mode]
def gpred(a, s, mode):
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); x = a['ms'] / (2 * np.pi * a['re'] ** 2) / SP
    c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * np.minimum(x, 1.) ** (s * gate(a, mode)); rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
def report(lab, s, mode):
    r = {p: round(hb(a, np.log10(a['go'] / gpred(a, s, mode))), 4) for p, a in P.items()}
    e = np.log10(dw['go'] / gpred(dw, s, mode)); r['dw_mean'] = round(float(np.mean(e)), 3)
    r['abc'] = np.round(B.FROZEN['abc'], 3).tolist(); r['s'] = round(float(s), 3); print(f'{lab:44s}', json.dumps(r), flush=True)
B.FROZEN['abc'] = ABC0; report('final formula', 0., 'none')
for mode in ('none', 'has gas', 'gas share'):
    for wd in (0., .1, 1.):
        def f(p):
            B.FROZEN['abc'] = np.asarray(p[:3]); r1 = np.log10(d['go'] / gpred(d, p[3], mode)) * d['hw'] / .15
            if wd == 0: return r1
            return np.concatenate([r1, np.sqrt(wd) * np.log10(dw['go'] / gpred(dw, p[3], mode)) * dw['hw'] / .4])
        fit = least_squares(f, list(ABC0) + [0.1]); B.FROZEN['abc'] = fit.x[:3]
        report(f'gate: {mode}; fit on disks' + (f' + dwarfs at {wd}' if wd else ''), fit.x[3], mode)
# reference: the same joint fits without the new term (s = 0), to see what joint fitting alone does
for wd in (.1, 1.):
    def f(p):
        B.FROZEN['abc'] = np.asarray(p); return np.concatenate([np.log10(d['go'] / gpred(d, 0., 'none')) * d['hw'] / .15, np.sqrt(wd) * np.log10(dw['go'] / gpred(dw, 0., 'none')) * dw['hw'] / .4])
    fit = least_squares(f, list(ABC0)); B.FROZEN['abc'] = fit.x; report(f's = 0; fit on disks + dwarfs at {wd}', 0., 'none')
B.FROZEN['abc'] = ABC0
