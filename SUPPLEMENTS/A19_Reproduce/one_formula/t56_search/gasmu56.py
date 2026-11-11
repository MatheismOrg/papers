"""SA73 T56 exploration (all data, no protocol): can the cold gas carry the bigger, looser halos of diffuse disks, so that the
stellar exponent B and with it the dwarf spheroidals stay put? M_ref = (M* + mu M_gas)/f_b + 10^A (M*/1e10)^B (R_e/3)^C,
c = c_DM14(M_ref) x min(x, 1)^(s w), w = 1 or the gas share. Fitted on disks only. Usage: python3 gasmu56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); cl = H.cluster_rows(); SP = 1e10 / (2 * np.pi * 9.)
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
     'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
d = P['disk']; dw = P['dwarfs']
def gpred(a, mu, s, mode):
    isd = a['kind'][0] == 'disk'; gas = (a['mp'] - a['ms']) if isd else 0.; g = np.clip(gas / a['mp'], 0, 1) if isd else np.zeros(len(a['ms']))
    mp = a['ms'] + mu * gas if isd else a['mp']
    m = G43.mass(th, a['ms'], mp, a['re'], a['z'], S.thresholds(a)); x = a['ms'] / (2 * np.pi * a['re'] ** 2) / SP
    w = {'none': 1., 'gas share': g, 'off': 0.}[mode]
    c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * np.minimum(x, 1.) ** (s * w); rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
for lab, mode, free_mu, free_s in (('final formula', 'off', False, False), ('mu only', 'off', True, False), ('mu + looseness (no gate)', 'none', True, True), ('mu + looseness (gas share)', 'gas share', True, True),
                                   ('mu + looseness (gas share), B fixed', 'gas share', True, True)):
    fixB = 'B fixed' in lab
    def unpack(p):
        abc = np.array([p[0], ABC0[1], p[1]]) if fixB else np.asarray(p[:3]); r = list(p[2 if fixB else 3:]); mu = abs(r.pop(0)) if free_mu else 1.; s = r.pop(0) if free_s else 0.; return abc, mu, s
    def f(p):
        abc, mu, s = unpack(p); B.FROZEN['abc'] = abc; return np.log10(d['go'] / gpred(d, mu, s, mode)) * d['hw']
    p0 = ([ABC0[0], ABC0[2]] if fixB else list(ABC0)) + ([1.5] if free_mu else []) + ([.1] if free_s else [])
    if len(p0) > 3 or fixB or free_mu: fit = least_squares(f, p0); abc, mu, s = unpack(fit.x)
    else: abc, mu, s = ABC0, 1., 0.
    B.FROZEN['abc'] = abc; r = {p: round(hb(a, np.log10(a['go'] / gpred(a, mu, s, mode))), 4) for p, a in P.items()}
    r['dw_mean'] = round(float(np.mean(np.log10(dw['go'] / gpred(dw, mu, s, mode)))), 3); r['abc'] = np.round(abc, 3).tolist(); r['mu'] = round(float(mu), 3); r['s'] = round(float(s), 3)
    print(f'{lab:38s}', json.dumps(r), flush=True)
B.FROZEN['abc'] = ABC0
