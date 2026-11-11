"""SA73 T56 exploration (disks only, all data, no protocol): concentration read from stellar surface density.
c = c_DM14(M_ref) x F(Sigma*), Sigma* = M*/(2 pi R_e^2), Sigma_p = 1e10/(2 pi 3^2) Msun/kpc^2 (the formula's own pivot).
Usage: python3 sig56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); d = H.arr(H.disk_rows(gs)); un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re'])
SP = 1e10 / (2 * np.pi * 9.); x = d['ms'] / (2 * np.pi * d['re'] ** 2) / SP; lms = np.log10(d['ms'])
print('log Sigma*/Sigma_p of the disks: percentiles 5, 25, 50, 75, 95:', np.round(np.percentile(np.log10(x[np.unique(gi, return_index=True)[1]]), [5, 25, 50, 75, 95]), 2))
FORMS = {'none': (0, lambda p: 1.), 'power (x^s)': (1, lambda p: x ** p[0]), 'boost (1+x)^s': (1, lambda p: (1 + x) ** p[0]), 'max(x,1)^s': (1, lambda p: np.maximum(x, 1.) ** p[0]),
         'min(x,1)^s': (1, lambda p: np.minimum(x, 1.) ** p[0]), 'power, B fixed': (1, lambda p: x ** p[0]), 'power, C = 0': (1, lambda p: x ** p[0]), 'two slopes': (2, lambda p: np.minimum(x, 1.) ** p[0] * np.maximum(x, 1.) ** p[1])}
def gpred(F):
    m = G43.mass(th, d['ms'], d['mp'], d['re'], d['z'], S.thresholds(d)); c = K.c_dm14(m, d['z']) * F
    rv = (H.G * m / (100 * H.hz(d['z']) ** 2)) ** (1 / 3)
    return d['gb'] + H.G * (1 - H.FB) * m * H.ff(d['r'] * c / rv) / H.ff(c) / d['r'] ** 2 * H.ACC, m
gmass = np.array([d['ms'][gi == i][0] for i in range(len(un))]); out = {}
for nm, (n, F) in FORMS.items():
    def unpack(p):
        if nm == 'power, B fixed': return np.array([p[0], ABC0[1], p[1]]), p[2:]
        if nm == 'power, C = 0': return np.array([p[0], p[1], 0.]), p[2:]
        return np.asarray(p[:3]), p[3:]
    def f(p):
        abc, q = unpack(p); B.FROZEN['abc'] = abc; return np.log10(d['go'] / gpred(F(q))[0]) * d['hw']
    p0 = ([ABC0[0], ABC0[2]] if nm == 'power, B fixed' else [ABC0[0], ABC0[1]] if nm == 'power, C = 0' else list(ABC0)) + [0.1] * n
    fit = least_squares(f, p0); abc, q = unpack(fit.x); B.FROZEN['abc'] = abc
    g, m = gpred(F(q)); e = np.log10(d['go'] / g); gm = np.bincount(gi, e) / ct; w = e - gm[gi]
    hbv = np.sqrt(np.mean(np.bincount(gi, e ** 2) / ct)); wi = np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    mh = lambda lm: float(np.median(np.log10(m[np.abs(lms - lm) < .25])))
    out[nm] = dict(disk=round(float(hbv), 4), within=round(float(wi), 4), abc=np.round(abc, 3).tolist(), s=np.round(q, 3).tolist())
    print(f'{nm:16s} disk {hbv:.4f} within {wi:.4f} means {np.sqrt(np.mean(gm**2)):.4f} | slope corr logM* {np.corrcoef(np.log10(gmass)[ok], sl[ok])[0,1]:+.2f} | A,B,C {np.round(abc, 3).tolist()} s {np.round(q, 3).tolist()} | log M_h at 8: {mh(8.):.2f} 9: {mh(9.):.2f} 10: {mh(10.):.2f} 10.6: {mh(10.6):.2f} 11: {mh(11.):.2f}', flush=True)
B.FROZEN['abc'] = ABC0; json.dump(out, open('res56_sig.json', 'w'), indent=1)
