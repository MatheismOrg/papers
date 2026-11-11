"""SA73 T56 exploration (disks only, all data, no protocol): does a concentration that depends on stellar mass, gas share or
surface brightness remove the within-disk pattern? c = c_DM14(M_ref) x F; A, B, C refitted with each F's exponents.
Usage: python3 ctilt56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); d = H.arr(H.disk_rows(gs)); un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re'])
lms = np.log10(d['ms']); lgas = np.log10(d['ms'] / d['mp']); lsb = np.log10(d['ms'] / (2 * np.pi * d['re'] ** 2) / 1e8)
FORMS = {'none': lambda p: 0., 'mass tilt': lambda p: p[0] * (lms - 10), 'gas share': lambda p: p[0] * lgas, 'surface density': lambda p: p[0] * lsb,
         'offset': lambda p: p[0] + 0 * lms, 'mass + gas': lambda p: p[0] * (lms - 10) + p[1] * lgas, 'offset + mass + gas': lambda p: p[0] + p[1] * (lms - 10) + p[2] * lgas}
NP = {'none': 0, 'mass tilt': 1, 'gas share': 1, 'surface density': 1, 'offset': 1, 'mass + gas': 2, 'offset + mass + gas': 3}
def gpred(dc):
    m = G43.mass(th, d['ms'], d['mp'], d['re'], d['z'], S.thresholds(d)); c = K.c_dm14(m, d['z']) * 10 ** dc
    rv = (H.G * m / (100 * H.hz(d['z']) ** 2)) ** (1 / 3)
    return d['gb'] + H.G * (1 - H.FB) * m * H.ff(d['r'] * c / rv) / H.ff(c) / d['r'] ** 2 * H.ACC, m
gmass = np.array([d['ms'][gi == i][0] for i in range(len(un))])
for nm, F in FORMS.items():
    n = NP[nm]
    def f(p):
        B.FROZEN['abc'] = np.asarray(p[:3]); return np.log10(d['go'] / gpred(F(p[3:]))[0]) * d['hw']
    fit = least_squares(f, list(ABC0) + [0.] * n, x_scale=[1, .3, .3] + [.1] * n); B.FROZEN['abc'] = fit.x[:3]
    g, m = gpred(F(fit.x[3:])); e = np.log10(d['go'] / g); gm = np.bincount(gi, e) / ct; w = e - gm[gi]
    hbv = np.sqrt(np.mean(np.bincount(gi, e ** 2) / ct)); wi = np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    stack = [round(float(np.mean(w[(lx >= np.log10(lo)) & (lx < np.log10(hi))])), 3) for lo, hi in ((.01, .3), (.3, .6), (.6, 1), (1, 2.5), (2.5, 7), (7, 50))]
    mh = lambda lm: float(np.median(np.log10(m[np.abs(lms - lm) < .2]))) if np.any(np.abs(lms - lm) < .2) else np.nan
    print(f'{nm:20s} disk {hbv:.4f} within {wi:.4f} means {np.sqrt(np.mean(gm**2)):.4f} | slope med {np.nanmedian(sl):+.3f} corr logM* {np.corrcoef(np.log10(gmass)[ok], sl[ok])[0,1]:+.2f} | stack {stack} | p {np.round(fit.x, 3).tolist()} | med log M_h at 10.6: {mh(10.6):.2f}, 10.85: {mh(10.85):.2f}', flush=True)
B.FROZEN['abc'] = ABC0
