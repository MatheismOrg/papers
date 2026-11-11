"""SA73 T56 exploration (disks only, all data, no protocol): cores where cold gas remains.
M_DM(<r) = M_NFW(<r) tanh(r / (eta R_e))^n (Read, Agertz and Collins 2016, eta = 1.75), with the core strength n set by the
gas share of the baryons, g = 1 - M*/M_p, in place of their star-formation time. A, B, C refitted on all disks.
Usage: python3 core56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
gs = H.disks(); d = H.arr(H.disk_rows(gs)); un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re'])
g = 1 - d['ms'] / d['mp']; lms = np.log10(d['ms']); x = d['ms'] / (2 * np.pi * d['re'] ** 2) / (1e10 / (2 * np.pi * 9.))
print('gas share of the disks: percentiles 5, 25, 50, 75, 95:', np.round(np.percentile(g[np.unique(gi, return_index=True)[1]], [5, 25, 50, 75, 95]), 2))
def gpred(n, eta):
    m = G43.mass(th, d['ms'], d['mp'], d['re'], d['z'], S.thresholds(d)); c = K.c_dm14(m, d['z']); rv = (H.G * m / (100 * H.hz(d['z']) ** 2)) ** (1 / 3)
    return d['gb'] + H.G * (1 - H.FB) * m * H.ff(d['r'] * c / rv) / H.ff(c) * np.tanh(d['r'] / (eta * d['re'])) ** n / d['r'] ** 2 * H.ACC, m, c, rv
def tdyn_n():
    _, m, c, rv = gpred(0., 1.75); rs = rv / c; t = 2 * np.pi * np.sqrt(rs ** 3 / (H.G * m * H.ff(1.) / H.ff(c))) * 0.9778; return np.tanh(0.04 * 14. / t)
FORMS = {'NFW': (0, lambda p: (0. * g, 1.75)), 'n = g': (0, lambda p: (g, 1.75)), 'n = g^2': (0, lambda p: (g ** 2, 1.75)), 'n = sqrt g': (0, lambda p: (np.sqrt(g), 1.75)),
         'n = 1 (all)': (0, lambda p: (1. + 0 * g, 1.75)), 'n = Read(t_dyn)': (0, lambda p: (tdyn_n(), 1.75)), 'n = Read x g': (0, lambda p: (tdyn_n() * g, 1.75)),
         'n = a g': (1, lambda p: (p[0] * g, 1.75)), 'n = g, eta free': (1, lambda p: (g, abs(p[0]) * 17.5 + 1e-3)), 'n = a g, eta free': (2, lambda p: (p[0] * g, abs(p[1]) * 17.5 + 1e-3)),
         'n = a g min(x,1)^-b': (2, lambda p: (p[0] * g * np.minimum(x, 1.) ** (-p[1]), 1.75))}
gmass = np.array([d['ms'][gi == i][0] for i in range(len(un))]); out = {}
for nm, (k, F) in FORMS.items():
    def f(p):
        B.FROZEN['abc'] = np.asarray(p[:3]); n, eta = F(p[3:]); return np.log10(d['go'] / gpred(n, eta)[0]) * d['hw']
    fit = least_squares(f, list(ABC0) + [1., .1][:k] if k else list(ABC0)); B.FROZEN['abc'] = fit.x[:3]; n, eta = F(fit.x[3:])
    gp, m, c, rv = gpred(n, eta); e = np.log10(d['go'] / gp); gm = np.bincount(gi, e) / ct; w = e - gm[gi]
    hbv = np.sqrt(np.mean(np.bincount(gi, e ** 2) / ct)); wi = np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    stack = [round(float(np.mean(w[(lx >= np.log10(lo)) & (lx < np.log10(hi))])), 3) for lo, hi in ((.01, .3), (.3, .6), (.6, 1), (1, 2.5), (2.5, 7), (7, 50))]
    extra = np.round(fit.x[3:], 3).tolist() if k else []
    if 'eta free' in nm: extra[-1] = round(float(eta), 3)
    out[nm] = dict(disk=round(float(hbv), 4), within=round(float(wi), 4), abc=np.round(fit.x[:3], 3).tolist(), extra=extra)
    print(f'{nm:20s} disk {hbv:.4f} within {wi:.4f} means {np.sqrt(np.mean(gm**2)):.4f} | slope med {np.nanmedian(sl):+.3f} corr logM* {np.corrcoef(np.log10(gmass)[ok], sl[ok])[0,1]:+.2f} | stack {stack} | A,B,C {np.round(fit.x[:3], 3).tolist()} {extra}', flush=True)
B.FROZEN['abc'] = ABC0; json.dump(out, open('res56_core.json', 'w'), indent=1)
