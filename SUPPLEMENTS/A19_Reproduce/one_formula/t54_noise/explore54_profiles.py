"""SA73 T54 exploration (all data, no protocol): do published feedback-made halo profiles remove the disks' within-galaxy
pattern (formula too high inside R_e in small, gas-rich disks)? Same halo masses and concentrations as the final formula,
A, B, C refitted on all disks for each profile:
  nfw      the final formula
  dc14     Di Cintio et al. (2014b): (alpha, beta, gamma) set by X = log10(M*/M_halo), -4.1 < X < -1.3 (NFW below, capped above),
           c = c_DM (1 + 3e-5 exp(3.4 (X + 4.5))), r_s from r_-2 = r_vir / c
  corenfw  Read, Agertz and Collins (2016): M(<r) = M_NFW(<r) tanh(r / (1.75 R_e))^n, n = tanh(0.04 t_SF / t_dyn), t_SF = 14 Gyr
Usage: python3 explore54_profiles.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119])
GLX, GLW = np.polynomial.legendre.leggauss(96)
def prof_int(x, al, be, ga):
    """int_0^x t^(2-ga) (1 + t^al)^(-(be-ga)/al) dt, vectorized (log-spaced Gauss-Legendre from 1e-7 x)"""
    lo = np.log(1e-7 * x); hi = np.log(x); lt = .5 * (hi - lo)[..., None] * (GLX + 1) + lo[..., None]; t = np.exp(lt)
    f = t ** (3 - ga[..., None]) * (1 + t ** al[..., None]) ** (-(be - ga)[..., None] / al[..., None])
    return .5 * (hi - lo) * (f @ GLW) + (1e-7 * x) ** (3 - ga) / (3 - ga)
def mdm(variant, m, c, rv, r, ms, re):
    nfw = m * H.ff(r * c / rv) / H.ff(c)
    if variant == 'nfw': return nfw
    if variant == 'dc14':
        X = np.log10(ms / m); low = X < -4.1; Xc = np.clip(X, -4.1, -1.3)
        al = 2.94 - np.log10(10 ** (-1.08 * (Xc + 2.33)) + 10 ** (2.29 * (Xc + 2.33))); be = 4.23 + 1.34 * Xc + 0.26 * Xc ** 2
        ga = -0.06 + np.log10(10 ** (-0.68 * (Xc + 2.56)) + 10 ** (Xc + 2.56))
        c2 = c * (1 + 3e-5 * np.exp(3.4 * (Xc + 4.5))); rs = rv / c2 / ((2 - ga) / (be - 2)) ** (1 / al)
        out = m * prof_int(r / rs, al, be, ga) / prof_int(rv / rs, al, be, ga)
        return np.where(low, nfw, out)
    if variant == 'corenfw':
        rs = rv / c; tdyn = 2 * np.pi * np.sqrt(rs ** 3 / (H.G * m * H.ff(1.) / H.ff(c))) * 0.9778          # Gyr (kpc/(km/s) = 0.9778 Gyr)
        n = np.tanh(0.04 * 14. / tdyn); return nfw * np.tanh(r / (1.75 * re)) ** n
def gpred(variant, a):
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * mdm(variant, m, c, rv, a['r'], a['ms'], a['re']) / a['r'] ** 2 * H.ACC
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
     'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
d = P['disk']; un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re'])
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
res = {}
for variant in ('nfw', 'dc14', 'corenfw'):
    def f(p):
        B.FROZEN['abc'] = np.asarray(p); return np.log10(d['go'] / gpred(variant, d)) * d['hw']
    fit = least_squares(f, ABC0, x_scale=[1, .3, .3]); B.FROZEN['abc'] = fit.x
    r = {'abc': np.round(fit.x, 4).tolist()}
    for p, a in P.items(): r[p] = round(hb(a, np.log10(a['go'] / gpred(variant, a))), 4)
    e = np.log10(d['go'] / gpred(variant, d)); gm = np.bincount(gi, e) / ct; w = e - gm[gi]
    r['disk_galaxy_means'] = round(float(np.sqrt(np.mean(gm ** 2))), 4); r['disk_within'] = round(float(np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))), 4)
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    lms = np.log10(np.array([d['ms'][gi == i][0] for i in range(len(un))]))
    r['disk_slope_median'] = round(float(np.nanmedian(sl)), 3); r['disk_slope_sd'] = round(float(np.nanstd(sl)), 3); r['disk_slope_vs_logMs'] = round(float(np.corrcoef(lms[ok], sl[ok])[0, 1]), 3)
    r['disk_stack'] = [round(float(np.mean(w[(lx >= np.log10(lo)) & (lx < np.log10(hi))])), 3) for lo, hi in ((.01, .3), (.3, .6), (.6, 1), (1, 1.5), (1.5, 2.5), (2.5, 4), (4, 7), (7, 50))]
    a = P['dwarfs']; e = np.log10(a['go'] / gpred(variant, a)); r['dwarf_mean'] = round(float(np.mean(e)), 3)
    Xm = np.c_[np.ones(len(e)), np.log10(a['ms']), np.log10(a['re'])]; r['dwarf_slope_logRe_at_fixed_Ms'] = round(float(np.linalg.lstsq(Xm, e, rcond=None)[0][2]), 3)
    res[variant] = r; print(variant, json.dumps(r), flush=True)
B.FROZEN['abc'] = ABC0
json.dump(res, open('res54_profiles.json', 'w'), indent=1)
