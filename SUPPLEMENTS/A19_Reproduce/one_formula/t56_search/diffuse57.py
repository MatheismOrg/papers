"""SA73 T57 look (all data, no protocol): bigger and looser halos only for diffuse, gas-rich galaxies, so that big spirals,
spheroidals and clusters keep the final formula. With g the gas share and L = -log10 min(Sigma*/Sigma_p, 1):
  M_ref -> M_ref x 10^(t g L),   c -> c x 10^(-s g L).   A, B, C, t, s fitted on all disks (or t, s alone with A, B, C frozen).
Usage: PYTHONPATH=<dmfit> python3 diffuse57.py"""
import sys, json, numpy as np
from scipy.optimize import least_squares
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; P44 = M.P; T = P44.T; Z = P44.Z; th = np.array([11.339]); ABC0 = np.array([11.7481, .8178, -.7119]); SP = 1e10 / (2 * np.pi * 9.)
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs)); LENS = {False: H.Lens(gs, False), True: H.Lens(gs, True)}
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows())}
un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re']); gmass = np.array([d['ms'][gi == i][0] for i in range(len(un))])
def gl(ms, mp, re): return np.clip(1 - ms / mp, 0, 1) * -np.log10(np.minimum(ms / (2 * np.pi * re ** 2) / SP, 1.))
def gpred(a, t, s):
    w = gl(a['ms'], a['mp'], a['re']) if a['kind'][0] == 'disk' else 0.
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)) * 10 ** (t * w); c = K.c_dm14(m, a['z']) * a.get('fac', 1.) * 10 ** (-s * w)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(st, t, s):
    w = gl(st.ms, st.mp, st.re); m = G43.mass(th, st.ms, st.mp, st.re, st.z, 0.) * 10 ** (t * w); c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m) * 10 ** (-s * w)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3); T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def hb(a, e):
    u, i = np.unique(a['name'], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(i, e ** 2) / np.bincount(i))))
from scipy.linalg import solve_triangular
def report(lab, t, s):
    r = dict(abc=np.round(B.FROZEN['abc'], 3).tolist(), t=round(float(t), 3), s=round(float(s), 3), disk=round(hb(d, np.log10(d['go'] / gpred(d, t, s))), 4), **{k: round(hb(a, np.log10(a['go'] / gpred(a, t, s))), 4) for k, a in ext.items()})
    for drop in (False, True):
        L = LENS[drop]; r['lens_' + ('strict' if drop else 'primary')] = round(float(sum(np.sum(solve_triangular(L.L[j], L.y[j] - esd(st, t, s)[L.masks[j]], lower=True) ** 2) for j, st in enumerate(L.stacks))), 1)
    e = np.log10(d['go'] / gpred(d, t, s)); gm = np.bincount(gi, e) / ct; w = e - gm[gi]; r['within_disk'] = round(float(np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))), 4)
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    r['slope_vs_mass'] = round(float(np.corrcoef(np.log10(gmass)[ok], sl[ok])[0, 1]), 2); r['inner_stack'] = round(float(np.mean(w[lx < np.log10(.6)])), 3); r['disk_means'] = round(float(np.sqrt(np.mean(gm ** 2))), 4)
    print(f'{lab:34s}', json.dumps(r), flush=True)
B.FROZEN['abc'] = ABC0; report('final formula', 0., 0.)
def f2(p): B.FROZEN['abc'] = ABC0; return np.log10(d['go'] / gpred(d, p[0], p[1])) * d['hw']
o = least_squares(f2, [.1, .1]); B.FROZEN['abc'] = ABC0; report('t, s only (A, B, C frozen)', *o.x)
def f1(p): B.FROZEN['abc'] = ABC0; return np.log10(d['go'] / gpred(d, 0., p[0])) * d['hw']
o = least_squares(f1, [.1]); B.FROZEN['abc'] = ABC0; report('s only (A, B, C frozen)', 0., o.x[0])
def f5(p): B.FROZEN['abc'] = np.asarray(p[:3]); return np.log10(d['go'] / gpred(d, p[3], p[4])) * d['hw']
o = least_squares(f5, list(ABC0) + [.1, .1]); B.FROZEN['abc'] = o.x[:3]; report('A, B, C, t, s on disks', o.x[3], o.x[4])
B.FROZEN['abc'] = ABC0
