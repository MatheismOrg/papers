"""SA73 T56 exploration (all data, no protocol): halo-mass scatter in the lensing stacks. A stack measures the mean excess
surface density of its lenses; if halo masses scatter log-normally about the formula's M_ref (the median) by sigma dex, the
stack sees the average over that spread. Final formula otherwise (T48: all-halo concentrations for the stacks).
Usage: python3 lscat56.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.linalg import solve_triangular
from scipy.stats import chi2
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; P44 = M.P; T = P44.T; Z = P44.Z; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
XN, XW = np.polynomial.hermite.hermgauss(9); XW = XW / np.sqrt(np.pi)
def esd(st, sig, mode='median'):
    m0 = S.mass(th, st.ms, st.mp, st.re, st.z, 0.) if hasattr(S, 'mass') else P44.G.mass(th, st.ms, st.mp, st.re, st.z, 0.)
    if mode == 'mean': m0 = m0 * 10 ** (-0.5 * np.log(10) * sig ** 2)          # formula's mass read as the mean, not the median
    tot = 0.
    T.F, T.dF = Z.Fc, Z.dFc
    try:
        for x, w in zip(XN, XW):
            m = m0 * 10 ** (sig * np.sqrt(2) * x); c = K.c_dm14(m, K.Z0) * P44.ratio_all_relaxed(m); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
            tot = tot + w * T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * tot, axis=2) + st.bary)
gs = H.disks(); out = {}
for mode in ('median', 'mean'):
    for sig in (0., .1, .15, .2, .3, .4, .5):
        row = {}
        for drop in (False, True):
            L = H.Lens(gs, drop); c2 = []; zs = []
            for j, st in enumerate(L.stacks):
                z = solve_triangular(L.L[j], L.y[j] - esd(st, sig, mode)[L.masks[j]], lower=True); c2.append(round(float(z @ z), 1)); zs.append(np.round(z, 1).tolist())
            row['strict' if drop else 'primary'] = (round(sum(c2), 1), c2)
            if not drop: row['z_outermost'] = [z[0] for z in zs]; row['z_innermost'] = [z[-1] for z in zs]
        out[f'{mode}|{sig}'] = row; print(mode, 'sigma', sig, 'primary', row['primary'], 'p %.4f' % chi2.sf(row['primary'][0], 29), '| strict', row['strict'][0], 'p %.3f' % chi2.sf(row['strict'][0], 25), '| z outermost', row['z_outermost'], 'innermost', row['z_innermost'], flush=True)
json.dump(out, open('res56_lscat.json', 'w'), indent=1)
