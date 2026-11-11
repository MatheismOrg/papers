"""All-data diagnostic (not held out): scatter of halo mass at fixed galaxy properties in the stacked lensing.
A stack measures the mean excess surface density of its lenses. If halo masses scatter log-normally around the rule's
M_ref with width sigma (dex), the stack sees the average over that spread, which is not the profile of the median halo.
Rotation curves and dispersions are single systems and keep the median. MOGF (T21) at the all-data coefficients.
Usage: python3 diag_scatter.py <sigma> [refit]"""
import sys, json
import numpy as np
from scipy.linalg import solve_triangular
import harness_T21 as F
Z = F.Z; H = Z.H; B = Z.B; K = H.K; T = Z.T
XN, XW = np.polynomial.hermite.hermgauss(11); XW = XW / np.sqrt(np.pi)
SIG = {'v': 0.}
def esd_sc(model, th, st):
    base = F.FLOOR[model][0]; m0 = F.mass(model, th, st.ms, st.mp, st.re, st.z); e = 0.
    nodes = [(0., 1.)] if SIG['v'] == 0 else zip(np.sqrt(2) * SIG['v'] * XN, XW)
    for dx, w in nodes:
        m = m0 * 10 ** dx; c = K.c_dm14(m, K.Z0); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
        T.F, T.dF = Z.Fc, Z.dFc
        try: e = e + w * T.esd_dm(m, c, rv, st.mp, st.R, Z.strength(base, th))
        finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
_esd = H.esd
def esd(model, th, st): return esd_sc(model, th, st) if model in F.FLOOR else _esd(model, th, st)
H.esd = esd
def score(model, th, a):
    l = H.losses(model, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
if __name__ == '__main__':
    SIG['v'] = float(sys.argv[1]); refit = len(sys.argv) > 2
    gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); out = {}
    for drop in (False, True):
        lens = H.Lens(gs, drop); th = np.array([11.182185, .173541])
        if refit: th = np.array(H.fit('MOGF', .1, d, ext, lens, [0, 1, 2, 3])['th'])
        z = [float(np.sum(lens.white('MOGF', th, [j]) ** 2)) for j in range(4)]
        ratio = [np.round(lens.y[j] / H.esd('MOGF', th, lens.stacks[j])[lens.masks[j]], 2).tolist() for j in range(4)]
        out['strict' if drop else 'primary'] = dict(sigma=SIG['v'], th=np.round(th, 4).tolist(), lens=round(sum(z), 2), by_group=np.round(z, 2).tolist(),
              disk=score('MOGF', th, d), **{k: score('MOGF', th, a) for k, a in ext.items()}, data_over_model=ratio)
        print(json.dumps(out['strict' if drop else 'primary']), flush=True)
