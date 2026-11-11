"""All-data diagnostic (not held out): the baryons missing from each galaxy, fb M_ref - Mp, placed in a hot halo with
density proportional to r^-1.5 inside r200 (a fraction f of them), on top of T19's MOGC. CLASH points are left out of the
hot term because their g_bar already includes the intracluster gas. Dynamics: g_hot = G M_hot(<r)/r^2; lensing: projected
numerically. The paper's rule is kept as fitted on disks without the hot halo."""
import numpy as np, sys
import harness_T19 as Z
H = Z.H; B = Z.B; K = H.K; T = Z.T
HOT = {'f': 0.5, 'p': 1.5}
def mhot_enc(r, m, mp, rv):
    tot = HOT['f'] * np.maximum(H.FB * m - mp, 0.)
    return tot * np.minimum(r / rv, 1.) ** (3 - HOT['p'])
def gp(th, a, skip_clash=True):
    g = H.gpred('MOGC', th, a)
    m = Z.mass('MOGC', th, a['ms'], a['mp'], a['re'], a['z']); rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    add = H.G * mhot_enc(a['r'], m, a['mp'], rv) / a['r'] ** 2 * H.ACC
    if skip_clash: add = np.where(np.char.startswith(a['kind'].astype(str), 'CLASH'), 0., add)
    return g + add
TN, TW = np.polynomial.legendre.leggauss(48); TMAX = 9.; TT = .5 * TMAX * (TN + 1); WT = .5 * TMAX * TW; CH, SH, TH = np.cosh(TT), np.sinh(TT), np.tanh(TT)
def es(th, st):
    base = H.esd('MOGC', th, st)
    m = Z.mass('MOGC', th, st.ms, st.mp, st.re, st.z); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    tot = HOT['f'] * np.maximum(H.FB * m - st.mp, 0.); q = 3 - HOT['p']
    R = st.R; R4 = R[..., None]; rv4 = rv[:, None, None, None]; tot4 = tot[:, None, None, None]
    Menc = lambda r: tot4 * np.minimum(r / rv4, 1.) ** q
    dM = lambda r: np.where(r < rv4, tot4 * q * r ** (q - 1) / rv4 ** q, 0.)
    r_ = R4 * CH; d = dM(r_)
    sigma = (d / CH) @ WT / (2 * np.pi * R)
    mcyl = Menc(R4)[..., 0] + (d * (1 - TH) * R4 * SH) @ WT
    e = (mcyl / (np.pi * R ** 2) - sigma) / 1e12
    return base + st.w @ np.sum(st.aw * e, axis=2)
from scipy.linalg import solve_triangular
if __name__ == '__main__':
    import diag_mosg as G, kids_split as KS
    d, ext, lens = G.d, G.ext, G.lens
    B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.182185, .173541])
    for f in [0., .5, 1.]:
        HOT['f'] = f
        z = [float(np.sum(solve_triangular(lens.L[j], lens.y[j] - es(th, lens.stacks[j])[lens.masks[j]], lower=True) ** 2)) for j in range(4)]
        print('hot f=%.1f | disk %.4f' % (f, G.sc(d, gp(th, d))), ' '.join('%s %.4f' % (n, G.sc(a, gp(th, a))) for n, a in ext.items()), '| lens %.1f' % sum(z), np.round(z, 1), flush=True)
