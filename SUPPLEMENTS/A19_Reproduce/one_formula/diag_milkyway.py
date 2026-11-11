"""Out-of-sample check (not a protocol run): the Milky Way's rotation curve predicted by the frozen T26 formula from the
Galaxy's published stars and gas, against the Gaia-based curve of Eilers et al. (2019, Table 1, 5 to 25 kpc).
Baryons: McMillan (2017) best-fitting model: bulge 9.23e9 Msun (point mass at these radii), thin disc Sigma0 896 Msun/pc^2,
Rd 2.50 kpc, z 0.3 kpc; thick disc 183, 3.02, 0.9; HI Sigma0 53.1, Rd 7, hole Rm 4, z 0.085; H2 2180, 1.5, Rm 12, z 0.045
(Sigma = Sigma0 exp(-Rm/R - R/Rd)). Disc forces by summing rings (potential of a ring offset by the scale height).
Halo: the formula's M_halo from M* = 5.43e10 (McMillan), M_p = M* + gas, R_e = 5.75 kpc (Lian et al. 2024); NFW with
Dutton and Maccio (2014) c, dark mass (1 - f_b) M_halo, as in every benchmark fit."""
import numpy as np
from scipy.special import ellipk
import harness_T26 as M
F = M.F; H = M.H; K = H.K; B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.1868])
G = 4.30091e-6                                   # kpc (km/s)^2 / Msun
def disc_v2(R, sig0, rd, rm, z0, rmax=60., n=6000):
    a = np.linspace(1e-3, rmax, n); da = a[1] - a[0]
    sig = sig0 * 1e6 * np.exp(-rm / a - a / rd)   # Msun/kpc^2
    m = 2 * np.pi * a * sig * da
    def phi(r):
        k2 = 4 * a * r / ((r + a) ** 2 + z0 ** 2)
        return np.sum(-2 * G * m / np.pi * ellipk(k2) / np.sqrt((r + a) ** 2 + z0 ** 2))
    h = 1e-3
    return np.array([r * (phi(r + h) - phi(r - h)) / (2 * h) for r in R])
d = np.loadtxt(H.S / 'ext/milkyway/eilers2019_table1.txt'); R, vobs, em, ep = d.T
discs = {'thin': (896., 2.50, 0., .3), 'thick': (183., 3.02, 0., .9), 'HI': (53.1, 7., 4., .085), 'H2': (2180., 1.5, 12., .045)}
v2b = G * 9.23e9 / R + sum(disc_v2(R, *p) for p in discs.values())
mstar = 5.43e10; a_ = np.linspace(1e-3, 60, 6000)
mgas = sum(2 * np.pi * np.sum(a_ * p[0] * 1e6 * np.exp(-p[2] / a_ - a_ / p[1])) * (a_[1] - a_[0]) for k, p in discs.items() if k in ('HI', 'H2'))
mp = mstar + mgas
def vdm2(m):
    c = K.c_dm14(m, 0.); rv = (H.G * m / (100 * H.hz(0.) ** 2)) ** (1 / 3)
    return G * (1 - H.FB) * m * H.ff(R * c / rv) / H.ff(c) / R, c, rv
mh = float(F.mass('MOGF0IE', th, np.array([mstar]), np.array([mp]), np.array([5.75]), np.array([0.]))[0])
v2d, c, rv = vdm2(mh); vmod = np.sqrt(v2b + v2d)
print('gas mass %.2e, M_p %.2e, formula M_halo %.3e (c %.2f, r200 %.0f kpc)' % (mgas, mp, mh, c, rv))
print('R      v_obs  v_model  v_bar  v_dm   ratio')
for r, vo, vm, vb, vd in zip(R, vobs, vmod, np.sqrt(v2b), np.sqrt(v2d)):
    if r in (5.27, 8.19, 10.26, 12.25, 15.22, 18.24, 20.27, 24.82): print('%5.2f %6.1f %7.1f %6.1f %6.1f %6.3f' % (r, vo, vm, vb, vd, vm / vo))
k = R < 21
print('mean model/observed, 5 to 21 kpc: %.3f (rms of log ratio %.3f dex in v, %.3f dex in g)' % (np.mean(vmod[k] / vobs[k]), np.sqrt(np.mean(np.log10(vmod[k] / vobs[k]) ** 2)), 2 * np.sqrt(np.mean(np.log10(vmod[k] / vobs[k]) ** 2))))
from scipy.optimize import brentq
best = brentq(lambda lm: np.mean(np.log10((v2b[k] + vdm2(10 ** lm)[0][k]) / vobs[k] ** 2)), 10.5, 13.5)
print('halo mass that matches the curve on average (same NFW and c(M)): %.3e' % 10 ** best)
