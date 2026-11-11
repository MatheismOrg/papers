"""SA73 T59: baryonic model for the Iorio et al. (2017) dwarf irregulars, on the scale of the SPARC benchmark.
Gas: razor-thin disk from the published intrinsic HI surface density, times 1.33 for helium (as SPARC), by a ring sum.
Stars: razor-thin exponential disk (Freeman 1970) of mass 0.5 L[3.6] (as SPARC), L[3.6] from the LVL flux of Dale et al.
(2009) at the distance of Iorio et al.; scale length R_* of Read et al. (2017).  Options are fixed by validate59.py on the
galaxies shared with SPARC, before any outside galaxy is scored."""
import numpy as np
from scipy.special import ellipk, ellipe, i0, i1, k0, k1
G = 4.30091e-6; ACC = 3.24078e-14

def sigma_profile(R, Sd, mode):
    """Returns f(a) in Msun/kpc^2 (HI only) and the outer radius of the model. Inside the first point: flat. Between
    points: linear. Beyond the last point: 'cut' = zero from R_last + half a step; 'exp' = exponential with the slope of
    a straight-line fit to ln(Sd) over the outer half of the points, to 3 R_last."""
    R = np.asarray(R, float); S = np.maximum(np.asarray(Sd, float), 0.) * 1e6; step = R[-1] - R[-2]
    if mode == 'cut':
        rout = R[-1] + step / 2
        return (lambda a: np.where(a <= rout, np.interp(a, R, S), 0.)), rout
    k = slice(len(R) // 2, None); ok = S[k] > 0
    b = np.polyfit(R[k][ok], np.log(S[k][ok]), 1)[0] if ok.sum() >= 2 else 0.
    b = min(b, -1. / (3 * R[-1]))                        # at least as steep as a scale length of 3 R_last
    rout = 3 * R[-1]
    return (lambda a: np.where(a <= R[-1], np.interp(a, R, S), S[-1] * np.exp(b * (a - R[-1])))), rout

def ring_v2(f, rout, Rq, n=50000):
    """V^2 at the radii Rq of a razor-thin disk with surface density f(a) (Msun/kpc^2) out to rout. For each radius the
    rings sit at (j + 1/2) d with d = R / ceil(R / (rout / n)), so that R is a node and the ring sum is a principal value."""
    out = []
    for R in np.atleast_1d(np.asarray(Rq, float)):
        d = R / np.ceil(R / (rout / n)); a = (np.arange(int(np.ceil(rout / d))) + .5) * d; m = 2 * np.pi * a * d * f(a); k2 = 4 * a * R / (a + R) ** 2
        out.append(G / np.pi * ((ellipk(k2) / (a + R) - ellipe(k2) / (a - R)) * m).sum())
    return np.array(out)

def hi_mass(f, rout, n=20000):
    a = (np.arange(n) + .5) * rout / n; return float((2 * np.pi * a * rout / n * f(a)).sum())

def gas(R, Sd, mode, he=1.33):
    f, rout = sigma_profile(R, Sd, mode); f0, r0 = sigma_profile(R, Sd, 'cut')
    return he * ring_v2(f, rout, R), hi_mass(f, rout), hi_mass(f0, r0)      # V_gas^2 with helium; HI mass of the model; HI mass to the last point

def stars_v2(R, M, Rd):
    y = np.asarray(R, float) / (2 * Rd); return 2 * G * M / Rd * y ** 2 * (i0(y) * k0(y) - i1(y) * k1(y))

def l36(f_jy, D_mpc):
    return 10 ** (-0.4 * (-2.5 * np.log10(f_jy / 280.9) - 5 * np.log10(D_mpc * 1e5) - 3.24))

if __name__ == '__main__':                                # self test against Freeman's exponential disk
    Rd, M = 1.3, 2e8; R = np.linspace(0.05, 12, 60); S = M / (2 * np.pi * Rd ** 2) * np.exp(-R / Rd) / 1e6
    f = lambda a: M / (2 * np.pi * Rd ** 2) * np.exp(-a / Rd); Rq = np.array([0.05, 0.1, 0.3, 1, 3, 6, 10.]); v2 = ring_v2(f, 40., Rq)
    print('ring sum against Freeman, fractional error in V^2 at', Rq.tolist(), np.round(v2 / stars_v2(Rq, M, Rd) - 1, 5).tolist(), 'mass', hi_mass(f, 40.) / M)
