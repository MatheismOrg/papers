"""Di Cintio et al. (2014, MNRAS 441, 2986) feedback-dependent halo profile, for dynamics and projected lensing.
rho = rho_s / [(r/r_s)^g (1 + (r/r_s)^a)^((b-g)/a)]   (eq. 1)
X = log10(M*/M_halo), valid -4.1 < X < -1.3 (eq. 3); NFW (1, 3, 1) below -4.1:
  a = 2.94 - log10[(10^(X+2.33))^-1.08 + (10^(X+2.33))^2.29]
  b = 4.23 + 1.34 X + 0.26 X^2
  g = -0.06 + log10[(10^(X+2.56))^-0.68 + 10^(X+2.56)]
Concentration (eq. 6): c_hydro / c_DMO = 1 + 0.00003 exp(3.4 (X + 4.5)), c = R_vir / r_-2;
r_-2 = [(2 - g)/(b - 2)]^(1/a) r_s (eq. 2). M_halo is the virial mass (Bryan and Norman overdensity here).
Enclosed moment: I(x) = int_0^x u^(2-g) (1+u^a)^(-(b-g)/a) du = x^(3-g)/(3-g) 2F1((b-g)/a, (3-g)/a; 1+(3-g)/a; -x^a).
Projected: D(X) = DeltaSigma/(rho_s r_s) = 4 I(X)/X^2 - 2X int_0^inf f(X cosh t) cosh t exp(-2t) dt."""
import numpy as np
from scipy.special import hyp2f1
def params(X):
    X = np.asarray(X, float); Xc = np.clip(X, -4.1, -1.3)
    a = 2.94 - np.log10((10 ** (Xc + 2.33)) ** -1.08 + (10 ** (Xc + 2.33)) ** 2.29)
    b = 4.23 + 1.34 * Xc + 0.26 * Xc ** 2
    g = -0.06 + np.log10((10 ** (Xc + 2.56)) ** -0.68 + 10 ** (Xc + 2.56))
    cf = 1 + 3e-5 * np.exp(3.4 * (Xc + 4.5))
    low = X <= -4.1
    return np.where(low, 1., a), np.where(low, 3., b), np.where(low, 1., g), np.where(low, 1., cf)
def blend(p, w):
    """Blend DC14 parameters toward NFW (1, 3, 1, 1) with weight w on DC14."""
    a, b, g, cf = p; return 1 + w * (a - 1), 3 + w * (b - 3), 1 + w * (g - 1), 1 + w * (cf - 1)
def I(x, a, b, g):
    x = np.asarray(x, float); s = 3 - g
    return x ** s / s * hyp2f1((b - g) / a, s / a, 1 + s / a, -x ** a)
def rs_over_rm2(a, b, g):
    return ((2 - g) / (b - 2)) ** (-1 / a)
TN, TW = np.polynomial.legendre.leggauss(48); TMAX = 12.; TT = .5 * TMAX * (TN + 1); WT = .5 * TMAX * TW
CH = np.cosh(TT); E2 = np.exp(-2 * TT)
def D(X, a, b, g):
    """Dimensionless DeltaSigma; X, a, b, g broadcastable arrays (no trailing t axis)."""
    X, a, b, g = np.broadcast_arrays(*[np.asarray(v, float) for v in (X, a, b, g)])
    x = X[..., None] * CH; aa, bb, gg = a[..., None], b[..., None], g[..., None]
    f = x ** -gg * (1 + x ** aa) ** (-(bb - gg) / aa)
    return 4 * I(X, a, b, g) / X ** 2 - 2 * X * ((f * CH * E2) @ WT)
def delta_vir(z, om0=.3175):
    e2 = om0 * (1 + z) ** 3 + 1 - om0; x = om0 * (1 + z) ** 3 / e2 - 1
    return 18 * np.pi ** 2 + 82 * x - 39 * x ** 2
def m200_to_mvir(m, c, z):
    """Virial mass (Bryan-Norman, critical) of an NFW halo with M200c = m and c200 = c."""
    ff = lambda t: np.log1p(t) - t / (1 + t); dv = delta_vir(np.asarray(z, float))
    lo, hi = np.full(np.shape(m), .5), np.full(np.shape(m), 3.)
    for _ in range(60):
        y = (lo + hi) / 2; dens = 200 * ff(c * y) / ff(c) / y ** 3      # mean density inside y r200, units of rho_crit
        lo = np.where(dens > dv, y, lo); hi = np.where(dens > dv, hi, y)
    y = (lo + hi) / 2
    return m * ff(c * y) / ff(c)
