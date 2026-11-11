"""Einasto halo with the published Dutton and Maccio (2014) relations, for dynamics and projected lensing.
rho/rho_-2 = exp(-(2/alpha)((r/r_-2)^alpha - 1)); c200 = r200/r_-2.
Einasto c200 (eqs 14-15): log10 c = a + b log10(M200/1e12 h^-1 Msun), b = -0.130 + 0.029 z, a = 0.459 + (0.977-0.459) exp(-0.490 z^1.303).
Shape (eq 17): alpha = 0.0095 nu^2 + 0.155. Peak height (eqs 18-19): log10 nu(M,0) = -0.11 + 0.146 m + 0.0138 m^2 + 0.00123 m^3,
m = log10(M200/1e12 h^-1 Msun); nu(M,z)/nu(M,0) = 0.033 + 0.79(1+z) + 0.176 exp(-1.356 z).
Projected excess surface density: D(X) = DeltaSigma/(rho_-2 r_-2) = 4 I(X)/X^2 - 2X int_0^inf f(X cosh t) cosh t exp(-2t) dt,
with f the dimensionless density and I(X) = int_0^X f(u) u^2 du (derivation in the T17 report). Tabulated in (alpha, log10 X)."""
import numpy as np
from scipy.special import gammainc, gammaln
from scipy.interpolate import RegularGridInterpolator
H_LITTLE = 0.671
def c_ein(m, z):
    z = np.asarray(z, float); a = 0.459 + (0.977 - 0.459) * np.exp(-0.490 * z ** 1.303); b = -0.130 + 0.029 * z
    return 10 ** (a + b * np.log10(np.asarray(m) * H_LITTLE / 1e12))
def alpha(m, z):
    t = np.log10(np.asarray(m) * H_LITTLE / 1e12); z = np.asarray(z, float)
    nu = 10 ** (-0.11 + 0.146 * t + 0.0138 * t ** 2 + 0.00123 * t ** 3) * (0.033 + 0.79 * (1 + z) + 0.176 * np.exp(-1.356 * z))
    return 0.0095 * nu ** 2 + 0.155
def log_I(x, a):
    """log of int_0^x exp(-(2/a)(u^a - 1)) u^2 du."""
    return 2 / a - np.log(a) + (3 / a) * np.log(a / 2) + gammaln(3 / a) + np.log(gammainc(3 / a, 2 / a * np.asarray(x) ** a))
def enc_ratio(x, c, a):
    """M(<x r_-2) / M(<c r_-2)."""
    return gammainc(3 / a, 2 / a * np.asarray(x) ** a) / gammainc(3 / a, 2 / a * np.asarray(c) ** a)
T = np.linspace(0, 40, 4001); WT = np.full(T.size, T[1] - T[0]); WT[0] = WT[-1] = WT[0] / 2     # trapezoid in t
def D_generic(X, f, logI):
    """Dimensionless DeltaSigma for density f(x) with log enclosed-moment function logI(x); X 1-d."""
    X = np.asarray(X, float)[:, None]; ct = np.cosh(T)[None, :]
    integ = (f(X * ct) * ct * np.exp(-2 * T)[None, :]) @ WT
    return 4 * np.exp(logI(X[:, 0])) / X[:, 0] ** 2 - 2 * X[:, 0] * integ
A_GRID = np.round(np.arange(0.14, 0.7001, 0.005), 4); LX_GRID = np.linspace(-3, 3, 481)
_TAB = None
def table():
    global _TAB
    if _TAB is None:
        X = 10 ** LX_GRID; rows = []
        for a in A_GRID:
            d = D_generic(X, lambda x: np.exp(-(2 / a) * (x ** a - 1)), lambda x: log_I(x, a)); rows.append(np.log(d))
        _TAB = RegularGridInterpolator((A_GRID, LX_GRID), np.array(rows), method='cubic')
    return _TAB
def D_ein(X, a):
    X, a = np.broadcast_arrays(np.asarray(X, float), np.asarray(a, float))
    if np.any(a < A_GRID[0]) or np.any(a > A_GRID[-1]) or np.any(X < 10 ** LX_GRID[0]) or np.any(X > 10 ** LX_GRID[-1]):
        raise ValueError('outside Einasto table')
    pts = np.stack([a.ravel(), np.log10(X).ravel()], axis=1)
    return np.exp(table()(pts)).reshape(X.shape)
