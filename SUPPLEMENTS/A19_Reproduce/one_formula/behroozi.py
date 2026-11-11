"""Behroozi, Wechsler and Conroy (2013, ApJ 770, 57; arXiv:1207.6105) median stellar mass - halo mass relation for central
galaxies, Section 5 best fit (read from ar5iv in the built-in browser on Gary's Mac, 6 November 2026):
log10 M*(Mh) = log10(eps M1) + f(log10(Mh/M1)) - f(0), f(x) = -log10(10^(alpha x) + 1) + delta (log10(1 + e^x))^gamma / (1 + exp(10^-x));
nu = exp(-4 a^2); log10 eps = -1.777 + (-0.006 (a-1) - 0.000 z) nu - 0.119 (a-1); log10 M1 = 11.514 + (-1.793 (a-1) - 0.251 z) nu;
alpha = -1.412 + 0.731 (a-1) nu; delta = 3.508 + (2.608 (a-1) - 0.043 z) nu; gamma = 0.316 + (1.319 (a-1) + 0.279 z) nu.
Halo mass: virial (Bryan and Norman); h = 0.7; Chabrier IMF. Inverted here by bisection in log Mh within [9, 16.5],
then converted to M200c with the NFW profile and the Dutton-Maccio concentration (dc14.m200_to_mvir, inverted)."""
import numpy as np
import harness_T15_hashed as H
import dc14 as DC
def logms(lmh, z):
    a = 1 / (1 + z); nu = np.exp(-4 * a * a)
    le = -1.777 + (-0.006 * (a - 1) - 0.000 * z) * nu - 0.119 * (a - 1)
    lm1 = 11.514 + (-1.793 * (a - 1) - 0.251 * z) * nu
    al = -1.412 + 0.731 * (a - 1) * nu; de = 3.508 + (2.608 * (a - 1) - 0.043 * z) * nu; ga = 0.316 + (1.319 * (a - 1) + 0.279 * z) * nu
    f = lambda x: -np.log10(10 ** (al * x) + 1) + de * np.log10(1 + np.exp(x)) ** ga / (1 + np.exp(10 ** -x))
    return le + lm1 + f(lmh - lm1) - f(0.)
CACHE = {}
def halo(ms, z):
    """M200c (Msun) of the central halo for stellar mass ms at redshift z."""
    out = np.empty(np.size(ms))
    for i, (m, zz) in enumerate(zip(np.ravel(ms), np.ravel(z))):
        k = (float(m), float(zz))
        if k not in CACHE:
            lo, hi = 9., 16.5; t = np.log10(m)
            for _ in range(60):
                mid = (lo + hi) / 2
                if logms(mid, zz) < t: lo = mid
                else: hi = mid
            mvir = 10 ** ((lo + hi) / 2)
            lo2, hi2 = np.log10(mvir) - .3, np.log10(mvir)        # invert M200c -> Mvir
            for _ in range(50):
                mid = (lo2 + hi2) / 2; m2 = 10 ** mid
                if DC.m200_to_mvir(m2, H.K.c_dm14(m2, zz), zz) < mvir: lo2 = mid
                else: hi2 = mid
            CACHE[k] = 10 ** ((lo2 + hi2) / 2)
        out[i] = CACHE[k]
    return out.reshape(np.shape(ms))
