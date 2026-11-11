"""SA73 T53: rebuild Tian et al. (2024) MaNGA BCG accelerations from their Table 1 (n, R_e, r_last/R_e, m, b, log M_bar,
<sigma>), to check the construction before swapping in the NSA Sersic profiles. g_obs: isotropic spherical Jeans with Abel
inversion (Tian et al. 2024, eqs. as published): rho sigma_r^2 = -(1/pi) int_r^inf d(I s^2)/dR dR / sqrt(R^2 - r^2), constant
M/L, s(R) = <sigma> (m R/R_e + b); g = -(1/nu) d(nu sigma_r^2)/dr, with R = r cosh t. g_bar: G M(<r)/r^2 with the Sersic
mass inside r, projected (2D) or deprojected (3D, numerical Abel). Usage: python3 -I jeans53.py <table1 txt> [nsa csv]"""
import sys, re, csv, json
import numpy as np
from scipy.special import gammainc
G = 4.30091e-6; ACC = 3.24078e-14
def b_n(n): return 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n ** 2)
T = np.linspace(0, 12, 24001)                                  # cosh(12) ~ 8e4
def sersic_terms(R, n, re):
    b = b_n(n); x = R / re; I = np.exp(-b * (x ** (1 / n) - 1)); u = -b / (n * re) * x ** (1 / n - 1)
    up = -b / (n * re ** 2) * (1 / n - 1) * x ** (1 / n - 2); return I, I * u, I * (u * u + up)
def gobs(r, n, re, sig0, m, b0, re_sig, mode='linear', r_last=None):
    R = r * np.cosh(T); I, I1, I2 = sersic_terms(R, n, re)
    s = sig0 * (m * R / re_sig + b0); s1 = np.full_like(R, sig0 * m / re_sig)
    if mode == 'flat':                                         # dispersion held at its last measured value beyond r_last
        out = R > r_last; s = np.where(out, sig0 * (m * r_last / re_sig + b0), s); s1 = np.where(out, 0., s1)
    h2 = I2 * s * s + 4 * I1 * s * s1 + 2 * I * s1 * s1
    num = np.trapezoid(h2 * np.cosh(T), T); den = np.trapezoid(I1, T)
    return -num / den                                          # (km/s)^2 / kpc
def f2d(x, n): return gammainc(2 * n, b_n(n) * x ** (1 / n))
def f3d(r, n, re):                                             # numerical deprojection: M3D(<r)/M_tot
    rr = np.geomspace(1e-4 * re, 200 * re, 4000); nu = np.array([-np.trapezoid(sersic_terms(q * np.cosh(T), n, re)[1], T) / np.pi for q in rr])
    m3 = np.concatenate([[0], np.cumsum(0.5 * (nu[1:] * rr[1:] ** 2 + nu[:-1] * rr[:-1] ** 2) * np.diff(rr))]) * 4 * np.pi
    return np.interp(r, rr, m3) / m3[-1]
rows = []
for l in open(sys.argv[1]):
    if not l.strip(): continue
    ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
    rows.append(dict(id=ident.replace('∗', '').replace('*', ''), z=v[0], n=v[1], re=v[2], xl=v[3], m=v[4], b=v[5], lgM=v[6], sig=v[8], lgb=v[10], lgo=v[12]))
d = {k: [] for k in ('o_lin', 'o_flat', 'b2', 'b3')}
for q in rows:
    r = q['xl'] * q['re']
    d['o_lin'].append(np.log10(gobs(r, q['n'], q['re'], q['sig'], q['m'], q['b'], q['re']) * ACC) - q['lgo'])
    d['o_flat'].append(np.log10(gobs(r, q['n'], q['re'], q['sig'], q['m'], q['b'], q['re'], 'flat', r) * ACC) - q['lgo'])
    d['b2'].append(np.log10(G * 10 ** q['lgM'] * f2d(q['xl'], q['n']) / r ** 2 * ACC) - q['lgb'])
    d['b3'].append(np.log10(G * 10 ** q['lgM'] * f3d(r, q['n'], q['re']) / r ** 2 * ACC) - q['lgb'])
for k, v in d.items():
    v = np.array(v); print(f'{k:7s} rebuilt minus Tian: median {np.median(v):+.3f}  sd {np.std(v):.3f}  max|dev from median| {np.max(np.abs(v - np.median(v))):.3f}')
json.dump({k: [round(float(x), 4) for x in v] for k, v in d.items()}, open('res53_jeans_check.json', 'w'))
