"""SA73 T53 exploration (not pre-registered, all data): MaNGA BCGs with the formula's halos contracted by their stars
(Cautun et al. 2020, eq. 11, strength s) and the stellar IMF factor (Posacki et al. 2015 trend, T24) kept or removed.
Usage: python3 -I f5_contract.py <dmfit dir>"""
import sys, re; sys.path.insert(0, sys.argv[1])
import numpy as np
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.P.P.S if hasattr(M.P, 'P') else None
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
I = __import__('harness_T24'); sig = I.sigmas()
a = H.arr([r for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG']); lre = np.log10(a['re'])
al = np.array([I.alpha(sig[n], I.DELTA['C']) for n in a['name']])
gp = H.gpred(M.MODEL, th, a); gdm = gp - a['gb']               # formula halo term (uncontracted)
mdm = gdm / H.ACC * a['r'] ** 2 / H.G / (1 - H.FB)             # M_NFW(<r) of the halo (total, as in the harness)
def F(eta): return 0.45 + 0.38 * (eta + 1.16) ** 0.53
def score(gb, s):
    eta = gb / H.ACC * a['r'] ** 2 / H.G / (H.FB * mdm); g = gb + gdm * F(eta) ** s; e = np.log10(a['go'] / g)
    return round(float(np.sqrt(np.mean(e ** 2))), 4), round(float(np.mean(e)), 3), round(float(np.corrcoef(lre, e)[0, 1]), 3), round(float(np.median(F(eta) ** s)), 2)
for lab, gb in (('IMF trend (Salpeter-like, current)', a['gb']), ('no IMF factor (Chabrier)', a['gb'] / al), ('Kroupa-like (x 10^0.035 on Chabrier)', a['gb'] / al * 10 ** .035)):
    for s in (0., .5, 1.):
        print(f'{lab:40s} s={s:3.1f}  rms, mean, corr(logRe), median contraction factor:', score(gb, s))
