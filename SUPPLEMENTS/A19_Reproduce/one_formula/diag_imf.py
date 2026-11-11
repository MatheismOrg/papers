"""All-data what-if (not held out, not part of any formula): the MaNGA BCG offset against a heavier IMF in massive
spheroids. Posacki et al. 2015 (ATLAS3D + SLACS): log alpha = 0.38 log(sigma_e/200 km/s) - 0.06, alpha relative to
Salpeter, rms 0.12 dex. If the MaNGA masses (constant M/L, Tian et al. 2024) are on a Chabrier-like scale, the stellar
g_bar is too low by alpha_Chab = alpha x 10^delta, delta = log(M_Salp/M_Chab) of 0.2 to 0.25. Each BCG's <sigma_los>
(Tian et al. 2024, Table 1) stands in for sigma_e. MOGF at the all-data coefficients."""
import re, numpy as np, harness_T21 as F
H = F.H; B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
dd = H.PKG / 'data/Clusters'; sig = {}
for l in open(dd / 'manga_bcg_table1_extracted.txt'):
    if not l.strip(): continue
    ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
    sig['MaNGA_' + ident.replace('∗', '').replace('*', '')] = v[8]
cl = H.cluster_rows(); a = H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']); th = np.array([11.182185, .173541])
s = np.array([sig[n] for n in a['name']]); print('sigma range', s.min(), s.max(), 'median', np.median(s))
def score(b):
    l = H.losses('MOGF', th, b); er = np.log10(b['go'] / H.gpred('MOGF', th, b)); return np.sqrt(np.mean(list(l.values()))), np.mean(er)
print('as published: rms %.4f mean %+.4f' % score(a))
for delta in (0., .215, .25):
    b = dict(a); b['gb'] = a['gb'] * 10 ** (.38 * np.log10(s / 200) - .06 + delta)
    print('alpha_Chab = Posacki x 10^%.3f (median alpha %.2f): rms %.4f mean %+.4f' % ((delta, np.median(10 ** (.38 * np.log10(s / 200) - .06 + delta))) + score(b)))
