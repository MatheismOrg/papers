"""SA73 T53: forecast F5 (FORECASTS_T52.txt, b674324f...) with Tian et al.'s own construction. Their accelerations are
reproduced by (i) g_bar = G M_bar f_2D(r*/R_e, n) / r*^2 (projected Sersic mass) and (ii) the isotropic Jeans-Abel g_obs of
their linear dispersion profile, both at one radius r* solved from their g_bar (jeans53.py, jeans53b.py: g_obs within
-0.020 +- 0.025 dex). F5 as written swaps in the NSA Sersic n, R_e and M* (DRPall v3_1_1) for g_bar and the halo; 'full'
also rebuilds g_obs with the new light profile. Ratios to Tian's construction, so the benchmark is recovered exactly with
Tian's values. Usage: python3 -I f5_full.py <dmfit dir> <nsa csv>"""
import sys, csv, re, json
dm, nsa_csv = sys.argv[1], sys.argv[2]
sys.argv = sys.argv[:1] + ['/dev/null']; exec(open('jeans53.py').read().split('rows = []')[0])
sys.path.insert(0, dm)
import numpy as np
from scipy.optimize import brentq
from colossus.cosmology import cosmology
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
nsa = {r['plateifu']: {k: float(v) for k, v in r.items() if k != 'plateifu'} for r in csv.DictReader(open(nsa_csv))}
tab = {}
for l in open(H.PKG / 'data/Clusters/manga_bcg_table1_extracted.txt'):
    if not l.strip(): continue
    ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
    tab[ident.replace('∗', '').replace('*', '')] = dict(z=v[0], n=v[1], re=v[2], xl=v[3], m=v[4], b=v[5], lgM=v[6], sig=v[8], lgb=v[10], lgo=v[12])
a = H.arr([r for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG']); ids = [n.replace('MaNGA_', '') for n in a['name']]
c70 = cosmology.setCosmology('flat70', dict(flat=True, H0=70., Om0=0.3, Ob0=0.049, sigma8=0.81, ns=0.96))
kpa = np.array([c70.angularDiameterDistance(nsa[i]['nsa_z']) / 0.7 * 1e3 * np.pi / 180 / 3600 for i in ids])
reS = np.array([nsa[i]['nsa_sersic_th50'] for i in ids]) * kpa; nS = np.array([nsa[i]['nsa_sersic_n'] for i in ids])
MS = np.array([nsa[i]['nsa_sersic_mass'] for i in ids]) / 0.7 ** 2
reT = np.array([tab[i]['re'] for i in ids]); nT = np.array([tab[i]['n'] for i in ids]); MT = 10 ** np.array([tab[i]['lgM'] for i in ids])
rs = np.array([brentq(lambda x, q=tab[i]: np.log10(G * 10 ** q['lgM'] * f2d(x, q['n']) / (x * q['re']) ** 2 * ACC) - q['lgb'], .01, 50) for i in ids]) * reT
fb_ratio = (MS * f2d(rs / reS, nS)) / (MT * f2d(rs / reT, nT))
go_ratio = np.array([gobs(rs[k], nS[k], reS[k], tab[i]['sig'], tab[i]['m'], tab[i]['b'], reT[k]) / gobs(rs[k], nT[k], reT[k], tab[i]['sig'], tab[i]['m'], tab[i]['b'], reT[k])
                     for k, i in enumerate(ids)])
def resid(gb, go, ms, re, r):
    b = dict(a); b['gb'] = gb; b['ms'] = ms; b['mp'] = ms; b['re'] = re; b['r'] = r; return np.log10(go / H.gpred(M.MODEL, th, b))
rng = np.random.default_rng(5353); IDX = rng.integers(0, 50, (20000, 50))
out = dict(r_star_over_benchmark_r=np.round(np.percentile(rs / a['r'], [0, 50, 100]), 3).tolist(),
           log_gbar_ratio=dict(median=round(float(np.median(np.log10(fb_ratio))), 3), sd=round(float(np.std(np.log10(fb_ratio))), 3)),
           log_gobs_ratio=dict(median=round(float(np.median(np.log10(go_ratio))), 3), sd=round(float(np.std(np.log10(go_ratio))), 3)))
def rep(lab, e):
    c = lambda x: round(float(np.corrcoef(np.log10(x), e)[0, 1]), 3)
    cb = np.array([np.corrcoef(np.log10(reS)[i], e[i])[0, 1] for i in IDX])
    d = dict(rms=round(float(np.sqrt(np.mean(e ** 2))), 4), mean=round(float(np.mean(e)), 3), corr_logRe_sersic=c(reS), corr_logRe_tian=c(reT),
             boot95_sersic=np.round(np.percentile(cb, [2.5, 97.5]), 4).tolist(), p_above_045=round(float(np.mean(cb > .45)), 3), p_below_030=round(float(np.mean(cb < .3)), 3),
             slope_sersic=round(float(np.polyfit(np.log10(reS), e, 1)[0]), 3)); out[lab] = d; print(lab, d, flush=True)
rep('B0 benchmark (Tian values, benchmark radius)', resid(a['gb'], a['go'], a['ms'], a['re'], a['r']))
rep('B1 benchmark values, halo at r*', resid(a['gb'], a['go'], a['ms'], a['re'], rs))
rep('F5 as written: NSA Sersic n, R_e, M* for g_bar and halo (g_obs as published)', resid(a['gb'] * fb_ratio, a['go'], MS, reS, a['r']))
rep('F5 as written, halo at r*', resid(a['gb'] * fb_ratio, a['go'], MS, reS, rs))
rep('F5 full: g_obs also rebuilt with the NSA Sersic light profile', resid(a['gb'] * fb_ratio, a['go'] * go_ratio, MS, reS, rs))
rep('F5 full, sizes and index only (Tian M*)', resid(a['gb'] * fb_ratio * MT / MS, a['go'] * go_ratio, MT, reS, rs))
for k in ('r_star_over_benchmark_r', 'log_gbar_ratio', 'log_gobs_ratio'): print(k, out[k])
json.dump(out, open('res53_f5_full.json', 'w'), indent=1)
# Robustness: the dispersion profile held flat beyond R_cut (instead of Tian's linear extrapolation), for both light profiles.
for cut in (1., 1.5, 2.):
    gr = np.array([gobs(rs[k], nS[k], reS[k], tab[i]['sig'], tab[i]['m'], tab[i]['b'], reT[k], 'flat', cut * reT[k]) /
                   gobs(rs[k], nT[k], reT[k], tab[i]['sig'], tab[i]['m'], tab[i]['b'], reT[k], 'flat', cut * reT[k]) for k, i in enumerate(ids)])
    rep(f'F5 full, dispersion flat beyond {cut} R_e', resid(a['gb'] * fb_ratio, a['go'] * gr, MS, reS, rs))
# Paired bootstrap of the error change, benchmark -> full rebuild
e0 = resid(a['gb'], a['go'], a['ms'], a['re'], a['r']); e1 = resid(a['gb'] * fb_ratio, a['go'] * go_ratio, MS, reS, rs)
dd = np.array([np.sqrt(np.mean(e1[i] ** 2)) - np.sqrt(np.mean(e0[i] ** 2)) for i in IDX])
out['paired_boot_rms_change'] = dict(median=round(float(np.median(dd)), 4), ci95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), p_not_better=round(float(np.mean(dd >= 0)), 4))
print('paired bootstrap rms change (full - benchmark):', out['paired_boot_rms_change'])
# MOND's relation and standard halos on the same two versions of the data
S = M.S; K = H.K; A0 = 1.2e-10
def nu(y): return 1 / (1 - np.exp(-np.sqrt(y)))
def g_lcdm(b):
    m = S.m_exp(b['ms'], b['z'], S.thresholds(b)); m = (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
    c = K.c_dm14(m, b['z']); rv = (H.G * m / (100 * H.hz(b['z']) ** 2)) ** (1 / 3)
    return b['gb'] + H.G * (1 - H.FB) * m * H.ff(b['r'] * c / rv) / H.ff(c) / b['r'] ** 2 * H.ACC
for lab, gb, go, ms, re_, r in (('benchmark', a['gb'], a['go'], a['ms'], a['re'], a['r']), ('full rebuild', a['gb'] * fb_ratio, a['go'] * go_ratio, MS, reS, rs)):
    b = dict(a); b['gb'] = gb; b['ms'] = ms; b['mp'] = ms; b['re'] = re_; b['r'] = r
    for mod, g in (('MOND', gb * nu(gb / A0)), ('standard halos', g_lcdm(b)), ('formula', H.gpred(M.MODEL, th, b))):
        e = np.log10(go / g); d = dict(rms=round(float(np.sqrt(np.mean(e ** 2))), 4), mean=round(float(np.mean(e)), 3), corr_logRe=round(float(np.corrcoef(np.log10(re_), e)[0, 1]), 3))
        out[f'compare|{lab}|{mod}'] = d; print('compare', lab, mod, d)
json.dump(out, open('res53_f5_full.json', 'w'), indent=1)
