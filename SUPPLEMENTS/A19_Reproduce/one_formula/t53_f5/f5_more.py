"""SA73 T53: test of forecast F5 (FORECASTS_T52.txt, b674324f...) with the NASA-Sloan Atlas values in the MaNGA DRPall
catalog (drpall-v3_1_1.fits, data.sdss.org DR17, sha256 9d34371c...; extracted CSV sha256 5d481c17...).
Usage: python3 -I f5_test.py <dmfit dir> <nsa csv>"""
import sys, csv, re, json
sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.special import gammainc
from colossus.cosmology import cosmology
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
I = __import__('harness_T24')
nsa = {r['plateifu']: {k: float(v) for k, v in r.items() if k != 'plateifu'} for r in csv.DictReader(open(sys.argv[2]))}
tab = {}
for l in open(H.PKG / 'data/Clusters/manga_bcg_table1_extracted.txt'):
    if not l.strip(): continue
    ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
    tab[ident.replace('∗', '').replace('*', '')] = dict(z=v[0], n=v[1], Re=v[2], lgM=v[6], sig=v[8], lgb=v[10], lgo=v[12])
raw = {r['name']: r for r in __import__('harness_T15_hashed').cluster_rows() if r['kind'] == 'MaNGA_BCG'}   # unpatched copy
a = H.arr([r for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG']); ids = [n.replace('MaNGA_', '') for n in a['name']]
c70 = cosmology.setCosmology('flat70', dict(flat=True, H0=70., Om0=0.3, Ob0=0.049, sigma8=0.81, ns=0.96))
def kpc_per_arcsec(z): return c70.angularDiameterDistance(z) / 0.7 * 1e3 * np.pi / 180 / 3600
def b_n(n): return 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n ** 2)
def f_enc(x, n):   # Prugniel-Simien deprojected Sersic: fraction of mass inside r = x R_e
    p = 1 - 0.6097 / n + 0.05463 / n ** 2; return gammainc(n * (3 - p), b_n(n) * x ** (1 / n))
alpha = np.array([I.alpha(tab[i]['sig'], I.DELTA['C']) for i in ids])
reT = np.array([tab[i]['Re'] for i in ids]); nT = np.array([tab[i]['n'] for i in ids]); MT = 10 ** np.array([tab[i]['lgM'] for i in ids])
s = np.array([kpc_per_arcsec(nsa[i]['nsa_z']) for i in ids])
reP = np.array([nsa[i]['nsa_elpetro_th50_r'] for i in ids]) * s; reS = np.array([nsa[i]['nsa_sersic_th50'] for i in ids]) * s
nS = np.array([nsa[i]['nsa_sersic_n'] for i in ids])
MP = np.array([nsa[i]['nsa_elpetro_mass'] for i in ids]) / 0.7 ** 2; MS = np.array([nsa[i]['nsa_sersic_mass'] for i in ids]) / 0.7 ** 2
r = a['r']; gb_raw = np.array([raw['MaNGA_' + i]['gb'] for i in ids])
out = {}
out['Tian_Re_vs_NSA_elpetro'] = dict(median_log_ratio=round(float(np.median(np.log10(reT / reP))), 3), sd=round(float(np.std(np.log10(reT / reP))), 3))
out['Tian_n_vs_NSA_sersic_n'] = round(float(np.max(np.abs(nT - np.round(nS, 1)))), 3)
out['NSA_sersic_over_elpetro_Re'] = dict(median_log=round(float(np.median(np.log10(reS / reP))), 3), sd=round(float(np.std(np.log10(reS / reP))), 3))
out['Tian_M_vs_NSA'] = dict(elpetro=round(float(np.median(np.log10(MT / MP))), 3), sersic=round(float(np.median(np.log10(MT / MS))), 3))
gT = H.G * MT * f_enc(r / reT, nT) / r ** 2 * H.ACC           # Tian's g_bar rebuilt from his own n, R_e, M*
out['check_Tian_gbar_rebuilt'] = dict(median_log_ratio=round(float(np.median(np.log10(gb_raw / gT))), 3), sd=round(float(np.std(np.log10(gb_raw / gT))), 3))
def resid(gb, ms, re):
    b = dict(a); b['gb'] = gb; b['ms'] = ms; b['mp'] = ms; b['re'] = re; return np.log10(a['go'] / H.gpred(M.MODEL, th, b))
def rep(lab, e, size):
    d = dict(rms=round(float(np.sqrt(np.mean(e ** 2))), 4), mean=round(float(np.mean(e)), 3), corr_logRe=round(float(np.corrcoef(np.log10(size), e)[0, 1]), 3),
             slope=round(float(np.polyfit(np.log10(size), e, 1)[0]), 3)); out[lab] = d; print(lab, d, flush=True)

e0 = resid(a['gb'], a['ms'], a['re'])
gT_shape = f_enc(r / reT, nT); fS = f_enc(r / reS, nS) / gT_shape
e5 = resid(a['gb'] * fS * MS / MT, MS, reS)
rng = np.random.default_rng(53)
def boot(e, x, nb=20000):
    idx = rng.integers(0, len(e), (nb, len(e))); c = [np.corrcoef(x[i], e[i])[0, 1] for i in idx]; return np.round(np.percentile(c, [2.5, 50, 97.5]), 3).tolist()
lr = lambda v: np.log10(v)
out2 = {}
out2['V5_corr_vs'] = dict(sersic=round(float(np.corrcoef(lr(reS), e5)[0, 1]), 3), tian=round(float(np.corrcoef(lr(reT), e5)[0, 1]), 3), elpetro=round(float(np.corrcoef(lr(reP), e5)[0, 1]), 3),
                          boot95_sersic=boot(e5, lr(reS)), frac_boot_above_045=None)
idx = rng.integers(0, 50, (20000, 50)); cc = np.array([np.corrcoef(lr(reS)[i], e5[i])[0, 1] for i in idx])
out2['V5_corr_vs']['frac_boot_above_045'] = round(float(np.mean(cc > .45)), 3); out2['V5_corr_vs']['frac_boot_below_030'] = round(float(np.mean(cc < .3)), 3)
out2['V0_boot95_tian'] = boot(e0, lr(reT))
out2['Tian_M_over_elpetro_sd'] = round(float(np.std(lr(MT / MP))), 4)
out2['corr_logMS_MT_vs_logReT'] = round(float(np.corrcoef(lr(MS / MT), lr(reT))[0, 1]), 3)
out2['corr_logMS_MP_vs_e0'] = round(float(np.corrcoef(lr(MS / MP), e0)[0, 1]), 3)
# partial correlations of e0
def partial(y, x, z):
    X = np.c_[np.ones_like(z), z]; ry = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]; rx = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]; return round(float(np.corrcoef(rx, ry)[0, 1]), 3)
out2['partial_e0_logReT_given_logMS_MT'] = partial(e0, lr(reT), lr(MS / MT))
out2['partial_e0_logMS_MT_given_logReT'] = partial(e0, lr(MS / MT), lr(reT))
for lab, e, size in [('V9 g_bar kept, halo from Sersic M* (Tian R_e)', resid(a['gb'], MS, reT), reT),
                     ('V10 g_bar kept, halo from Sersic M* and Sersic R_e', resid(a['gb'], MS, reS), reS),
                     ('V11 g_bar kept, halo from elpetro M* (Tian R_e)', resid(a['gb'], MP, reT), reT),
                     ('V12 g_bar kept, halo from Tian M* x 10^0.1 (Tian R_e)', resid(a['gb'], MT * 10 ** .1, reT), reT)]:
    d = dict(rms=round(float(np.sqrt(np.mean(e ** 2))), 4), mean=round(float(np.mean(e)), 3), corr_logRe=round(float(np.corrcoef(lr(size), e)[0, 1]), 3),
             corr_logReS=round(float(np.corrcoef(lr(reS), e)[0, 1]), 3), slope=round(float(np.polyfit(lr(size), e, 1)[0]), 3)); out2[lab] = d
for k, v in out2.items(): print(k, v, flush=True)
json.dump(out2, open('res53_f5_more.json', 'w'), indent=1)
