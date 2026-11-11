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
e0 = resid(a['gb'], a['ms'], a['re']); rep('V0 benchmark, vs Tian R_e', e0, a['re'])
rep('V1 benchmark residuals, vs NSA Sersic R_e', e0, reS)
rep('V1b benchmark residuals, vs NSA elpetro R_e', e0, reP)
rep('V2 F5: g_bar from NSA Sersic (n, R_e, M*) x IMF, halo from NSA Sersic M*, R_e', resid(H.G * MS * f_enc(r / reS, nS) / r ** 2 * H.ACC * alpha, MS, reS), reS)
rep('V3 sizes only: NSA Sersic n, R_e with Tian M*', resid(H.G * MT * f_enc(r / reS, nS) / r ** 2 * H.ACC * alpha, MT, reS), reS)
rep('V4 masses only: NSA Sersic M* with Tian n, R_e', resid(gT * (MS / MT) * alpha, MS, reT), reT)
for k in ('Tian_Re_vs_NSA_elpetro', 'Tian_n_vs_NSA_sersic_n', 'NSA_sersic_over_elpetro_Re', 'Tian_M_vs_NSA', 'check_Tian_gbar_rebuilt'): print(k, out[k])
json.dump(out, open('res53_f5.json', 'w'), indent=1)
# Differential rebuilds: keep Tian's g_bar (patched with the IMF factor) and change only what the NSA values imply.
gT_shape = f_enc(r / reT, nT)
fS = f_enc(r / reS, nS) / gT_shape                       # mass inside r if the Sersic size and index held, relative to Tian's
rep('V5 differential: Sersic size and Sersic mass', resid(a['gb'] * fS * MS / MT, MS, reS), reS)
rep('V6 differential: Sersic size only (Tian M*)', resid(a['gb'] * fS, MT, reS), reS)
rep('V7 differential: Sersic mass only (Tian R_e)', resid(a['gb'] * MS / MT, MS, reT), reT)
rep('V8 differential: elpetro mass only (Tian R_e)', resid(a['gb'] * MP / MT, MP, reT), reT)
print('median enclosed-fraction ratio Sersic/Tian at r', round(float(np.median(fS)), 3), 'sd log', round(float(np.std(np.log10(fS))), 3))
print('corr(log M_S/M_T, e0)', round(float(np.corrcoef(np.log10(MS / MT), e0)[0, 1]), 3), ' corr(log R_S/R_P, e0)', round(float(np.corrcoef(np.log10(reS / reP), e0)[0, 1]), 3),
      ' corr(n, e0)', round(float(np.corrcoef(nT, e0)[0, 1]), 3))
json.dump(out, open('res53_f5.json', 'w'), indent=1)
