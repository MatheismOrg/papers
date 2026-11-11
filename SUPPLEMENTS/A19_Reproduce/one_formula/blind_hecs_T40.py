"""T40 (PROTOCOL_T40.txt): blind test of the frozen T37 formula on the HeCS-omnibus clusters (Sohn et al. 2020, Tables 2, 3):
predicted halo mass from the BCG's stellar mass against the caustic M200. Nothing fitted.
Usage: python3 -I blind_hecs_T40.py <t2.csv> <t3.csv> <out.json>"""
import sys, json, csv, re
sys.path.insert(0, '/tmp/claude-0/-home-claude/0a1f618d-0405-5db3-a0ab-84f1f9545369/scratchpad/dmfit')
import numpy as np
import harness_T37 as S
from scipy.optimize import brentq
H = S.H; K = H.K; B = S.B; F = S.F
ABC = np.array([11.7481, .8178, -.7119]); B.FROZEN['abc'] = ABC; TH37 = np.array([11.1768]); TH26 = np.array([11.1868])
man = [(r['ms'], r['re']) for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG']
man = {(a, b) for a, b in man}; P_RE = np.polyfit(np.log10([a for a, _ in man]), np.log10([b for _, b in man]), 1)
t2 = {r['ID']: r for r in csv.DictReader(open(sys.argv[1]))}; t3 = {r['ID']: r for r in csv.DictReader(open(sys.argv[2]))}
rows = []
for k, r in t3.items():
    if k not in t2: continue
    lm, m200 = float(r['Mass']), float(t2[k]['M200'])
    if lm <= 0 or m200 <= 0: continue
    rows.append(dict(ID=k, lm=lm, z=float(t2[k]['zspec']), m200=m200 * 1e14, r200=float(t2[k]['R200']) * 1e3, mx=float(t2[k]['Mx'])))
def norm(s):
    s = re.sub(r'[^A-Z0-9+\-]', '', s.upper()); m = re.match(r'^(A|ABELL)0*(\d+)$', s)
    return 'A' + m.group(2) if m else s
CLASH = ['A209', 'A2261', 'A383', 'A611', 'MS2137', 'RXJ1347', 'RXJ1532', 'RXJ2129', 'RXJ2248', 'MACSJ0329', 'MACSJ0416', 'MACSJ0429',
         'MACSJ0647', 'MACSJ0717', 'MACSJ0744', 'MACSJ1115', 'MACSJ1149', 'MACSJ1206', 'MACSJ1720', 'MACSJ1931']
def in_clash(i):
    n = norm(i)
    return any(n == c or (n.startswith(c) and not n[len(c)].isdigit()) for c in CLASH)
def evaluate(rows, mstar_fac=1., m200_fac=1.):
    ms = 10 ** np.array([x['lm'] for x in rows]) * mstar_fac; z = np.array([x['z'] for x in rows]); mp = ms.copy()
    re_ = 10 ** np.polyval(P_RE, np.log10(ms)); m200 = np.array([x['m200'] for x in rows]) * m200_fac
    r200 = np.array([x['r200'] for x in rows]) * m200_fac ** (1 / 3); mx = np.array([x['mx'] for x in rows])
    rule = mp / H.FB + 10 ** ABC[0] * (ms / 1e10) ** ABC[1] * (re_ / 3) ** ABC[2]
    preds = {'T37_selected_1e14': S.mass(TH37, ms, mp, re_, z, np.full(len(ms), 1e14)),
             'T37_no_selection': S.mass(TH37, ms, mp, re_, z, np.zeros(len(ms))),
             'T26_mean_inversion': F.mass('MOGF0IE', TH26, ms, mp, re_, z),
             'rule_only': (rule ** 8 + 1e8 ** 8) ** (1 / 8),
             'T37_selected_10^13.8': S.mass(TH37, ms, mp, re_, z, np.full(len(ms), 10 ** 13.8)),
             'T37_selected_10^14.2': S.mass(TH37, ms, mp, re_, z, np.full(len(ms), 10 ** 14.2))}
    obs = np.log10(m200); out = {}
    def nfw_in(m, zz, r):
        c = K.c_dm14(m, zz); rv = (H.G * m / (100 * H.hz(zz) ** 2)) ** (1 / 3); return m * H.ff(r * c / rv) / H.ff(c), c, rv
    for lab, p in preds.items():
        d = np.log10(p) - obs; min_, c, rv = nfw_in(p, z, r200); dg = np.log10(min_ / m200)
        m500 = []
        for m, cc, rr, zz in zip(p, c, rv, z):
            rho = lambda r: m * H.ff(r * cc / rr) / H.ff(cc) / (4 / 3 * np.pi * r ** 3)
            rc = m / (4 / 3 * np.pi * rr ** 3) / 200                      # critical density: M200 = 200 rho_c (4/3) pi r200^3
            r5 = brentq(lambda r: rho(r) - 500 * rc, 1e-3 * rr, rr); m500.append(m * H.ff(r5 * cc / rr) / H.ff(cc))
        k = mx > 0; dx = np.log10(np.array(m500)[k] / (mx[k] * 1e14))
        out[lab] = dict(N=int(len(d)), mean=round(float(d.mean()), 3), median=round(float(np.median(d)), 3), rms=round(float(np.sqrt(np.mean(d ** 2))), 3),
                        corr_logMstar=round(float(np.corrcoef(d, np.log10(ms))[0, 1]), 3), frac_within_0p3=round(float(np.mean(np.abs(d) < .3)), 3),
                        accel_R200_mean=round(float(dg.mean()), 3), accel_R200_rms=round(float(np.sqrt(np.mean(dg ** 2))), 3),
                        mcxc_N=int(k.sum()), mcxc_mean=round(float(dx.mean()), 3), mcxc_rms=round(float(np.sqrt(np.mean(dx ** 2))), 3))
    d0 = np.median(obs) - obs; out['constant_median_obs'] = dict(N=int(len(d0)), rms=round(float(np.sqrt(np.mean(d0 ** 2))), 3))
    out['rule_term_fraction_of_T37_max'] = round(float(np.max(rule / preds['T37_selected_1e14'])), 3)
    out['obs_logM200_range'] = [round(float(obs.min()), 2), round(float(np.median(obs)), 2), round(float(obs.max()), 2)]
    out['logMstar_range'] = [round(float(np.log10(ms).min()), 2), round(float(np.median(np.log10(ms))), 2), round(float(np.log10(ms).max()), 2)]
    return out
res = dict(N_t2=len(t2), N_t3=len(t3), N_used=len(rows), clash_overlap=[x['ID'] for x in rows if in_clash(x['ID'])],
           P_RE=[round(float(v), 4) for v in P_RE], all=evaluate(rows), without_clash=evaluate([x for x in rows if not in_clash(x['ID'])]),
           h_converted=evaluate(rows, (0.70 / 0.671) ** 2, 0.70 / 0.671))
a = res['all']; t = a['T37_selected_1e14']
res['passes'] = bool(abs(t['mean']) <= .10 and t['rms'] <= .30 and t['rms'] < a['T37_no_selection']['rms'] and t['rms'] < a['T26_mean_inversion']['rms'])
json.dump(res, open(sys.argv[3], 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('h_converted',)}, indent=1)); print('T40 passes:', res['passes'])
