"""SA74 T31: fill SUPPLEMENT_README_T31.tmpl.md from the result files (multihist_T30.json and multihist_skip50_T30.json for 4^4 to 12^4,
a1_results_T31.json for 16^4 and 20^4, au_results_T31.json for a = 0.02, fss/fss_results_T30.json for autocorrelation on 4^4 to 10^4,
the 8^4 rough chain) and write SUPPLEMENT_README.md; prints the numbers."""
import json, math, numpy as np
N = '/home/claude/sa74/regge/r2/night/'; MAC = '/home/claude/sa74/NIGHT/mac_pull/'
M = json.load(open(N + 'multihist_T30.json')); M2 = json.load(open(N + 'multihist_skip50_T30.json'))
A = json.load(open(N + 'a1_results_T31.json')); U = json.load(open(N + 'au_results_T31.json')); F = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results_T30.json'))
Ls = ['4', '5', '6', '8', '10', '12']; FK = A['fixed_k']; P = A['fixed_k_power']; PB = A['fixed_k_power_lowbias']
def pw(p): return '%+.2f(%d)' % (p['power'], round(100 * p['se']))
def fk(k): return ', '.join('%.1f(%d)' % (FK[k][L]['chi'], round(10 * FK[k][L]['se'])) for L in Ls + ['16', '20'])
ch = {r['chain']: r for r in A['chains']}
c12 = [r for r in A['chains'] if r['L'] == 12 and r['finished']]; rr = {r['chain']: r for r in c12}
tau = {}
for key, rows in F['per_chain'].items():
    L, k = key.split('_')
    if 4 <= int(L) <= 10 and 0.2425 <= float(k) <= 0.2575:
        for r in rows:
            if r['n'] >= 800: tau.setdefault(int(L), []).append(r['n'] / max(r['ess'], 1e-9))
tau[12] = [r['tau'] for r in c12 if r['chain'].startswith('A1_') and 0.2425 <= r['k'] <= 0.2575]
for L in (16, 20): tau[L] = [r['tau'] for r in A['chains'] if r['L'] == L and r['finished'] and 0.2425 <= r['k'] <= 0.2575]
def Rs(r): return '%+.4f(%d)' % (r['R'], round(1e4 * r['R_se']))
def cs(r): return '%.1f(%d)' % (r['chi'], round(10 * r['chi_se']))
def cross(L):
    by = {}
    for r in A['chains']:
        if r['L'] == L and r['finished']: by.setdefault(r['k'], []).append(r['R'])
    p = sorted((k, float(np.mean(v))) for k, v in by.items())
    for (k1, r1), (k2, r2) in zip(p, p[1:]):
        if r1 < 0 <= r2: return k1 + (k2 - k1) * (-r1) / (r2 - r1)
    (k1, r1), (k2, r2) = p[-2], p[-1]; return k1 + (k2 - k1) * (-r1) / (r2 - r1)
h = lambda L: (0.16 + 1.11 / L) * L ** 0.64; hp = math.log(h(20) / h(4)) / math.log(5)
sc = A['start_check']
def scs(key): s = sc[key]; return '%+.4f (%+.1f standard errors)' % (s['dR'], s['z'])
a = np.loadtxt(MAC + 'Q_L8_k0.250_rough2.dat', comments='#'); m = a[a[:, 1] == 1]; b = np.array_split(np.arange(len(m)), 6)
f = np.loadtxt(MAC + 'N_L8_k0.250.dat', comments='#'); f = f[f[:, 1] == 1]
au = {(r['L'], round(r['k'], 3)): r for r in U['chains']}; x8 = au[(8, 0.46)]; cu = {int(L): v for L, v in U['crossing'].items()}
s16 = ['A1_L16_k0.248', 'A2_L16_k0.248_r2', 'A2_L16_k0.244_r1', 'A2_L16_k0.244_r2', 'A3_L16_k0.244_f1', 'A3_L16_k0.244_f2', 'A2_L16_k0.250', 'A1_L16_k0.252', 'A2_L16_k0.252_r2']
V = dict(f244_list=fk('0.2440'), f248_list=fk('0.2480'), f252_list=fk('0.2520'),
         p244=pw(P['0.2440']), p248=pw(P['0.2480']), p252=pw(P['0.2520']), p244b=pw(PB['0.2440']), p248b=pw(PB['0.2480']), p252b=pw(PB['0.2520']),
         hp='%.2f' % hp, hr='%.2f' % (h(20) / h(4)), g1='%.1f' % ((hp - P['0.2440']['power']) / P['0.2440']['se']), g2='%.1f' % ((hp - PB['0.2440']['power']) / PB['0.2440']['se']),
         crit='%.1f' % 5 ** 0.94,
         peak_list=', '.join('%.1f(%d)' % (M['per_L'][L]['peak_chi'], round(10 * M['per_L'][L]['peak_chi_boot_sd'])) for L in Ls),
         pw=pw(M['peak_power']), pw2=pw(M2['peak_power']), n12=str(len([r for r in c12 if r['chain'].startswith('A1_')])),
         dq12='%.3f' % max(abs(r['meanq_halves'][1] - r['meanq_halves'][0]) for r in c12 if r['chain'].startswith('A1_')),
         dR12='%.4f' % max(abs(r['R_halves'][1] - r['R_halves'][0]) for r in c12 if r['chain'].startswith('A1_')),
         d12_list=', '.join('%+.4f' % (rr['A1_L12_k%.3f_r1' % k]['R'] - rr['A1_L12_k%.3f_r2' % k]['R']) for k in (0.244, 0.246, 0.248, 0.250, 0.252, 0.254)),
         c14a='%.1f' % ch['Q_L14_k0.2475']['chi'], c14b='%.1f' % ch['Q_L14_k0.2525']['chi'],
         chains16='; '.join('%s: R %s, chi %s, tau %d' % (c, Rs(ch[c]), cs(ch[c]), round(ch[c]['tau'])) for c in s16),
         p16_244='%.1f(%d)' % (FK['0.2440']['16']['chi'], round(10 * FK['0.2440']['16']['se'])), l16_244=' and '.join('%.2f' % x for x in FK['0.2440']['16']['lineage_chi']),
         sc248=scs('A1_L16_k0.248 - A2_L16_k0.248_r2'), sc252=scs('A1_L16_k0.252 - A2_L16_k0.252_r2'), sc244=scs('A2_L16_k0.244_r1 - A2_L16_k0.244_r2'),
         scf1='%.4f (%.1f standard errors)' % (abs(sc['A2_L16_k0.244_r1 - A3_L16_k0.244_f1']['dR']), abs(sc['A2_L16_k0.244_r1 - A3_L16_k0.244_f1']['z'])),
         scf2='%.4f (%.1f standard errors)' % (abs(sc['A2_L16_k0.244_r2 - A3_L16_k0.244_f2']['dR']), abs(sc['A2_L16_k0.244_r2 - A3_L16_k0.244_f2']['z'])),
         six_r1=' '.join('%+.4f' % v for v in ch['A2_L16_k0.244_r1']['R_sixths']), six_r2=' '.join('%+.4f' % v for v in ch['A2_L16_k0.244_r2']['R_sixths']),
         chains20='; '.join('k = %.3f: R %s, chi %s, tau %d, R by sixths %s' % (ch[c]['k'], Rs(ch[c]), cs(ch[c]), round(ch[c]['tau']), ' '.join('%+.4f' % v for v in ch[c]['R_sixths']))
                             for c in ('A5_L20_k0.244', 'A5_L20_k0.248', 'A5_L20_k0.252')),
         cr16='%.4f' % cross(16), cr20='%.4f' % cross(20),
         cr_lo='%.4f' % min([M['per_L'][L]['crossing'] for L in Ls] + [cross(16), cross(20)]), cr_hi='%.4f' % max([M['per_L'][L]['crossing'] for L in Ls] + [cross(16), cross(20)]),
         cr12='%.4f' % M['per_L']['12']['crossing'],
         tau_list=', '.join('%d' % round(float(np.median(tau[L]))) for L in (4, 5, 6, 8, 10, 12, 16, 20)),
         bias16='%d%% to %d%%' % (round(100 * min(ch[c]['chi_low_bias'] for c in s16)), round(100 * max(ch[c]['chi_low_bias'] for c in s16))),
         rough2_250='R by sixths of its 3,000 measured sweeps ' + ', '.join('%+.3f' % m[i, 4].mean() for i in b) + ', last half %+.4f against %+.4f for the flat start' % (m[len(m) // 2:, 4].mean(), f[:, 4].mean()),
         au46=', '.join('%.2f' % au[(L, 0.46)]['chi'] for L in (4, 6, 8, 10, 12)), au_pow46='%+.2f' % U['chi_power']['0.46']['power'],
         au_R10=', '.join('%s' % Rs(au[(10, k)]) for k in (0.45, 0.46, 0.47)), au_R12=', '.join('%s' % Rs(au[(12, k)]) for k in (0.45, 0.46, 0.47)),
         au_c10=', '.join('%.2f' % au[(10, k)]['chi'] for k in (0.45, 0.46, 0.47)), au_c12=', '.join('%.2f' % au[(12, k)]['chi'] for k in (0.45, 0.46, 0.47)),
         au_cr4='%.3f' % cu[4], au_cr6='%.3f' % cu[6], au_cr8='%.3f' % (0.46 - x8['R'] / x8['slope']), au_cr10='%.3f' % cu[10], au_cr12='%.3f' % cu[12], au_R8=Rs(x8))
s = open(N + 'SUPPLEMENT_README_T31.tmpl.md').read()
for k, v in V.items(): s = s.replace('<<' + k + '>>', v)
assert '<<' not in s, s[s.index('<<'):s.index('<<') + 40]
open(N + 'SUPPLEMENT_README.md', 'w').write(s)
for k, v in V.items(): print(f'{k}: {v}')
