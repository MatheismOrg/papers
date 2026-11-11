"""SA74 T30: fill SUPPLEMENT_README_T30.tmpl.md from the result files (multihist_T30.json, multihist_skip50_T30.json, a1_results_T30.json,
au_results.json, fss/fss_results.json, the 8^4 rough chain) and write SUPPLEMENT_README.md; prints the numbers for the reply."""
import json, sys, numpy as np
N = '/home/claude/sa74/regge/r2/night/'; MAC = '/home/claude/sa74/NIGHT/mac_pull/'
M = json.load(open(N + 'multihist_T30.json')); M2 = json.load(open(N + 'multihist_skip50_T30.json'))
A = json.load(open(N + 'a1_results_T30.json')); U = json.load(open(N + 'au_results.json')); F = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results_T30.json'))
Ls = ['4', '5', '6', '8', '10', '12']
def pw(p): return '%+.2f(%d)' % (p['power'], round(100 * p['se']))
def fk(k, with16=True):
    v = ['%.1f(%d)' % (M['per_L'][L]['fixed_k'][k]['chi'], round(10 * M['per_L'][L]['fixed_k'][k]['boot_sd'])) for L in Ls]
    if with16: x = A['fixed_k'][k]['16']; v.append('%.1f(%d)' % (x['chi'], round(10 * x['se'])))
    return ', '.join(v)
c12 = [r for r in A['chains'] if r['L'] == 12 and r['finished']]
rr = {r['chain']: r for r in c12}
r244 = sorted((r['chain'], r['R'], r['R_se']) for r in c12 if abs(r['k'] - 0.244) < 1e-9)
tau = {}
for key, rows in F['per_chain'].items():
    L, k = key.split('_')
    if 4 <= int(L) <= 10 and 0.2425 <= float(k) <= 0.2575:
        for r in rows:
            if r['n'] >= 800: tau.setdefault(int(L), []).append(r['n'] / max(r['ess'], 1e-9))
tau[12] = [r['tau'] for r in c12 if 0.2425 <= r['k'] <= 0.2575]
c16 = {r['k']: r for r in A['chains'] if r['L'] == 16 and r['finished']}; tau[16] = [r['tau'] for r in c16.values()]
Lt = Ls + ['16']
a = np.loadtxt(MAC + 'Q_L8_k0.250_rough2.dat', comments='#'); m = a[a[:, 1] == 1]; b = np.array_split(np.arange(len(m)), 6)
f = np.loadtxt(MAC + 'N_L8_k0.250.dat', comments='#'); f = f[f[:, 1] == 1]
au = {(r['L'], round(r['k'], 3)): r for r in U['chains']}; x8 = au[(8, 0.46)]
cr16 = 0.248 + 0.004 * (0 - c16[0.248]['R']) / (c16[0.252]['R'] - c16[0.248]['R'])
V = dict(f248_list=fk('0.2480'), f252_list=fk('0.2520'), f244_list=fk('0.2440', False),
         p248=pw(A['fixed_k_power']['0.2480']), p252=pw(A['fixed_k_power']['0.2520']), p244=pw(A['fixed_k_power']['0.2440']),
         peak_list=', '.join('%.1f(%d)' % (M['per_L'][L]['peak_chi'], round(10 * M['per_L'][L]['peak_chi_boot_sd'])) for L in Ls),
         pw=pw(M['peak_power']), pw2=pw(M2['peak_power']), n12=str(len(c12)),
         r12_244=' and '.join('%+.4f(%d)' % (R, round(1e4 * se)) for _, R, se in r244),
         dq12='%.3f' % max(abs(r['meanq_halves'][1] - r['meanq_halves'][0]) for r in c12),
         dR12='%.4f' % max(abs(r['R_halves'][1] - r['R_halves'][0]) for r in c12),
         c14a='%.1f' % next(r['chi'] for r in A['chains'] if r['L'] == 14 and abs(r['k'] - 0.2475) < 1e-9),
         c14b='%.1f' % next(r['chi'] for r in A['chains'] if r['L'] == 14 and abs(r['k'] - 0.2525) < 1e-9),
         cr_lo='%.4f' % min([M['per_L'][L]['crossing'] for L in Ls] + [cr16]), cr_hi='%.4f' % max([M['per_L'][L]['crossing'] for L in Ls] + [cr16]), cr12='%.4f' % M['per_L']['12']['crossing'], cr16='%.4f' % cr16,
         R16a='%+.4f(%d)' % (c16[0.248]['R'], round(1e4 * c16[0.248]['R_se'])), R16b='%+.4f(%d)' % (c16[0.252]['R'], round(1e4 * c16[0.252]['R_se'])),
         c16a='%.1f(%d)' % (c16[0.248]['chi'], round(10 * c16[0.248]['chi_se'])), c16b='%.1f(%d)' % (c16[0.252]['chi'], round(10 * c16[0.252]['chi_se'])), tau16='%d and %d' % (round(c16[0.248]['tau']), round(c16[0.252]['tau'])),
         bias16='%d%% and %d%%' % (round(100 * c16[0.248]['chi_low_bias']), round(100 * c16[0.252]['chi_low_bias'])), bias16max='%d%%' % round(100 * max(c16[0.248]['chi_low_bias'], c16[0.252]['chi_low_bias'])),
         r16a='%.2f' % (c16[0.248]['chi'] / M['per_L']['4']['fixed_k']['0.2480']['chi']), r16b='%.2f' % (c16[0.252]['chi'] / M['per_L']['4']['fixed_k']['0.2520']['chi']),
         d12_list=', '.join('%+.4f' % (rr['A1_L12_k%.3f_r1' % k]['R'] - rr['A1_L12_k%.3f_r2' % k]['R']) for k in (0.244, 0.246, 0.248, 0.250, 0.252, 0.254)),
         tau_list=', '.join('%d' % round(float(np.median(tau[int(L)]))) for L in Lt),
         rough2_250='R by sixths of its 3,000 measured sweeps ' + ', '.join('%+.3f' % m[i, 4].mean() for i in b) + ', last half %+.4f against %+.4f for the flat start' % (m[len(m) // 2:, 4].mean(), f[:, 4].mean()),
         au46=', '.join('%.2f' % au[(L, 0.46)]['chi'] for L in (4, 6, 8)), au_pow46='%+.2f' % U['chi_power']['0.46']['power'] + ' over 3⁴ to 8⁴',
         au_cr4='%.3f' % U['crossing']['4'], au_cr6='%.3f' % U['crossing']['6'], au_cr8='%.3f' % (0.46 - x8['R'] / x8['slope']), au_R8='%+.4f(%d)' % (x8['R'], round(1e4 * x8['R_se'])),
         a1_status=sys.argv[1] if len(sys.argv) > 1 else 'status at T29 in NIGHT_LOG.md')
s = open(N + 'SUPPLEMENT_README_T30.tmpl.md').read()
for k, v in V.items(): s = s.replace('<<' + k + '>>', v)
assert '<<' not in s, s[s.index('<<'):s.index('<<') + 40]
open(N + 'SUPPLEMENT_README.md', 'w').write(s)
for k, v in V.items(): print(f'{k}: {v}')
