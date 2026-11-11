"""PROTOCOL_T65, Test 2: the benchmark's own dwarf-spheroidal points (one acceleration per galaxy at 4/3 of the half-light
radius) for the 20 dwarfs with a quenching time in Weisz et al. (2015), under the readings of lib65.py. t_SF = 13.8 Gyr
minus 10^tau_90 years. Error of a reading = rms over galaxies of log10(g_obs / g_pred); paired bootstrap over galaxies.
Also the amplitude A of a history reading's predicted shift in the residuals of the formula: least squares through the
origin of (log g_obs - log g_G) on (log g_H - log g_G); A = 0 if gas-free dwarfs keep their cusps, 1 if the history reading
is right. P: T64's table at 1 cm^2/g (t64 res64_sidm_points.json).
Usage: PYTHONPATH=../dmfit:../t59 python3 points65.py sfh65_inputs.json res65_alldata_fit.json ../t64/res64_sidm_points.json"""
import sys, json, numpy as np
import lib65 as L
H = L.H; wz = json.load(open(sys.argv[1]))['weisz2015']; fit = json.load(open(sys.argv[2])); TK = (fit['t'], fit['kappa'])
sidm = {n: o for p, n, o in json.load(open(sys.argv[3])) if p == 'dwarfs'}
d = H.arr(H.dwarf_rows()); rows = {}; KINDS = ('G', 'H1', 'H2', 'P')
for wn, nm in L.WEISZ.items():
    k = d['name'] == nm; assert k.sum() == 1; a = {key: np.asarray(v)[k] for key, v in d.items()}; tsf = max(L.TU - 10 ** wz[wn]['tau90'] / 1e9, 0.05); r = a['r']
    g = {kk: float((a['gb'] + H.G * L.reading(a, r, kk, tsf=tsf, tk=TK)[0] / r ** 2 * H.ACC)[0]) for kk in ('G', 'H1', 'H2')}; go = float(a['go'][0])
    e = {kk: float(np.log10(go / g[kk])) for kk in g}; e['P'] = float(sidm[nm]['by_s0']['1.0']['er'][0]); a1 = a; h = L.halo(a1)
    rows[nm] = dict(weisz=wn, tau90=wz[wn]['tau90'], t_SF=tsf, err=e, shift_H1=float(np.log10(g['H1'] / g['G'])), shift_H2=float(np.log10(g['H2'] / g['G'])), n_H1=float(L.n_read(h, tsf)[0]), w_hist=float(L.w_hist(a1, tsf)[0]), logMh=float(np.log10(h[0][0])), ms=float(a['ms'][0]))
    print('%-18s tau90 %5.2f t_SF %5.2f | log obs/pred: %s | shift H1 %+.3f H2 %+.3f' % (nm, wz[wn]['tau90'], tsf, ' '.join('%s %+.2f' % (kk, e[kk]) for kk in KINDS), rows[nm]['shift_H1'], rows[nm]['shift_H2']), flush=True)
names = list(rows); rng = np.random.default_rng(652); idx = rng.integers(0, len(names), (10000, len(names)))
E = {kk: np.array([rows[n]['err'][kk] for n in names]) for kk in KINDS}; out = dict(galaxies=rows, n=len(names), rms={kk: round(float(np.sqrt(np.mean(v ** 2))), 4) for kk, v in E.items()}, mean={kk: round(float(np.mean(v)), 4) for kk, v in E.items()})
for kk in ('H1', 'H2', 'P'):
    dd = np.sqrt((E[kk][idx] ** 2).mean(1)) - np.sqrt((E['G'][idx] ** 2).mean(1))
    out[kk + '_minus_G'] = dict(diff=round(float(np.sqrt(np.mean(E[kk] ** 2)) - np.sqrt(np.mean(E['G'] ** 2))), 4), pct=round(float(100 * (np.sqrt(np.mean(E[kk] ** 2) / np.mean(E['G'] ** 2)) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(),
                                n_better=int(np.sum(np.abs(E[kk]) < np.abs(E['G']))), n_worse=int(np.sum(np.abs(E[kk]) > np.abs(E['G']))))
for kk in ('H1', 'H2'):
    s = np.array([rows[n]['shift_' + kk] for n in names]); A = float(np.sum(E['G'] * s) / np.sum(s ** 2)); Ab = np.sum(E['G'][idx] * s[idx], 1) / np.sum(s[idx] ** 2, 1)
    out['amplitude_' + kk] = dict(A=round(A, 3), boot95=np.round(np.percentile(Ab, [2.5, 97.5]), 3).tolist(), sum_shift_sq=round(float(np.sum(s ** 2)), 4), largest_shift=round(float(np.min(s)), 3))
ts = np.array([rows[n]['t_SF'] for n in names]); from scipy.stats import spearmanr
rho, p = spearmanr(ts, E['G']); perm = np.array([spearmanr(rng.permutation(ts), E['G'])[0] for _ in range(5000)])
out['spearman_residual_vs_tSF'] = dict(rho=round(float(rho), 3), p_two_sided_permutation=round(float(np.mean(np.abs(perm) >= abs(rho))), 4))
long_ = ts >= 6.9; out['mean_residual_G'] = dict(long_history=dict(n=int(long_.sum()), mean=round(float(E['G'][long_].mean()), 3), sem=round(float(E['G'][long_].std(ddof=1) / np.sqrt(long_.sum())), 3)), short_history=dict(n=int((~long_).sum()), mean=round(float(E['G'][~long_].mean()), 3), sem=round(float(E['G'][~long_].std(ddof=1) / np.sqrt((~long_).sum())), 3)))
print(json.dumps({k: v for k, v in out.items() if k != 'galaxies'}, indent=0)); json.dump(out, open('res65_points.json', 'w'), indent=1)
