"""PROTOCOL_T65, Test 1: the dark-matter density at 150 pc of the sixteen dwarfs of Read, Walker and Steger (2019), measured
by them from stellar kinematics (8 dwarf spheroidals) and HI rotation curves (8 dwarf irregulars), against the frozen
formula's halos under five readings of the core (lib65.py: N, G, H1, H2, P). A19's own inputs for every galaxy; nothing is
fitted. Error of a reading = rms over galaxies of log10(rho_obs / rho_pred); paired bootstrap over galaxies (10,000 draws).
Usage: PYTHONPATH=../dmfit:../t59 python3 dens65.py ../t59 sfh65_inputs.json res65_alldata_fit.json"""
import sys, json, numpy as np
import lib65 as L
T59 = sys.argv[1]; sf = json.load(open(sys.argv[2]))['read2019']; fit = json.load(open(sys.argv[3])); TK = (fit['t'], fit['kappa'])
R150 = np.array([0.15]); KINDS = ('N', 'G', 'H1', 'H2', 'P'); rows = {}
for lab, v in sf.items():
    a, dsph = L.galaxy(T59, lab); tsf = L.TU - v['t_trunc']; obs = v['rho150'] * 1e8; sig = 0.5 * (np.log10((v['rho150'] + v['hi']) / v['rho150']) + np.log10(v['rho150'] / (v['rho150'] - v['lo'])))
    pred = {k: float(L.reading(a, R150, k, tsf=tsf, tk=TK)[1][0]) for k in ('N', 'G', 'H1', 'H2')}; mm, rho, flag, r1 = L.sidm(a, R150, dsph); pred['P'] = float(rho[0])
    a1 = {k: np.asarray(x)[:1] for k, x in a.items()}; h = L.halo(a1)
    rows[lab] = dict(type=v['type'], t_trunc=v['t_trunc'], t_SF=tsf, obs=obs, sig_dex=float(sig), pred=pred, err={k: float(np.log10(obs / pred[k])) for k in KINDS}, logMh=float(np.log10(h[0][0])), c=float(h[1][0]), Re=float(a1['re'][0]),
                     n_H1=float(L.n_read(h, tsf)[0]), w_gas=float(L.w_gas(a1)[0]), w_hist=float(L.w_hist(a1, tsf)[0]), sidm_flag=flag, sidm_r1=float(r1))
    print('%-9s %s t_SF %5.2f  obs %.2e (+-%.2f dex) | pred N %.2e G %.2e H1 %.2e H2 %.2e P %.2e | log obs/pred: %s | P %s r1 %.2f' % (lab, v['type'], tsf, obs, sig, *[pred[k] for k in KINDS], ' '.join('%s %+.2f' % (k, rows[lab]['err'][k]) for k in KINDS), flag, r1), flush=True)
rng = np.random.default_rng(651); out = dict(galaxies=rows, frozen=dict(t_gas=L.T_G, kappa_gas=L.K_G, t_hist=TK[0], kappa_hist=TK[1], kappa_read=L.KR, eta=L.ETA, s0=1.0))
for grp, sel in (('dSph', [k for k, v in rows.items() if v['type'] == 'dSph']), ('dIrr', [k for k, v in rows.items() if v['type'] == 'dIrr']), ('all', list(rows))):
    e = {k: np.array([rows[g]['err'][k] for g in sel]) for k in KINDS}; sg = np.array([rows[g]['sig_dex'] for g in sel]); idx = rng.integers(0, len(sel), (10000, len(sel)))
    o = dict(n=len(sel), rms={k: round(float(np.sqrt(np.mean(v ** 2))), 4) for k, v in e.items()}, mean={k: round(float(np.mean(v)), 4) for k, v in e.items()},
             chi2_obs_plus_0p3={k: round(float(np.sum(v ** 2 / (sg ** 2 + 0.3 ** 2))), 2) for k, v in e.items()})
    for k in ('N', 'H1', 'H2', 'P'):
        dd = np.sqrt((e[k][idx] ** 2).mean(1)) - np.sqrt((e['G'][idx] ** 2).mean(1))
        o[k + '_minus_G'] = dict(diff=round(float(np.sqrt(np.mean(e[k] ** 2)) - np.sqrt(np.mean(e['G'] ** 2))), 4), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_better=int(np.sum(np.abs(e[k]) < np.abs(e['G']))), n_worse=int(np.sum(np.abs(e[k]) > np.abs(e['G']))))
    dd = np.sqrt((e['P'][idx] ** 2).mean(1)) - np.sqrt((e['H1'][idx] ** 2).mean(1)); o['P_minus_H1'] = dict(diff=round(float(np.sqrt(np.mean(e['P'] ** 2)) - np.sqrt(np.mean(e['H1'] ** 2))), 4), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist())
    out[grp] = o; print(grp, json.dumps(o), flush=True)
json.dump(out, open('res65_density.json', 'w'), indent=1)
