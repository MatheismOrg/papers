"""SA74 plan A1 (Mac runner 2, sequential edge order, fixed 03:27 PST before its runs) and the 14^4 chains of plan Q1: the
curvature-squared ensemble at a = 0.005 on 12^4, 14^4 and 16^4. For each chain: the sweeps planned and done, the average
curvature R of the measurement sweeps with a batch-means error, R and the mean squared edge in the two halves of the
measurement (the drift check), the curvature fluctuation chi = 4 var(sum dA)/<sum V> with a five-block error, tau = n/ESS of R,
and the longest squared edge. Then chi at fixed couplings across all grids: 4^4 to 12^4 from the multiple-histogram curves
(multihist_T29.json, plan A1 chains included, with bootstrap errors), 14^4 and 16^4 from their chains at the couplings run, and
the fitted power of chi in L at each fixed coupling. Writes a1_results.json. T30: a single chain's chi error is the larger of the five-block
error and the error from the effective sample size of the squared deviations (CHISE=blocks gives T29's five-block error); chi_low_bias
estimates the low bias of a single chain's variance, 2/ESS of the curvature sum."""
import json, os, re, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import batches, ess
MAC = Path('/home/claude/sa74/NIGHT/mac_pull'); N = Path('/home/claude/sa74/regge/r2/night')
plan = {}
for line in open(MAC / 'queue2.txt'):
    if line.startswith('#') or not line.strip(): continue
    f = line.split(); plan[f[0]] = dict(L=int(f[1]), k=float(f[2]), therm=int(f[3]), meas=int(f[4]))
plan['Q_L14_k0.2475'] = dict(L=14, k=0.2475, therm=800, meas=1600); plan['Q_L14_k0.2525'] = dict(L=14, k=0.2525, therm=800, meas=1600)
rows = []
A1ONLY = [x for x in os.environ.get('A1ONLY', '').split(',') if x]   # T30: rebuild an earlier state (only these plan A1 chains)
for tag, p in plan.items():
    if A1ONLY and tag.startswith('A1_') and tag not in A1ONLY: continue
    if tag.startswith('A2_'): continue   # plan A2 has its own analysis
    fn = MAC / f'{tag}.dat'
    if not fn.exists(): continue
    a = np.loadtxt(fn, comments='#')
    if a.ndim != 2: continue
    th = a[a[:, 1] == 0]; m = a[a[:, 1] == 1]
    r = dict(chain=tag, L=p['L'], k=p['k'], therm_planned=p['therm'], meas_planned=p['meas'], therm_done=int(len(th)), meas_done=int(len(m)),
             finished=bool(len(m) >= p['meas']))
    if len(m) >= 400:
        x = m[:, 4]; h = len(m) // 2; bl = np.array_split(np.arange(len(m)), 5)
        cb = [4 * m[b, 3].var() / m[b, 2].mean() for b in bl]
        y = (m[:, 3] - m[:, 3].mean()) ** 2; ey = ess(y); ex = ess(m[:, 3])   # T30: error of chi from the effective sample size of the squared deviations
        se_ess = 4 / m[:, 2].mean() * y.std(ddof=1) / np.sqrt(ey); se_blk = np.std(cb, ddof=1) / np.sqrt(5)
        r.update(R=float(x.mean()), R_se=float(max(b['se'] for b in batches(x))), R_halves=[float(x[:h].mean()), float(x[h:].mean())],
                 meanq_halves=[float(m[:h, 5].mean()), float(m[h:, 5].mean())], chi=float(4 * m[:, 3].var() / m[:, 2].mean()),
                 chi_se=float(max(se_blk, se_ess) if os.environ.get('CHISE', 'max') == 'max' else se_blk), chi_se_blocks=float(se_blk), chi_se_ess=float(se_ess),
                 chi_low_bias=float(2 / ex), tau=float(len(x) / ess(x)), ess=float(ess(x)), maxq=float(m[:, 9].mean()))
    rows.append(r)
MH = json.load(open(N / os.environ.get('MHFILE', 'multihist_T29.json')))
fixed = {}
for L, v in MH['per_L'].items():
    for k, f in v['fixed_k'].items(): fixed.setdefault(k, {})[int(L)] = dict(chi=f['chi'], se=f['boot_sd'], source='reweighted')
for r in rows:
    if r['L'] >= 14 and r.get('chi') is not None and r['finished']:
        fixed.setdefault('%.4f' % r['k'], {})[r['L']] = dict(chi=r['chi'], se=r['chi_se'], source='chain ' + r['chain'])
power = {}
for k, d in sorted(fixed.items()):
    Ls = sorted(d)
    if len(Ls) < 3: continue
    x = np.log(Ls); y = np.log([d[L]['chi'] for L in Ls]); w = np.array([(d[L]['chi'] / max(d[L]['se'], 1e-9)) ** 2 for L in Ls])
    A = np.vstack([x, np.ones_like(x)]).T; cov = np.linalg.inv(A.T @ (w[:, None] * A)); c = cov @ A.T @ (w * y)
    chi2 = float(((y - A @ c) ** 2 * w).sum()); dof = len(x) - 2
    power[k] = dict(L=Ls, chi=[d[L]['chi'] for L in Ls], power=float(c[0]), se=float(np.sqrt(cov[0, 0] * max(1.0, chi2 / max(dof, 1)))), chi2=chi2, dof=dof)
json.dump(dict(description=__doc__, chains=rows, fixed_k=fixed, fixed_k_power=power), open(N / ('a1_results%s.json' % os.environ.get('OUTSUFFIX', '')), 'w'), indent=1)
for r in rows:
    s = f"{r['chain']:20s} th {r['therm_done']:5d}/{r['therm_planned']} me {r['meas_done']:5d}/{r['meas_planned']}{' done' if r['finished'] else '     '}"
    if 'R' in r: s += f" R {r['R']:+.4f}({1e4 * r['R_se']:.0f}) halves {r['R_halves'][0]:+.4f} {r['R_halves'][1]:+.4f} q {r['meanq_halves'][0]:.3f} {r['meanq_halves'][1]:.3f} chi {r['chi']:.2f}({100 * r['chi_se']:.0f}) tau {r['tau']:.0f} maxq {r['maxq']:.0f}"
    print(s)
for k, p in power.items(): print(f"k={k}: chi " + ', '.join(f"{L}^4 {c:.2f}" for L, c in zip(p['L'], p['chi'])) + f"; power {p['power']:+.2f}({100 * p['se']:.0f}), chi2 {p['chi2']:.1f}/{p['dof']}")
