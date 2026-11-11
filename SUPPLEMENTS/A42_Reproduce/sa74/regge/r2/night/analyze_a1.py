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
for line in open(MAC / 'queue3.txt'):   # T31: plan A5 (20^4, Mac runner 3, multi-core engine; same columns as queue2.txt)
    if line.startswith('#') or not line.strip(): continue
    f = line.split(); plan[f[0]] = dict(L=int(f[1]), k=float(f[2]), therm=int(f[3]), meas=int(f[4]))
plan['Q_L14_k0.2475'] = dict(L=14, k=0.2475, therm=800, meas=1600); plan['Q_L14_k0.2525'] = dict(L=14, k=0.2525, therm=800, meas=1600)
# T31: plan A3's forks share their parent's history up to the fork, so a fork and its parent count as one lineage
LINEAGE = {'A3_L16_k0.244_f1': 'A2_L16_k0.244_r1', 'A3_L16_k0.244_f2': 'A2_L16_k0.244_r2'}
rows = []
A1ONLY = [x for x in os.environ.get('A1ONLY', '').split(',') if x]   # T30: rebuild an earlier state (only these plan A1 chains)
for tag, p in plan.items():
    if A1ONLY and tag.startswith('A1_') and tag not in A1ONLY: continue
    if tag.startswith(('A2_a', 'A4_')): continue   # the a = 0.02 chains go to analyze_au.py
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
                 chi_low_bias=float(2 / ex), tau=float(len(x) / ess(x)), ess=float(ess(x)), maxq=float(m[:, 9].mean()),
                 R_sixths=[float(v.mean()) for v in np.array_split(x, 6)])
    rows.append(r)
MH = json.load(open(N / os.environ.get('MHFILE', 'multihist_T29.json')))
fixed = {}
for L, v in MH['per_L'].items():
    for k, f in v['fixed_k'].items(): fixed.setdefault(k, {})[int(L)] = dict(chi=f['chi'], se=f['boot_sd'], source='reweighted')
# T31: on 14^4 and larger, the finished chains at one coupling are pooled: the fluctuation from all their sweeps about the common mean
# (so differences between chains count), its error the larger of the combined single-chain errors and the spread between chains
groups = {}
for r in rows:
    if r['L'] >= 14 and r.get('chi') is not None and r['finished']: groups.setdefault((r['L'], '%.4f' % r['k']), []).append(r)
for (L, k), rs in groups.items():
    xs, Vs, lin = [], [], []
    for r in rs:
        a = np.loadtxt(MAC / f"{r['chain']}.dat", comments='#'); m = a[a[:, 1] == 1]; xs.append(m[:, 3]); Vs.append(m[:, 2]); lin.append(LINEAGE.get(r['chain'], r['chain']))
    x = np.concatenate(xs); chi = float(4 * ((x - x.mean()) ** 2).mean() / np.concatenate(Vs).mean())
    # T31: lineages (a chain with its forks) are the independent units; a lineage's error is its parent chain's (a fork adds sweeps
    # but is not counted as reducing the error), combined over lineages with weights by sweeps; then the spread between lineages
    gids = sorted(set(lin)); par = {r['chain']: r for r in rs}
    ng = np.array([sum(len(v) for v, l in zip(xs, lin) if l == g) for g in gids], float)
    se = float(np.sqrt(((ng / ng.sum()) ** 2 * np.array([par[g]['chi_se'] for g in gids]) ** 2).sum()))
    chig = []
    for g in gids:
        xg = np.concatenate([v for v, l in zip(xs, lin) if l == g]); Vg = np.concatenate([v for v, l in zip(Vs, lin) if l == g])
        chig.append(float(4 * ((xg - xg.mean()) ** 2).mean() / Vg.mean()))
    if len(gids) >= 2: se = max(se, float(np.std(chig, ddof=1) / np.sqrt(len(gids))))
    fixed.setdefault(k, {})[L] = dict(chi=chi, se=se, lineages=len(gids), lineage_chi=chig,
                                      source=('chain ' if len(rs) == 1 else 'chains ') + ', '.join(r['chain'] for r in rs))
# T31: start check on 16^4: at k = 0.248 and 0.252 the plan A1 chain (flat start, 3,000 sweeps at its coupling) against plan A2's second
# chain (started from the other coupling's final configuration, 2,000 sweeps at its coupling); at 0.244 the two A2 lineages
byc = {r['chain']: r for r in rows if r.get('R') is not None}
start_check = {}
for a_, b_ in (('A1_L16_k0.248', 'A2_L16_k0.248_r2'), ('A1_L16_k0.252', 'A2_L16_k0.252_r2'), ('A2_L16_k0.244_r1', 'A2_L16_k0.244_r2'),
               ('A2_L16_k0.244_r1', 'A3_L16_k0.244_f1'), ('A2_L16_k0.244_r2', 'A3_L16_k0.244_f2')):
    if a_ in byc and b_ in byc:
        ra, rb = byc[a_], byc[b_]; d = ra['R'] - rb['R']; s = float(np.hypot(ra['R_se'], rb['R_se']))
        dq = ra['meanq_halves'][0] / 2 + ra['meanq_halves'][1] / 2 - rb['meanq_halves'][0] / 2 - rb['meanq_halves'][1] / 2
        start_check[a_ + ' - ' + b_] = dict(dR=float(d), z=float(d / s), dmeanq=float(dq), chi=[ra['chi'], rb['chi']], chi_se=[ra['chi_se'], rb['chi_se']])
def fit_power(fx):
    out = {}
    for k, d in sorted(fx.items()):
        Ls = sorted(d)
        if len(Ls) < 3: continue
        x = np.log(Ls); y = np.log([d[L]['chi'] for L in Ls]); w = np.array([(d[L]['chi'] / max(d[L]['se'], 1e-9)) ** 2 for L in Ls])
        A = np.vstack([x, np.ones_like(x)]).T; cov = np.linalg.inv(A.T @ (w[:, None] * A)); c = cov @ A.T @ (w * y)
        chi2 = float(((y - A @ c) ** 2 * w).sum()); dof = len(x) - 2
        out[k] = dict(L=Ls, chi=[d[L]['chi'] for L in Ls], power=float(c[0]), se=float(np.sqrt(cov[0, 0] * max(1.0, chi2 / max(dof, 1)))), chi2=chi2, dof=dof)
    return out
power = fit_power(fixed)
# T31: sensitivity to the low bias of a variance from few independent sweeps: values on 14^4 and larger divided by (1 - 2/ESS), with ESS of
# the curvature sum added over the lineages' parent chains (forks not counted, which makes the correction larger)
fixed_lb = {k: {L: dict(v) for L, v in d.items()} for k, d in fixed.items()}
for (L, k), rs in groups.items():
    b = 2 / sum(2 / r['chi_low_bias'] for r in rs if r['chain'] not in LINEAGE)
    fixed_lb[k][L]['chi'] = fixed[k][L]['chi'] / (1 - b); fixed_lb[k][L]['low_bias'] = b
power_lb = fit_power(fixed_lb)
json.dump(dict(description=__doc__, chains=rows, fixed_k=fixed, fixed_k_power=power, fixed_k_power_lowbias=power_lb, start_check=start_check), open(N / ('a1_results%s.json' % os.environ.get('OUTSUFFIX', '')), 'w'), indent=1)
for r in rows:
    s = f"{r['chain']:20s} th {r['therm_done']:5d}/{r['therm_planned']} me {r['meas_done']:5d}/{r['meas_planned']}{' done' if r['finished'] else '     '}"
    if 'R' in r: s += f" R {r['R']:+.4f}({1e4 * r['R_se']:.0f}) halves {r['R_halves'][0]:+.4f} {r['R_halves'][1]:+.4f} q {r['meanq_halves'][0]:.3f} {r['meanq_halves'][1]:.3f} chi {r['chi']:.2f}({100 * r['chi_se']:.0f}) tau {r['tau']:.0f} maxq {r['maxq']:.0f}"
    print(s)
for k, p in power.items(): print(f"k={k}: chi " + ', '.join(f"{L}^4 {c:.2f}" for L, c in zip(p['L'], p['chi'])) + f"; power {p['power']:+.2f}({100 * p['se']:.0f}), chi2 {p['chi2']:.1f}/{p['dof']}")
