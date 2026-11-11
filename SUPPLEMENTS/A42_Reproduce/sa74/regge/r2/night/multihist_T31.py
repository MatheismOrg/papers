"""SA74 night: multiple-histogram reweighting in k (Ferrenberg and Swendsen, 1989) of all a = 0.005 flat-start chains on
each grid L <= LMAX. Each chain j at coupling k_j samples exp(k_j x) times a k-independent weight, with x = sum dA (column 3).
Chains are combined with their statistical inefficiencies g_j = n_j / ESS_j (ESS of x), solving
   exp(f_l) = sum_s exp(k_l x_s) / g_j(s) / sum_m (n_m / g_m) exp(k_m x_s - f_m)
by iteration; then any coupling k between the chains gets weights w_s(k) proportional to the same expression with k in place
of k_l. Reports chi(k) = 4 var_k(x) / <sum V>_k and R(k), their largest values, and the fitted power of the largest chi in L, on grids from
LMIN = 4 (Hamber 1993 left out smaller grids, since the curvature-squared term couples next-nearest neighbors) to LMAX = 10.
SKIPFRAC=0.5 repeats it on the second halves of the chains (multihist_skip50.json).
Every THIN-th sweep is used (default 10). T29: A1INC=1 adds the plan A1 chains (runner 2) and drops the earlier 12^4 and larger chains; KPAD > 0 takes the largest chi only within KPAD of the couplings run; chi at the couplings KFIX is reported with bootstrap errors. Errors: the largest chi is recomputed on NBOOT (default 40) bootstrap resamples (chains drawn with replacement within each k,
then five blocks of each chain drawn with replacement)."""
import json, os, sys, numpy as np
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import ess
src = open('/home/claude/sa74/regge/r2/night/reweight.py').read().split('grid = ')[0]
exec(src)   # builds `chains` (flat starts, a = 0.005, measurement sweeps, SKIPFRAC applied)
LMAX = int(os.environ.get('LMAX', '10')); LMIN = int(os.environ.get('LMIN', '4')); KLO, KHI = 0.2415, 0.258; THIN = int(os.environ.get('THIN', '10')); NBOOT = int(os.environ.get('NBOOT', '40'))
A1INC = int(os.environ.get('A1INC', '0')); KPAD = float(os.environ.get('KPAD', '0'))
if A1INC:   # T29: plan A1 chains (runner 2, sequential edge order; same target distribution) replace the earlier, drifting 12^4 and larger chains
    import re as _re
    chains = [c for c in chains if not (c['name'].startswith('mac_') and c['L'] >= 12)]
    A1PLAN = {l.split()[0]: int(l.split()[4]) for l in open(MAC / 'queue2.txt') if l.strip() and not l.startswith('#')}   # measurement sweeps planned per chain
    A1ONLY = [x for x in os.environ.get('A1ONLY', '').split(',') if x]   # T30: rebuild an earlier state (only these plan A1 chains)
    for fn in sorted(MAC.glob('A1_L*_k*.dat')):
        if A1ONLY and fn.stem not in A1ONLY: continue
        mm = _re.match(r'A1_L(\d+)_k([0-9.]+)(_r\d+)?$', fn.stem)
        if not mm: continue
        a = np.loadtxt(fn, comments='#'); a = a[a[:, 1] == 1]
        if len(a) < max(MINROWS, A1PLAN.get(fn.stem, 10 ** 9)): continue   # finished chains only
        a = a[int(SKIP * len(a)):]
        chains.append(dict(L=int(mm.group(1)), k=round(float(mm.group(2)), 4), name='mac_' + fn.stem, V=a[:, 2], DA=a[:, 3], R=a[:, 4]))
for c in chains:   # keep every THIN-th sweep (autocorrelation times are 90 to 170 sweeps), ESS from the thinned series
    for key in ('DA', 'V', 'R'): c[key] = c[key][::THIN]
grid = np.round(np.arange(KLO, KHI + 1e-9, 0.00025), 5)
def lse(a, axis=None):
    m = np.max(a, axis=axis, keepdims=True); return (m + np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True))).squeeze()
def solve(cs):
    x = np.concatenate([c['DA'] for c in cs]); V = np.concatenate([c['V'] for c in cs]); R = np.concatenate([c['R'] for c in cs])
    g = np.concatenate([np.full(len(c['DA']), len(c['DA']) / max(ess(c['DA']), 1.0)) for c in cs])
    kj = np.array([c['k'] for c in cs]); Nj = np.array([len(c['DA']) / (len(c['DA']) / max(ess(c['DA']), 1.0)) for c in cs])
    xc = x - x.mean(); f = np.zeros(len(cs))
    for _ in range(800):
        logden = lse(np.log(Nj)[None, :] + np.outer(xc, kj) - f[None, :], axis=1)   # per sample
        fn = np.array([lse(k * xc - np.log(g) - logden) for k in kj]); fn -= fn[0]
        if np.max(np.abs(fn - f)) < 1e-8: f = fn; break
        f = fn
    logden = lse(np.log(Nj)[None, :] + np.outer(xc, kj) - f[None, :], axis=1)
    out = []
    for k in grid:
        lw = k * xc - np.log(g) - logden; w = np.exp(lw - lw.max()); w /= w.sum()
        mx = (w * x).sum(); vx = (w * (x - mx) ** 2).sum(); mV = (w * V).sum()
        out.append((float(k), float(4 * vx / mV), float((w * R).sum()), float(1.0 / (w ** 2).sum())))
    return out
res = {}
rng = np.random.default_rng(7)
for L in sorted(set(c['L'] for c in chains)):
    if L > LMAX or L < LMIN: continue
    cs = [c for c in chains if c['L'] == L and KLO - 0.004 <= c['k'] <= KHI + 0.004]
    if len(set(c['k'] for c in cs)) < 2: continue
    curve = solve(cs)
    klo, khi = (min(c['k'] for c in cs) - KPAD, max(c['k'] for c in cs) + KPAD) if KPAD > 0 else (KLO, KHI)   # KPAD > 0: largest chi only within KPAD of the couplings run
    inwin = lambda k: klo - 1e-9 <= k <= khi + 1e-9
    pk = max([r for r in curve if inwin(r[0])], key=lambda r: r[1])
    bys = {}
    for c in cs: bys.setdefault(c['k'], []).append(c)
    KFIX = [float(v) for v in os.environ.get('KFIX', '0.244,0.2475,0.248,0.252,0.2525').split(',')]; ifix = [int(np.argmin(np.abs(grid - kf))) for kf in KFIX]
    boot = []; bootfix = []
    for b in range(NBOOT):
        bs = []
        for cl in bys.values():
            for i in rng.integers(0, len(cl), len(cl)):   # chains with replacement within each k, then 5 blocks of each with replacement
                c0 = cl[i]; blk = np.array_split(np.arange(len(c0['DA'])), 5); pick = np.concatenate([blk[j] for j in rng.integers(0, 5, 5)])
                bs.append(dict(k=c0['k'], DA=c0['DA'][pick], V=c0['V'][pick], R=c0['R'][pick]))
        cb = solve(bs); boot.append(max(r[1] for r in cb if inwin(r[0]))); bootfix.append([cb[i][1] for i in ifix])
    cross = None
    for (k1, _, r1, _), (k2, _, r2, _) in zip(curve, curve[1:]):
        if r1 < 0 <= r2: cross = k1 + (k2 - k1) * (-r1) / (r2 - r1); break
    fixed = {('%.4f' % KFIX[j]): dict(chi=curve[i][1], boot_sd=float(np.std([b[j] for b in bootfix])), R=curve[i][2], weight_ess=curve[i][3]) for j, i in enumerate(ifix)}
    res[L] = dict(chains=len(cs), ks=sorted(bys), window=[klo, khi], fixed_k=fixed, peak_k=pk[0], peak_chi=pk[1], peak_chi_boot_sd=float(np.std(boot)), crossing=cross,
                  min_eff_samples=float(min(r[3] for r in curve)), curve=curve)
    print(f"L={L:2d} chains={len(cs):2d} largest chi {pk[1]:.2f} +- {np.std(boot):.2f} at k={pk[0]:.5f}; R crosses 0 at {cross if cross is None else round(cross, 5)}; min weight-ESS {min(r[3] for r in curve):.0f}")
Ls = np.array(sorted(res)); P = np.array([res[L]['peak_chi'] for L in Ls]); E = np.array([res[L]['peak_chi_boot_sd'] for L in Ls])
x = np.log(Ls); y = np.log(P); w = (P / E) ** 2; A = np.vstack([x, np.ones_like(x)]).T
cov = np.linalg.inv(A.T @ (w[:, None] * A)); c = cov @ A.T @ (w * y); chi2 = float(((y - A @ c) ** 2 * w).sum()); dof = len(x) - 2
se = float(np.sqrt(cov[0, 0] * max(1.0, chi2 / max(dof, 1))))
print(f"largest chi ~ L^{c[0]:+.2f} +- {se:.2f} (chi2 {chi2:.1f} / {dof}) over L = {list(Ls)}")
fn = '/home/claude/sa74/regge/r2/night/multihist%s%s.json' % ('' if SKIP == 0 else '_skip%02d' % int(100 * SKIP), os.environ.get('SUFFIX', ''))
json.dump(dict(description=__doc__, k_range=[KLO, KHI], L_range=[LMIN, LMAX], skipfrac=SKIP, per_L={str(L): v for L, v in res.items()},
               peak_power=dict(power=float(c[0]), se=se, chi2=chi2, dof=dof, L=[int(v) for v in Ls])),
          open(fn, 'w'), indent=1)
