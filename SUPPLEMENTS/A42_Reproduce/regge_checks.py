"""The numbers Section 8.10 of A42 quotes for its larger grids, each with its standard error, computed from the result files the analysis
scripts write (a1_results_T37.json, multihist_T30.json, multihist_skip50_T30.json, au_results.json, fss_results.json).

The same function is run twice. make_paper_numbers.py ran it on the result files of the original chains and wrote PAPER_NUMBERS.json;
verify.py --regge runs it on the result files of the regenerated chains and compares the two, number by number.

A regenerated chain is a new sample of the same Markov chain, so its numbers differ from the paper's by their statistical errors.
Every comparison is therefore a z value, z = (regenerated - paper) / sqrt(se_regenerated^2 + se_paper^2), and a check passes when
|z| <= 3 for a single number. A check that compares m numbers at once passes when the largest |z| stays below zmax(m), the bound
with the same chance of a false alarm (0.27%) for the whole check. No check uses a hand-set range.

Where the errors come from:
  fluctuation at a fixed coupling   4^4 to 12^4: bootstrap of the reweighted curve (multihist.py); 16^4 to 32^4: analyze_a1.py (lineages)
  fitted powers                     weighted least squares in log L, the error scaled by sqrt(chi2/dof) when that exceeds 1
  zero crossing of the curvature    4^4 to 12^4: bootstrap of the reweighted curve; 16^4, 20^4 and a = 0.02: the batch-means errors of the
                                    two bracketing curvatures carried through the linear interpolation
  autocorrelation time per grid     the median over chains near the crossing; its error in the logarithm is 1.2533 sd(ln tau)/sqrt(n),
                                    the scatter between chains (two independent 12^4 chains at one coupling differ by a factor 1.6 on average)
  a = 0.02 fluctuation              the larger of the five-block and effective-sample-size errors (analyze_au.py)"""
import json, math
from pathlib import Path
from statistics import NormalDist, median

KS = ('0.2440', '0.2480', '0.2520')


def zmax(m):
    """Largest |z| allowed among m numbers so that the whole check has the false-alarm chance of one number at 3 standard errors."""
    alpha = 2 * (1 - NormalDist().cdf(3.0))
    return NormalDist().inv_cdf(1 - alpha / (2 * m))


def interp_cross(k1, r1, s1, k2, r2, s2):
    """Zero of the line through (k1, r1) and (k2, r2), with the error the two curvature errors give it."""
    c = k1 + (k2 - k1) * (-r1) / (r2 - r1)
    return c, abs(k2 - k1) / (r2 - r1) ** 2 * math.hypot(r2 * s1, r1 * s2)


def mean_se(rows):
    """Average curvature of the finished chains at one coupling and its error (errors of independent chains combined)."""
    return sum(r['R'] for r in rows) / len(rows), math.sqrt(sum(r['R_se'] ** 2 for r in rows)) / len(rows)


def cross_from_chains(rows):
    """Zero crossing on one grid from chain averages per coupling: interpolation between the bracketing couplings, or a short
    extrapolation from the last two when every average is still negative (20^4)."""
    by = {}
    for r in rows: by.setdefault(r['k'], []).append(r)
    q = sorted((k,) + mean_se(v) for k, v in by.items())
    for (k1, r1, s1), (k2, r2, s2) in zip(q, q[1:]):
        if r1 < 0 <= r2: return interp_cross(k1, r1, s1, k2, r2, s2)
    if len(q) >= 2:
        (k1, r1, s1), (k2, r2, s2) = q[-2], q[-1]
        return interp_cross(k1, r1, s1, k2, r2, s2)
    return None


def med_ln(v):
    lv = [math.log(x) for x in v]; n = len(lv); m = sum(lv) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in lv) / (n - 1)) if n > 1 else float('nan')
    return dict(median=float(median(v)), n=n, se_ln=1.2533 * sd / math.sqrt(n))


def numbers(night, fss):
    night, fss = Path(night), Path(fss)
    A = json.load(open(night / 'a1_results_T37.json')); M = json.load(open(night / 'multihist_T30.json'))
    Mb = json.load(open(night / 'multihist_skip50_T30.json')); U = json.load(open(night / 'au_results.json')); F = json.load(open(fss / 'fss_results.json'))
    out = {}
    # the fluctuation at three fixed couplings on every grid, and its fitted power in L (with and without the low-bias correction)
    out['fixed_k'] = {k: {str(L): dict(v=A['fixed_k'][k][str(L)]['chi'], se=A['fixed_k'][k][str(L)]['se']) for L in A['fixed_k_power'][k]['L']} for k in KS}
    out['power'] = {k: dict(v=A['fixed_k_power'][k]['power'], se=A['fixed_k_power'][k]['se'], L=A['fixed_k_power'][k]['L']) for k in KS}
    out['power_lowbias'] = {k: dict(v=A['fixed_k_power_lowbias'][k]['power'], se=A['fixed_k_power_lowbias'][k]['se']) for k in KS}
    h = lambda L: (0.16 + 1.11 / L) * L ** 0.64   # Hamber's (1993) fitted curve at k = 0.244
    out['hamber_power'] = math.log(h(32) / h(4)) / math.log(8)   # T37: over the range now covered, 4^4 to 32^4
    # the largest fluctuation near the crossing on 4^4 to 12^4, and its power (all sweeps; second halves)
    out['peak'] = {L: dict(v=v['peak_chi'], se=v['peak_chi_boot_sd']) for L, v in M['per_L'].items()}
    out['peak_power'] = dict(v=M['peak_power']['power'], se=M['peak_power']['se'])
    out['peak_power_second_halves'] = dict(v=Mb['peak_power']['power'], se=Mb['peak_power']['se'])
    # zero crossing of the average curvature on every grid
    cr = {L: dict(v=v['crossing'], se=v.get('crossing_boot_sd')) for L, v in M['per_L'].items()}
    fin = [r for r in A['chains'] if r['finished'] and r.get('R') is not None]
    for L in (16, 20, 24):
        c = cross_from_chains([r for r in fin if r['L'] == L])
        if c: cr[str(L)] = dict(v=c[0], se=c[1])
    out['crossing'] = cr
    # 16^4: chains started from another coupling's final configuration against flat starts
    out['start_check_z'] = {k: A['start_check'][k]['z'] for k in ('A1_L16_k0.248 - A2_L16_k0.248_r2', 'A1_L16_k0.252 - A2_L16_k0.252_r2')}
    # autocorrelation times near the crossing: the median over chains on each grid
    tau = {}
    for key, rows in F['per_chain'].items():
        L, k = key.split('_')
        if 4 <= int(L) <= 10 and 0.2425 <= float(k) <= 0.2575:
            tau.setdefault(int(L), []).extend(r['n'] / max(r['ess'], 1e-9) for r in rows if r['n'] >= 800)
    near = [r for r in A['chains'] if r['finished'] and r.get('tau') and 0.2425 <= r['k'] <= 0.2575]
    tau[12] = [r['tau'] for r in near if r['L'] == 12 and r['chain'].startswith('A1_')]
    for L in (16, 20, 24, 32): tau[L] = [r['tau'] for r in near if r['L'] == L]
    out['tau'] = {str(L): med_ln(v) for L, v in sorted(tau.items()) if v}
    # T37: a grid with one chain (32^4) has no scatter of its own; it takes the scatter between chains pooled over the grids with three or more
    lv = [[math.log(x) for x in v] for v in tau.values() if len(v) >= 3]
    pooled = math.sqrt(sum(sum((x - sum(l) / len(l)) ** 2 for x in l) for l in lv) / sum(len(l) - 1 for l in lv))
    for d in out['tau'].values():
        if d['n'] == 1: d['se_ln'] = 1.2533 * pooled
    # a = 0.02: zero crossings and the fluctuation at k = 0.46
    au = {}
    for r in U['chains']: au.setdefault(r['L'], []).append(r)
    out['au_crossing'] = {}
    for L, v in U['crossing'].items():
        if v is None: continue
        c = cross_from_chains(au[int(L)])
        out['au_crossing'][L] = dict(v=v, se=c[1], v_check=c[0])
    out['au_chi_046'] = {str(r['L']): dict(v=r['chi'], se=r['chi_se']) for r in sorted(U['chains'], key=lambda r: r['L']) if abs(r['k'] - 0.46) < 1e-9 and r['L'] >= 4}
    return out


def zlist(new, paper, log=False):
    """z values of matching entries {key: {v, se}} (log: {median, se_ln}, compared in the logarithm)."""
    z = {}
    for key, p in paper.items():
        if key not in new: z[key] = float('inf'); continue
        n = new[key]
        d, s = ((math.log(n['median']) - math.log(p['median']), math.hypot(n['se_ln'], p['se_ln'])) if log
                else (n['v'] - p['v'], math.hypot(n['se'] or 0.0, p['se'] or 0.0)))
        z[key] = d / s if s > 0 and s == s else (0.0 if abs(d) < 1e-12 else float('inf'))   # no error on either side: the numbers must agree
    return z
