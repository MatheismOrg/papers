"""SA74 night: flat against rough starts at a = 0.005 near the zero crossing (hysteresis test).
Chains: Mac 8^4 (N2b/Q1 flat starts, their continuations, Q1/Q3 rough starts from the a = 0, k = 0.1 preparations and their
continuations) and cloud HY1 (4^4, 6^4). For each (L, k, start) the curvature R over the last half of all measurement sweeps
(continuations appended to their parent chain), with a batch-means error; then the flat minus rough gap. Writes
hysteresis.json and hysteresis_8.png (R against sweep on 8^4, 50-sweep running means)."""
import json, re, sys, glob, numpy as np
from pathlib import Path
sys.path.insert(0, '/home/claude/sa74/regge/r2/audit_a005'); from analyze_r2 import batches, ess
MAC = Path('/home/claude/sa74/NIGHT/mac_pull'); RUNS = Path('/home/claude/sa74/regge/r2/audit_a005/runs')
def rows(fn):
    a = np.loadtxt(fn, comments='#'); return a if a.ndim == 2 else None
series = {}   # (L, k, start, chain id) -> full R series (therm + meas), and measurement-only series
def add(key, a):
    t = a[a[:, 1] == 0][:, 4]; m = a[a[:, 1] == 1][:, 4]
    if key in series: series[key]['full'] = np.concatenate([series[key]['full'], a[:, 4]]); series[key]['meas'] = np.concatenate([series[key]['meas'], m])
    else: series[key] = dict(full=a[:, 4].copy(), meas=m.copy(), ntherm=len(t))
# Mac 8^4 and larger: parents first, then continuations
fl = sorted(MAC.glob('[NQ]_L*_k*.dat'))
for fn in fl:
    m = re.match(r'[NQ]_L(\d+)_k([0-9.]+)(_rough2?|_r\d)?$', fn.stem)   # _rN: replica flat starts (plans Q4, Q5)
    if not m: continue
    st = {None: 'flat', '_rough': 'rough', '_rough2': 'rough2'}.get(m.group(3), 'flat')
    a = rows(fn)
    if a is not None: add((int(m.group(1)), round(float(m.group(2)), 4), st, fn.stem), a)
for fn in fl:
    m = re.match(r'Q_L(\d+)_k([0-9.]+)(_rough)?_cont$', fn.stem)
    if not m: continue
    L, k, st = int(m.group(1)), round(float(m.group(2)), 4), 'rough' if m.group(3) else 'flat'
    parent = [key for key in series if key[:3] == (L, k, st) and not key[3].endswith('_cont')]
    a = rows(fn)
    if a is None: continue
    if parent: add(parent[0], a)
    elif (L, k) == (8, 0.26): add((8, 0.26, 'flat', 'FSS_L8_k0.26_flat'), a)   # continues the cloud chain
# cloud chains: the 8^4 k = 0.26 parent, and HY1
for fn in [RUNS / 'FSS_L8_k0.26_flat.dat'] + sorted(RUNS.glob('HY1_L*_k*_*.dat')):
    if not fn.exists() or '# done' not in fn.read_text()[-20000:]: continue
    m = re.match(r'(?:HY1_|FSS_)L(\d+)_k([0-9.]+)_(flat|rough)(?:_s(\d))?', fn.stem)
    key = (int(m.group(1)), round(float(m.group(2)), 4), m.group(3), fn.stem)
    a = rows(fn)
    if key in series:   # parent loaded after its continuation: put the parent first
        cont = series.pop(key); add(key, a); add(key, np.column_stack([np.zeros(len(cont['full'])), np.ones(len(cont['full'])), np.zeros(len(cont['full'])), np.zeros(len(cont['full'])), cont['full']]))
    else: add(key, a)
out = []
for (L, k, st, cid), s in sorted(series.items()):
    m = s['meas']
    if len(m) < 200: continue
    h = m[len(m) // 2:]
    out.append(dict(L=L, k=k, start=st, chain=cid, n_meas=len(m), R_all=float(m.mean()), R_last_half=float(h.mean()),
                    se_last_half=float(max(b['se'] for b in batches(h))) if len(h) >= 32 else float('nan'), ess_last_half=float(ess(h))))
gaps = []
for L, k in sorted(set((r['L'], r['k']) for r in out)):
    F = [r for r in out if (r['L'], r['k']) == (L, k) and r['start'] == 'flat']; Rg = [r for r in out if (r['L'], r['k']) == (L, k) and r['start'] != 'flat']
    if F and Rg:
        f = np.mean([r['R_last_half'] for r in F]); g = np.mean([r['R_last_half'] for r in Rg])
        ef = np.sqrt(sum(r['se_last_half'] ** 2 for r in F)) / len(F); eg = np.sqrt(sum(r['se_last_half'] ** 2 for r in Rg)) / len(Rg)
        gaps.append(dict(L=L, k=k, flat=float(f), rough=float(g), gap=float(g - f), gap_se=float(np.hypot(ef, eg)), n_flat=len(F), n_rough=len(Rg)))
json.dump(dict(description=__doc__, chains=out, gaps=gaps), open('/home/claude/sa74/regge/r2/night/hysteresis.json', 'w'), indent=1)
for r in out: print(f"{r['L']:>2} {r['k']:.4f} {r['start']:<6} {r['chain']:<26} n={r['n_meas']:5d} R(all)={r['R_all']:+.4f} R(last half)={r['R_last_half']:+.4f}±{r['se_last_half']:.4f} ess={r['ess_last_half']:.0f}")
for g in gaps: print(f"GAP L={g['L']} k={g['k']}: flat {g['flat']:+.4f} rough {g['rough']:+.4f} gap {g['gap']:+.4f}±{g['gap_se']:.4f} ({g['n_flat']}+{g['n_rough']} chains)")
# figure: 8^4 time series
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2 = '#0b0b0b', '#52514e'; COL = {'flat': '#2a78d6', 'rough': '#eb6834', 'rough2': '#1baf7a'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.6})
ks = sorted(set(k for (L, k, st, c) in series if L == 8 and any(kk == k and LL == 8 and stt != 'flat' for (LL, kk, stt, cc) in series)))
if ks:
    fig, axs = plt.subplots(len(ks), 1, figsize=(6.4, 2.3 * len(ks) + 0.6), dpi=200, sharex=True, squeeze=False)
    for ax, k in zip(axs[:, 0], ks):
        seen = set()
        for (L, kk, st, cid), s in sorted(series.items()):
            if L != 8 or kk != k: continue
            y = np.convolve(s['full'], np.ones(50) / 50, mode='valid')
            lab = {'flat': 'flat starts', 'rough': 'rough start', 'rough2': 'rough start, second preparation'}[st]
            ax.plot(np.arange(len(y)) + 25, y, color=COL[st], lw=1.2, label=None if lab in seen else lab); seen.add(lab)
        ax.axhline(0, color=INK2, lw=0.8); ax.set_ylim(-0.06, 0.12); ax.set_ylabel('R'); ax.set_title(f'8⁴, k = {k:g}', fontsize=10, color=INK, loc='left')
        ax.legend(frameon=False, fontsize=8.5, loc='upper right')
    axs[-1, 0].set_xlabel('sweep (50-sweep running mean)')
    fig.suptitle('Flat against rough starts, a = 0.005', fontsize=11.5, color=INK)
    fig.tight_layout(); fig.savefig('/home/claude/sa74/regge/r2/night/hysteresis_8.png'); print('wrote hysteresis_8.png')
