"""SA74 night figures: (1) curvature fluctuation chi against grid size L at fixed k (log-log), with the growth a critical point
would force, chi ~ L^(2/nu - 4) (Hamber 2015): slope 0.94 for nu = 0.405 (delta = 0.62, a = 0.005) and slope 2 for nu = 1/3;
(2) average curvature R(k) for each grid size near the zero crossing. Reads fss/fss_results.json (run fss/analyze_fss.py first)."""
import json, numpy as np, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R2 = '/home/claude/sa74/regge/r2'
S = json.load(open(f'{R2}/fss/fss_results.json'))['summary']
INK, INK2 = '#0b0b0b', '#52514e'
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']   # reference palette, fixed order
MK = ['o', 's', 'D', '^', 'v', 'P', 'X', '*']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.6})
KS = [float(x) for x in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['0.245', '0.2475', '0.25', '0.2525'])]
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=200)
for i, k in enumerate(KS):
    pts = sorted((s['L'], s['chi'], s['chi_se']) for s in S if abs(s['k'] - k) < 1e-6)
    if not pts: continue
    L = np.array([p[0] for p in pts]); c = np.array([p[1] for p in pts]); e = np.array([p[2] for p in pts])
    ax.errorbar(L * (1 + 0.012 * (i - 1.5)), c, yerr=e, color=SER[i], marker=MK[i], ms=6, lw=1.4, capsize=2, label=f'k = {k:g}')
Lr = np.array([3, 16.5])
for p, txt in ((0.94, 'critical growth, ν = 0.405 (a = 0.005)'), (2.0, 'critical growth, ν = 1/3 (a = 0)')):
    ax.plot(Lr, 4.5 * (Lr / 3) ** p, color=INK2, lw=1.0, ls='--' if p < 1 else ':', label=txt)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.5, 60)
ax.set_xticks([3, 4, 5, 6, 8, 10, 12, 14, 16]); ax.set_xticklabels(['3', '4', '5', '6', '8', '10', '12', '14', '16'])
ax.set_xlabel('grid size L (L⁴ points)'); ax.set_ylabel('curvature fluctuation χ')
ax.set_title('Curvature fluctuation against grid size, a = 0.005', fontsize=11.5, color=INK)
ax.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig(f'{R2}/night/chi_vs_L.png'); print('wrote chi_vs_L.png')
# (2) R(k) by grid size near the crossing
fig, ax = plt.subplots(figsize=(6.4, 4.6), dpi=200)
Ls = sorted(set(s['L'] for s in S if s['L'] >= 4))
for i, Lv in enumerate(Ls):
    pts = sorted((s['k'], s['R'], s['R_se']) for s in S if s['L'] == Lv and 0.235 <= s['k'] <= 0.265)
    if len(pts) < 2: continue
    ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts], color=SER[i % 8], marker=MK[i % 8], ms=5, lw=1.2, capsize=2, label=f'{Lv}⁴')
ax.axhline(0, color=INK2, lw=0.8); ax.axvline(0.2443, color=INK2, lw=1.0, ls=':')
ax.text(0.2446, ax.get_ylim()[0] * 0.95, 'Hamber k$_c$ = 0.2443', fontsize=8.5, color=INK2, rotation=90, va='bottom')
ax.set_xlabel('k (inverse Newton coupling, lattice units)'); ax.set_ylabel('average curvature R')
ax.set_title('Average curvature near the zero crossing, a = 0.005', fontsize=11.5, color=INK)
ax.legend(frameon=False, fontsize=9, ncol=2, loc='upper left', title='grid', title_fontsize=9)
fig.tight_layout(); fig.savefig(f'{R2}/night/R_vs_k.png'); print('wrote R_vs_k.png')
