"""SA74 T31 figure: the curvature fluctuation chi at fixed couplings against grid size, a = 0.005, at k = 0.244 (Hamber's k_c), 0.248, 0.252.
Grids 4^4 to 12^4 from the multiple-histogram curves (multihist_T30.json via a1_results_T31.json, bootstrap errors; filled markers, lines);
16^4 (two independent chains per coupling pooled, plus two forks at 0.244) and 20^4 (one chain per coupling) from chains (hollow markers);
14^4 single chains at 0.2475 and 0.2525 plotted with the nearer coupling. Markers are nudged sideways a little so the three couplings do not
overlap. Reference lines: critical growth L^0.94 (delta = 0.62) from the mean 4^4 value, and Hamber's (1993) fit (0.16 + 1.11/L) L^0.64 at
k = 0.244 scaled to the 4^4 value at that coupling. Writes chi_fixedk.png (T30 version: fig_fixedk_T30.py, chi_fixedk_T30.png)."""
import json, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
N = '/home/claude/sa74/regge/r2/night/'
A = json.load(open(N + os.environ.get('A1FILE', 'a1_results_T31.json'))); FK = A['fixed_k']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e1'
SER = {'0.2480': '#2a78d6', '0.2520': '#eb6834', '0.2440': '#1baf7a'}   # 0.248 and 0.252 keep their T29/T30 colors; 0.244 takes the next slot
MRK = {'0.2440': 's', '0.2480': 'o', '0.2520': '^'}; NUDGE = {'0.2440': 0.985, '0.2480': 1.0, '0.2520': 1.015}
NEAR = {'0.2440': ('0.2440',), '0.2480': ('0.2475', '0.2480'), '0.2520': ('0.2520', '0.2525')}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
for k in ('0.2440', '0.2480', '0.2520'):
    d = FK[k]; Ls = sorted(int(L) for L, v in d.items() if v['source'] == 'reweighted'); c = SER[k]; nx = NUDGE[k]
    ax.errorbar([L * nx for L in Ls], [d[str(L)]['chi'] for L in Ls], yerr=[d[str(L)]['se'] for L in Ls], fmt='-' + MRK[k], color=c, lw=2, ms=7,
                mec='white', mew=1.5, capsize=0, elinewidth=1.3, zorder=4)
    ch = sorted({int(L): v for kk in NEAR[k] for L, v in FK.get(kk, {}).items() if v['source'].startswith('chain')}.items())
    for L, v in ch:
        ax.errorbar([L * nx], [v['chi']], yerr=[v['se']], fmt=MRK[k], color=c, mfc='white', mew=1.8, ms=7, capsize=0, elinewidth=1.3, zorder=5)
x0 = np.mean([FK[k]['4']['chi'] for k in ('0.2440', '0.2480', '0.2520')]); Lr = np.linspace(4, 20, 60)
ax.plot(Lr, x0 * (Lr / 4) ** 0.94, color=INK2, lw=1.2, ls='--', zorder=2)
ax.text(19.5, x0 * (19.5 / 4) ** 0.94 * 1.05, 'critical point: L$^{0.94}$', color=INK2, fontsize=9, ha='right', va='bottom')
h = lambda L: (0.16 + 1.11 / L) * L ** 0.64; x4 = FK['0.2440']['4']['chi']
ax.plot(Lr, x4 * h(Lr) / h(4), color=INK2, lw=1.2, ls=':', zorder=2)
ax.text(19.5, x4 * h(19.5) / h(4) * 1.05, "Hamber's 1993 fit at k = 0.244, scaled", color=INK2, fontsize=9, ha='right', va='bottom')
ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
ticks = [4, 5, 6, 8, 10, 12, 14, 16, 20]
ax.set_xticks(ticks); ax.set_xticklabels([f'{t}⁴' for t in ticks]); ax.set_yticks([2, 3, 5, 10, 20]); ax.set_yticklabels(['2', '3', '5', '10', '20'])
ax.set_xlim(3.7, 23); ax.set_ylim(2, 26)
ax.set_xlabel('grid'); ax.set_ylabel('curvature fluctuation χ')
ax.set_title('Curvature fluctuation at fixed coupling against grid size, a = 0.005', fontsize=11, color=INK)
hd = [Line2D([], [], color=SER[k], lw=2, marker=MRK[k], ms=7, mec='white', mew=1.5) for k in ('0.2440', '0.2480', '0.2520')] + \
     [Line2D([], [], color=INK2, lw=0, marker='o', ms=7, mfc='white', mew=1.8)]
ax.legend(hd, ["k = 0.244 (Hamber's k$_c$)", 'k = 0.248', 'k = 0.252', 'chains: 14⁴ (k = 0.2475, 0.2525), 16⁴ pooled, 20⁴'], frameon=False, fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig(N + 'chi_fixedk.png'); print('chi_fixedk.png')
