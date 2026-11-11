"""SA74 T29 figure: the curvature fluctuation chi at fixed couplings against grid size, a = 0.005. Grids 4^4 to 12^4 from the
multiple-histogram curves (multihist_T29.json, plan A1 chains included, bootstrap errors); 14^4 and 16^4 from their finished
single chains at the couplings run (hollow markers, five-block errors), plotted with the nearer of the two couplings. Reference
lines from the 4^4 value: critical growth L^0.94 (delta = 0.62) and the shape of Hamber's (1993) fitted curve (0.16 + 1.11/L) L^0.64
at k = 0.244. Writes chi_fixedk.png (a1_results.json supplies the values)."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
N = '/home/claude/sa74/regge/r2/night/'
import os
A = json.load(open(N + os.environ.get('A1FILE', 'a1_results_T30.json'))); FK = A['fixed_k']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e1'
SER = ['#2a78d6', '#eb6834']   # reference palette, fixed order (as fig_reweight.py)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
fig, ax = plt.subplots(figsize=(6.8, 4.6), dpi=200)
Lmax = 4
for i, (k, near) in enumerate((('0.2480', ('0.2475', '0.2480')), ('0.2520', ('0.2520', '0.2525')))):
    d = FK[k]; Ls = sorted(int(L) for L, v in d.items() if v['source'] == 'reweighted')
    y = [d[str(L)]['chi'] for L in Ls]; e = [d[str(L)]['se'] for L in Ls]
    ax.errorbar(Ls, y, yerr=e, fmt='-o', color=SER[i], lw=2, ms=8, mec='white', mew=2, capsize=0, elinewidth=1.4, label=f'k = {float(k):.3f}', zorder=4)
    single = sorted({int(L): v for kk in near for L, v in FK.get(kk, {}).items() if v['source'].startswith('chain')}.items())
    for L, v in single:
        ax.errorbar([L], [v['chi']], yerr=[v['se']], fmt='o', color=SER[i], mfc='white', mew=2, ms=8, capsize=0, elinewidth=1.4, zorder=5)
        Lmax = max(Lmax, L)
    Lmax = max(Lmax, max(Ls))
x0 = np.mean([FK['0.2480']['4']['chi'], FK['0.2520']['4']['chi']]); Lr = np.linspace(4, max(Lmax, 12), 50)
ax.plot(Lr, x0 * (Lr / 4) ** 0.94, color=INK2, lw=1.2, ls='--', zorder=2)
ax.text(Lr[-1] * 0.99, x0 * (Lr[-1] / 4) ** 0.94, 'critical point: L$^{0.94}$', color=INK2, fontsize=9, ha='right', va='bottom')
h = lambda L: (0.16 + 1.11 / L) * L ** 0.64
ax.plot(Lr, x0 * h(Lr) / h(4), color=INK2, lw=1.2, ls=':', zorder=2)
ax.text(Lr[-1] * 0.99, x0 * h(Lr[-1]) / h(4) * 1.04, "shape of Hamber's 1993 fit", color=INK2, fontsize=9, ha='right', va='bottom')
ax.set_xscale('log'); ax.set_yscale('log'); ax.minorticks_off()
ticks = [4, 5, 6, 8, 10, 12] + ([14] if Lmax >= 14 else []) + ([16] if Lmax >= 16 else [])
ax.set_xticks(ticks); ax.set_xticklabels([f'{t}⁴' for t in ticks]); ax.set_yticks([2, 3, 5, 10, 20]); ax.set_yticklabels(['2', '3', '5', '10', '20'])
ax.set_xlim(3.7, max(Lmax, 12) * 1.25); ax.set_ylim(2, 22)
ax.set_xlabel('grid'); ax.set_ylabel('curvature fluctuation χ')
ax.set_title('Curvature fluctuation at fixed coupling against grid size, a = 0.005', fontsize=11, color=INK)
from matplotlib.lines import Line2D
hd = [Line2D([], [], color=SER[0], lw=2, marker='o', ms=7, mec='white', mew=1.5), Line2D([], [], color=SER[1], lw=2, marker='o', ms=7, mec='white', mew=1.5),
      Line2D([], [], color=INK2, lw=0, marker='o', ms=7, mfc='white', mew=1.8)]
ax.legend(hd, ['k = 0.248', 'k = 0.252', 'single chains (14⁴ at k = 0.2475, 0.2525; 16⁴)'], frameon=False, fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig(N + 'chi_fixedk.png'); print('chi_fixedk.png')
