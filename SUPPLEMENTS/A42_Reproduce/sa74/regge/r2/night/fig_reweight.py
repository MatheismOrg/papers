"""SA74 night figure: the curvature fluctuation chi(k) on each grid from 3^4 to 10^4, combined from all flat-start chains on that
grid by multiple-histogram reweighting (multihist.py), with the largest value between k = 0.2415 and 0.258 marked. A critical point
would make that value grow as L^0.94 (delta = 0.62); the inset compares the largest values, with bootstrap errors, with that."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R2 = '/home/claude/sa74/regge/r2'
M = json.load(open(f'{R2}/night/multihist.json')); D = M['per_L']
INK, INK2 = '#0b0b0b', '#52514e'
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']   # reference palette, fixed order
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.6})
fig, ax = plt.subplots(figsize=(6.8, 4.8), dpi=200)
Ls = [L for L in sorted(int(x) for x in D) if L <= 10]
lo, hi = M['k_range']
peaks = []
for i, L in enumerate(Ls):
    c = [r for r in D[str(L)]['curve'] if lo - 1e-9 <= r[0] <= hi + 1e-9]
    k = np.array([r[0] for r in c]); chi = np.array([r[1] for r in c])
    ax.plot(k, chi, color=SER[i], lw=2.0, label=f'{L}⁴', solid_capstyle='round')
    j = int(np.argmax(chi)); peaks.append((L, chi[j], D[str(L)]['peak_chi_boot_sd']))
    ax.plot(k[j], chi[j], 'o', ms=8, color=SER[i], mec='white', mew=2, zorder=5)
    for kk in D[str(L)]['ks']:
        if lo <= kk <= hi: ax.plot([kk], [0.3 + 0.25 * i], '|', color=SER[i], ms=7, mew=1.6)
ax.text(lo + 0.0003, 1.95, 'simulated couplings (ticks)', fontsize=8, color=INK2)
ax.axvline(0.2443, color=INK2, lw=1.0, ls=':')
ax.text(0.24445, 6.75, 'Hamber k$_c$ = 0.2443', fontsize=8.5, color=INK2, rotation=90, va='bottom')
ax.set_xlim(lo, hi); ax.set_ylim(0, 12.5)
ax.set_xlabel('k (inverse Newton coupling, lattice units)'); ax.set_ylabel('curvature fluctuation χ (reweighted)')
ax.set_title('Curvature fluctuation between the simulated couplings, a = 0.005', fontsize=11.5, color=INK)
ax.legend(frameon=False, fontsize=9, ncol=3, loc='upper left', title='grid', title_fontsize=9)
pp = M['peak_power']
ax.text(0.2462, 9.0, "largest χ grows as L$^{%+.2f \\pm %.2f}$\n(a critical point: L$^{0.94}$)" % (pp['power'], pp['se']), fontsize=8.5, color=INK, va='top')
ins = fig.add_axes([0.70, 0.60, 0.25, 0.22])
Lp = np.array([p[0] for p in peaks]); P = np.array([p[1] for p in peaks]); E = np.array([p[2] for p in peaks])
ins.errorbar(Lp, P, yerr=E, fmt='o', color=INK, ms=3.5, lw=1, capsize=1.5)
Lr = np.array([3, 10]); ins.plot(Lr, P[0] * (Lr / 3) ** 0.94, color=INK2, lw=1.0, ls='--')
ins.text(4.6, P[0] * (4.6 / 3) ** 0.94 * 1.25, 'L$^{0.94}$', fontsize=8, color=INK2)
ins.set_xscale('log'); ins.set_yscale('log'); ins.minorticks_off()
ins.set_xticks([3, 4, 6, 8, 10]); ins.set_xticklabels(['3', '4', '6', '8', '10'], fontsize=7)
ins.set_yticks([3, 5, 10, 20]); ins.set_yticklabels(['3', '5', '10', '20'], fontsize=7); ins.set_ylim(2.5, 25)
ins.set_title('largest χ against L', fontsize=8, color=INK2); ins.grid(True, color='#e6e5e1', lw=0.5)
fig.savefig(f'{R2}/night/chi_reweighted.png', bbox_inches='tight'); print('wrote chi_reweighted.png', [(a, round(b, 2), round(c, 2)) for a, b, c in peaks])
