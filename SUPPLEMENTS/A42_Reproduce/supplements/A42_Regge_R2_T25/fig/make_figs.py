"""Figures for the T25 curvature-squared rerun. Flat starts blue solid, rough starts orange dashed (validated pair)."""
import json, sys, glob, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R2 = Path(__file__).resolve().parent.parent
AUD = R2 / 'audit_a005'
BLUE, ORANGE, INK, INK2 = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.6})
sys.path.insert(0, str(AUD)); from analyze_r2 import batches

def traj(planfile, keys, out, ylim=(-0.5, 1.5)):
    plan = json.loads((AUD / planfile).read_text())
    fig, axs = plt.subplots(2, 2, figsize=(10, 6.4), sharex=True)
    for ax, (L, k) in zip(axs.flat, keys):
        jobs = [j for j in plan['main'] if j['L'] == L and abs(j['k'] - k) < 1e-9]
        nth = jobs[0]['thermal']
        ax.axvspan(0, nth, color='#efeeea', zorder=0)
        for j in jobs:
            a = np.loadtxt(AUD / 'runs' / f"{j['tag']}.dat"); x = a[:, 4]; n = len(x) // 100
            m = x[:n * 100].reshape(n, 100).mean(1); t = np.arange(n) * 100 + 50
            flat = j['start'] == 'flat'
            ax.plot(t, m, color=BLUE if flat else ORANGE, ls='-' if flat else '--', lw=1.6,
                    label=('flat start' if flat else 'rough start') if j['tag'].endswith('s0') else None)
        ax.set_ylim(*ylim); ax.set_title(f'L = {L}, k = {k}', fontsize=11, color=INK)
        ax.text(nth - 150, ylim[0] + 0.62 * (ylim[1] - ylim[0]), 'equilibration', ha='right', va='center', rotation=90, fontsize=9, color=INK2)
    for ax in axs[1]: ax.set_xlabel('sweep')
    for ax in axs[:, 0]: ax.set_ylabel('average curvature R')
    axs[0, 0].legend(frameon=False, loc='upper right')
    fig.tight_layout(); fig.savefig(out, dpi=200, bbox_inches='tight'); plt.close(fig)

def scan(out):
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.axvspan(0.0, 2 * np.sqrt(0.005), color='#eef4fc', zorder=0)
    ax.text(0.071, 0.27, 'action bounded below\n(k < 2√(aλ) = 0.141)', ha='center', va='top', fontsize=9, color=INK2)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.axvline(0.2443, color=INK2, lw=1.0, ls=':')
    ax.text(0.2463, 0.29, 'Hamber k$_c$ = 0.2443', fontsize=9, color=INK2, rotation=90, va='top')
    for st, col, mk, ls in [('flat', BLUE, 'o', '-'), ('rough', ORANGE, 's', '--')]:
        ks, ms, es = [], [], []
        for fn in sorted(glob.glob(str(R2 / 'scan' / f'L3_k*_{st}.dat'))):
            k = float(re.search(r'L3_k([0-9.]+)_', fn).group(1)); a = np.loadtxt(fn); x = a[a[:, 1] == 1][:, 4]
            ks.append(k); ms.append(x.mean()); es.append(max(b['se'] for b in batches(x)))
        ax.errorbar(ks, ms, yerr=es, color=col, marker=mk, ms=6, lw=1.4, ls=ls, capsize=2, label=f'3⁴, {st} start (scan)')
    for st, col in [('flat', BLUE), ('rough', ORANGE)]:
        ks, ms, es = [], [], []
        for fn in sorted(glob.glob(str(R2 / 'scan4' / f'L4_k*_{st}.dat'))):
            k = float(re.search(r'L4_k([0-9.]+)_', fn).group(1)); a = np.loadtxt(fn); x = a[a[:, 1] == 1][:, 4]
            ks.append(k); ms.append(x.mean()); es.append(max(b['se'] for b in batches(x)))
        if ks: ax.errorbar(ks, ms, yerr=es, color=col, marker='^', ms=7, ls='none', capsize=2, label=f'4⁴, {st} start (scan)')
    g = []
    for name in ('results_A.json', 'results_B.json'):
        if (AUD / name).exists(): g += json.loads((AUD / name).read_text())['groups']
    for L, mk in [(3, 'D'), (4, '^')]:
        pts = sorted((x['k'], np.mean([x['by_start']['flat']['mean_R'], x['by_start']['rough']['mean_R']])) for x in g if x['L'] == L)
        if pts: ax.plot(*zip(*pts), ls='none', marker=mk, ms=8, mfc='none', mec=INK, mew=1.4, label=f'{L}⁴ reruns, both starts')
    ax.set_xlabel('k (inverse Newton constant, lattice units)'); ax.set_ylabel('average curvature R')
    ax.set_xlim(0.0, 0.31); ax.set_ylim(-0.38, 0.3)
    ax.legend(frameon=False, loc='lower right', fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200, bbox_inches='tight'); plt.close(fig)

if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'A': traj('plan.json', [(3, 0.05), (3, 0.065), (4, 0.05), (4, 0.065)], R2 / 'fig' / 'regge_r2_audit.png')
    if what == 'B': traj('plan_B.json', [(3, 0.20), (3, 0.24), (4, 0.20), (4, 0.24)], R2 / 'fig' / 'regge_r2_nearkc.png', ylim=(-0.4, 1.6))
    if what == 'scan': scan(R2 / 'fig' / 'regge_r2_scan.png')
    print('ok', what)
