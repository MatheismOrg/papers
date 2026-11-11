"""SA74 night: one pass over everything at the end of the night. Copies the staged Mac files into mac_pull, reruns the size
analysis (fss/analyze_fss.py), the flat-against-rough comparison (hysteresis.py), the distribution and memory statistics
(shape_stats.py) and the figures (fig_night.py), then writes night_summary.json with the numbers the morning summary and the
paper paragraph use: crossings by size, the curvature fluctuation against size at each k (chain means, chain-to-chain errors),
fitted growth exponents, autocorrelation times, and the flat/rough gaps; then the later analyses (reweight.py, multihist.py,
fig_reweight.py, defects.py, analyze_au.py, analyze_as5.py), whose results it adds to night_summary.json. T29: multihist.py
with plan A1's finished 12^4 chains (multihist_T29.json, multihist_skip50_T29.json) and analyze_a1.py (a1_results.json)."""
import json, subprocess, shutil, glob, os, numpy as np
from pathlib import Path
N = Path('/home/claude/sa74/regge/r2/night'); FSS = Path('/home/claude/sa74/regge/r2/fss'); MAC = Path('/home/claude/sa74/NIGHT/mac_pull')
UP = Path('/mnt/user-data/uploads/gr-qm-tests/sa74_night_regge')
for f in glob.glob(str(UP / '*.dat')) + glob.glob(str(UP / '*.log')): shutil.copy(f, MAC)
env = dict(os.environ, MINROWS=os.environ.get('MINROWS', '800'))
for cmd in (['python3', str(FSS / 'analyze_fss.py')], ['python3', str(N / 'hysteresis.py')], ['python3', str(N / 'shape_stats.py')]):
    subprocess.run(cmd, check=True, capture_output=True, env=env, cwd=str(N))
S = json.load(open(FSS / 'fss_results.json'))
out = dict(crossings=S['crossings'], chi_vs_L={}, R_vs_L={}, growth={})
for k in (0.245, 0.2475, 0.25, 0.2525, 0.255):
    pts = sorted((s['L'], s['chi'], s['chi_se'], s['chains'], s['R'], s['R_se']) for s in S['summary'] if abs(s['k'] - k) < 1e-6)
    out['chi_vs_L'][str(k)] = [dict(L=p[0], chi=p[1], chi_se=p[2], chains=p[3]) for p in pts]
    out['R_vs_L'][str(k)] = [dict(L=p[0], R=p[4], R_se=p[5], chains=p[3]) for p in pts]
    if len(pts) >= 3:
        L = np.log([p[0] for p in pts]); y = np.log([p[1] for p in pts]); w = 1 / np.maximum(np.array([p[2] / p[1] for p in pts]), 0.05) ** 2
        A = np.vstack([L, np.ones_like(L)]).T; cov = np.linalg.inv(A.T @ (A * w[:, None])); b = cov @ (A.T @ (w * y))
        out['growth'][str(k)] = dict(slope=float(b[0]), slope_se=float(np.sqrt(cov[0, 0])), sizes=[p[0] for p in pts],
                                     expected_nu_0405=0.94, expected_nu_third=2.0)
H = json.load(open(N / 'hysteresis.json')); out['hysteresis'] = H['gaps']
T = json.load(open(N / 'shape_stats.json')); out['tau'] = [dict(L=r['L'], k=r['k'], tau=r['tau'], chain=r['chain']) for r in T]
json.dump(out, open(N / 'night_summary.json', 'w'), indent=1)
subprocess.run(['python3', str(N / 'fig_night.py'), '0.245,0.2475,0.25,0.2525'], check=True, capture_output=True, cwd=str(N))
# later additions (04:00-04:30): reweighting, multiple histograms and their figure; long edges in the rough starts; a = 0.02; small-a mixing
os.makedirs(MAC / 'cfg', exist_ok=True)
for f in glob.glob(str(UP / 'cfg' / '*.cfg')): shutil.copy(f, MAC / 'cfg')
for f in ('queue.txt', 'queue2.txt', 'night_queue.log', 'night_queue2.log'):   # the runners' final plan files and logs
    if (UP / f).exists(): os.makedirs(MAC / 'runner', exist_ok=True); shutil.copy(UP / f, MAC / 'runner' / f); shutil.copy(UP / f, MAC / f)
# T29: multihist.py with plan A1's finished chains (A1INC=1) on 4^4 to 12^4, largest chi within 0.002 of the couplings run (KPAD); the T28
# results stay frozen in multihist_T28.json and multihist_skip50_T28.json (and fig_reweight.py's figure, from multihist.json, is T28's).
t29 = dict(env, A1INC='1', LMAX='12', KPAD='0.002', SUFFIX='_T29')
for cmd, e in ((['python3', str(N / 'reweight.py')], env), (['python3', str(N / 'multihist.py')], t29), (['python3', str(N / 'multihist.py')], dict(t29, SKIPFRAC='0.5')),
               (['python3', str(N / 'analyze_a1.py')], env), (['python3', str(N / 'defects.py')], env), (['python3', str(N / 'analyze_au.py')], env),
               (['python3', str(N / 'analyze_as5.py')], env)):
    subprocess.run(cmd, check=True, capture_output=True, env=e, cwd=str(N))
M = json.load(open(N / 'multihist_T29.json'))
out['multihist'] = dict(k_range=M['k_range'], peak_power=M['peak_power'], peaks={L: dict(chi=v['peak_chi'], sd=v['peak_chi_boot_sd'], k=v['peak_k'], crossing=v['crossing']) for L, v in M['per_L'].items()})
out['defects'] = json.load(open(N / 'defects.json'))['configs']; out['a002_a003'] = {k: v for k, v in json.load(open(N / 'as5_results.json')).items() if k != 'description'}
out['a02'] = {k: v for k, v in json.load(open(N / 'au_results.json')).items() if k != 'description'}
out['a1'] = {k: v for k, v in json.load(open(N / 'a1_results.json')).items() if k != 'description'}
json.dump(out, open(N / 'night_summary.json', 'w'), indent=1)
print('largest chi (multiple histograms):', {L: round(v['chi'], 2) for L, v in out['multihist']['peaks'].items()}, 'power', round(M['peak_power']['power'], 3), '+-', round(M['peak_power']['se'], 3))
for k, g in out['growth'].items(): print(f"k={k}: chi ~ L^{g['slope']:+.2f} +- {g['slope_se']:.2f} over L={g['sizes']}")
for c in out['crossings']: print('crossing', c['L'], c.get('k0'))
for g in out['hysteresis']: print(f"gap L={g['L']} k={g['k']}: {g['gap']:+.4f} +- {g['gap_se']:.4f} ({g['n_flat']}+{g['n_rough']})")
