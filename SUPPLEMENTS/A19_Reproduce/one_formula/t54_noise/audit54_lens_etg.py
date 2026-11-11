"""SA73 T54: is the formula distinguishable from a perfect one in the lensing stacks (chi-square against its own covariance,
nothing fitted to lensing) and in the 16 early types (noise floor from Lelli et al. 2017's published errors)?
Usage: python3 -I audit54_lens_etg.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.stats import chi2
from scipy.linalg import solve_triangular
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
gs = H.disks(); out = {}
import compare46 as C
for drop in (False, True):
    L = H.Lens(gs, drop); tot = 0; n = 0; rows = []
    for j, st in enumerate(L.stacks):
        for nm, ef in (('formula', lambda s: H.esd(M.MODEL, th, s)), ('MOND', C.esd_mond), ('standard', C.esd_lcdm)):
            z = solve_triangular(L.L[j], L.y[j] - ef(st)[L.masks[j]], lower=True); c2 = float(z @ z)
            rows.append((j + 1, nm, len(z), round(c2, 2), round(float(chi2.sf(c2, len(z))), 4), np.round(z, 2).tolist() if nm == 'formula' else None,
                         np.round(L.y[j], 2).tolist() if nm == 'formula' else None))
    for nm in ('formula', 'MOND', 'standard'):
        c2 = sum(r[3] for r in rows if r[1] == nm); n = sum(r[2] for r in rows if r[1] == nm)
        out[f'lens|{"strict" if drop else "primary"}|{nm}'] = dict(chi2=round(c2, 2), n=n, p=float(chi2.sf(c2, n)), by_stack=[(r[0], r[2], r[3], r[4]) for r in rows if r[1] == nm])
        print('lensing', 'strict' if drop else 'primary', nm, 'chi2', round(c2, 2), 'N', n, 'p', '%.2e' % chi2.sf(c2, n), [(r[2], r[3]) for r in rows if r[1] == nm])
    for r in rows:
        if r[1] == 'formula': print('   stack', r[0], 'whitened residuals', r[5], '| data', r[6])
# held-out lensing chi-squares of the benchmark, same N
for tag, c2, n in (('held out primary', 64.7, out['lens|primary|formula']['n']), ('held out strict', 39.4, out['lens|strict|formula']['n'])):
    print(tag, c2, 'N', n, 'p', '%.2e' % chi2.sf(c2, n)); out[tag] = dict(chi2=c2, n=n, p=float(chi2.sf(c2, n)))
# early types
meta = {}
for line in (H.AST / 'inputs/data/ETG/ETG_Lelli2017.mrt').read_text().splitlines():
    if line and not line.startswith('#'):
        t = line.split(); meta[t[0]] = np.array(t[1:], float)
rows = H.etg_rows(); a = H.arr(rows)
eo = np.array([meta[r['name']][12 + 2 * (r['kind'] == 'ETG_outer')] for r in rows]); eb = np.array([meta[r['name']][16 + 2 * (r['kind'] == 'ETG_outer')] for r in rows])
eD = np.array([meta[r['name']][1] / meta[r['name']][0] for r in rows])
g0 = H.gpred(M.MODEL, th, a); e_obs = np.log10(a['go'] / g0); hs = np.unique(a['name'])
def rms(e): return float(np.sqrt(np.mean([np.mean(e[a['name'] == h] ** 2) for h in hs])))
rng = np.random.default_rng(542); sims = []; gi = np.array([list(hs).index(n) for n in a['name']])
for it in range(2000):
    zb = rng.standard_normal(len(hs))[gi]                     # mass-to-light error shared by a galaxy's two radii
    db = eb * zb; do = eo * rng.standard_normal(len(eo)); b = dict(a); b['gb'] = a['gb'] * 10 ** db; b['ms'] = a['ms'] * 10 ** db; b['mp'] = a['mp'] * 10 ** db
    sims.append(rms(np.log10(g0 * 10 ** do / H.gpred(M.MODEL, th, b))))
sims = np.array(sims)
out['ETG'] = dict(observed_all_data=round(rms(e_obs), 4), held_out=0.185, noise=round(float(np.median(sims)), 4), p16_84=np.round(np.percentile(sims, [16, 84]), 4).tolist(),
                  p2_98=np.round(np.percentile(sims, [2.5, 97.5]), 4).tolist(), p_noise_at_least_held_out=round(float(np.mean(sims >= 0.185)), 4),
                  med_err_gobs=float(np.median(eo)), med_err_gbar=float(np.median(eb)), med_frac_err_D=round(float(np.median(eD)), 3),
                  mean_resid_inner=round(float(np.mean(e_obs[a['kind'] == 'ETG_inner'])), 3), mean_resid_outer=round(float(np.mean(e_obs[a['kind'] == 'ETG_outer'])), 3))
print('ETG', out['ETG'])
json.dump(out, open('res54_lens_etg.json', 'w'), indent=1)
