"""SA73 T54: what is the lowest lensing chi-square any isolated-halo profile could reach? For each KiDS stack, the best
profile that is non-negative and does not rise with radius (true of any centrally concentrated mass) is found by
non-negative least squares in the whitened space; its null distribution comes from data redrawn about the formula's own
prediction with the published covariance. Also the three accounts on the points left when each stack's outermost bin and
the one negative point are set aside. Usage: python3 lensfloor54.py <dmfit dir>"""
import sys, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.optimize import nnls
from scipy.stats import chi2
from scipy.linalg import solve_triangular
import final_T62 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
import compare46 as C
gs = H.disks(); rng = np.random.default_rng(545); out = {}
def floor(Lc, y):
    n = len(y); A = solve_triangular(Lc, np.tril(np.ones((n, n))), lower=True); b = solve_triangular(Lc, y, lower=True); d, rn = nnls(A, b); return rn ** 2
for drop in (False, True):
    L = H.Lens(gs, drop); tag = 'strict' if drop else 'primary'; fl = []; c2 = []; sim_fl = np.zeros(4000); sim_c2 = np.zeros(4000); n = 0; mono = True
    for j, st in enumerate(L.stacks):
        e = H.esd(M.MODEL, th, st)[L.masks[j]]; mono &= bool(np.all(np.diff(e) >= -1e-9) and e[0] >= 0)
        z = solve_triangular(L.L[j], L.y[j] - e, lower=True); c2.append(float(z @ z)); fl.append(floor(L.L[j], L.y[j])); n += len(e)
        for it in range(4000):
            zz = rng.standard_normal(len(e)); sim_c2[it] += zz @ zz; sim_fl[it] += floor(L.L[j], e + L.L[j] @ zz)
    out[tag] = dict(n=n, formula_chi2=round(sum(c2), 2), floor_chi2=round(sum(fl), 2), floor_by_stack=np.round(fl, 2).tolist(), formula_by_stack=np.round(c2, 2).tolist(),
                    formula_profile_monotone=mono, floor_null_median=round(float(np.median(sim_fl)), 2), floor_null_95=np.round(np.percentile(sim_fl, [2.5, 97.5]), 2).tolist(),
                    p_floor=float(np.mean(sim_fl >= sum(fl))), excess_over_floor=round(sum(c2) - sum(fl), 2),
                    excess_null_median=round(float(np.median(sim_c2 - sim_fl)), 2), excess_null_95=np.round(np.percentile(sim_c2 - sim_fl, [2.5, 97.5]), 2).tolist(),
                    p_excess=float(np.mean(sim_c2 - sim_fl >= sum(c2) - sum(fl))))
    print(tag, json.dumps(out[tag]), flush=True)
# clean set: strict mask without stack 1's negative point
L = H.Lens(gs, True)
for nm, ef in (('formula', lambda s: H.esd(M.MODEL, th, s)), ('MOND', C.esd_mond), ('standard', C.esd_lcdm)):
    tot = 0; n = 0
    for j, st in enumerate(L.stacks):
        k = np.flatnonzero(L.y[j] > 0); Cj = L.L[j] @ L.L[j].T; Lk = np.linalg.cholesky(Cj[np.ix_(k, k)])
        z = solve_triangular(Lk, (L.y[j] - ef(st)[L.masks[j]])[k], lower=True); tot += float(z @ z); n += len(k)
    out['clean|' + nm] = dict(chi2=round(tot, 3), n=n, p=float(chi2.sf(tot, n))); print('clean', nm, out['clean|' + nm])
json.dump(out, open('res62_lensfloor.json', 'w'), indent=1)
