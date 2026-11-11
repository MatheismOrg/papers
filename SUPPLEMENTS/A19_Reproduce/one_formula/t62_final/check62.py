"""SA73 T62: secondary numbers of Sections 5.4 and 5.8 for the final formula with its core (final_T62): the KiDS early-late
split, the lowest-mass stack's negative point, the held-out lensing shares, and the core's variable w = g L for the blue
SDSS centrals of Mandelbaum et al. (2016) and the Milky Way (zero means the core and its mass term do not touch them).
Usage: python3 check62.py <dmfit dir>"""
import sys, csv, json; sys.path.insert(0, sys.argv[1])
import numpy as np
from scipy.stats import chi2
from scipy.linalg import solve_triangular
import final_T62 as M8
import kids_split as KS
M = M8.P; H = M.H; B = M.B; T61 = M8.T; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); out = {}
gs = H.disks(); SS = KS.setup(gs); out['kids_split'] = KS.check(M.MODEL, th, SS); print('KiDS early-late split, final formula:', out['kids_split'])
T61.U['t'], T61.U['k'] = 0., 0.; out['kids_split_without_core'] = KS.check(M.MODEL, th, SS); T61.U['t'], T61.U['k'] = 0.1501, 0.6229; print('   without the core:', out['kids_split_without_core'])
L = H.Lens(gs, True); st = L.stacks[0]; e = H.esd(M.MODEL, th, st)[L.masks[0]]; y = L.y[0]; C = L.L[0] @ L.L[0].T; j = int(np.argmin(y)); sig = float(np.sqrt(C[j, j]))
zall = solve_triangular(L.L[0], y - e, lower=True); k = np.array([i for i in range(len(y)) if i != j]); zk = solve_triangular(np.linalg.cholesky(C[np.ix_(k, k)]), (y - e)[k], lower=True)
out['negative_point'] = dict(value=float(y[j]), sigma=sig, model=float(e[j]), z=float((y[j] - e[j]) / sig), cost=float(zall @ zall - zk @ zk), stack_chi2=float(zall @ zall)); print('negative point:', out['negative_point'])
for tag, c2, n in (('held out primary', 64.52, 29), ('held out strict', 39.37, 25)): out[tag] = dict(chi2=c2, n=n, p=float(chi2.sf(c2, n))); print(tag, out[tag])
(p1, p0), (q1, q0) = H.pop_coeff(gs); rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv'))); out['blue'] = []
for r in rows:
    if r['color'] != 'blue': continue
    lm = float(r['logMs_eff_VAGC']); ms = 10 ** lm; re = 10 ** (p0 + p1 * lm); mp = ms + 1.33 * 10 ** (q0 + q1 * lm); w = float(T61.gl(np.array([ms]), np.array([mp]), np.array([re]), 0.)[0])
    out['blue'].append(dict(logMs=lm, w=round(w, 4), halo_mass_factor=round(10 ** (0.1501 * w), 4)))
print('blue centrals:', out['blue'])
mstar = 5.43e10; out['milky_way_w'] = float(T61.gl(np.array([mstar]), np.array([mstar * 1.2]), np.array([5.75]), 0.)[0]); print('Milky Way w:', out['milky_way_w'])
json.dump(out, open('res62_check.json', 'w'), indent=1, default=str)
