"""SA73 T54 (10 Nov 2026): noise floor for the 149 SPARC disks, as noise50 did for dwarfs and clusters. Truth = the final
formula (T48, frozen A, B, C) at the catalog values. Each draw perturbs what the analyst would have measured, from the
published errors: distance (e_D, SPARC Table 1), inclination (e_Inc), luminosity (e_L[3.6]), rotation velocities (errV, per
point), and the adopted errors of Lelli et al. (2017, Section 4.1): 0.11 dex in the stellar mass-to-light ratio of disk and
bulge, 10 percent in the HI flux. Then the benchmark's host-balanced error is computed on the perturbed data.
Usage: python3 -I noise54_disks.py <dmfit dir> [ndraw]"""
import sys, os, json; sys.path.insert(0, sys.argv[1])
import numpy as np
import final_T62 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
ND = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
t1 = {}
lines = open(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt').read().split('\n'); last = max(i for i, l in enumerate(lines) if l.startswith('-----'))
for l in lines[last + 1:]:
    p = l.split()
    if len(p) >= 18: t1[p[0]] = dict(D=float(p[2]), eD=float(p[3]), fD=int(p[4]), Inc=float(p[5]), eInc=float(p[6]), L=float(p[7]) * 1e9, eL=float(p[8]) * 1e9,
                                     Re=float(p[9]), SB=float(p[10]), MHI=float(p[13]) * 1e9, Vf=float(p[15]), Q=int(p[17]))
gal = []
for g, r in sorted(t1.items()):
    if not (r['Q'] < 3 and r['Inc'] >= 30): continue
    d = np.loadtxt(H.PKG / 'data/SPARC/Rotmod_LTG' / (g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= 0.1) & (vb2 > 0)
    if k.sum(): gal.append(dict(name=g, R=R[k], V=V[k], eV=eV[k], Vg=Vg[k], Vd=Vd[k], Vb=Vb[k], **r))
assert len(gal) == 149 and sum(len(g['R']) for g in gal) == 2700
gi = np.concatenate([[i] * len(g['R']) for i, g in enumerate(gal)]); cat = lambda k: np.concatenate([g[k] for g in gal]); per = lambda k: np.array([g[k] for g in gal], float)
R, V, eV, Vg, Vd, Vb = (cat(k) for k in ('R', 'V', 'eV', 'Vg', 'Vd', 'Vb'))
Dg, eD, Inc, eInc, L, eL, Re, MHI = (per(k) for k in ('D', 'eD', 'Inc', 'eInc', 'L', 'eL', 'Re', 'MHI'))
name = np.array([gal[i]['name'] for i in gi]); n = len(gal); ct = np.bincount(gi)
def build(fD, fi, fUd, fUb, fHI, fL, fV):
    """analyst's data set: perturbed observed acceleration (as a factor on the true one) and perturbed inputs of the formula"""
    gb = (Vg * np.abs(Vg) * fHI[gi] + 0.5 * fUd[gi] * fL[gi] * Vd * np.abs(Vd) + 0.7 * fUb[gi] * fL[gi] * Vb * np.abs(Vb)) / R * H.ACC
    ms = 0.5 * L * fUd * fL * fD ** 2; mp = ms + 1.33 * MHI * fHI * fD ** 2
    a = dict(r=R * fD[gi], gb=np.maximum(gb, 1e-16), ms=ms[gi], mp=mp[gi], re=(Re * fD)[gi], z=np.zeros(len(R)), name=name, kind=np.array(['disk'] * len(R)))
    a['hw'] = 1 / np.sqrt(n * ct[gi]); a['go'] = np.ones(len(R))
    return a, fV ** 2 * fi[gi] / fD[gi]
one = np.ones(n); a0, _ = build(one, one, one, one, one, one, np.ones(len(R)))
g0 = H.gpred(M.MODEL, th, a0); go_cat = V ** 2 / R * H.ACC; e_obs = np.log10(go_cat / g0)
def hb(e): return float(np.sqrt(np.mean(np.bincount(gi, e ** 2) / ct)))
def gmean(e): return np.bincount(gi, e) / ct
print('observed (all data, frozen coefficients): host-balanced rms', round(hb(e_obs), 4), ' point rms', round(float(np.sqrt(np.mean(e_obs ** 2))), 4),
      ' rms of galaxy means', round(float(np.std(gmean(e_obs))), 4), flush=True)
print('median errors: D', round(float(np.median(eD / Dg)), 3), 'Inc deg', np.median(eInc), 'errV/V', round(float(np.median(eV / V)), 3), 'L', round(float(np.median(eL / L)), 3))

import compare46 as C
props = {'log M*': np.log10(0.5 * L), 'log R_e': np.log10(Re), 'log gas/stars': np.log10(np.maximum(1.33 * MHI / (0.5 * L), 1e-4)), 'log SB': np.log10(per('SB')), 'log Vflat': np.log10(np.maximum(per('Vf'), 1.))}
lx = np.log10(R / Re[gi]); rng = np.random.default_rng(7)
def within(e): return e - gmean(e)[gi]
def slopes(e):
    w = within(e); sl = np.full(n, np.nan)
    for i in range(n):
        k = gi == i
        if k.sum() >= 5 and np.ptp(lx[k]) > .3: sl[i] = np.polyfit(lx[k], w[k], 1)[0]
    return sl
models = {'formula': g0, 'MOND': C.g_mond(a0), 'standard halos': C.g_lcdm(a0)}
for nm, g in models.items():
    e = np.log10(go_cat / g); w = within(e); sl = slopes(e); ok = np.isfinite(sl)
    print(f'{nm:15s} host-bal {hb(e):.4f}  galaxy means {np.sqrt(np.mean(gmean(e)**2)):.4f}  within {hb(w):.4f}  median slope {np.nanmedian(sl):+.3f}  sd slopes {np.nanstd(sl):.3f}  n {ok.sum()}')
    print('    stacked within-galaxy residual by r/R_e:', ' '.join(f'[{lo:.1f},{hi:.1f}) {np.mean(w[(lx>=np.log10(lo))&(lx<np.log10(hi))]):+.3f}' for lo, hi in ((.01,.3),(.3,.6),(.6,1),(1,1.5),(1.5,2.5),(2.5,4),(4,7),(7,50))))
    print('    slope vs property:', {k: round(float(np.corrcoef(v[ok], sl[ok])[0, 1]), 2) for k, v in props.items()})
# noise: sd of per-galaxy slopes
sds = []
for it in range(200):
    fV = np.clip(1 + eV / V * rng.standard_normal(len(R)), .05, None); sl = slopes(np.log10(fV ** 2)); sds.append(np.nanstd(sl))
print('noise (velocity errors) sd of per-galaxy slopes: median', round(float(np.median(sds)), 3), '95%', np.round(np.percentile(sds, [2.5, 97.5]), 3))
