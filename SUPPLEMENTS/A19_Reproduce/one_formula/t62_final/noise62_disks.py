"""SA73 T62 (noise54_disks.py with the final formula): noise floor for the 149 SPARC disks. Truth = the final
formula with its core (final_T62; frozen A, B, C, t, kappa) at the catalog values. Each draw perturbs what the analyst would have measured, from the
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
rng = np.random.default_rng(54)
def draw(which):
    z = lambda k: rng.standard_normal(n) if (which == 'all' or which == k) else np.zeros(n)
    fD = np.clip(1 + eD / Dg * z('D'), .3, None); i2 = np.clip(Inc + eInc * z('Inc'), 10., 90.); fi = (np.sin(np.radians(Inc)) / np.sin(np.radians(i2))) ** 2
    fUd = 10 ** (0.11 * z('ML')); fUb = 10 ** (0.11 * z('ML')); fHI = np.clip(1 + .1 * z('HI'), .3, None); fL = np.clip(1 + eL / L * z('L'), .3, None)
    fV = np.clip(1 + eV / V * rng.standard_normal(len(R)), .05, None) if which in ('all', 'V') else np.ones(len(R))
    a, fo = build(fD, fi, fUd, fUb, fHI, fL, fV); return np.log10(g0 * fo / H.gpred(M.MODEL, th, a)), a
out = {}
for which in ('V', 'Inc', 'D', 'ML', 'HI', 'L', 'all'):
    nd = ND if which == 'all' else max(100, ND // 5); s = []; pt = []; gm = []
    for it in range(nd):
        e, a = draw(which); s.append(hb(e)); pt.append(np.sqrt(np.mean(e ** 2))); gm.append(np.std(gmean(e)))
    s = np.array(s); out[which] = dict(median=round(float(np.median(s)), 5), p16_84=np.round(np.percentile(s, [16, 84]), 5).tolist(), p2_98=np.round(np.percentile(s, [2.5, 97.5]), 5).tolist(),
                                       point_rms=round(float(np.median(pt)), 4), rms_galaxy_means=round(float(np.median(gm)), 4))
    if which == 'all': out['p_noise_at_least_observed'] = dict(vs_all_data=round(float(np.mean(s >= hb(e_obs))), 4), vs_held_out_0153=round(float(np.mean(s >= 0.1531)), 4), MOND_0155=round(float(np.mean(s >= 0.1546)), 4), standard_0166=round(float(np.mean(s >= 0.1655)), 4))
    print(which, out[which], flush=True)
print(out['p_noise_at_least_observed'])
out['observed'] = dict(host_balanced=round(hb(e_obs), 4), point=round(float(np.sqrt(np.mean(e_obs ** 2))), 4), rms_galaxy_means=round(float(np.std(gmean(e_obs))), 4))
json.dump(out, open('res62_noise_disks.json', 'w'), indent=1)
