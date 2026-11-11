"""SA73 T62 (noise60_dropped.py with the final formula): what the benchmark's cut on velocity errors does by itself. Truth = the final
formula with its core (final_T62; the keys named T48 below are that formula) at the catalog values for every SPARC point; each draw redraws distance, inclination, luminosity,
mass-to-light ratio, HI flux (as noise54_disks.py) and every velocity from its published error, then applies the cut to
the redrawn velocity. Scored: the formula, MOND's relation and standard halos on the points the cut drops and on those it
keeps. Usage: PYTHONPATH=<dmfit> python3 noise60_dropped.py [ndraw]"""
import sys, os, json, numpy as np
import final_T62 as M8
import compare46 as C
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
ND = int(sys.argv[1]) if len(sys.argv) > 1 else 500
t1 = {}
lines = open(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt').read().split('\n'); last = max(i for i, l in enumerate(lines) if l.startswith('-----'))
for l in lines[last + 1:]:
    p = l.split()
    if len(p) >= 18: t1[p[0]] = dict(D=float(p[2]), eD=float(p[3]), Inc=float(p[5]), eInc=float(p[6]), L=float(p[7]) * 1e9, eL=float(p[8]) * 1e9, Re=float(p[9]), MHI=float(p[13]) * 1e9, Q=int(p[17]))
gal = []
for g, r in sorted(t1.items()):
    d = np.loadtxt(H.PKG / 'data/SPARC/Rotmod_LTG' / (g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (vb2 > 0)
    if k.sum(): gal.append(dict(name=g, R=R[k], V=V[k], eV=eV[k], Vg=Vg[k], Vd=Vd[k], Vb=Vb[k], **r))
gi = np.concatenate([[i] * len(g['R']) for i, g in enumerate(gal)]); cat = lambda k: np.concatenate([g[k] for g in gal]); per = lambda k: np.array([g[k] for g in gal], float)
R, V, eV, Vg, Vd, Vb = (cat(k) for k in ('R', 'V', 'eV', 'Vg', 'Vd', 'Vb'))
Dg, eD, Inc, eInc, L, eL, Re, MHI = (per(k) for k in ('D', 'eD', 'Inc', 'eInc', 'L', 'eL', 'Re', 'MHI'))
name = np.array([gal[i]['name'] for i in gi]); n = len(gal)
def build(fD, fi, fUd, fUb, fHI, fL):
    gb = (Vg * np.abs(Vg) * fHI[gi] + 0.5 * fUd[gi] * fL[gi] * Vd * np.abs(Vd) + 0.7 * fUb[gi] * fL[gi] * Vb * np.abs(Vb)) / R * H.ACC
    ms = 0.5 * L * fUd * fL * fD ** 2; mp = ms + 1.33 * MHI * fHI * fD ** 2
    return dict(r=R * fD[gi], gb=np.maximum(gb, 1e-16), ms=ms[gi], mp=mp[gi], re=(Re * fD)[gi], z=np.zeros(len(R)), name=name, kind=np.array(['disk'] * len(R)), hw=np.ones(len(R)), go=np.ones(len(R)))
one = np.ones(n); a0 = build(one, one, one, one, one, one); g0 = H.gpred(M.MODEL, th, a0); Vt = np.sqrt(g0 * R / H.ACC)
def hb(e, k):
    c = np.bincount(gi[k], minlength=n); s = np.bincount(gi[k], e[k] ** 2, minlength=n); m = np.bincount(gi[k], e[k], minlength=n); ok = c > 0
    return float(np.sqrt(np.mean(s[ok] / c[ok]))), float(np.mean(m[ok] / c[ok]))
def score(a, go, k): return {nm: hb(np.log10(go / gf), k) for nm, gf in (('T48', H.gpred(M.MODEL, th, a)), ('MOND', C.g_mond(a)), ('STD', C.g_lcdm(a)))}
go_cat = V ** 2 / R * H.ACC; kd = eV / V > 0.1; obs = dict(dropped=score(a0, go_cat, kd), kept=score(a0, go_cat, ~kd), n_dropped=int(kd.sum()), n_kept=int((~kd).sum()))
print('observed', json.dumps(obs), flush=True)
rng = np.random.default_rng(60); res = {s: {m: [] for m in ('T48', 'MOND', 'STD', 'T48_mean', 'MOND_mean', 'T48_minus_MOND', 'n')} for s in ('dropped', 'kept')}
for it in range(ND):
    z = lambda: rng.standard_normal(n)
    fD = np.clip(1 + eD / Dg * z(), .3, None); i2 = np.clip(Inc + eInc * z(), 10., 90.); fi = (np.sin(np.radians(Inc)) / np.sin(np.radians(i2))) ** 2
    fUd = 10 ** (0.11 * z()); fUb = 10 ** (0.11 * z()); fHI = np.clip(1 + .1 * z(), .3, None); fL = np.clip(1 + eL / L * z(), .3, None)
    Vo = Vt + eV * rng.standard_normal(len(R)); pos = Vo > 0; Vs = np.where(pos, Vo, 1.)
    a = build(fD, fi, fUd, fUb, fHI, fL); go = g0 * (Vs / Vt) ** 2 * fi[gi] / fD[gi]
    for s, k in (('dropped', pos & (eV / Vs > 0.1)), ('kept', pos & (eV / Vs <= 0.1))):
        sc = score(a, go, k); r = res[s]
        for m in ('T48', 'MOND', 'STD'): r[m].append(sc[m][0])
        r['T48_mean'].append(sc['T48'][1]); r['MOND_mean'].append(sc['MOND'][1]); r['T48_minus_MOND'].append(sc['T48'][0] - sc['MOND'][0]); r['n'].append(int(k.sum()))
out = dict(observed=obs, draws=ND)
for s in res:
    out[s] = {m: dict(median=round(float(np.median(v)), 5), p2_98=np.round(np.percentile(v, [2.5, 97.5]), 5).tolist()) for m, v in res[s].items()}
    o = obs[s]; out[s]['share_T48_at_least_observed'] = round(float(np.mean(np.array(res[s]['T48']) >= o['T48'][0])), 4)
    out[s]['share_at_least_MOND_observed'] = round(float(np.mean(np.array(res[s]['T48']) >= o['MOND'][0])), 4); out[s]['share_at_least_STD_observed'] = round(float(np.mean(np.array(res[s]['T48']) >= o['STD'][0])), 4)
    out[s]['share_diff_at_least_observed'] = round(float(np.mean(np.array(res[s]['T48_minus_MOND']) >= o['T48'][0] - o['MOND'][0])), 4)
    print(s, json.dumps(out[s]), flush=True)
json.dump(out, open('res62_noise_dropped.json', 'w'), indent=1)
