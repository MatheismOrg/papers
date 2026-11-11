"""FORECAST_T59: the T57 term (halos of diffuse gas-rich disks bigger and looser) and the final formula T48, all
coefficients frozen, on dwarf irregulars of Iorio et al. (2017) that are not in SPARC. Inputs built by model59.py; the two
photometric factors KL, KR come from validate59.py (seven galaxies shared with SPARC, input ratios only).
Usage: PYTHONPATH=<dmfit> python3 test59.py"""
import os, sys, json, numpy as np
import model59 as M
import harness_T57d as M7
import harness_T48 as M8
import compare46 as C
H = M7.H; B = M7.B; D = H.D
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); M7.U['t'] = 0.361; M7.U['s'] = 0.146
KL, KR = 1.172, 1.424
CLEAN = ['DDO 52', 'DDO 210', 'UGC 8508', 'CVnIdwA']
gal = json.load(open('lt59_inputs.json'))

def rows(names, kl, kr, mode='exp', vkey='Vc'):
    out = []
    for nm in names:
        g = gal[nm]
        if g['f36'] is None or g['inc'] < 30: continue
        R = np.array(g['R']); V = np.array(g[vkey]); eV = np.array(g['eVc'])
        Ms = 0.5 * M.l36(g['f36'][0], g['D']) / kl; Re = 1.678 * g['read']['Rstar'] * g['D'] / g['read']['D'] / kr
        vg2, mhi, mcut = M.gas(R, g['Sd'], mode); vb2 = vg2 + M.stars_v2(R, Ms, Re / 1.678)
        k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= 0.1) & (vb2 > 0)
        Mp = Ms + 1.33 * (mhi if mode == 'exp' else mcut)
        out += [dict(name=nm, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=Re, z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * H.ACC, V[k] ** 2 / R[k] * H.ACC)]
    return out
def sparc_rows(keep):
    t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); out = []
    for g, r in sorted(t1.items()):
        if not keep(g, r): continue
        d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
        R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
        vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= 0.1) & (vb2 > 0)
        Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']
        out += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
    return out
def per_galaxy(rw):
    a = H.arr(rw); un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi)
    e = dict(T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), T57=np.log10(a['go'] / H.gpred('MOGF0IESBD', th, a)),
             MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
    ms = {k: np.bincount(gi, v ** 2) / ct for k, v in e.items()}; mean = {k: np.bincount(gi, v) / ct for k, v in e.items()}
    w = np.array([M7.gl(a['ms'][gi == i][:1], a['mp'][gi == i][:1], a['re'][gi == i][:1], 0.)[0] for i in range(len(un))])
    info = {n: dict(points=int(ct[i]), Ms=float(a['ms'][gi == i][0]), Mp=float(a['mp'][gi == i][0]), Re=float(a['re'][gi == i][0]), gL=round(float(w[i]), 3),
                    **{k + '_rms': round(float(np.sqrt(ms[k][i])), 4) for k in ms}, **{k + '_mean': round(float(mean[k][i]), 4) for k in ms}) for i, n in enumerate(un)}
    return un, ms, mean, ct, info
rng = np.random.default_rng(59)
def summary(un, ms, mean, ct, sel=None):
    sel = np.ones(len(un), bool) if sel is None else sel; r = dict(n=int(sel.sum()), points=int(ct[sel].sum()))
    if r['n'] == 0: return r
    for k in ms:
        r[k] = round(float(np.sqrt(ms[k][sel].mean())), 4); r[k + '_mean'] = round(float(mean[k][sel].mean()), 4)
        r[k + '_within'] = round(float(np.sqrt(np.mean(np.maximum(ms[k][sel] - mean[k][sel] ** 2, 0)))), 4)
    x, y = ms['T48'][sel], ms['T57'][sel]; idx = rng.integers(0, len(x), (10000, len(x))); dd = np.sqrt(y[idx].mean(1)) - np.sqrt(x[idx].mean(1))
    r.update(pct=round(float(100 * (np.sqrt(y.mean() / x.mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), p_not_better=round(float(np.mean(dd >= 0)), 4),
             n_better=int(np.sum(y < x)), n_worse=int(np.sum(y > x)))
    for o in ('MOND', 'STD'):
        z = ms[o][sel]; d2 = np.sqrt(x[idx].mean(1)) - np.sqrt(z[idx].mean(1)); r['T48_minus_' + o + '_boot95'] = np.round(np.percentile(d2, [2.5, 97.5]), 4).tolist()
    return r
test = [n for n, g in gal.items() if g['sparc'] is None]; shared = [n for n, g in gal.items() if g['sparc'] is not None]
out = dict(KL=KL, KR=KR, no_lvl_flux=[n for n in test if gal[n]['f36'] is None])
def run(lab, names, kl, kr, mode='exp', vkey='Vc', clean=True):
    un, ms, mean, ct, info = per_galaxy(rows(names, kl, kr, mode, vkey)); o = dict(all=summary(un, ms, mean, ct), galaxies=info)
    if clean: o['clean'] = summary(un, ms, mean, ct, np.isin(un, CLEAN))
    print(lab, json.dumps({k: v for k, v in o.items() if k != 'galaxies'})); out[lab] = o; return un, ms, mean, ct
P = run('test_primary', test, KL, KR)
run('test_uncalibrated', test, 1., 1.); run('test_cut_gas', test, KL, KR, 'cut'); run('test_no_drift_correction', test, KL, KR, 'exp', 'Vr')
run('shared_pipeline', shared, KL, KR, clean=False); run('shared_pipeline_uncalibrated', shared, 1., 1., clean=False)
sp = {gal[n]['sparc'] for n in shared}; un, ms, mean, ct, info = per_galaxy(sparc_rows(lambda g, r: g in sp)); out['shared_sparc_inputs'] = dict(all=summary(un, ms, mean, ct), galaxies=info); print('shared_sparc_inputs', json.dumps(out['shared_sparc_inputs']['all']))
u8, m8, me8, c8, i8 = per_galaxy(sparc_rows(lambda g, r: not (r['Q'] < 3 and r['Inc'] >= 30)))            # the 17 galaxies of FORECAST_T58
un = np.concatenate([P[0], u8]); ms = {k: np.concatenate([P[1][k], m8[k]]) for k in m8}; mean = {k: np.concatenate([P[2][k], me8[k]]) for k in m8}; ct = np.concatenate([P[3], c8])
out['T58_alone'] = summary(u8, m8, me8, c8); out['union_with_T58'] = summary(un, ms, mean, ct); print('T58_alone', json.dumps(out['T58_alone'])); print('union_with_T58', json.dumps(out['union_with_T58']))
json.dump(out, open('res59_test.json', 'w'), indent=1)
for n, v in out['test_primary']['galaxies'].items(): print('%-9s' % n, v)
