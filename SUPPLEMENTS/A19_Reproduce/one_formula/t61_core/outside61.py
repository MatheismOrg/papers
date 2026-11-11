"""PROTOCOL_T61, outside checks (d) and (e): T61 with A, B, C, t, kappa frozen at the all-data values, on the eight dwarf
irregulars of Iorio et al. (2017) outside SPARC (inputs as FORECAST_T59) and on the 26 SPARC galaxies outside the benchmark.
Usage: PYTHONPATH=<dmfit>:<t59 folder> python3 outside61.py <t59 folder>"""
import sys, os, json, numpy as np
import harness_T61 as M
import harness_T57d as M7
import harness_T48 as M8
import compare46 as C
import model59 as L
H = M.H; B = M.B; T59 = sys.argv[1]
ABC = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); TK = (0.1501, 0.6229); KL, KR = 1.172, 1.424
gal = json.load(open(os.path.join(T59, 'lt59_inputs.json')))
def rows(names, cut):
    out = []
    for nm in names:
        g = gal[nm]
        if g['f36'] is None or g['inc'] < 30: continue
        R = np.array(g['R']); V = np.array(g['Vc']); eV = np.array(g['eVc'])
        Ms = 0.5 * L.l36(g['f36'][0], g['D']) / KL; Re = 1.678 * g['read']['Rstar'] * g['D'] / g['read']['D'] / KR
        vg2, mhi, mcut = L.gas(R, g['Sd'], 'exp'); vb2 = vg2 + L.stars_v2(R, Ms, Re / 1.678)
        k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= cut) & (vb2 > 0); Mp = Ms + 1.33 * mhi
        out += [dict(name=nm, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=Re, z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * H.ACC, V[k] ** 2 / R[k] * H.ACC)]
    return H.arr(out)
rng = np.random.default_rng(611)
def score(a, lab):
    B.FROZEN['abc'] = ABC; M.U['t'], M.U['k'] = TK; M7.U['t'], M7.U['s'] = 0.361, 0.146
    e = dict(T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)), T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), T57=np.log10(a['go'] / H.gpred('MOGF0IESBD', th, a)),
             MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
    un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi); f = {m: np.bincount(gi, v ** 2) / ct for m, v in e.items()}; idx = rng.integers(0, len(un), (10000, len(un)))
    o = dict(galaxies=len(un), points=int(len(gi)), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()}, per_galaxy={n: {m: round(float(np.sqrt(f[m][i])), 4) for m in f} for i, n in enumerate(un)} if len(un) <= 30 else None)
    for ref in ('T48', 'MOND'):
        dd = np.sqrt(f['T61'][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
        o['T61_minus_' + ref] = dict(pct=round(float(100 * (np.sqrt(f['T61'].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f['T61'] < f[ref])), n_higher=int(np.sum(f['T61'] > f[ref])))
    print(lab, json.dumps({k: v for k, v in o.items() if k != 'per_galaxy'}), flush=True); return o
test = [n for n, g in gal.items() if g['sparc'] is None]; shared = [n for n, g in gal.items() if g['sparc'] is not None]; out = dict(frozen=dict(t=TK[0], kappa=TK[1]))
out['d_iorio_outside_all_points'] = score(rows(test, 9.), '(d) Iorio outside, every point'); out['iorio_outside_benchmark_cut'] = score(rows(test, .1), 'Iorio outside, 10 percent cut')
out['iorio_shared_all_points'] = score(rows(shared, 9.), 'Iorio shared with SPARC, every point')
bench = {g['name'] for g in H.disks()}; A = M.ALL; k = ~np.isin(A['name'], list(bench)); a = M.pick(A, k); print('SPARC galaxies outside the benchmark:', len(np.unique(a['name'])))
out['e_sparc_outside_benchmark_all_points'] = score(a, '(e) SPARC outside the benchmark, every point')
out['sparc_outside_benchmark_cut'] = score(M.pick(a, a['fe'] <= .1), 'SPARC outside the benchmark, 10 percent cut'); out['sparc_outside_benchmark_dropped'] = score(M.pick(a, a['fe'] > .1), 'SPARC outside the benchmark, dropped points')
json.dump(out, open('res61_outside.json', 'w'), indent=1)
