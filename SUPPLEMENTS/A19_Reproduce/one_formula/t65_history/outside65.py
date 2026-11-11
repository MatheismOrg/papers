"""PROTOCOL_T65, Test 3, outside checks (d) and (e) for the gate "any cold gas" (H), every number frozen at its all-data
value (res65_alldata_fit.json): Iorio et al.'s eight dwarf irregulars outside SPARC and the 26 SPARC galaxies outside the
benchmark, every point, against T48, the gas gate (T61) and MOND's relation.
Usage: PYTHONPATH=../dmfit:../t59 python3 outside65.py ../t59 res65_alldata_fit.json"""
import sys, os, json, numpy as np
import lib65 as L
import harness_T48 as M8
import compare46 as C
H = L.H; B = L.B; M5 = L.M5; M61 = L.M61; M = L.M; T59 = sys.argv[1]; FIT = json.load(open(sys.argv[2])); th = L.TH
gal = json.load(open(os.path.join(T59, 'lt59_inputs.json'))); rng = np.random.default_rng(654)
def stack(names):
    rows = [L.lt_rows(T59, n) for n in names if gal[n]['f36'] is not None and gal[n]['inc'] >= 30]; return {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}
def score(a, lab):
    B.FROZEN['abc'] = L.ABC; M61.U['t'], M61.U['k'] = L.T_G, L.K_G; M5.U['t'], M5.U['k'] = FIT['t'], FIT['kappa']
    e = dict(H=np.log10(a['go'] / H.gpred(M5.MODEL, th, a)), T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)), T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), MOND=np.log10(a['go'] / C.g_mond(a)))
    un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi); f = {m: np.bincount(gi, v ** 2) / ct for m, v in e.items()}; idx = rng.integers(0, len(un), (10000, len(un)))
    o = dict(galaxies=len(un), points=int(len(gi)), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()})
    for ref in ('T48', 'MOND', 'T61'):
        dd = np.sqrt(f['H'][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1)); o['H_minus_' + ref] = dict(pct=round(float(100 * (np.sqrt(f['H'].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f['H'] < f[ref])), n_higher=int(np.sum(f['H'] > f[ref])))
    print(lab, json.dumps(o), flush=True); return o
out = dict(frozen=dict(t=FIT['t'], kappa=FIT['kappa']))
out['d_iorio_outside_all_points'] = score(stack([n for n, g in gal.items() if g['sparc'] is None]), '(d) Iorio outside, every point')
bench = {g['name'] for g in H.disks()}; A = M.ALL; a = M.pick(A, ~np.isin(A['name'], list(bench))); out['e_sparc_outside_benchmark_all_points'] = score(a, '(e) SPARC outside the benchmark, every point')
json.dump(out, open('res65_outside.json', 'w'), indent=1)
