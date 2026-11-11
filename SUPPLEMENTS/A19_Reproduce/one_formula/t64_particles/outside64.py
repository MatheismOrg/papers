"""PROTOCOL_T64, outside checks (d) and (e) for the fuzzy and warm models, every number frozen at its all-data value
(A, B, C = 11.7481, 0.8178, -0.7119; M_t 11.339; fuzzy and warm numbers from res64_alldata_fit.json; the gated core's
t, kappa = 0.1501, 0.6229): the eight dwarf irregulars of Iorio et al. (2017) outside SPARC (inputs as FORECAST_T59), every
point, and the 26 SPARC galaxies outside the benchmark, every point. The self-interacting model's outside scores come from
sidm64.py's table at the all-data cross-section (sidm64_heldout.py).
Usage: PYTHONPATH=../dmfit:../t59 python3 outside64.py ../t59 res64_alldata_fit.json"""
import sys, os, json, numpy as np
import harness_T64 as M
import harness_T48 as M8
import compare46 as C
import model59 as L
H = M.H; B = M.B; M61 = M.M61; T59 = sys.argv[1]; FIT = json.load(open(sys.argv[2]))
ABC = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); KL, KR = 1.172, 1.424
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
rng = np.random.default_rng(641)
def score(a, lab):
    B.FROZEN['abc'] = ABC; M61.U['t'], M61.U['k'] = 0.1501, 0.6229; M.V['lm22'] = FIT['MOGF0IESBF']['par']; M.V['lhm'] = FIT['MOGF0IESBW']['par']
    e = dict(F=np.log10(a['go'] / H.gpred('MOGF0IESBF', th, a)), W=np.log10(a['go'] / H.gpred('MOGF0IESBW', th, a)), T61=np.log10(a['go'] / H.gpred('MOGF0IESBC', th, a)),
             T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
    un, gi = np.unique(a['name'], return_inverse=True); ct = np.bincount(gi); f = {m: np.bincount(gi, v ** 2) / ct for m, v in e.items()}; idx = rng.integers(0, len(un), (10000, len(un)))
    o = dict(galaxies=len(un), points=int(len(gi)), **{m: round(float(np.sqrt(v.mean())), 4) for m, v in f.items()})
    for mod in ('F', 'W'):
        for ref in ('T48', 'MOND', 'T61'):
            dd = np.sqrt(f[mod][idx].mean(1)) - np.sqrt(f[ref][idx].mean(1))
            o[mod + '_minus_' + ref] = dict(pct=round(float(100 * (np.sqrt(f[mod].mean() / f[ref].mean()) - 1)), 2), boot95=np.round(np.percentile(dd, [2.5, 97.5]), 4).tolist(), n_lower=int(np.sum(f[mod] < f[ref])), n_higher=int(np.sum(f[mod] > f[ref])))
    print(lab, json.dumps(o), flush=True); return o
test = [n for n, g in gal.items() if g['sparc'] is None]; out = dict(frozen=dict(lm22=FIT['MOGF0IESBF']['par'], lhm=FIT['MOGF0IESBW']['par']))
out['d_iorio_outside_all_points'] = score(rows(test, 9.), '(d) Iorio outside, every point')
bench = {g['name'] for g in H.disks()}; A = M.ALL; a = M.pick(A, ~np.isin(A['name'], list(bench))); print('SPARC galaxies outside the benchmark:', len(np.unique(a['name'])))
out['e_sparc_outside_benchmark_all_points'] = score(a, '(e) SPARC outside the benchmark, every point')
json.dump(out, open('res64_outside.json', 'w'), indent=1)
