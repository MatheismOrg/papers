"""SA73 T54: the noise floors of noise50_dwarfs.py and noise50_clusters.py again, keeping the full distributions, to place the
formula's held-out errors and those of MOND's relation and standard halos (Table 8) within them; and the early types with
their distance errors also moving the formula's inputs. Truth = the final formula. Usage: python3 audit54_rest.py <dmfit dir>"""
import sys, csv, json; sys.path.insert(0, sys.argv[1])
import numpy as np
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
out = {}
def place(tag, sims, vals):
    sims = np.asarray(sims); d = dict(noise_median=round(float(np.median(sims)), 4), noise68=np.round(np.percentile(sims, [16, 84]), 4).tolist(), noise95=np.round(np.percentile(sims, [2.5, 97.5]), 4).tolist())
    for k, v in vals.items(): d[k] = dict(error=v, p_noise_at_least=round(float(np.mean(sims >= v)), 4))
    out[tag] = d; print(tag, json.dumps(d), flush=True)
# dwarfs (as noise50_dwarfs.py)
rows = H.dwarf_rows(); a = H.arr(rows)
rr, _ = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(H.PKG / 'blind_test')); R = {x['key']: x for x in rr}
def f(x, k, d=np.nan):
    try: return float(x[k])
    except Exception: return d
es = []; er = []; ed = []; em = []
for n in a['name']:
    x = R[n]; s = f(x, 'vlos_sigma'); es.append(.5 * (f(x, 'vlos_sigma_em', 0.) + f(x, 'vlos_sigma_ep', 0.)) / s)
    er.append(.5 * (f(x, 'rhalf_em', 0.) + f(x, 'rhalf_ep', 0.)) / f(x, 'rhalf')); ed.append(.5 * (f(x, 'distance_modulus_em', .1) + f(x, 'distance_modulus_ep', .1)))
    em.append(.5 * (f(x, 'apparent_magnitude_v_em', .2) + f(x, 'apparent_magnitude_v_ep', .2)))
es, er, ed, em = map(lambda v: np.nan_to_num(np.array(v), nan=0.), (es, er, ed, em))
g0 = H.gpred(M.MODEL, th, a); e_obs = np.log10(a['go'] / g0)
def slope(e, b): X = np.c_[np.ones(len(e)), np.log10(b['ms']), np.log10(b['re'])]; return float(np.linalg.lstsq(X, e, rcond=None)[0][2])
rng = np.random.default_rng(50); rm = []; sl = []
for it in range(2000):
    fs = np.clip(1 + es * rng.standard_normal(len(es)), .05, None); dd = 10 ** (0.2 * ed * rng.standard_normal(len(es)))
    fr = np.clip(1 + er * rng.standard_normal(len(es)), .1, None) * dd; fm = 10 ** (-0.4 * em * rng.standard_normal(len(es))) * dd ** 2
    b = dict(a); b['re'] = a['re'] * fr; b['r'] = a['r'] * fr; b['ms'] = a['ms'] * fm; b['mp'] = a['mp'] * fm; b['gb'] = a['gb'] * fm / fr ** 2
    e = np.log10(g0 * fs ** 2 / fr / H.gpred(M.MODEL, th, b)); rm.append(np.sqrt(np.mean(e ** 2))); sl.append(slope(e, b))
place('dwarfs', rm, {'formula': 0.370, 'MOND': 0.633, 'standard': 0.440})
s_obs = slope(e_obs, a); sl = np.array(sl); out['dwarfs_slope_logRe_fixed_Ms'] = dict(observed=round(s_obs, 3), noise_median=round(float(np.median(sl)), 3), noise95=np.round(np.percentile(sl, [2.5, 97.5]), 3).tolist(),
                                                                                    p_two_sided=round(float(2 * min(np.mean(sl <= s_obs), np.mean(sl >= s_obs))), 4))
print('dwarf slope', out['dwarfs_slope_logRe_fixed_Ms'])
# MaNGA and CLASH (as noise50_clusters.py)
D = H.PKG / 'data/Clusters'; E = {}
for fn in ('manga_bcg_standardized.csv', 'clash_standardized.csv'):
    for r in csv.DictReader(open(D / fn)): E[(r['canonical_cluster_name'] if 'clash' in fn else r['system'], round(float(r['radius_kpc']), 1))] = (float(r['loggbar_err_dex']), float(r['loggobs_err_dex']))
cl = H.cluster_rows(); rng = np.random.default_rng(51)
COMP = {'MaNGA': {'formula': 0.116, 'MOND': 0.186, 'standard': 0.116}, 'CLASH_centre': {'formula': 0.102, 'MOND': 0.460, 'standard': 0.102}, 'CLASH_outer': {'formula': 0.100, 'MOND': 0.526, 'standard': 0.101}}
for pop, kind in (('MaNGA', 'MaNGA_BCG'), ('CLASH_centre', 'CLASH_centre'), ('CLASH_outer', 'CLASH_outer')):
    rows = [r for r in cl if r['kind'] == kind]; a = H.arr(rows)
    def err(r):
        k = (r['name'], round(float(r['r']), 1))
        if k in E: return E[k]
        c = [v for (n, rad), v in E.items() if (n == r['name'] or n.endswith(r['name']) or r['name'].endswith(n)) and abs(rad - r['r']) < .06 * r['r'] + .2]
        return c[0] if c else (np.nan, np.nan)
    eb, eo = np.array([err(r) for r in rows]).T; eb = np.nan_to_num(eb, nan=np.nanmedian(eb)); eo = np.nan_to_num(eo, nan=np.nanmedian(eo))
    g0 = H.gpred(M.MODEL, th, a); hs = np.unique(a['name']); e_obs = np.log10(a['go'] / g0)
    def rms(e): return np.sqrt(np.mean([np.mean(e[a['name'] == h] ** 2) for h in hs]))
    sims = []; means = []
    for it in range(2000):
        db = eb * rng.standard_normal(len(eb)); do = eo * rng.standard_normal(len(eo)); b = dict(a); b['gb'] = a['gb'] * 10 ** db
        if pop == 'MaNGA': b['ms'] = a['ms'] * 10 ** db; b['mp'] = a['mp'] * 10 ** db
        e = np.log10(g0 * 10 ** do / H.gpred(M.MODEL, th, b)); sims.append(rms(e)); means.append(np.mean([np.mean(e[a['name'] == h]) for h in hs]))
    place(pop, sims, COMP[pop])
    mo = float(np.mean([np.mean(e_obs[a['name'] == h]) for h in hs])); means = np.array(means)
    out[pop]['mean_offset'] = dict(observed=round(mo, 3), noise95=np.round(np.percentile(means, [2.5, 97.5]), 3).tolist(), p_two_sided=round(float(2 * min(np.mean(means <= mo), np.mean(means >= mo))), 4))
    print('   mean offset', out[pop]['mean_offset'])
# early types, with distance errors on the formula's inputs as well
meta = {}
for line in (H.AST / 'inputs/data/ETG/ETG_Lelli2017.mrt').read_text().splitlines():
    if line and not line.startswith('#'):
        t = line.split(); meta[t[0]] = np.array(t[1:], float)
rows = H.etg_rows(); a = H.arr(rows); hs = np.unique(a['name']); gi = np.array([list(hs).index(n) for n in a['name']])
eo = np.array([meta[r['name']][12 + 2 * (r['kind'] == 'ETG_outer')] for r in rows]); eb = np.array([meta[r['name']][16 + 2 * (r['kind'] == 'ETG_outer')] for r in rows])
eD = np.array([meta[h][1] / meta[h][0] for h in hs]); g0 = H.gpred(M.MODEL, th, a)
def rms(e): return float(np.sqrt(np.mean([np.mean(e[a['name'] == h] ** 2) for h in hs])))
rng = np.random.default_rng(542)
for tag, useD in (('ETG', False), ('ETG_with_distance_on_inputs', True)):
    sims = []
    for it in range(2000):
        zb = rng.standard_normal(len(hs))[gi]; db = eb * zb; do = eo * rng.standard_normal(len(eo)); fD = np.clip(1 + eD * rng.standard_normal(len(hs)), .3, None)[gi] if useD else 1.
        b = dict(a); b['gb'] = a['gb'] * 10 ** db; b['ms'] = a['ms'] * 10 ** db * fD ** 2; b['mp'] = a['mp'] * 10 ** db * fD ** 2; b['re'] = a['re'] * fD; b['r'] = a['r'] * fD
        sims.append(rms(np.log10(g0 * 10 ** do / H.gpred(M.MODEL, th, b))))
    place(tag, sims, {'formula': 0.185, 'MOND': 0.182, 'standard': 0.182})
json.dump(out, open('res54_rest.json', 'w'), indent=1)
