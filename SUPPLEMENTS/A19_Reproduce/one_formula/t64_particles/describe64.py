"""T64 descriptive, all data, no marks (PROTOCOL_T64, "reported without marks"); A, B, C frozen at the all-disk values.
 1. The fuzzy model with the soliton joined at 2, 3 and 4 core radii: the disks' weighted loss over the mass grid.
 2. The value each population alone would choose for each particle number, and its error there against the cold end.
 3. The warm number as a thermal relic mass.
 4. Where the formula without a core (T48) errs at small radius, in halos of like mass: gas-rich diffuse disks (gate w
    above 0.3), other disks, and dwarf spheroidals. Mean of log10(g_obs / g_pred): negative = the formula predicts too much.
Usage: PYTHONPATH=../dmfit python3 describe64.py"""
import json, numpy as np
import harness_T64 as M
import harness_T48 as M8
H = M.H; B = M.B; M61 = M.M61; S = M.S
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
gs = H.disks(); cl = H.cluster_rows(); names = [g['name'] for g in gs]
A = M.pick(M.ALL, np.isin(M.ALL['name'], names)); nw = M.weights(A)
pops = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
def hb(er, nm): return float(np.sqrt(np.mean([np.mean(er[nm == n] ** 2) for n in np.unique(nm)])))
out = {}
# 1. transition radius
g = np.arange(-1., 6.001, 0.125); tr = {}
for rt in (2., 3., 4.):
    M.RT = rt; L = []
    for x in g: M.V['lm22'] = float(x); e = nw * np.log10(A['go'] / H.gpred('MOGF0IESBF', TH, A)); L.append(float(e @ e))
    j = int(np.argmin(L)); tr[str(rt)] = dict(best_lm22=float(g[j]), loss=L[j], loss_cold=L[-1], loss_at_0=L[list(np.round(g, 3)).index(0.)], loss_at_1=L[list(np.round(g, 3)).index(1.)], loss_at_2=L[list(np.round(g, 3)).index(2.)])
    print('join at %.0f r_c: lowest loss %.5f at lm22 = %.3f; cold end %.5f; at 1e-22, 1e-21, 1e-20 eV: %.5f %.5f %.5f' % (rt, L[j], g[j], L[-1], tr[str(rt)]['loss_at_0'], tr[str(rt)]['loss_at_1'], tr[str(rt)]['loss_at_2']))
M.RT = 3.; out['fuzzy_join_radius'] = tr
# 2. each population alone
pref = {}
for model, key, grid, cold in (('MOGF0IESBF', 'lm22', np.arange(-1., 6.001, 0.125), 6.), ('MOGF0IESBW', 'lhm', np.arange(3., 12.501, 0.125), 3.)):
    for p, a in dict(disks_weighted=A, **pops).items():
        v = []
        for x in grid:
            M.V[key] = float(x); e = np.log10(a['go'] / H.gpred(model, TH, a)); v.append(float(np.sqrt((nw * e) @ (nw * e))) if p == 'disks_weighted' else hb(e, a['name']))
        j = int(np.argmin(v)); M.V[key] = cold; pref[model + ' ' + p] = dict(best=float(grid[j]), err=round(v[j], 4), err_cold=round(v[list(np.round(grid, 3)).index(cold)], 4))
        print('%s %-15s alone chooses %s = %6.3f: %.4f against %.4f at the cold end' % (model, p, key, grid[j], v[j], pref[model + ' ' + p]['err_cold']))
out['each_population_alone'] = pref; M.V['lm22'] = 6.; M.V['lhm'] = 3.
# 3. warm number as a thermal relic mass (Schneider et al. 2012, Table 1 and eq. 5: M_hm = 1.3e10 / h Msun at 1 keV, proportional to m^-3.33)
h = H.D.H0 / 100.
def kev(lhm): return float((1.3e10 / h / 10 ** lhm) ** (1 / 3.33))
fit = json.load(open('res64_alldata_fit.json'))['MOGF0IESBW']['par']; ho = json.load(open('res64_heldout.json'))['fits']['lhm']
out['warm_relic'] = dict(h=h, lhm_all_data=fit, keV_all_data=kev(fit), lhm_partitions=[min(ho), max(ho)], keV_partitions=[kev(max(ho)), kev(min(ho))])
print('warm: all-data log M_hm %.3f = %.2f keV; partitions %.2f to %.2f = %.2f to %.2f keV (h = %.4f)' % (fit, kev(fit), min(ho), max(ho), kev(max(ho)), kev(min(ho)), h))
# 4. inner residuals of T48 in halos of like mass
def inner(a, kind):
    m, c, rv = M.halo(TH, a); e = np.log10(a['go'] / H.gpred(M8.MODEL, TH, a)); w = M61.gl(a['ms'], a['mp'], a['re'], S.thresholds(a)); return m, e, w, a['r'] / a['re']
rows = {}
mA, eA, wA, xA = inner(A, 'disk'); d = pops['dwarfs']; mD, eD, wD, xD = inner(d, 'dwarfs')
for lab, lo, hi in (('halo 10^9.7 to 10^10.3', 9.7, 10.3), ('halo 10^10.3 to 10^11', 10.3, 11.), ('halo below 10^9.7', 0., 9.7)):
    r = {}
    for nm, k in (('gas-rich diffuse disks, r below 1.75 R_e', (wA > .3) & (xA < 1.75)), ('other disks, r below 1.75 R_e', (wA <= .3) & (xA < 1.75)), ('gas-rich diffuse disks, r above 1.75 R_e', (wA > .3) & (xA >= 1.75))):
        kk = k & (np.log10(mA) >= lo) & (np.log10(mA) < hi)
        if kk.sum() >= 5:
            gm = [np.mean(eA[kk & (A['name'] == n)]) for n in np.unique(A['name'][kk])]; r[nm] = dict(points=int(kk.sum()), galaxies=len(gm), mean=round(float(np.mean(gm)), 3), sem=round(float(np.std(gm, ddof=1) / np.sqrt(len(gm))), 3) if len(gm) > 1 else None)
    kk = (np.log10(mD) >= lo) & (np.log10(mD) < hi)
    if kk.sum() >= 3: r['dwarf spheroidals (one point each, near the half-light radius)'] = dict(points=int(kk.sum()), galaxies=int(kk.sum()), mean=round(float(np.mean(eD[kk])), 3), sem=round(float(np.std(eD[kk], ddof=1) / np.sqrt(kk.sum())), 3), median_r_over_re=round(float(np.median(xD[kk])), 2))
    rows[lab] = r; print(lab); [print('   %-62s %s' % (k_, json.dumps(v_))) for k_, v_ in r.items()]
out['inner_residuals_T48'] = rows
json.dump(out, open('res64_describe.json', 'w'), indent=1)
