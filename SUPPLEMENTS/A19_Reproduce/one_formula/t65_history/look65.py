"""T65, after the protocol's runs: an all-data look for a possible next round (not a test, no protocol; A, B, C frozen).
The published history recipe (Read, Agertz and Collins 2016: n = tanh(0.04 t_SF / t_dyn), r_c = 1.75 R_e, t_SF = 13.8 Gyr
for a disk still forming stars, no mass term, nothing fitted) on every rotation-curve point of the 149 benchmark disks:
(a) in every disk; (b) only in diffuse disks (L above 0); (c) as (b) with the strength also scaled by L capped at 1.
Usage: PYTHONPATH=../dmfit:../t59 python3 look65.py"""
import json, numpy as np
import lib65 as L
import harness_T48 as M8
import compare46 as C
H = L.H; M = L.M; M5 = L.M5; M61 = L.M61; B = L.B
names = [g['name'] for g in H.disks()]; A = M.pick(M.ALL, np.isin(M.ALL['name'], names)); nw = M.weights(A); B.FROZEN['abc'] = L.ABC
def hb(e, k=None):
    k = np.ones(len(e), bool) if k is None else k; return float(np.sqrt(np.mean([np.mean(e[k & (A['name'] == n)] ** 2) for n in np.unique(A['name'][k])])))
def row(g): e = np.log10(A['go'] / g); return dict(loss=round(float((nw * e) @ (nw * e)), 5), kept=round(hb(e, A['fe'] <= .1), 4), dropped=round(hb(e, A['fe'] > .1), 4), all=round(hb(e), 4))
h = L.halo(A); n0 = L.n_read(h, L.TU); Ld = M5.diffuse(A['ms'], A['re']); rc = L.ETA * A['re']
def g_core(n): return A['gb'] + H.G * L.core_m(h, A['r'], n, rc) / A['r'] ** 2 * H.ACC
out = {}; M61.U['t'], M61.U['k'] = L.T_G, L.K_G; fit = json.load(open('res65_alldata_fit.json')); M5.U['t'], M5.U['k'] = fit['t'], fit['kappa']
out['no core (T48)'] = row(H.gpred(M8.MODEL, L.TH, A)); out['formula, gas gate (T61)'] = row(H.gpred('MOGF0IESBC', L.TH, A)); out['formula, any-gas gate (T65)'] = row(H.gpred(M5.MODEL, L.TH, A)); out['MOND'] = row(C.g_mond(A))
out['(a) published recipe in every disk'] = row(g_core(n0)); out['(b) published recipe in diffuse disks only'] = row(g_core(np.where(Ld > 0, n0, 0.)))
out['(c) published recipe times min(L, 1)'] = row(g_core(n0 * np.minimum(Ld, 1.)))
un, i = np.unique(A['name'], return_index=True); out['n_published_percentiles_5_50_95'] = np.round(np.percentile(n0[i], [5, 50, 95]), 2).tolist()
for k, v in out.items(): print('%-44s %s' % (k, json.dumps(v)))
json.dump(out, open('res65_look.json', 'w'), indent=1)
