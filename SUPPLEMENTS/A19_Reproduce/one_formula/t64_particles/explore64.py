"""T64 all-data look (not a test; A, B, C frozen at the all-disk values): what the 149 benchmark disks' points, with T61's
weights, ask of a fuzzy particle mass and of a warm half-mode mass, and what each value does to every population.
Usage: PYTHONPATH=../dmfit python3 explore64.py"""
import json, numpy as np
import harness_T64 as M
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B; M61 = M.M61
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
gs = H.disks(); cl = H.cluster_rows(); names = [g['name'] for g in gs]
A = M.pick(M.ALL, np.isin(M.ALL['name'], names)); nw = M.weights(A)
pops = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
        'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
def hb(er, nm, k=None):
    k = np.ones(len(er), bool) if k is None else k; return float(np.sqrt(np.mean([np.mean(er[k & (nm == n)] ** 2) for n in np.unique(nm[k])])))
def row(model):
    e = np.log10(A['go'] / H.gpred(model, TH, A)); o = dict(loss=float((nw * e) @ (nw * e)), kept=hb(e, A['name'], A['fe'] <= .1), dropped=hb(e, A['name'], A['fe'] > .1), all=hb(e, A['name']))
    for p, a in pops.items(): o[p] = hb(np.log10(a['go'] / H.gpred(model, TH, a)), a['name'])
    return o
out = {}; M61.U['t'], M61.U['k'] = 0.1501, 0.6229
out['T48'] = row(M8.MODEL); out['T61'] = row('MOGF0IESBC')
e = np.log10(A['go'] / C.g_mond(A)); out['MOND'] = dict(loss=float((nw * e) @ (nw * e)), kept=hb(e, A['name'], A['fe'] <= .1), dropped=hb(e, A['name'], A['fe'] > .1), all=hb(e, A['name']))
cols = ['loss', 'kept', 'dropped', 'all', 'dwarfs', 'ETG', 'MaNGA', 'CLASH_centre', 'CLASH_outer']
print('%-22s' % 'model', ' '.join('%12s' % c for c in cols))
for k in ('T48', 'T61', 'MOND'): print('%-22s' % k, ' '.join('%12.4f' % out[k].get(c, np.nan) for c in cols))
for model, key, grid in (('MOGF0IESBF', 'lm22', np.arange(-1., 6.01, 0.25)), ('MOGF0IESBW', 'lhm', np.arange(3., 12.51, 0.5))):
    for x in grid:
        M.V[key] = float(x); o = row(model); out['%s %s=%.2f' % (model, key, x)] = o
        print('%-22s' % ('%s %.2f' % (key, x)), ' '.join('%12.4f' % o[c] for c in cols), flush=True)
json.dump(out, open('res64_explore.json', 'w'), indent=1)
