"""SA73 T61 all-data look (before PROTOCOL_T61; not a test): harness_T61's formula with A, B, C frozen at the all-data rule
and (t, kappa) fitted as the protocol will fit them, on all points of the 149 benchmark disks; every scale against the
final formula (T48), MOND's relation and standard halos. Usage: PYTHONPATH=<dmfit> python3 describe61.py"""
import json, numpy as np
import harness_T61 as M
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B
ABC = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); B.FROZEN['abc'] = ABC
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs)); LENS = {False: H.Lens(gs, False), True: H.Lens(gs, True)}
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
names = np.unique(d['name']).tolist(); (t, kap), loss, att, bh, npts = M.fit_disks('MOGF0IESBC', names, ABC)
print('all-data fit on', len(names), 'disks,', npts, 'points: t %.4f kappa %.4f' % (t, kap), 'bound hits', bh, 'attempts', [(np.round(x['th'], 4).tolist(), round(x['loss'], 6)) for x in att])
A = M.pick(M.ALL, np.isin(M.ALL['name'], names)); kept = A['fe'] <= 0.1; nw = M.weights(A)
def hb(e, k):
    un, gi = np.unique(A['name'][k], return_inverse=True); return round(float(np.sqrt(np.mean(np.bincount(gi, e[k] ** 2) / np.bincount(gi)))), 4)
def scores(model):
    B.FROZEN['abc'] = ABC
    if model == 'T61': M.U['t'], M.U['k'] = t, kap; g = H.gpred('MOGF0IESBC', th, A); mod = 'MOGF0IESBC'
    elif model == 'T48': g = H.gpred(M8.MODEL, th, A); mod = M8.MODEL
    elif model == 'MOND': g = C.g_mond(A); mod = None
    else: g = C.g_lcdm(A); mod = None
    e = np.log10(A['go'] / g); r = dict(kept=hb(e, kept), dropped=hb(e, ~kept), all=hb(e, np.ones(len(e), bool)), weighted_all=round(float(np.sqrt(np.sum((nw * e) ** 2))), 4), slow_dropped=hb(e, ~kept & (A['go'] * A['r'] / H.ACC < 900)))
    if mod:
        for p, a in ext.items(): r[p] = round(float(np.sqrt(np.mean(list(H.losses(mod, th, a).values())))), 4)
        for drop in (False, True): r['lens_' + ('strict' if drop else 'primary')] = round(float(sum(np.sum(LENS[drop].white(mod, th, [j]) ** 2) for j in range(4))), 2)
    return r
out = dict(t=t, kappa=kap, n_disks=len(names), n_points=npts, n_kept=int(kept.sum()), n_dropped=int((~kept).sum()), n_gal_dropped=int(len(np.unique(A['name'][~kept]))))
for m in ('T48', 'T61', 'MOND', 'STD'): out[m] = scores(m); print(m, out[m], flush=True)
json.dump(out, open('res61_describe.json', 'w'), indent=1)
