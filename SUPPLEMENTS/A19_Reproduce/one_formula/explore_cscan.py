"""Exploratory all-data scan (no protocol; a look before any hashed round): T43's formula with every halo concentration
multiplied by 10^dc (rule's disk fit and all populations, dynamics and lensing), A, B, C refitted on all disks, then M_t.
Usage: python3 explore_cscan.py <dc>"""
import sys, json, warnings
import numpy as np
import harness_T43 as G
S = G.S; H = S.H; K = H.K; B = S.B
DC = float(sys.argv[1]); _c = K.c_dm14
K.c_dm14 = lambda m, z: _c(m, z) * 10 ** DC
_mref = H.mref
def mref(model, th, ms, mp, re):
    m, c = _mref(model, th, ms, mp, re); return (m, c * 10 ** DC) if model == 'ORIG' else (m, c)
H.mref = mref
def score(model, th, a):
    l = H.losses(model, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
out = dict(dc=DC)
for drop in (False, True):
    lens = H.Lens(gs, drop)
    abc = np.array(H.fit('ORIG', None, d, ext, lens, [0, 1, 2])['th']); B.FROZEN['abc'] = abc
    f = H.fit('MOGF0IESG', .1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
    out['strict' if drop else 'primary'] = dict(abc=np.round(abc, 4).tolist(), logMt=round(float(th[0]), 3), disk=score('MOGF0IESG', th, d),
        **{k: score('MOGF0IESG', th, a) for k, a in ext.items()}, lens=round(float(np.sum(lens.white('MOGF0IESG', th, [0, 1, 2, 3]) ** 2)), 2),
        lens_by_group=[round(float(np.sum(lens.white('MOGF0IESG', th, [j]) ** 2)), 1) for j in range(4)])
print(json.dumps(out), flush=True)
