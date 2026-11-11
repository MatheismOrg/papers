"""All-data check for T27: fit the final formula's free coefficient without MaNGA and CLASH (primary mask), then score
every population, MaNGA and CLASH as predictions. Usage: python3 describe27.py <harness module> <model>"""
import sys, json, importlib
import numpy as np
M = importlib.import_module(sys.argv[1]); model = sys.argv[2]
H = __import__('harness_T15_hashed'); B = __import__('harness_T15b')
def score(th, a):
    l = H.losses(model, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
lens = H.Lens(gs, False); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
train = {k: v for k, v in ext.items() if k not in ('MaNGA', 'CLASH')}
f = H.fit(model, .1, d, train, lens, [0, 1, 2, 3]); th = np.array(f['th'])
print(json.dumps(dict(module=sys.argv[1], model=model, trained_on='disks, dwarfs, ETG, lensing', th=np.round(th, 4).tolist(), bound_hits=f['bound_hits'],
      disk=score(th, d), **{k: score(th, a) for k, a in ext.items()}, lens=round(float(np.sum(lens.white(model, th, [0, 1, 2, 3]) ** 2)), 2))), flush=True)
