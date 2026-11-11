"""All-data refit (primary mask) of one model from a harness module: coefficients, all-data scores, lensing chi2.
Usage: python3 describe24.py <module> <model>"""
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
f = H.fit(model, .1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
er = np.log10(ext['MaNGA']['go'] / H.gpred(model, th, ext['MaNGA']))
print(json.dumps(dict(module=sys.argv[1], model=model, th=np.round(th, 4).tolist(), bound_hits=f['bound_hits'], disk=score(th, d),
      **{k: score(th, a) for k, a in ext.items()}, MaNGA_mean=round(float(np.mean(er)), 4),
      lens=round(float(np.sum(lens.white(model, th, [0, 1, 2, 3]) ** 2)), 2))), flush=True)
