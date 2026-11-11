"""T45 sanity check (all data): harness_T45 reproduces explore45's calibrated M_t and scores."""
import numpy as np, json
import harness_T45 as M
H = M.H; B = M.B
assert H is __import__('harness_T15_hashed')
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
f = H.fit(M.MODEL, .1, None, None, None, [0, 1, 2, 3]); th = np.array(f['th']); print(json.dumps(f))
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
sc = {'disk': H.losses(M.MODEL, th, d)}; sc.update({k: H.losses(M.MODEL, th, a) for k, a in ext.items()})
out = {k: round(float(np.sqrt(np.mean(list(v.values())))), 4) for k, v in sc.items()}
out['lens'] = round(float(np.sum(H.Lens(gs, False).white(M.MODEL, th, [0, 1, 2, 3]) ** 2)), 2)
out['lens_strict'] = round(float(np.sum(H.Lens(gs, True).white(M.MODEL, th, [0, 1, 2, 3]) ** 2)), 2)
dev = M.red_dev(th[0]); out['red_rms'] = round(float(np.sqrt(np.mean(np.square(dev)))), 4); out['red_dev'] = np.round(dev, 3).tolist()
print(json.dumps(out))
