"""T64 all-data fit (not a test): the fuzzy and warm numbers on all points of the 149 benchmark disks, A, B, C frozen."""
import json, numpy as np
import harness_T64 as M
H = M.H; names = [g['name'] for g in H.disks()]; abc = np.array([11.7481, .8178, -.7119]); out = {}
for model in M.MODELS:
    x, loss, att, bh, npts, cold = M.fit_disks(model, names, abc)
    out[model] = dict(par=x, loss=loss, loss_cold=cold, bound_hits=bh, points=npts, attempts=[(round(a['th'][0], 4), round(a['loss'], 6)) for a in att])
    print(model, json.dumps(out[model]))
json.dump(out, open('res64_alldata_fit.json', 'w'), indent=1)
