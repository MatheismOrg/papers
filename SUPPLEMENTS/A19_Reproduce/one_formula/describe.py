"""All-data descriptive fits of the T15 formulas (for coefficients and dwarf upper limits; not validation scores)."""
import sys, json
import numpy as np
import harness_T15_hashed as H

def limit_rows():
    bt = H.PKG / 'blind_test'; rules = json.load(open(bt / 'frozen_disk_rules.json'))['frozen_on_all_disks']
    rr, src = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(bt))
    aa = H.V4.arrays(rr, rules, ml=2., major=False); pp = H.V3.predict(rr, rules, ml=2., major=False)
    out = [dict(name=x['key'], kind='limit', r=pp['r'][i], gb=pp['gb'][i], go=10 ** aa['lo'][i], ms=pp['Ms'][i], mp=pp['Ms'][i],
                re=pp['Re'][i], z=0.) for i, x in enumerate(rr) if aa['ul'][i]]
    assert len(out) == 10; return H.arr(out)

def main(configs, out):
    gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    lim = limit_rows(); res = {}
    for drop in (False, True):
        lens = H.Lens(gs, drop)
        for cfg in configs:
            model, w = cfg.split('@'); w = None if w == 'disk' else float(w)
            f = H.fit(model, w, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
            gl = H.gpred(model, th, lim)
            res[f"{cfg}|{'strict' if drop else 'primary'}"] = dict(
                coefficients=dict(zip(H.MODELS[model][0], [round(x, 6) for x in th])), bound_hits=f['bound_hits'],
                limits_exceeded=int(np.sum(gl > lim['go'])), lensing_chi2_all=round(float(np.sum(lens.white(model, th, [0, 1, 2, 3]) ** 2)), 3))
            print(cfg, 'strict' if drop else 'primary', res[f"{cfg}|{'strict' if drop else 'primary'}"], flush=True)
    json.dump(res, open(out, 'w'), indent=1)

if __name__ == '__main__':
    main(sys.argv[2:], sys.argv[1])
