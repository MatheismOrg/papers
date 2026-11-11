"""All-data descriptive fits of the T17b switch formula (coefficients, dwarf upper limits; not validation scores)."""
import json, sys
import numpy as np
import harness_T17b as Y
H = Y.H; B = Y.B
def limit_rows(sub=True):
    bt = H.PKG / 'blind_test'; rules = json.load(open(bt / 'frozen_disk_rules.json'))['frozen_on_all_disks']
    rr, src = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(bt))
    aa = H.V4.arrays(rr, rules, ml=2., major=False); pp = H.V3.predict(rr, rules, ml=2., major=False)
    out = [dict(name=x['key'], kind='limit', r=pp['r'][i], gb=pp['gb'][i], go=10 ** aa['lo'][i], ms=pp['Ms'][i], mp=pp['Ms'][i],
                re=pp['Re'][i], z=0., fac=float(pp['fac'][i]) if sub else 1.) for i, x in enumerate(rr) if aa['ul'][i]]
    assert len(out) == 10; return H.arr(out)
def score(model, th, a):
    l = H.losses(model, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
if __name__ == '__main__':
    out = sys.argv[1]; configs = sys.argv[2:]
    gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    lim, lim0 = limit_rows(True), limit_rows(False); res = {}
    for drop in (False, True):
        lens = H.Lens(gs, drop)
        B.FROZEN['abc'] = np.array(H.fit('ORIG', None, d, ext, lens, [0, 1, 2])['th'])
        for cfg in configs:
            model, w = cfg.split('@'); f = H.fit(model, float(w), d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
            key = f"{cfg}|{'strict' if drop else 'primary'}"
            res[key] = dict(frozen_abc=[round(x, 6) for x in B.FROZEN['abc']], coefficients=dict(zip(H.MODELS[model][0], [round(x, 6) for x in th])),
                            bound_hits=f['bound_hits'], limits_exceeded=int(np.sum(H.gpred(model, th, lim) > lim['go'])),
                            limits_exceeded_without_subhalo=int(np.sum(H.gpred(model, th, lim0) > lim0['go'])),
                            all_data_scores=dict(disk=score(model, th, d), **{k: score(model, th, a) for k, a in ext.items()}),
                            lensing_chi2_all=round(float(np.sum(lens.white(model, th, [0, 1, 2, 3]) ** 2)), 3),
                            switch_on_share={k: round(float(np.mean(Y.switch(model, th, a['ms'], a['re']) > .5)), 3) for k, a in [('disk', d)] + list(ext.items())})
            print(key, res[key], flush=True)
    json.dump(res, open(out, 'w'), indent=1)
