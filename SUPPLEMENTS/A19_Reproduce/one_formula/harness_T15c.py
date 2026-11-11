"""T15c: the T15/T15b formulas with the redshift-dependent Dutton and Maccio (2014) concentration at each halo's own z.
Imports harness_T15b (and through it harness_T15_hashed) unchanged. Usage: python3 harness_T15c.py <primary|strict> <out.json> cfg..."""
import sys
import numpy as np
import harness_T15b as B
H = B.H
BASE = {'ORIGZ': 'ORIG', 'KNEESFZ': 'KNEESF', 'KNEESZ': 'KNEES', 'KNEESEZ': 'KNEESE'}
for z, b in BASE.items(): H.MODELS[z] = H.MODELS[b]
_gpred, _esd, _mref = H.gpred, H.esd, H.mref
def dc_of(model, th, ms, re):
    """The trait shift in log c already applied by the base model (none for these four)."""
    return 0.
def cz(m, z): return H.K.c_dm14(m, z)
def gpred(model, th, a):
    if model not in BASE: return _gpred(model, th, a)
    b = BASE[model]; m, _ = _mref(b, th, a['ms'], a['mp'], a['re']); c = cz(m, a['z'])
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    g = a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
    if b == 'KNEESE':
        q = 10 ** th[7] * rv; g = g + H.G * th[6] * (1 - H.FB) * m * a['r'] / (a['r'] ** 2 + q ** 2) ** 1.5 * H.ACC
    return g
def esd(model, th, st):
    if model not in BASE: return _esd(model, th, st)
    b = BASE[model]; m, _ = _mref(b, th, st.ms, st.mp, st.re); c = cz(m, H.K.Z0)
    rc = H.K.rho_crit(H.K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c
    rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * H.K.g_wb(st.R / rs[:, None, None]) / 1e12
    if b == 'KNEESE':
        q = 10 ** th[7] * rv; e = e + th[6] * (1 - H.FB) * m[:, None, None] * st.R ** 2 / (np.pi * (st.R ** 2 + q[:, None, None] ** 2) ** 2) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
import json, csv, time, hashlib
from pathlib import Path
def run(mask, out, configs):
    """As harness_T15b.run, plus ORIGZ@disk: the paper's rule fitted on training disks only (the frozen A, B, C),
    evaluated with the redshift-dependent concentration."""
    gs = H.disks(); folds = list(csv.DictReader(open(H.PKG / 'data/SPARC/A19_Disk_Folds.csv')))
    maps = json.load(open(H.AST / 'T56_fold_maps.json'))['external_fold_maps']; cl = H.cluster_rows()
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    recs = []; t0 = time.monotonic()
    for asg in sorted(maps):
        for fold in range(5):
            held = {r['galaxy'] for r in folds if r[asg] == str(fold)}
            tr = [g for g in gs if g['name'] not in held]; te = [g for g in gs if g['name'] in held]
            dtr, dte = H.arr(H.disk_rows(tr)), H.arr(H.disk_rows(te)); etr, ete = {}, {}
            for p, a in ext.items():
                k = np.array([maps[asg][p][n] == fold for n in a['name']]); etr[p] = H.sub(a, ~k); ete[p] = H.sub(a, k)
            lens = H.Lens(tr, mask == 'strict')
            B.FROZEN['abc'] = np.array(H.fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])
            for j in range(4):
                groups = [i for i in range(4) if i != j]
                for cfg in configs:
                    model, w = cfg.split('@')
                    if cfg == 'ORIGZ@disk': f = dict(th=B.FROZEN['abc'].tolist(), attempts=[], bound_hits=[])
                    else: f = H.fit(model, float(w), dtr, etr, lens, groups)
                    th = np.array(f['th'])
                    if model == 'KNEESFZ': f['frozen_abc'] = B.FROZEN['abc'].tolist()
                    l = {'disk': H.losses(model, th, dte)}
                    for p, a in ete.items(): l[p] = H.losses(model, th, a)
                    for kind in ['CLASH_centre', 'CLASH_outer', 'ETG_inner', 'ETG_outer']:
                        a = ete['CLASH' if kind.startswith('CLASH') else 'ETG']; l[kind] = H.losses(model, th, H.sub(a, a['kind'] == kind))
                    z = lens.white(model, th, [j])
                    recs.append(dict(config=cfg, assignment=asg, fold=fold, lens=j, fit=f, losses=l, test_lensing_chi2=float(z @ z)))
            print(mask, asg, fold, round(time.monotonic() - t0, 1), flush=True)
    res = dict(mask=mask, configs=configs, summary=H.summarize(recs), records=recs,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    json.dump(res, open(out, 'w')); print(json.dumps(res['summary'], indent=1))
if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3:])
