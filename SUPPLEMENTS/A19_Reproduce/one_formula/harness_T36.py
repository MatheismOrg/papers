"""T36: the frozen T26 formula with halo contraction that grows with stellar mass. Contraction as Cautun et al. (2020)
eq. 11, contraction branch only (as T19): M_DM(<r) = (1 - fb) M_NFW(<r) max(F(eta), 1)^s, F = 0.45 + 0.38 (eta + 1.16)^0.53,
eta = M_bar(<r)/(fb M_NFW(<r)), with the strength depending on the galaxy's stellar mass:
  s(M*) = 1 / (1 + (M_s/M*)^2),
1 (the published strength) for galaxies well above M_s, 0 well below. The paper's rule A, B, C and log10 M_s are refitted on
each fold's training disks with this contraction and without the switch (as T18); then log10 M_t is fitted jointly as in
T26. Baryons, gate, floor and IMF factors as T26 (harness_T26, imported first). Dynamics: M_bar(<r) = g_bar r^2/G; lensing:
lens baryons as a point mass (harness_T18.esd_dm).
MOGF0IEC: 5 coefficients (A, B, C, log10 M_s on disks; log10 M_t joint).
Usage: python3 harness_T36.py run <mask> <parts.jsonl> <budget s> [k n]   (resumable, parts split as k of n)
       python3 harness_T36.py merge <mask> <res.json> <parts.jsonl> ...
       python3 harness_T36.py describe                                       (all-data check)"""
import sys, json, csv, time, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import harness_T26 as M
F = M.F; Z = F.Z; H = M.H; B = Z.B; K = H.K; T = Z.T
LMS = {'v': 10.5}
MODEL = 'MOGF0IEC'
H.MODELS[MODEL] = H.MODELS['MOGF0']
def s_of(ms): return 1. / (1. + (10 ** LMS['v'] / np.asarray(ms, float)) ** 2)
def mass(th, ms, mp, re, z): return F.mass('MOGF0', th, ms, mp, re, z)
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model != MODEL: return _gpred(model, th, a)
    m = mass(th, a['ms'], a['mp'], a['re'], a['z']); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * Z.Fc(eta) ** s_of(a['ms']) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MODEL: return _esd(model, th, st)
    m = mass(th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, s_of(st.ms)[:, None, None, None])
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
OFF = np.array([12.5])                              # switch off for the disk fit (as T18: without the switch)
def disk_fit(dtr):
    def fun(p):
        B.FROZEN['abc'] = p[:3]; LMS['v'] = p[3]
        return dtr['hw'] * np.log10(dtr['go'] / gpred(MODEL, OFF, dtr)) / .15
    lo, hi, st = [9, -.5, -2.5, 8.], [14.5, 2, 1.5, 13.], [11.7481, .8178, -.7119, 10.5]
    best = None
    for s0 in (10.5, 9.5, 11.5):
        o = least_squares(fun, st[:3] + [s0], bounds=(lo, hi), max_nfev=400, ftol=1e-9, xtol=1e-9, gtol=1e-9)
        if best is None or o.cost < best.cost: best = o
    B.FROZEN['abc'] = best.x[:3]; LMS['v'] = float(best.x[3])
    return dict(abc=best.x[:3].tolist(), log_Ms=LMS['v'], success=bool(best.success), cost=float(best.cost),
                log_Ms_at_bound=bool(min(best.x[3] - 8., 13. - best.x[3]) < 1e-5))
def data():
    gs = H.disks(); cl = H.cluster_rows()
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    return gs, ext
def run(mask, out, budget, kk=0, nn=1):
    t0 = time.monotonic(); gs, ext = data()
    folds = list(csv.DictReader(open(H.PKG / 'data/SPARC/A19_Disk_Folds.csv'))); maps = json.load(open(H.AST / 'T56_fold_maps.json'))['external_fold_maps']
    done = set()
    if Path(out).exists(): done = {(json.loads(l)['assignment'], json.loads(l)['fold']) for l in open(out)}
    i = -1
    for asg in sorted(maps):
        for fold in range(5):
            i += 1
            if i % nn != kk or (asg, fold) in done: continue
            if time.monotonic() - t0 > budget: return
            held = {r['galaxy'] for r in folds if r[asg] == str(fold)}
            tr = [g for g in gs if g['name'] not in held]; te = [g for g in gs if g['name'] in held]
            dtr, dte = H.arr(H.disk_rows(tr)), H.arr(H.disk_rows(te)); etr, ete = {}, {}
            for p, a in ext.items():
                sel = np.array([maps[asg][p][n] == fold for n in a['name']]); etr[p] = H.sub(a, ~sel); ete[p] = H.sub(a, sel)
            lens = H.Lens(tr, mask == 'strict'); dfit = disk_fit(dtr); recs = []
            for j in range(4):
                groups = [g for g in range(4) if g != j]
                f = H.fit(MODEL, .1, dtr, etr, lens, groups); th = np.array(f['th']); f['disk_fit'] = dfit
                l = {'disk': H.losses(MODEL, th, dte)}
                for p, a in ete.items(): l[p] = H.losses(MODEL, th, a)
                for kind in ['CLASH_centre', 'CLASH_outer', 'ETG_inner', 'ETG_outer']:
                    a = ete['CLASH' if kind.startswith('CLASH') else 'ETG']; l[kind] = H.losses(MODEL, th, H.sub(a, a['kind'] == kind))
                z = lens.white(MODEL, th, [j])
                recs.append(dict(config=MODEL + '@0.1', assignment=asg, fold=fold, lens=j, fit=f, losses=l, test_lensing_chi2=float(z @ z)))
            with open(out, 'a') as fh: fh.write(json.dumps(dict(assignment=asg, fold=fold, records=recs)) + '\n')
            print(mask, asg, fold, round(time.monotonic() - t0, 1), flush=True)
def merge(mask, res, *outs):
    recs = [r for out in outs for l in open(out) for r in json.loads(l)['records']]
    assert len(recs) == 60 and len({(r['assignment'], r['fold'], r['lens']) for r in recs}) == 60, len(recs)
    json.dump(dict(mask=mask, configs=[MODEL + '@0.1'], summary=H.summarize(recs), records=recs,
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), open(res, 'w'))
    print(json.dumps(H.summarize(recs), indent=1))
def describe():
    gs, ext = data(); d = H.arr(H.disk_rows(gs)); dfit = disk_fit(d)
    def score(th, a):
        l = H.losses(MODEL, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
    out = dict(disk_fit=dfit)
    for drop in (False, True):
        lens = H.Lens(gs, drop); f = H.fit(MODEL, .1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
        out['strict' if drop else 'primary'] = dict(th=np.round(th, 4).tolist(), disk=score(th, d), **{k: score(th, a) for k, a in ext.items()},
                                                    lens=round(float(np.sum(lens.white(MODEL, th, [0, 1, 2, 3]) ** 2)), 2))
    print(json.dumps(out), flush=True); return out
if __name__ == '__main__':
    if sys.argv[1] == 'run': run(sys.argv[2], sys.argv[3], float(sys.argv[4]), *(int(x) for x in sys.argv[5:7]))
    elif sys.argv[1] == 'merge': merge(sys.argv[2], sys.argv[3], *sys.argv[4:])
    else: describe()
