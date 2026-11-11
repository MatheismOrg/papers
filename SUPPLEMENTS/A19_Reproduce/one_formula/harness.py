"""A19 cross-scale test of few-coefficient connected halo formulas (session SA73, turn 15).
Scored on the SO73 T54-T63 benchmark: 149 SPARC disks, 75 measured LVDB dwarfs, 16 early types (32 points),
50 MaNGA central galaxies, 20 CLASH clusters (centre and outer points) and the four KiDS Fig-9 lensing stacks,
with the same 60 partitions (three disk fold assignments x five folds x four held-out lensing groups) and the
same host-balanced scores. Own code; data loaders follow the A19 package and SO73 cross_scale.py.
Usage: python3 harness.py run <mask: primary|strict> <out.json> [model@weight ...]"""
import sys, json, csv, zipfile, io, time, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular

S = Path('/tmp/claude-0/-home-claude/0a1f618d-0405-5db3-a0ab-84f1f9545369/scratchpad')
PKG = S / 'pkg/A19_Data_Package'; AST = S / 'so73_t63'
sys.path.insert(0, str(PKG / 'code'))
import A19_Disk_Fit_SPARC as D
import A19_KiDS_v9 as K
import A19_Combined_v4 as V4
import A19_New_Dwarfs_v3 as V3
import A19_Clusters_v5 as C5

G, ACC, h, FB = D.G, D.ACC, D.h, D.FB
def hz(z): return D.H0 / 1000 * np.sqrt(.3175 * (1 + np.asarray(z, float)) ** 3 + .6825)     # km/s/kpc
def ff(t): return np.log1p(t) - t / (1 + t)

# ---------------- data ----------------
def arr(rows):
    a = {k: np.array([r[k] for r in rows], float) for k in ['r', 'gb', 'go', 'ms', 'mp', 're', 'z']}
    a['name'] = np.array([r['name'] for r in rows]); a['kind'] = np.array([r['kind'] for r in rows])
    un, inv, ct = np.unique(a['name'], return_inverse=True, return_counts=True)
    a['hw'] = 1 / np.sqrt(len(un) * ct[inv])        # sum of squares = mean over hosts of each host's mean square
    return a
def sub(a, k):
    b = {n: v[k] for n, v in a.items() if n != 'hw'}
    un, inv, ct = np.unique(b['name'], return_inverse=True, return_counts=True); b['hw'] = 1 / np.sqrt(len(un) * ct[inv]); return b
def disks():
    gs = D.load(D.table1(PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'), PKG / 'data/SPARC/Rotmod_LTG')
    assert len(gs) == 149 and sum(len(g['r']) for g in gs) == 2700
    return gs
def disk_rows(gs):
    return [dict(name=g['name'], kind='disk', r=r, gb=gb, go=go, ms=g['Ms'], mp=g['Mp'], re=g['Re'], z=0.)
            for g in gs for r, gb, go in zip(g['r'], g['gbar'], g['gobs'])]
def dwarf_rows():
    bt = PKG / 'blind_test'; rules = json.load(open(bt / 'frozen_disk_rules.json'))['frozen_on_all_disks']
    rr, src = V4.load_rows(str(PKG / 'data/LVDB/dwarf_all.csv'), str(bt))
    aa = V4.arrays(rr, rules, ml=2., major=False); pp = V3.predict(rr, rules, ml=2., major=False)
    out = [dict(name=x['key'], kind='dwarfs', r=pp['r'][i], gb=pp['gb'][i], go=10 ** aa['lo'][i], ms=pp['Ms'][i],
                mp=pp['Ms'][i], re=pp['Re'][i], z=0.) for i, x in enumerate(rr) if not aa['ul'][i]]
    assert len(out) == 75; return out
def etg_rows():
    """Lelli et al. 2017 SPARC early types, two radii each, as in SO73 cross_scale.etg()."""
    meta = {}
    for line in (AST / 'inputs/data/ETG/ETG_Lelli2017.mrt').read_text().splitlines():
        if line and not line.startswith('#'):
            t = line.split(); meta[t[0]] = np.array(t[1:], float)
    rows = []
    with zipfile.ZipFile(AST / 'inputs/data/ETG/Rotmod_ETG.zip') as zf:
        for n in sorted(zf.namelist()):
            if not n.endswith('_rotmod.dat'): continue
            name = Path(n).name.replace('_rotmod.dat', ''); m = meta[name]; a = np.loadtxt(io.BytesIO(zf.read(n)))
            for j in range(2):
                rows.append(dict(name=name, kind='ETG_inner' if j == 0 else 'ETG_outer', r=a[j, 0], gb=10 ** m[15 + 2 * j],
                                 go=10 ** m[11 + 2 * j], ms=.8 * m[5] * 1e9, mp=.8 * m[5] * 1e9, re=m[7], z=0.))
    assert len(rows) == 32; return rows
def cluster_rows():
    rr, _ = C5.load(str(PKG / 'data/Clusters'))
    return [dict(name=x['host'], kind=x['group'], r=x['r'], gb=x['gbar'], go=x['gobs'], ms=x['Ms'], mp=x['Ms'], re=x['Re'], z=x['z']) for x in rr]

# ---------------- lensing (KiDS Fig-9 stacks, as in the A19 KiDS v9 model and SO73 run_joint) ----------------
def pop_coeff(gs):
    lm = np.log10([g['Ms'] for g in gs]); lr = np.log10([g['Re'] for g in gs])
    hi = np.array([(g['Mp'] - g['Ms']) / 1.33 for g in gs]); ok = hi > 0
    return np.polyfit(lm, lr, 1), np.polyfit(lm[ok], np.log10(hi[ok]), 1)
class Stack:
    def __init__(self, rows, coeff):
        rows = [(lm, n, kind == 'early') for kind, rr in rows for lm, n, _ in rr]
        lm = np.array([r[0] for r in rows]); n = np.array([r[1] for r in rows]); early = np.array([r[2] for r in rows])
        self.ms = 10 ** lm
        self.mp = self.ms + np.where(early, 0., 1.33 * 10 ** np.polyval(coeff[1], lm))
        self.re = 10 ** (np.polyval(coeff[0], lm) + np.where(early, np.interp(lm, K.VDW_LM, K.VDW_D), 0.))
        self.z = np.full(len(lm), K.Z0)
        mg = K.m_gal(lm); self.w = n * mg; self.w /= self.w.sum()
        edges = np.sqrt(K.GMSUN * mg[:, None] / K.GEDGES[None, :]) / K.MPC
        rout, rin = edges[:, :-1], edges[:, 1:]; lo, hi = np.log(rin), np.log(rout)
        self.R = np.exp((hi - lo)[:, :, None] * .5 * (K.GL_X[None, None, :] + 1) + lo[:, :, None])
        self.aw = .5 * (hi - lo)[:, :, None] * K.GL_W[None, None, :] * self.R ** 2 / ((rout ** 2 - rin ** 2) / 2)[:, :, None]
        self.bary = self.mp[:, None] / np.pi * np.log(rout / rin) * 2 / (rout ** 2 - rin ** 2) / 1e12
class Lens:
    def __init__(self, gs, drop):
        kd = PKG / 'data/Lensing/KiDS'
        pop = K.populations(json.load(open(kd / 'brouwer2021_fig7_digitized.json')), 'primary')
        C, _ = K.load_cov(kd / 'Fig-9_Massbins_covmatrix_compact.json')
        self.coeff = pop_coeff(gs)
        self.stacks = [Stack(pop[f'fig9_bin{i + 1}'], self.coeff) for i in range(4)]
        self.masks = [np.flatnonzero(K.inner_mask(m))[int(drop):] for m in K.FIG9_MGAL]
        self.y = [K.load_profile(kd / f'Fig-9_RAR-KiDS-isolated_Massbin-{i + 1}.txt')['esd'][k] for i, k in enumerate(self.masks)]
        self.L = [np.linalg.cholesky(C[np.ix_(k + i * 15, k + i * 15)]) for i, k in enumerate(self.masks)]
    def white(self, model, th, groups):
        return np.concatenate([solve_triangular(self.L[j], self.y[j] - esd(model, th, self.stacks[j])[self.masks[j]], lower=True) for j in groups])

# ---------------- the formulas ----------------
# Every model: one NFW halo per system, dark mass (1 - fb) M_ref inside r200 (r200 at the system's redshift), untruncated,
# the same density for dynamics and for projected lensing. log10 c = 0.905 - 0.101 log10(M_ref h / 1e12) + dc.
MODELS = {
    # name: (parameter names, lower bounds, upper bounds, starts)
    'ORIG':   (['A', 'B', 'C'], [9, -.5, -2.5], [14.5, 2, 1.5], [[11.7481, .8178, -.7119]]),
    'KNEE':   (['A', 'B', 'C', 'A2', 'B2'], [9, -.5, -2.5, 9, 0], [14.5, 2, 1.5, 16, 5],
               [[11.7481, .8178, -.7119, 12.5, 2.5], [11.7481, .8178, -.7119, 13.0, 1.5], [11.7481, .8178, -.7119, 12.0, 3.5]]),
    'KNEES':  (['A', 'B', 'C', 'A2', 'B2', 'C2'], [9, -.5, -2.5, 9, 0, -2], [14.5, 2, 1.5, 16, 5, 3],
               [[11.7481, .8178, -.7119, 12.5, 2.5, 0.], [11.7481, .8178, -.7119, 12.5, 2.0, 1.], [11.7481, .8178, -.7119, 12.0, 3.0, .5]]),
    'KNEEC':  (['A', 'B', 'C', 'A2', 'B2', 'v', 's'], [9, -.5, -2.5, 9, 0, -.5, -1], [14.5, 2, 1.5, 16, 5, .5, 1],
               [[11.7481, .8178, -.7119, 12.5, 2.5, 0., 0.], [11.7481, .8178, -.7119, 13.0, 1.5, -.1, .2], [11.7481, .8178, -.7119, 12.0, 3.5, .1, -.2]]),
    'KNEESC': (['A', 'B', 'C', 'A2', 'B2', 'C2', 'v', 's'], [9, -.5, -2.5, 9, 0, -2, -.5, -1], [14.5, 2, 1.5, 16, 5, 3, .5, 1],
               [[11.7481, .8178, -.7119, 12.5, 2.5, 0., 0., 0.], [11.7481, .8178, -.7119, 12.5, 2.0, 1., -.1, .2], [11.7481, .8178, -.7119, 12.0, 3.0, .5, .1, -.2]]),
    'KNEEE':  (['A', 'B', 'C', 'A2', 'B2', 'eta', 'logq'], [9, -.5, -2.5, 9, 0, 0, -1], [14.5, 2, 1.5, 16, 5, 10, 1],
               [[11.7481, .8178, -.7119, 12.5, 2.5, .1, -.3], [11.7481, .8178, -.7119, 12.5, 2.5, 1., -.5], [11.7481, .8178, -.7119, 12.5, 2.5, 3., 0.]]),
}
def mref(model, th, ms, mp, re):
    m = mp / FB + 10 ** th[0] * (ms / 1e10) ** th[1] * (re / 3) ** th[2]
    if model in ('KNEE', 'KNEEC', 'KNEEE'): m = m + 10 ** th[3] * (ms / 1e11) ** th[4]
    if model in ('KNEES', 'KNEESC'): m = m + 10 ** th[3] * (ms / 1e11) ** th[4] * (re / 10) ** th[5]
    dc = 0.
    if model == 'KNEEC': dc = th[5] + th[6] * np.log10(ms / (2 * np.pi * re ** 2) / 1e8)
    if model == 'KNEESC': dc = th[6] + th[7] * np.log10(ms / (2 * np.pi * re ** 2) / 1e8)
    return m, 10 ** (.905 - .101 * np.log10(m * h / 1e12) + dc)
def gpred(model, th, a):
    m, c = mref(model, th, a['ms'], a['mp'], a['re'])
    rv = (G * m / (100 * hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    g = a['gb'] + G * (1 - FB) * m * ff(x) / ff(c) / a['r'] ** 2 * ACC
    if model == 'KNEEE':                                   # Plummer envelope, mass eta (1-fb) M_ref, scale q r200
        q = 10 ** th[6] * rv; g = g + G * th[5] * (1 - FB) * m * a['r'] / (a['r'] ** 2 + q ** 2) ** 1.5 * ACC
    return g
def esd(model, th, st):
    m, c = mref(model, th, st.ms, st.mp, st.re)
    rc = K.rho_crit(K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c
    rhos = 200 / 3 * c ** 3 / ff(c) * rc
    e = (1 - FB) * rs[:, None, None] * rhos[:, None, None] * K.g_wb(st.R / rs[:, None, None]) / 1e12
    if model == 'KNEEE':                                   # projected Plummer: M R^2 / (pi (R^2 + a^2)^2) per unit area
        q = 10 ** th[6] * rv; e = e + th[5] * (1 - FB) * m[:, None, None] * st.R ** 2 / (np.pi * (st.R ** 2 + q[:, None, None] ** 2) ** 2) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)

# ---------------- fitting and scoring ----------------
SIG = {'dwarfs': .4, 'ETG': .2, 'MaNGA': .2, 'CLASH': .2}     # SO73 working scales; disks 0.15 dex
def fit(model, weight, dtr, etr, lens, groups):
    names, lo, hi, starts = MODELS[model]; nb = sum(len(lens.masks[j]) for j in groups)
    def fun(th):
        rr = [dtr['hw'] * np.log10(dtr['go'] / gpred(model, th, dtr)) / .15]
        if weight is not None:
            rr.append(np.sqrt(.1 / nb) * lens.white(model, th, groups))
            rr += [np.sqrt(weight / 4) * a['hw'] * np.log10(a['go'] / gpred(model, th, a)) / SIG[k] for k, a in etr.items()]
        return np.concatenate(rr)
    att = []
    for st in starts:
        o = least_squares(fun, st, bounds=(lo, hi), max_nfev=400, ftol=1e-9, xtol=1e-9, gtol=1e-9)
        att.append(dict(th=o.x.tolist(), loss=float(o.fun @ o.fun), success=bool(o.success), nfev=int(o.nfev)))
    good = [x for x in att if x['success']] or att
    best = min(good, key=lambda x: x['loss']); th = np.array(best['th'])
    bh = [n for n, t, l, u in zip(names, th, lo, hi) if min(abs(t - l), abs(t - u)) < 1e-5]
    return dict(th=th.tolist(), loss=best['loss'], attempts=att, bound_hits=bh)
def losses(model, th, a):
    er = np.log10(a['go'] / gpred(model, th, a)); return {n: float(np.mean(er[a['name'] == n] ** 2)) for n in np.unique(a['name'])}

def run(mask, out, configs):
    gs = disks(); folds = list(csv.DictReader(open(PKG / 'data/SPARC/A19_Disk_Folds.csv')))
    maps = json.load(open(AST / 'T56_fold_maps.json'))['external_fold_maps']
    cl = cluster_rows()
    ext = {'dwarfs': arr(dwarf_rows()), 'ETG': arr(etg_rows()), 'MaNGA': arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    recs = []; t0 = time.monotonic()
    for asg in sorted(maps):
        for fold in range(5):
            held = {r['galaxy'] for r in folds if r[asg] == str(fold)}
            tr = [g for g in gs if g['name'] not in held]; te = [g for g in gs if g['name'] in held]
            dtr, dte = arr(disk_rows(tr)), arr(disk_rows(te)); etr, ete = {}, {}
            for p, a in ext.items():
                k = np.array([maps[asg][p][n] == fold for n in a['name']]); etr[p] = sub(a, ~k); ete[p] = sub(a, k)
                assert not set(etr[p]['name']) & set(ete[p]['name'])
            lens = Lens(tr, mask == 'strict'); cache = {}
            for j in range(4):
                groups = [i for i in range(4) if i != j]
                for cfg in configs:
                    model, w = cfg.split('@'); w = None if w == 'disk' else float(w)
                    if w is None and model in cache: f = cache[model]
                    else:
                        f = fit(model, w, dtr, etr, lens, groups)
                        if w is None: cache[model] = f
                    th = np.array(f['th']); l = {'disk': losses(model, th, dte)}
                    for p, a in ete.items(): l[p] = losses(model, th, a)
                    for kind in ['CLASH_centre', 'CLASH_outer', 'ETG_inner', 'ETG_outer']:
                        a = ete['CLASH' if kind.startswith('CLASH') else 'ETG']; l[kind] = losses(model, th, sub(a, a['kind'] == kind))
                    z = lens.white(model, th, [j])
                    recs.append(dict(config=cfg, assignment=asg, fold=fold, lens=j, fit=f, losses=l, test_lensing_chi2=float(z @ z)))
            print(mask, asg, fold, round(time.monotonic() - t0, 1), flush=True)
            json.dump(dict(mask=mask, records=recs), open(out, 'w'))
    res = dict(mask=mask, configs=configs, summary=summarize(recs), records=recs,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    json.dump(res, open(out, 'w')); print(json.dumps(res['summary'], indent=1))

def summarize(recs):
    out = {}
    for cfg in sorted(set(r['config'] for r in recs)):
        rr = [r for r in recs if r['config'] == cfg]; pooled = {}
        for r in rr:
            for p, dd in r['losses'].items():
                for host, mse in dd.items(): pooled.setdefault(p, {}).setdefault(host, []).append(mse)
        sc = {}
        for p, dd in pooled.items():
            assert all(len(v) == 12 for v in dd.values()), (cfg, p)
            sc[p] = round(float(np.sqrt(np.mean([np.mean(v) for v in dd.values()]))), 6)
        lc = [float(np.mean([r['test_lensing_chi2'] for r in rr if r['lens'] == j])) for j in range(4)]
        out[cfg] = dict(scores=sc, lensing_chi2=round(sum(lc), 4), lensing_by_group=[round(x, 3) for x in lc],
                        bound_hit_partitions=sum(bool(r['fit']['bound_hits']) for r in rr),
                        failed_attempts=sum(not a['success'] for r in rr for a in r['fit']['attempts']))
    return out

if __name__ == '__main__':
    if sys.argv[1] == 'run': run(sys.argv[2], sys.argv[3], sys.argv[4:])
