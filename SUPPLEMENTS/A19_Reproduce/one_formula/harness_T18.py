"""T18: T17b's switch formula (MOSS) + halo contraction from Cautun et al. (2020, eq. 11) at strength s, in dynamics and lensing:
  M_DM(<r) = (1 - fb) M_NFW(<r) F(eta)^s,  F = 0.45 + 0.38 (eta + 1.16)^0.53,  eta = M_bar(<r) / (fb M_NFW(<r)),
with M_NFW the uncontracted total (dark-matter-only) NFW mass. Dynamics: M_bar(<r) = g_bar r^2 / G at the data radius.
Lensing: the galaxy's baryons M_p as a point mass (bins lie at 4 or more lens effective radii), projected numerically.
The paper's rule A, B, C (and s, for CAUD) are refitted on each fold's training disks, with contraction and without the switch;
then R_t is fitted as in T17b. CAUD: s fitted on disks in [0, 1.5]. CAU5: s = 0.5. CAU1: s = 1, the published strength.
Imports harness_T17b (and through it T17, T16, T15c, T15b and the hashed T15 harness) unchanged.
Usage: python3 harness_T18.py <mask> <out> cfg..."""
import sys, json, csv, time, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import harness_T17b as Y
H = Y.H; B = Y.B; K = H.K
S_NOW = {'s': 0.}
CFG = {'CAUD': None, 'CAU5': 0.5, 'CAU1': 1.0}
for k in CFG: H.MODELS[k] = H.MODELS['MOSS']
def F(eta): return 0.45 + 0.38 * (eta + 1.16) ** 0.53
def dF(eta): return 0.2014 * (eta + 1.16) ** (-0.47)
def dyn(a, m, s):
    c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    eta = a['gb'] * a['r'] ** 2 / (H.G * H.ACC) / (H.FB * mt)
    return a['gb'] + H.G * (1 - H.FB) * mt * F(eta) ** s / a['r'] ** 2 * H.ACC
def rule_mass(a):
    abc = B.FROZEN['abc']; return a['mp'] / H.FB + 10 ** abc[0] * (a['ms'] / 1e10) ** abc[1] * (a['re'] / 3) ** abc[2]
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a):
    if model not in CFG: return _gpred(model, th, a)
    return dyn(a, Y.mass('MOSS', th, a['ms'], a['mp'], a['re'], a['z']), S_NOW['s'])
TN, TW = np.polynomial.legendre.leggauss(48); TMAX = 9.; TT = .5 * TMAX * (TN + 1); WT = .5 * TMAX * TW
CH, SH, TH = np.cosh(TT), np.sinh(TT), np.tanh(TT)
def esd_dm(m, c, rv, mp, R, s):
    """Projected excess surface density (Msun/pc^2) of the contracted dark matter; m, c, rv, mp per sample, R (n, bins, nodes) in Mpc."""
    k = (c / rv)[:, None, None, None]; M4 = m[:, None, None, None]; P4 = mp[:, None, None, None]; fc = H.ff(c)[:, None, None, None]
    def parts(r):
        x = r * k; mn = M4 * H.ff(x) / fc; dmn = M4 * k * x / (1 + x) ** 2 / fc
        eta = P4 / (H.FB * mn); f = F(eta)
        return (1 - H.FB) * mn * f ** s, (1 - H.FB) * dmn * (f ** s - s * f ** (s - 1) * dF(eta) * eta)
    R4 = R[..., None]
    md, _ = parts(R4)
    _, dmd = parts(R4 * CH)
    sigma = (dmd / CH) @ WT / (2 * np.pi * R)
    mcyl = md[..., 0] + (dmd * (1 - TH) * R4 * SH) @ WT
    return (mcyl / (np.pi * R ** 2) - sigma) / 1e12
def esd(model, th, st):
    if model not in CFG: return _esd(model, th, st)
    m = Y.mass('MOSS', th, st.ms, st.mp, st.re, st.z); c = K.c_dm14(m, K.Z0)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    return st.w @ (np.sum(st.aw * esd_dm(m, c, rv, st.mp, st.R, S_NOW['s']), axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
def disk_fit(dtr, s_fixed):
    """The paper's rule (and s if not fixed) on training disks, with contraction, without the switch."""
    def fun(p):
        B.FROZEN['abc'] = p[:3]; s = p[3] if s_fixed is None else s_fixed
        return dtr['hw'] * np.log10(dtr['go'] / dyn(dtr, rule_mass(dtr), s)) / .15
    lo, hi, st = [9, -.5, -2.5], [14.5, 2, 1.5], [11.7481, .8178, -.7119]
    if s_fixed is None: lo, hi, st = lo + [0.], hi + [1.5], st + [.3]
    o = least_squares(fun, st, bounds=(lo, hi), max_nfev=400, ftol=1e-9, xtol=1e-9, gtol=1e-9)
    B.FROZEN['abc'] = o.x[:3]; S_NOW['s'] = float(o.x[3]) if s_fixed is None else s_fixed
    return dict(abc=o.x[:3].tolist(), s=S_NOW['s'], success=bool(o.success), s_at_bound=bool(s_fixed is None and min(o.x[3], 1.5 - o.x[3]) < 1e-5))
def run(mask, out, configs):
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
            for cfg in configs:
                model, w = cfg.split('@'); dfit = disk_fit(dtr, CFG[model])
                for j in range(4):
                    groups = [i for i in range(4) if i != j]
                    f = H.fit(model, float(w), dtr, etr, lens, groups); th = np.array(f['th']); f['disk_fit'] = dfit
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
