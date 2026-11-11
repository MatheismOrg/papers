"""SA73 T56 all-data look (no protocol yet): harness_T56's formula fitted on all disks, against the final formula (T48).
Usage: python3 describe56.py   (from the dmfit folder)"""
import sys, csv, json, numpy as np
from scipy.linalg import solve_triangular
import harness_T56 as M6
import harness_T48 as M8
P = M6.P; H = M6.H; B = M6.B; S = M6.S; K = H.K
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
(p1, p0), (q1, q0) = H.pop_coeff(gs); rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
LENS0 = H.Lens(gs, False)
def score(model, th, a): return round(float(np.sqrt(np.mean(list(H.losses(model, th, a).values())))), 4)
def look(model, lab):
    if model != M6.MODEL: B.FROZEN['abc'] = np.array(H.fit('ORIG', None, d, ext, LENS0, [0, 1, 2])['th'])
    f = H.fit(model, .1, d, ext, LENS0, [0, 1, 2]); th = np.array(f['th'])
    r = dict(abc=np.round(B.FROZEN['abc'], 4).tolist(), mu=round(M6.U['mu'], 3) if model == M6.MODEL else 1., logMt=round(float(th[0]), 3), red_chi2=round(f['loss'], 1), disk=score(model, th, d), **{k: score(model, th, a) for k, a in ext.items()})
    for drop in (False, True):
        L = H.Lens(gs, drop); c2 = [round(float(np.sum(L.white(model, th, [j]) ** 2)), 1) for j in range(4)]; r['lens_' + ('strict' if drop else 'primary')] = (round(sum(c2), 1), c2)
    dw = ext['dwarfs']; e = np.log10(dw['go'] / H.gpred(model, th, dw)); X = np.c_[np.ones(len(e)), np.log10(dw['ms']), np.log10(dw['re'])]
    r['dwarf_slope'] = round(float(np.linalg.lstsq(X, e, rcond=None)[0][2]), 3); r['dwarf_mean'] = round(float(np.mean(e)), 3)
    mfun = (lambda *a: M6.mass(th, *a)) if model == M6.MODEL else (lambda *a: P.P.G.mass(th, *a))
    c2 = {'red': 0., 'blue': 0.}; devs = []
    for w in rows:
        lm, z, obs, ep, em = float(w['logMs_eff_VAGC']), float(w['z_eff']), float(w['logM200m_hinv']), float(w['err_plus']), float(w['err_minus']); red = w['color'] == 'red'
        ms = 10 ** lm; re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
        m = float(mfun(np.array([ms]), np.array([mp]), np.array([re]), np.array([z]), 0.)[0]); dev = L8.to_m200m(m, K.c_dm14(m, z), z) - obs
        c2[w['color']] += (dev / (ep if dev > 0 else em)) ** 2; devs.append((w['color'][0], lm, round(dev, 2)))
    r['mandelbaum_chi2'] = {k: round(v, 1) for k, v in c2.items()}; r['blue_dev'] = [x[2] for x in devs if x[0] == 'b']; r['red_dev'] = [x[2] for x in devs if x[0] == 'r']
    mstar = 5.43e10; a_ = np.linspace(1e-3, 60, 6000); discs = {'HI': (53.1, 7., 4., .085), 'H2': (2180., 1.5, 12., .045)}
    mgas = sum(2 * np.pi * np.sum(a_ * p[0] * 1e6 * np.exp(-p[2] / a_ - a_ / p[1])) * (a_[1] - a_[0]) for p in discs.values())
    r['MW_log_M200c'] = round(float(np.log10(mfun(np.array([mstar]), np.array([mstar + mgas]), np.array([5.75]), np.array([0.]), 0.)[0])), 3)
    print(lab, json.dumps(r), flush=True); return r
out = {'T48': look(M8.MODEL, 'T48 final '), 'T56': look(M6.MODEL, 'T56 gas mu')}
json.dump(out, open('res56_describe.json', 'w'), indent=1)
