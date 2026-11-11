"""SA73 T45 exploration (all data, no protocol, not a test): (1) set the gate threshold M_t outside the benchmark, on the
SDSS red central halo masses of Mandelbaum et al. (2016); (2) a mixture form of the halo, q M_E + (1 - q)(M_p/f_b + rule),
in place of T44's sum M_p/f_b + rule + q M_E. Usage: python3 explore45.py"""
import sys, csv, json, time
import numpy as np
import harness_T44 as M
S = M.S; H = M.H; K = H.K; T = M.T; Z = M.Z; B = S.B
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
ABC = np.array([11.7481, .8178, -.7119]); B.FROZEN['abc'] = ABC.copy()
MIX = 'MXGF0IESR'; H.MODELS[MIX] = H.MODELS['MOGF0']
def mass_mix(th, ms, mp, re, z, thr):
    abc = B.FROZEN['abc']; thr = np.broadcast_to(np.asarray(thr, float), np.shape(ms))
    q = np.where(thr > 0, 1., Z.switch('MOG', th, ms, mp, re))
    m = q * S.m_exp(ms, z, thr) + (1 - q) * (mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2])
    return (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
def mass_add(th, ms, mp, re, z, thr): return M.G.mass(th, ms, mp, re, z, thr)
_gp, _es = H.gpred, H.esd
def gpred(model, th, a):
    if model != MIX: return _gp(model, th, a)
    m = mass_mix(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model != MIX: return _es(model, th, st)
    m = mass_mix(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * M.ratio_all_relaxed(m)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
FORMS = {'add': ('MOGF0IESR', mass_add), 'mix': (MIX, mass_mix)}
gs = H.disks(); (p1, p0), (q1, q0) = H.pop_coeff(gs)
MB = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
def mand_pred(massf, lmt, row):
    lm, z = float(row['logMs_eff_VAGC']), float(row['z_eff']); red = row['color'] == 'red'; ms = 10 ** lm
    re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
    m = float(massf(np.array([lmt]), np.array([ms]), np.array([mp]), np.array([re]), np.array([z]), 0.)[0])
    return L8.to_m200m(m, K.c_dm14(m, z), z)
def mand_chi2(massf, lmt, color='red'):
    c2, dev = 0., []
    for r in MB:
        if r['color'] != color: continue
        p = mand_pred(massf, lmt, r); d = p - float(r['logM200m_hinv']); e = float(r['err_plus'] if d > 0 else r['err_minus'])
        c2 += (d / e) ** 2; dev.append(round(d, 3))
    return c2, dev
def calib(massf):
    grid = np.round(np.arange(10.0, 12.5 + 1e-9, 0.01), 3); c2 = np.array([mand_chi2(massf, g)[0] for g in grid])
    k = int(np.argmin(c2)); g0 = grid[k]
    fine = np.round(np.arange(g0 - .01, g0 + .01 + 1e-9, .001), 4); c2f = np.array([mand_chi2(massf, g)[0] for g in fine])
    j = int(np.argmin(c2f)); ok = c2f <= c2f[j] + 1
    lo = fine[ok].min(); hi = fine[ok].max()
    # 1-sigma range on the coarse grid
    okc = c2 <= c2[k] + 1
    return float(fine[j]), float(c2f[j]), [float(grid[okc].min()), float(grid[okc].max())]
cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
lens = H.Lens(gs, False); lens_s = H.Lens(gs, True)
def scores(model, th):
    sc = {'disk': H.losses(model, th, d)}; sc.update({k: H.losses(model, th, a) for k, a in ext.items()})
    out = {k: round(float(np.sqrt(np.mean(list(v.values())))), 4) for k, v in sc.items()}
    out['lens'] = round(float(np.sum(lens.white(model, th, [0, 1, 2, 3]) ** 2)), 2)
    out['lens_strict'] = round(float(np.sum(lens_s.white(model, th, [0, 1, 2, 3]) ** 2)), 2)
    return out
res = {}
for name, (model, massf) in FORMS.items():
    t0 = time.time()
    lmt, c2, rng = calib(massf)
    c2_44, dev44 = mand_chi2(massf, 11.1182); _, devc = mand_chi2(massf, lmt); _, devb = mand_chi2(massf, lmt, 'blue')
    fitb = H.fit(model, .1, d, {k: a for k, a in ext.items() if not k.startswith('CLASH')} | {'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}, lens, [0, 1, 2, 3])
    thb = np.array(fitb['th'])
    res[name] = dict(Mt_mandelbaum=lmt, chi2_red=round(c2, 2), Mt_1sigma=rng, red_dev_at_cal=devc, blue_dev=devb,
                     chi2_red_at_benchmark_Mt=round(c2_44, 2), red_dev_at_11p118=dev44,
                     Mt_benchmark_fit=round(float(thb[0]), 4), bound_hits=fitb['bound_hits'],
                     red_dev_at_benchmark_fit=mand_chi2(massf, float(thb[0]))[1],
                     scores_cal=scores(model, np.array([lmt])), scores_benchmark_fit=scores(model, thb), sec=round(time.time() - t0))
    print(name, json.dumps(res[name]), flush=True)
json.dump(res, open('res45_explore.json', 'w'), indent=1)
