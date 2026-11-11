"""SA73 T46 comparison (Gary 11:14: "what's the dex in all scales now v before and best current physics"): the same
benchmark, scoring and lensing covariance for two published frameworks with nothing fitted, on all data (equal to held out
when nothing is fitted): MOND's radial acceleration relation (McGaugh et al. 2016, a0 = 1.2e-10 m/s^2; lensing from its
phantom mass, projected) and standard LCDM halos from the published Moster relation (expected halo M_E with the Tinker HMF
and each sample's selection, NFW, Dutton-Maccio c, Moline subhalo factor, 1e8 floor, all-halo c for the lensing stacks).
Each is scored with the original stellar g_bar and with the spheroid IMF trend applied (as the final formula uses)."""
import json, numpy as np
from scipy.linalg import solve_triangular
H = __import__('harness_T15_hashed')
gs = H.disks()
def pops(etg_rows, cluster_rows):
    cl = cluster_rows()
    return {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(etg_rows()),
            'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']),
            'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
RAW = pops(H.etg_rows, H.cluster_rows)                 # before any IMF patch
import harness_T45 as M                                # patches the spheroid rows with the IMF trend
S = M.S; K = H.K; T = M.P.T; Z = M.P.Z
IMF = pops(H.etg_rows, H.cluster_rows)
A0 = 1.2e-10
def nu(y): return 1 / (1 - np.exp(-np.sqrt(y)))
def g_mond(a): return a['gb'] * nu(a['gb'] / A0)
def g_lcdm(a):
    m = S.m_exp(a['ms'], a['z'], S.thresholds(a)); m = (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
    c = K.c_dm14(m, a['z']) * a.get('fac', 1.); rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(a['r'] * c / rv) / H.ff(c) / a['r'] ** 2 * H.ACC
def score(gfun, P):
    out = {}
    for p, a in P.items():
        er = np.log10(a['go'] / gfun(a)); hosts = np.unique(a['name'])
        out[p] = round(float(np.sqrt(np.mean([np.mean(er[a['name'] == n] ** 2) for n in hosts]))), 4)
        out[p + '_mean'] = round(float(np.mean([np.mean(er[a['name'] == n]) for n in hosts])), 4)
    return out
CH, SH, TH, WT = T.CH, T.SH, T.TH, T.WT
def esd_mond(st):
    mp = st.mp[:, None, None, None]; R4 = st.R[..., None]
    def parts(r):
        y = K.GMSUN * mp / (r * K.MPC) ** 2 / A0; sy = np.sqrt(y); e = np.exp(-sy)
        mph = mp * (nu(y) - 1); dnu = -e / (2 * sy * (1 - e) ** 2)
        return mph, mp * dnu * (-2 * y / r)
    md, _ = parts(R4); _, dmd = parts(R4 * CH)
    sigma = (dmd / CH) @ WT / (2 * np.pi * st.R); mcyl = md[..., 0] + (dmd * (1 - TH) * R4 * SH) @ WT
    e = (mcyl / (np.pi * st.R ** 2) - sigma) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def esd_lcdm(st):
    m = S.m_exp(st.ms, st.z, 0.); m = (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
    c = K.c_dm14(m, K.Z0) * M.P.ratio_all_relaxed(m); rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
def lens_chi2(efun, drop):
    L = H.Lens(gs, drop); out = []
    for j, st in enumerate(L.stacks):
        e = efun(st)[L.masks[j]]; z = solve_triangular(L.L[j], L.y[j] - e, lower=True); out.append(round(float(z @ z), 2))
    return out
res = {}
for name, gfun, efun in [('MOND_RAR', g_mond, esd_mond), ('LCDM_Moster_ME', g_lcdm, esd_lcdm)]:
    for tag, P in (('stdIMF', RAW), ('spheroidIMF', IMF)):
        r = score(gfun, P)
        lp, ls = lens_chi2(efun, False), lens_chi2(efun, True)
        r.update(lens=round(sum(lp), 2), lens_by_stack=lp, lens_strict=round(sum(ls), 2), lens_strict_by_stack=ls)
        res[f'{name}_{tag}'] = r; print(name, tag, json.dumps(r), flush=True)
json.dump(res, open('res46_compare.json', 'w'), indent=1)
