"""SA73 T46 exploration (all data): standard halos (Moster expected halo with selection, Dutton-Maccio c, no subhalo
factor) contracted by their baryons with the published relation of Cautun et al. (2020, eq. 11), strength s:
M_DM(<r) = (1 - fb) M_NFW(<r) F(eta)^s, F = 0.45 + 0.38 (eta + 1.16)^0.53, eta = M_bar(<r) / (fb M_NFW(<r)),
with M_bar(<r) = g_bar r^2 / G at each measured radius. Also the formula's halos (rule at all-disk A, B, C, not refitted)."""
import json, numpy as np
import harness_T45 as M
H = M.H; S = M.S; K = H.K; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()),
     'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']),
     'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
def F(eta): return 0.45 + 0.38 * (eta + 1.16) ** 0.53
def g(a, model, s):
    thr = S.thresholds(a)
    if model == 'std':
        m = S.m_exp(a['ms'], a['z'], thr); m = (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8); c = K.c_dm14(m, a['z'])
    else:
        m = M.P.G.mass(TH, a['ms'], a['mp'], a['re'], a['z'], thr); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mn = m * H.ff(a['r'] * c / rv) / H.ff(c)
    mbar = a['gb'] / H.ACC * a['r'] ** 2 / H.G; eta = mbar / (H.FB * mn)
    return a['gb'] + H.G * (1 - H.FB) * mn * F(eta) ** s / a['r'] ** 2 * H.ACC
out = {}
for model in ('std', 'formula'):
    for s in (0., .25, .5, 1.):
        r = {}
        for p, a in P.items():
            er = np.log10(a['go'] / g(a, model, s)); hosts = np.unique(a['name'])
            r[p] = round(float(np.sqrt(np.mean([np.mean(er[a['name'] == h] ** 2) for h in hosts]))), 4)
        out[f'{model}|{s}'] = r; print(model, s, r, flush=True)
json.dump(out, open('res46_contract.json', 'w'), indent=1)
