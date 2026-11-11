"""SA73 T46 exploration (all data): standard halos (Moster expected halo with selection) whose concentration comes from
the galaxy's size, R_e = 0.02 (c/10)^-0.7 R_vir (Jiang et al. 2019, MNRAS 488, 4801, zoom simulations), read here with
200c definitions, for galaxy-scale halos inside the simulations' range (1e9.5 to 1e12.5 Msun); Dutton-Maccio c elsewhere
(dwarfs below the range, cluster-scale halos). R_e taken as the projected radius, and as 4/3 of it (3-D)."""
import json, numpy as np
import harness_T45 as M
H = M.H; S = M.S; K = H.K; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()),
     'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']),
     'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
def g(a, model, mode, f3d):
    thr = S.thresholds(a)
    if model == 'std': m = S.m_exp(a['ms'], a['z'], thr); m = (m ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8); fac = 1.
    else: m = M.P.G.mass(TH, a['ms'], a['mp'], a['re'], a['z'], thr); fac = a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); c = K.c_dm14(m, a['z']) * fac
    if mode == 'size':
        cj = 10 * (f3d * a['re'] / (0.02 * rv)) ** (-1 / 0.7)
        ok = (m >= 10 ** 9.5) & (m <= 10 ** 12.5) & (np.asarray(thr) == 0)
        c = np.where(ok, cj, c)
    mn = m * H.ff(a['r'] * c / rv) / H.ff(c)
    return a['gb'] + H.G * (1 - H.FB) * mn / a['r'] ** 2 * H.ACC, c
out = {}
for model in ('std', 'formula'):
    for mode, f3d in (('dm14', 1.), ('size', 1.), ('size', 4 / 3)):
        r = {}; cs = []
        for p, a in P.items():
            gp, c = g(a, model, mode, f3d); er = np.log10(a['go'] / gp); hosts = np.unique(a['name'])
            r[p] = round(float(np.sqrt(np.mean([np.mean(er[a['name'] == h] ** 2) for h in hosts]))), 4)
            r[p + '_mean'] = round(float(np.mean([np.mean(er[a['name'] == h]) for h in hosts])), 3)
            if p == 'disk': cs = c
        key = f'{model}|{mode}|{f3d:.2f}'; out[key] = r
        print(key, {k: v for k, v in r.items() if not k.endswith('_mean')}, 'disk mean', r['disk_mean'], 'disk c median', round(float(np.median(cs)), 1), flush=True)
json.dump(out, open('res46_sizec.json', 'w'), indent=1)
