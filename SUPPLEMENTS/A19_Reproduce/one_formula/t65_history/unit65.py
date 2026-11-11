"""T65 unit checks (inputs and limits only; no score of any reading is computed here)."""
import sys, json, numpy as np
import lib65 as L
import harness_T48 as M8
H = L.H; M5 = L.M5; M61 = L.M61; B = L.B; S = L.S
B.FROZEN['abc'] = L.ABC; TH = L.TH
gs = H.disks(); cl = H.cluster_rows()
pops = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
M5.U['t'] = 0.; M5.U['k'] = 0.
for p, a in pops.items():
    thr = S.thresholds(a); wg = M61.gl(a['ms'], a['mp'], a['re'], thr); wh = M5.hl(a['ms'], a['mp'], a['re'], thr); un, i = np.unique(a['name'], return_index=True)
    print('%-7s hosts %3d | t = kappa = 0 reproduces T48: max |dlog g| %.1e | hosts with gas gate above 0: %3d, with history-harness gate above 0: %3d | hosts with M_p above M*: %3d' % (
        p, len(un), np.max(np.abs(np.log10(H.gpred(M5.MODEL, TH, a) / H.gpred(M8.MODEL, TH, a)))), int(np.sum(wg[i] > 0)), int(np.sum(wh[i] > 0)), int(np.sum(a['mp'][i] > a['ms'][i] * (1 + 1e-9)))))
a = pops['disk']; un, i = np.unique(a['name'], return_index=True); g = 1 - a['ms'][i] / a['mp'][i]; Ld = M5.diffuse(a['ms'][i], a['re'][i])
print('149 disks: gas fraction g percentiles (5, 25, 50, 75, 95):', np.round(np.percentile(g, [5, 25, 50, 75, 95]), 2), '| diffuseness L:', np.round(np.percentile(Ld, [5, 25, 50, 75, 95]), 2), '| L above 0 in', int(np.sum(Ld > 0)), '| correlation of g and L: %.2f' % np.corrcoef(g, Ld)[0, 1], '| among L above 0: %.2f' % np.corrcoef(g[Ld > 0], Ld[Ld > 0])[0, 1])
d = pops['dwarfs']; print('75 dwarf spheroidals: diffuseness L percentiles (5, 50, 95):', np.round(np.percentile(M5.diffuse(d['ms'], d['re']), [5, 50, 95]), 2))
# the cored profile: limits and the density formula against a numerical derivative of the mass
k = d['name'] == 'fornax_1'; a1 = {key: np.asarray(v)[k] for key, v in d.items()}; h = L.halo(a1); r = np.array([0.05, 0.15, 0.5, 0.927, 3.])
print('Fornax halo: log M %.3f c %.2f r_200 %.2f kpc, t_dyn %.3f Gyr, G = %.5e, 1 - f_b = %.4f' % (np.log10(h[0][0]), h[1][0], h[2][0], L.tdyn(h)[0], H.G, 1 - H.FB))
for n in (0., 0.3, 1.):
    rc = L.ETA * a1['re']; e = 1e-4; num = (L.core_m(h, r * (1 + e), n, rc) - L.core_m(h, r * (1 - e), n, rc)) / (4 * np.pi * r ** 2 * 2 * e * r)
    print('  n = %.1f: density formula / numerical derivative of mass: %s | mass / NFW mass: %s' % (n, np.round(L.core_rho(h, r, n, rc) / num, 5), np.round(L.core_m(h, r, n, rc) / L.nfw_m(h, r), 4)))
print('  reading G equals N for a gas-free dwarf:', bool(np.allclose(L.reading(a1, r, 'G')[1], L.reading(a1, r, 'N')[1])))
sf = json.load(open('sfh65_inputs.json'))
print('published history numbers n = tanh(0.04 t_SF / t_dyn) for the eight dwarf spheroidals of Read et al.:')
for lab, v in sf['read2019'].items():
    if v['type'] != 'dSph': continue
    kk = d['name'] == L.DSPH[lab]; a2 = {key: np.asarray(x)[kk] for key, x in d.items()}; h2 = L.halo(a2); tsf = L.TU - v['t_trunc']
    print('  %-9s log M_h %.2f (Read et al. %.2f) R_e %.3f t_SF %5.2f t_dyn %.2f Gyr n %.2f | history gate w = s L = %.2f' % (lab, np.log10(h2[0][0]), np.log10(v['M200_1e9'] * 1e9), a2['re'][0], tsf, L.tdyn(h2)[0], L.n_read(h2, tsf)[0], L.w_hist(a2, tsf)[0]))
print('Weisz et al. (2015) dwarfs in the benchmark:', sum(1 for kq, nm in L.WEISZ.items() if nm in set(d['name']) and kq in sf['weisz2015']), 'of', len(L.WEISZ))
