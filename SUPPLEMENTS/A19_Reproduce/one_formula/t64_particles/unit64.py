"""T64 unit checks (no scores): cold limits reproduce T48; the soliton's numbers; the warm factor."""
import numpy as np
import harness_T64 as M
import harness_T48 as M8
H = M.H; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
gs = H.disks(); cl = H.cluster_rows()
pops = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
M.V['lm22'] = 3.0; M.V['lhm'] = 3.0
for p, a in pops.items():
    g8 = H.gpred(M8.MODEL, TH, a); gf = H.gpred('MOGF0IESBF', TH, a); gw = H.gpred('MOGF0IESBW', TH, a)
    m, c, rv = M.halo(TH, a)
    print('%-7s n %5d hosts %3d  cold limits: fuzzy max |dlog g| %.2e  warm %.2e   log M_h %.2f to %.2f  r %.3g to %.4g kpc' % (p, len(a['r']), len(np.unique(a['name'])), np.max(np.abs(np.log10(gf / g8))), np.max(np.abs(np.log10(gw / g8))), np.log10(m.min()), np.log10(m.max()), a['r'].min(), a['r'].max()))
print('all points:', len(M.ALL['r']), 'galaxies', len(np.unique(M.ALL['name'])))
for mh in (1e9, 1e10, 1e11, 1e12, 1e15):
    rc, rhoc = M.soliton(mh, 0., 1.0); print('M_h %.0e m22 1: r_c %.3f kpc, rho_c %.3g Msun/pc^3, M_sol(3 r_c) %.3g, M(r_c) %.3g, share inside 3 r_c %.3f' % (mh, rc, rhoc / 1e9, rhoc * rc ** 3 * M.fsol(3.), rhoc * rc ** 3 * M.fsol(1.), M.fsol(3.) / 11.58))
print('zeta(0) %.1f; half density at r_c: %.4f; warm factor at M = M_hm, 10 M_hm, 100 M_hm: %s' % (M.zeta(0.), (1 + 0.091) ** -8, np.round(M.wdm_fac(np.array([1., 10., 100.]), 1.), 3)))
