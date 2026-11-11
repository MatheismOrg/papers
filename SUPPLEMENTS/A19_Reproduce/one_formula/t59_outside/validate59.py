"""SA73 T59: the inputs built for the Iorio et al. (2017) galaxies, checked on the six galaxies shared with SPARC.
No formula enters. Usage: python3 validate59.py"""
import os, sys, json, numpy as np
import model59 as M
S = '/tmp/claude-0/-home-claude/0a1f618d-0405-5db3-a0ab-84f1f9545369/scratchpad'
sys.path.insert(0, S + '/pkg/A19_Data_Package/code'); import A19_Disk_Fit_SPARC as D
PK = S + '/pkg/A19_Data_Package/data/SPARC/'; t1 = D.table1(PK + 'SPARC_Lelli2016c.mrt'); gal = json.load(open('lt59_inputs.json'))
print('independent names found in SPARC:', [g['lvl'] for g in gal.values() if g['sparc'] is None and g['lvl'] in t1])
out = {}
for nm, g in gal.items():
    sp = g['sparc']
    if sp is None: continue
    r = t1[sp]; d = np.loadtxt(PK + 'Rotmod_LTG/' + sp + '_rotmod.dat', comments='#', ndmin=2); Rs, Vs, eVs, Vg, Vd = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]
    R = np.array(g['R']); f = g['D'] / r['D']                    # Iorio distance over SPARC distance
    o = dict(D_iorio=g['D'], D_sparc=r['D'], inc_iorio=g['inc'], inc_sparc=r['Inc'], Q=r['Q'])
    o['L36_lvl_over_sparc'] = M.l36(g['f36'][0], r['D']) / r['L36']
    o['Re_1.678Rstar_over_sparc'] = 1.678 * g['read']['Rstar'] * r['D'] / g['read']['D'] / r['Reff']
    for mode in ('cut', 'exp'):
        v2, mm, mo = M.gas(R, g['Sd'], mode)
        o['MHI_%s_over_sparc' % mode] = (mm if mode == 'exp' else mo) / f ** 2 / r['MHI']
        k = (R / f >= Rs[0]) & (R / f <= Rs[-1]); vs2 = np.interp(R[k] / f, Rs, Vg * np.abs(Vg)) * f; kk = vs2 > 0
        o['Vgas2_%s_over_sparc_median' % mode] = float(np.median(v2[k][kk] / vs2[kk])); o['Vgas2_%s_outer_third' % mode] = float(np.median((v2[k][kk] / vs2[kk])[-max(1, kk.sum() // 3):]))
    o['MHI_read_over_sparc'] = g['read']['Mgas_1e7'] * 1e7 / 1.33 * (r['D'] / g['read']['D']) ** 2 / r['MHI']
    k = (R / f >= Rs[0]) & (R / f <= Rs[-1]); o['Vc2_iorio_over_sparc_median'] = float(np.median(np.array(g['Vc'])[k] ** 2 / np.interp(R[k] / f, Rs, Vs) ** 2))
    o['n_iorio'] = len(R); o['n_iorio_err10'] = int(np.sum(np.array(g['eVc']) / np.array(g['Vc']) <= 0.1)); o['n_sparc_err10'] = int(np.sum(eVs / Vs <= 0.1))
    out[nm] = o
keys = [k for k in next(iter(out.values()))]
for k in keys: print('%-30s' % k, ' '.join('%8.3f' % out[n][k] for n in out), '| median %.3f' % np.median([out[n][k] for n in out]))
print(' ' * 30, ' '.join('%8s' % n.replace(' ', '') for n in out))
json.dump(out, open('res59_validate.json', 'w'), indent=1)
