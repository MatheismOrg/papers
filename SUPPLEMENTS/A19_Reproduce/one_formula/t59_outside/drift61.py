"""SA73 T61 look (no formula): on the seven galaxies shared by SPARC and Iorio et al. (2017), SPARC's rotation speeds at the
points the benchmark drops and keeps, against Iorio et al.'s circular speeds (corrected for asymmetric drift) and their
rotation speeds before the correction, at the same physical radius. Usage: python3 drift61.py"""
import sys, json, numpy as np
S = '/tmp/claude-0/-home-claude/0a1f618d-0405-5db3-a0ab-84f1f9545369/scratchpad'
sys.path.insert(0, S + '/pkg/A19_Data_Package/code'); import A19_Disk_Fit_SPARC as D
PK = S + '/pkg/A19_Data_Package/data/SPARC/'; t1 = D.table1(PK + 'SPARC_Lelli2016c.mrt'); gal = json.load(open('lt59_inputs.json')); rows = []; out = {}
for nm, g in gal.items():
    sp = g['sparc']
    if sp is None: continue
    r = t1[sp]; d = np.loadtxt(PK + 'Rotmod_LTG/' + sp + '_rotmod.dat', comments='#', ndmin=2); Rs, Vs, eVs = d[:, 0], d[:, 1], d[:, 2]
    f = g['D'] / r['D']; R = np.array(g['R']) / f; Vc = np.array(g['Vc']); Vr = np.array(g['Vr'])
    # inclination: bring SPARC's speeds to Iorio's inclination so that only the kinematic treatment differs
    si = np.sin(np.radians(r['Inc'])) / np.sin(np.radians(g['inc']))
    k = (Rs >= R[0]) & (Rs <= R[-1]) & (Vs > 0); vc = np.interp(Rs[k], R, Vc); vr = np.interp(Rs[k], R, Vr); drop = (eVs / Vs)[k] > 0.1
    for a, b, c_, dd, e in zip(Vs[k] * si, vc, vr, drop, Rs[k] / r['Reff']): rows.append((nm, a, b, c_, dd, e))
    o = {}
    for lab, kk in (('dropped', drop), ('kept', ~drop)):
        if kk.sum(): o[lab] = dict(n=int(kk.sum()), sparc_over_vc=round(float(np.median((Vs[k] * si)[kk] / vc[kk])), 3), vr_over_vc=round(float(np.median(vr[kk] / vc[kk])), 3), sparc_over_vr=round(float(np.median((Vs[k] * si)[kk] / vr[kk])), 3))
    out[nm] = o; print('%-9s' % nm, o)
A = np.array([(a, b, c_, dd, e) for _, a, b, c_, dd, e in rows], float)
for lab, kk in (('dropped', A[:, 3] > 0), ('kept', A[:, 3] == 0), ('dropped, V_sparc < 30', (A[:, 3] > 0) & (A[:, 0] < 30)), ('all inside R_e', A[:, 4] < 1), ('all V_sparc < 30', A[:, 0] < 30)):
    x = A[kk]; o = dict(n=int(kk.sum()), sparc_over_vc=round(float(np.median(x[:, 0] / x[:, 1])), 3), vr_over_vc=round(float(np.median(x[:, 2] / x[:, 1])), 3), sparc_over_vr=round(float(np.median(x[:, 0] / x[:, 2])), 3),
                        dex_in_g_sparc_vs_vc=round(float(2 * np.mean(np.log10(x[:, 0] / x[:, 1]))), 3), dex_in_g_drift=round(float(2 * np.mean(np.log10(x[:, 2] / x[:, 1]))), 3))
    out['pooled ' + lab] = o; print('pooled', lab, o)
json.dump(out, open('res61_drift.json', 'w'), indent=1)
