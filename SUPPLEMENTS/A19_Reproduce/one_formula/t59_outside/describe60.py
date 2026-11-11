"""SA73 T60, description after Part A (no forecast): the dropped SPARC points by size of the velocity error, by radius, and
with a robust score. Usage: PYTHONPATH=<dmfit> python3 describe60.py"""
import os, json, numpy as np
src = open('test59.py').read().split("\ntest = [")[0]; exec(compile(src, 'test59_defs', 'exec'))
t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); rw = []; ex = []
for g, r in sorted(t1.items()):
    d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (vb2 > 0)
    Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']
    rw += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
    ex += [(f, a / r['Reff'], v, i / len(R)) for f, a, v, i in zip((eV / V)[k], R[k], V[k], np.arange(len(R))[k])]
a = H.arr(rw); ex = np.array(ex); fe, x, V, rank = ex.T
e = dict(T48=np.log10(a['go'] / H.gpred(M8.MODEL, th, a)), T57=np.log10(a['go'] / H.gpred('MOGF0IESBD', th, a)), MOND=np.log10(a['go'] / C.g_mond(a)), STD=np.log10(a['go'] / C.g_lcdm(a)))
def hb(v, k, clip=None):
    un, gi = np.unique(a['name'][k], return_inverse=True); w = v[k] if clip is None else np.clip(v[k], -clip, clip)
    return round(float(np.sqrt(np.mean(np.bincount(gi, w ** 2) / np.bincount(gi)))), 4), round(float(np.mean(np.bincount(gi, v[k]) / np.bincount(gi))), 4), len(un), int(k.sum())
out = {}
for lab, k in (('kept (error <= 10%)', fe <= .1), ('dropped (error > 10%)', fe > .1), ('dropped, 10-20%', (fe > .1) & (fe <= .2)), ('dropped, 20-40%', (fe > .2) & (fe <= .4)), ('dropped, > 40%', fe > .4),
               ('dropped, inside R_e', (fe > .1) & (x < 1)), ('dropped, outside R_e', (fe > .1) & (x >= 1)), ('dropped, V < 30 km/s', (fe > .1) & (V < 30)), ('dropped, V >= 30 km/s', (fe > .1) & (V >= 30)), ('all points', fe > -1)):
    out[lab] = {m: hb(v, k) for m, v in e.items()}; out[lab]['robust_clip_0.5'] = {m: hb(v, k, .5)[0] for m, v in e.items()}
    print('%-24s gal %3d pts %4d |' % (lab, out[lab]['T48'][2], out[lab]['T48'][3]), ' '.join('%s %.3f (%+.3f)' % (m, out[lab][m][0], out[lab][m][1]) for m in e), '| clipped', ' '.join('%.3f' % out[lab]['robust_clip_0.5'][m] for m in e))
k = fe > .1; print('dropped points: median V %.1f km/s, median radius %.2f R_e, share inside R_e %.2f, share in first third of the curve %.2f; kept: median V %.1f' % (np.median(V[k]), np.median(x[k]), np.mean(x[k] < 1), np.mean(rank[k] < 1 / 3), np.median(V[~k])))
out['dropped_description'] = dict(median_V=float(np.median(V[k])), median_r_over_Re=float(np.median(x[k])), share_inside_Re=float(np.mean(x[k] < 1)), share_first_third=float(np.mean(rank[k] < 1 / 3)), median_V_kept=float(np.median(V[~k])))
json.dump(out, open('res60_describe.json', 'w'), indent=1)
