"""FORECAST_T60 Part A: T57 against T48, all frozen, on the SPARC points the benchmark drops (velocity error above 10
percent). Usage: PYTHONPATH=<dmfit> python3 dropped60.py"""
import os, json, numpy as np
src = open('test59.py').read().split("\ntest = [")[0]; exec(compile(src, 'test59_defs', 'exec'))
t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); rw = []; inb = {}
for g, r in sorted(t1.items()):
    d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) > 0.1) & (vb2 > 0)
    Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']; inb[g] = bool(r['Q'] < 3 and r['Inc'] >= 30)
    rw += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
un, ms, mean, ct, info = per_galaxy(rw); w = np.array([info[n]['gL'] for n in un]); b = np.array([inb[n] for n in un])
out = dict(all=summary(un, ms, mean, ct), benchmark_galaxies=summary(un, ms, mean, ct, b), excluded_galaxies=summary(un, ms, mean, ct, ~b),
           strong_term=summary(un, ms, mean, ct, w > np.median(w[w > 0])), weak_or_no_term=summary(un, ms, mean, ct, w <= np.median(w[w > 0])), galaxies=info)
for k in ('all', 'benchmark_galaxies', 'excluded_galaxies', 'strong_term', 'weak_or_no_term'): print(k, json.dumps(out[k]))
json.dump(out, open('res60_dropped.json', 'w'), indent=1)
