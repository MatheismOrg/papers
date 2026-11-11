"""SA73 T59, after the fact (not part of FORECAST_T59): the same test with looser cuts on the velocity error, since the
benchmark's 10 percent cut leaves only three of the eight outside galaxies. Usage: PYTHONPATH=<dmfit> python3 explore59_cuts.py"""
import json, numpy as np
src = open('test59.py').read().split("\ntest = [")[0]; exec(compile(src, 'test59_defs', 'exec'))
test = [n for n, g in gal.items() if g['sparc'] is None]; out = {}
def rows_cut(names, cut):
    o = []
    for nm in names:
        g = gal[nm]
        if g['f36'] is None or g['inc'] < 30: continue
        R = np.array(g['R']); V = np.array(g['Vc']); eV = np.array(g['eVc'])
        Ms = 0.5 * M.l36(g['f36'][0], g['D']) / KL; Re = 1.678 * g['read']['Rstar'] * g['D'] / g['read']['D'] / KR
        vg2, mhi, mcut = M.gas(R, g['Sd'], 'exp'); vb2 = vg2 + M.stars_v2(R, Ms, Re / 1.678)
        k = (R > 0) & (V > 0) & (eV / np.where(V > 0, V, 1) <= cut) & (vb2 > 0); Mp = Ms + 1.33 * mhi
        o += [dict(name=nm, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=Re, z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * H.ACC, V[k] ** 2 / R[k] * H.ACC)]
    return o
for cut in (0.1, 0.2, 0.3, 9.):
    un, ms, mean, ct, info = per_galaxy(rows_cut(test, cut)); s = summary(un, ms, mean, ct); c = summary(un, ms, mean, ct, np.isin(un, CLEAN))
    out['cut_%g' % cut] = dict(all=s, clean=c, galaxies=info)
    print('cut', cut, {k: s[k] for k in ('n', 'points', 'T48', 'T57', 'MOND', 'STD', 'pct', 'boot95', 'p_not_better', 'n_better', 'n_worse', 'T48_within', 'T57_within')})
    print('   clean', {k: c[k] for k in ('n', 'points', 'T48', 'T57', 'MOND', 'STD', 'pct', 'boot95', 'n_better', 'n_worse')})
for n, v in out['cut_9']['galaxies'].items(): print('%-9s' % n, {k: v[k] for k in ('points', 'gL', 'T48_rms', 'T57_rms', 'MOND_rms', 'STD_rms', 'T48_mean', 'T57_mean')}, 'Ms %.2e Mp %.2e Re %.2f' % (v['Ms'], v['Mp'], v['Re']))
json.dump(out, open('res59_explore_cuts.json', 'w'), indent=1)
