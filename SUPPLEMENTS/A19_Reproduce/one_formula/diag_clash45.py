"""SA73 T45 exploration: the formula's halo mass for each CLASH cluster (M_E with S: M200c > 5e14) against Umetsu et al.
(2016) Table 2 (M200c in 1e14 h70^-1 Msun, converted to h = 0.671), and the CLASH outer residuals."""
import csv, numpy as np
import harness_T44 as M
H = M.H; S = M.S
B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
cl = [r for r in H.cluster_rows() if r['kind'] != 'MaNGA_BCG']; a = H.arr(cl); th = np.array([11.1182])
er = np.log10(a['go'] / H.gpred('MOGF0IESR', th, a))
U = {r['host']: r for r in csv.DictReader(open(H.PKG / 'data/Clusters/umetsu2016_table2.csv'))}
h70 = 0.671 / 0.7
out = []
for n in np.unique(a['name']):
    k = a['name'] == n; ms, z = a['ms'][k][0], a['z'][k][0]
    me = float(np.atleast_1d(S.m_exp(np.array([ms]), np.array([z]), np.array([5e14])))[0])
    u = U.get(n); mu = float(u['M200c_1e14_Msun_h70inv']) * 1e14 / h70 if u else np.nan
    ro = er[k & (a['kind'] == 'CLASH_outer')]; rc = er[k & (a['kind'] == 'CLASH_centre')]
    out.append((n, z, np.log10(ms), np.log10(me), np.log10(mu), float(np.mean(ro)), float(np.mean(rc)), u['route'] if u else '-'))
for o in sorted(out, key=lambda t: t[4]):
    print(f'{o[0]:18s} z {o[1]:.3f} logM*BCG {o[2]:.2f}  logM_E {o[3]:.2f}  logM_Umetsu {o[4]:.2f}  diff {o[4]-o[3]:+.2f}  res outer {o[5]:+.3f} centre {o[6]:+.3f} {o[7]}')
d = np.array([o[4] - o[3] for o in out]); print('mean diff', np.nanmean(d), 'median', np.nanmedian(d))
ro = np.array([o[5] for o in out]); print('corr(diff, outer res)', np.corrcoef(d, ro)[0, 1])
