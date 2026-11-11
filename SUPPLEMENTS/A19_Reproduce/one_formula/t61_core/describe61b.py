"""SA73 T61, after the held-out run (description, frozen all-data values t 0.1501, kappa 0.6229): all 175 SPARC galaxies,
the points the benchmark drops and all points, and the scatter within disks on the benchmark's points.
Usage: PYTHONPATH=<dmfit> python3 describe61b.py"""
import json, numpy as np
import harness_T61 as M
import harness_T57d as M7
import harness_T48 as M8
import compare46 as C
H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339]); M.U['t'], M.U['k'] = 0.1501, 0.6229; M7.U['t'], M7.U['s'] = 0.361, 0.146
A = M.ALL; V = np.sqrt(A['go'] * A['r'] / H.ACC); bench = np.isin(A['name'], [g['name'] for g in H.disks()])
e = dict(T48=np.log10(A['go'] / H.gpred(M8.MODEL, th, A)), T61=np.log10(A['go'] / H.gpred('MOGF0IESBC', th, A)), T57=np.log10(A['go'] / H.gpred('MOGF0IESBD', th, A)), MOND=np.log10(A['go'] / C.g_mond(A)), STD=np.log10(A['go'] / C.g_lcdm(A)))
def hb(v, k):
    un, gi = np.unique(A['name'][k], return_inverse=True); ct = np.bincount(gi); ms = np.bincount(gi, v[k] ** 2) / ct; mn = np.bincount(gi, v[k]) / ct
    return dict(rms=round(float(np.sqrt(ms.mean())), 4), mean=round(float(mn.mean()), 4), within=round(float(np.sqrt(np.mean(np.maximum(ms - mn ** 2, 0)))), 4), galaxies=len(un), points=int(k.sum()))
out = {}
for lab, k in (('benchmark points (149)', bench & (A['fe'] <= .1)), ('dropped, all 175', A['fe'] > .1), ('all points, all 175', A['fe'] > -1), ('dropped, slower than 30 km/s', (A['fe'] > .1) & (V < 30)), ('dropped, 30 km/s or faster', (A['fe'] > .1) & (V >= 30)),
               ('kept, all 175', A['fe'] <= .1)):
    out[lab] = {m: hb(v, k) for m, v in e.items()}; print('%-30s' % lab, 'gal %d pts %d |' % (out[lab]['T48']['galaxies'], out[lab]['T48']['points']), ' '.join('%s %.4f (%+.3f; w %.4f)' % (m, r['rms'], r['mean'], r['within']) for m, r in out[lab].items()))
w = M.gl(A['ms'], A['mp'], A['re'], 0.); un, ix = np.unique(A['name'], return_index=True); wg = w[ix]; nn = np.tanh(0.6229 * wg)
out['core_strength'] = dict(galaxies_with_core=int(np.sum(wg > 0)), n_median_nonzero=round(float(np.median(nn[wg > 0])), 3), n_above_half=int(np.sum(nn > .5)), n_max=round(float(nn.max()), 3)); print(out['core_strength'])
json.dump(out, open('res61_describe_b.json', 'w'), indent=1)
