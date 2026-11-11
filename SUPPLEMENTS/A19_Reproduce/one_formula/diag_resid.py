"""All-data residual pattern of a model (not held out): host-mean log10(g_obs/g_pred) by population and kind, against
host properties. Usage: python3 diag_resid.py <harness module> <model> <th...>"""
import sys, importlib
import numpy as np
M = importlib.import_module(sys.argv[1]); H = M.H if hasattr(M, 'H') else M.Z.H
B = __import__('harness_T15b')
model = sys.argv[2]; th = np.array([float(x) for x in sys.argv[3:]])
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
gs = H.disks(); cl = H.cluster_rows()
pops = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()),
        'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
for p, a in pops.items():
    er = np.log10(a['go'] / H.gpred(model, th, a))
    for kind in np.unique(a['kind']):
        k = a['kind'] == kind; hosts = np.unique(a['name'][k]); hm = np.array([np.mean(er[k & (a['name'] == n)]) for n in hosts])
        ms = np.array([a['ms'][k & (a['name'] == n)][0] for n in hosts]); re = np.array([a['re'][k & (a['name'] == n)][0] for n in hosts])
        lo = ms < np.median(ms); sm = np.log10(ms / (2 * np.pi * re ** 2))
        hi_s = sm > np.median(sm)
        print(f'{p:7s} {kind:13s} n={len(hosts):3d} mean {np.mean(hm):+.3f} rms {np.sqrt(np.mean(hm ** 2)):.3f} | low M* {np.mean(hm[lo]):+.3f} high M* {np.mean(hm[~lo]):+.3f}'
              f' | low Sigma* {np.mean(hm[~hi_s]):+.3f} high Sigma* {np.mean(hm[hi_s]):+.3f} | corr(logM*) {np.corrcoef(np.log10(ms), hm)[0, 1]:+.2f} corr(logSigma*) {np.corrcoef(sm, hm)[0, 1]:+.2f}', flush=True)
