"""SA73 T61 all-data looks, second pass: the trade between the benchmark's kept points and the dropped ones as the core
strength kappa grows (t refitted on the kept points at each kappa), and fits of (t, kappa) on all points of the benchmark
galaxies with each point weighted by 1 / (sigma_V^2 + sigma_0^2), sigma_V = (2 / ln 10) eV / V, sigma_0 = 0.131 dex (the
galaxy-to-galaxy scatter of Section 5.5). Usage: PYTHONPATH=<dmfit> python3 explore61b.py"""
import json, numpy as np
src = open('explore61.py').read().split("\nout = {}")[0]; exec(compile(src, 'explore61_defs', 'exec'))
from scipy.optimize import least_squares, minimize_scalar
t1e = {}; fe = []
for g, r in sorted(t1.items()):
    d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
    R, Vv, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (Vv > 0) & (vb2 > 0); fe += list((eV / Vv)[k])
fe = np.array(fe); sig = 2 / np.log(10) * fe; S0 = 0.131; pw = 1 / (sig ** 2 + S0 ** 2)
def whb(e, k):
    """error-weighted, host-balanced rms: per galaxy sum(w e^2)/sum(w), averaged over galaxies"""
    un, gi = np.unique(a['name'][k], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(gi, pw[k] * e[k] ** 2) / np.bincount(gi, pw[k]))))
def rep(lab, g):
    e = np.log10(a['go'] / g); r = {k: round(hb(e, v), 4) for k, v in SETS.items()}; r['w_all149'] = round(whb(e, SETS['all149']), 4); r['w_all175'] = round(whb(e, SETS['all175']), 4)
    print('%-40s' % lab, ' '.join('%s %.4f' % kv for kv in r.items()), flush=True); return r
out = {}
out['T48'] = rep('T48', gp((0, 0, 0, 1.75))); out['MOND'] = rep('MOND', C.g_mond(a)); out['STD'] = rep('standard halos', C.g_lcdm(a)); out['T57'] = rep('T57', gp((.361, .146, 0, 1.75)))
kb_ = SETS['kept149']; wt = hw(kb_)
for kap in (0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5):
    o = minimize_scalar(lambda t: float(np.sum((wt * np.log10(a['go'] / gp((t, 0, kap, 1.75))))[kb_] ** 2)), bounds=(-1, 2), method='bounded')
    out['scan kappa %.1f' % kap] = dict(t=float(o.x), **rep('kappa %.1f, t %.3f (t on kept)' % (kap, o.x), gp((o.x, 0, kap, 1.75))))
ka = SETS['all149']; un, gi = np.unique(a['name'][ka], return_inverse=True); nw = np.zeros(len(fe)); nw[ka] = np.sqrt(pw[ka] / np.bincount(gi, pw[ka])[gi] / len(un))
for free, start, lab in (([0, 2], [(.2, .5), (.1, 1.)], 'weighted (t, kappa) on all149'), ([2], [(.5,), (1.,)], 'weighted (kappa) on all149'), ([0, 1], [(.3, .15), (.6, .3)], 'weighted T57 form (t, s) on all149')):
    base = np.array([0., 0., 0., 1.75]); lo = np.array([-2., -2., 0., .3]); hi = np.array([2., 2., 20., 6.]); best = None
    def fun(q):
        p = base.copy(); p[free] = q; return (nw * np.log10(a['go'] / gp(p)))[ka]
    for s0 in start:
        o = least_squares(fun, np.array(s0, float), bounds=(lo[free], hi[free]), ftol=1e-10, xtol=1e-10)
        if best is None or o.cost < best.cost: best = o
    p = base.copy(); p[free] = best.x; out[lab] = dict(p=p.tolist(), **rep('%s -> t %.3f s %.3f kap %.3f' % (lab, p[0], p[1], p[2]), gp(p)))
json.dump(out, open('res61_explore_b.json', 'w'), indent=1)
