"""SA73 T61 all-data looks (no protocol; A, B, C frozen at 11.7481, 0.8178, -0.7119): cores for diffuse gas-rich disks.
M_DM(<r) = M_NFW(<r; M_ref 10^(t w), c 10^(-s w)) x tanh(r / (eta R_e))^n, w = g L (as T57), n = tanh(kappa w) (Read, Agertz
and Collins 2016 use tanh(r/r_c)^n with r_c = 1.75 R_1/2). Fits on the benchmark's kept points, or on all points of the
benchmark galaxies; scored on kept, dropped and all points. Usage: PYTHONPATH=<dmfit> python3 explore61.py"""
import os, sys, json, numpy as np
from scipy.optimize import least_squares
import harness_T57d as M7
import harness_T48 as M8
import compare46 as C
H = M7.H; B = M7.B; D = H.D; S = M7.S; K = H.K; G43 = M7.G43
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
t1 = D.table1(H.PKG / 'data/SPARC/SPARC_Lelli2016c.mrt'); rw = []; flag = []
for g, r in sorted(t1.items()):
    d = np.loadtxt(os.path.join(H.PKG / 'data/SPARC/Rotmod_LTG', g + '_rotmod.dat'), comments='#', ndmin=2)
    R, V, eV, Vg, Vd, Vb = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]
    vb2 = Vg * np.abs(Vg) + 0.5 * Vd * np.abs(Vd) + 0.7 * Vb * np.abs(Vb); k = (R > 0) & (V > 0) & (vb2 > 0)
    Ms = 0.5 * r['L36']; Mp = Ms + 1.33 * r['MHI']; bench = bool(r['Q'] < 3 and r['Inc'] >= 30)
    rw += [dict(name=g, kind='disk', r=a, gb=b, go=c, ms=Ms, mp=Mp, re=r['Reff'], z=0.) for a, b, c in zip(R[k], vb2[k] / R[k] * D.ACC, V[k] ** 2 / R[k] * D.ACC)]
    flag += [(bench, f <= 0.1, v) for f, v in zip((eV / V)[k], V[k])]
a = H.arr(rw); flag = np.array(flag, float); bench = flag[:, 0] > 0; kept = flag[:, 1] > 0; V = flag[:, 2]
kb = bench & kept; assert kb.sum() == 2700 and len(np.unique(a['name'][kb])) == 149
w = M7.gl(a['ms'], a['mp'], a['re'], 0.)
def gp(p):
    t, s, kap, eta = p
    m = G43.mass(th, a['ms'], a['mp'], a['re'], a['z'], 0.) * 10 ** (t * w); c = K.c_dm14(m, a['z']) * 10 ** (-s * w)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); mt = m * H.ff(a['r'] * c / rv) / H.ff(c)
    if kap > 0: mt = mt * np.tanh(a['r'] / (eta * a['re'])) ** np.tanh(kap * w)
    return a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
def hb(e, k):
    un, gi = np.unique(a['name'][k], return_inverse=True); return float(np.sqrt(np.mean(np.bincount(gi, e[k] ** 2) / np.bincount(gi))))
def hw(k):
    un, gi = np.unique(a['name'][k], return_inverse=True); ct = np.bincount(gi); x = np.zeros(len(k)); x[k] = 1 / np.sqrt(len(un) * ct[gi]); return x
SETS = dict(kept149=kb, dropped118=bench & ~kept, all149=bench, dropped_all=~kept, all175=np.ones(len(V), bool), slow_dropped=(~kept) & (V < 30))
def report(lab, g):
    e = np.log10(a['go'] / g); r = {k: round(hb(e, v), 4) for k, v in SETS.items()}; print('%-44s' % lab, ' '.join('%s %.4f' % kv for kv in r.items()), flush=True); return r
out = {}
out['T48'] = report('T48', gp((0, 0, 0, 1.75))); out['T57 (t .361, s .146)'] = report('T57 (t .361, s .146)', gp((.361, .146, 0, 1.75)))
out['MOND'] = report('MOND', C.g_mond(a)); out['STD'] = report('standard halos', C.g_lcdm(a))
def fit(free, k, start, lab):
    """free: indices of (t, s, kappa, eta) to fit; k: points used; host-balanced log residuals as in the benchmark"""
    wt = hw(k); base = np.array([0., 0., 0., 1.75]); lo = np.array([-2., -2., 0., .3]); hi = np.array([2., 2., 20., 6.])
    def fun(q):
        p = base.copy(); p[free] = q; return (wt * np.log10(a['go'] / gp(p)))[k]
    best = None
    for s0 in start:
        o = least_squares(fun, np.array(s0, float), bounds=(lo[free], hi[free]), ftol=1e-10, xtol=1e-10)
        if best is None or o.cost < best.cost: best = o
    p = base.copy(); p[free] = best.x; r = report('%s -> t %.3f s %.3f kap %.3f eta %.2f' % (lab, *p), gp(p)); r['p'] = p.tolist(); out[lab] = r; return p
for setname in ('kept149', 'all149'):
    k = SETS[setname]
    fit([0, 1], k, [(.3, .15), (0, 0)], 'looser (t, s) fit on ' + setname)
    fit([2], k, [(1.,), (3.,)], 'core (kappa) fit on ' + setname)
    fit([0, 2], k, [(.2, 1.), (0, 3.)], 'core + mass (t, kappa) fit on ' + setname)
    fit([0, 1, 2], k, [(.3, .1, 1.), (.2, 0, 3.)], 'core + looser (t, s, kappa) fit on ' + setname)
    fit([0, 2, 3], k, [(.2, 1., 1.75), (.2, 2., 1.)], 'core + mass, eta free (t, kappa, eta) fit on ' + setname)
json.dump(out, open('res61_explore.json', 'w'), indent=1)
