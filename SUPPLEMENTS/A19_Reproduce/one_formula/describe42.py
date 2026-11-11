"""T42 (PROTOCOL_T42.txt): the T37 formula with the concentration relation chosen by A19_CM (dm14, diemer19, ludlow16):
A, B, C refitted on all disks, then M_t; all-data scores (primary mask). Usage: A19_CM=<model> python3 describe42.py"""
import os, json, warnings
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import harness_T37 as S
H = S.H; K = H.K; B = S.B
warnings.filterwarnings('ignore')
CM = os.environ.get('A19_CM', 'dm14')
if CM != 'dm14':
    from colossus.halo import concentration
    LG, ZG = np.arange(6., 16.5 + 1e-9, .05), np.arange(0., 1.5 + 1e-9, .05)
    tab = np.array([np.log10(concentration.concentration(10 ** LG * S.HH, '200c', z, model=CM)) for z in ZG]).T
    itp = RegularGridInterpolator((LG, ZG), tab, bounds_error=False, fill_value=None)
    def cnew(m, z):
        m, z = np.broadcast_arrays(np.asarray(m, float), np.asarray(z, float))
        return 10 ** itp(np.c_[np.log10(m).ravel(), z.ravel()]).reshape(m.shape)
    K.c_dm14 = cnew
    _mref = H.mref
    def mref(model, th, ms, mp, re):
        m, c = _mref(model, th, ms, mp, re)
        return (m, cnew(m, 0.)) if model == 'ORIG' else (m, c)
    H.mref = mref
def score(model, th, a):
    l = H.losses(model, th, a); return round(float(np.sqrt(np.mean(list(l.values())))), 4)
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs))
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
lens = H.Lens(gs, False)
abc = np.array(H.fit('ORIG', None, d, ext, lens, [0, 1, 2])['th']); B.FROZEN['abc'] = abc
f = H.fit('MOGF0IES', .1, d, ext, lens, [0, 1, 2, 3]); th = np.array(f['th'])
mstar = 5.43e10; a_ = np.linspace(1e-3, 60, 6000); discs = {'HI': (53.1, 7., 4., .085), 'H2': (2180., 1.5, 12., .045)}
mp = mstar + sum(2 * np.pi * np.sum(a_ * p[0] * 1e6 * np.exp(-p[2] / a_ - a_ / p[1])) * (a_[1] - a_[0]) for p in discs.values())
mw = float(S.mass(th, np.array([mstar]), np.array([mp]), np.array([5.75]), np.array([0.]), 0.)[0])
print(json.dumps(dict(cm=CM, abc=np.round(abc, 4).tolist(), logMt=round(float(th[0]), 4), bound_hits=f['bound_hits'], disk=score('MOGF0IES', th, d),
      **{k: score('MOGF0IES', th, a) for k, a in ext.items()}, lens=round(float(np.sum(lens.white('MOGF0IES', th, [0, 1, 2, 3]) ** 2)), 2),
      MW_halo=float('%.3g' % mw), MW_c=round(float(K.c_dm14(np.array([mw]), np.array([0.]))[0]), 2))), flush=True)
