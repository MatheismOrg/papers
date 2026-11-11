"""T16: KNEESFZ + Moline subhalo concentration for dwarf satellites + (1+z)^D2 evolution of the high-mass branch.
Imports harness_T15c (and through it T15b and the hashed T15 harness) unchanged. Usage: python3 harness_T16.py <mask> <out> cfg..."""
import sys, json, csv, time, hashlib
from pathlib import Path
import numpy as np
import harness_T15c as Cc
B = Cc.B; H = B.H

_arr = H.arr
def arr(rows):
    a = _arr(rows); a['fac'] = np.array([r.get('fac', 1.) for r in rows], float); return a
def sub(a, k):
    b = {n: v[k] for n, v in a.items() if n != 'hw'}
    un, inv, ct = np.unique(b['name'], return_inverse=True, return_counts=True); b['hw'] = 1 / np.sqrt(len(un) * ct[inv]); return b
def dwarf_rows():
    bt = H.PKG / 'blind_test'; rules = json.load(open(bt / 'frozen_disk_rules.json'))['frozen_on_all_disks']
    rr, src = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(bt))
    aa = H.V4.arrays(rr, rules, ml=2., major=False); pp = H.V3.predict(rr, rules, ml=2., major=False)
    out = [dict(name=x['key'], kind='dwarfs', r=pp['r'][i], gb=pp['gb'][i], go=10 ** aa['lo'][i], ms=pp['Ms'][i], mp=pp['Ms'][i],
                re=pp['Re'][i], z=0., fac=float(pp['fac'][i])) for i, x in enumerate(rr) if not aa['ul'][i]]
    assert len(out) == 75; return out
H.arr, H.sub = arr, sub

NEW = {  # name: (base for the rule and branch, subhalo factor?, z-evolution?, joint?)
    'SUBH': (True, False, False), 'ZEVO': (False, True, False), 'BEST1': (True, True, False), 'BEST1J': (True, True, True)}
H.MODELS['SUBH'] = H.MODELS['KNEESF']
H.MODELS['ZEVO'] = (['A2', 'B2', 'C2', 'D2'], [9, 0, -2, -3], [16, 5, 3, 3], [[12.06, 2.63, 2.5, 0.], [12.06, 2.63, 2.5, 1.], [12.3, 2.3, 2.0, -1.]])
H.MODELS['BEST1'] = H.MODELS['ZEVO']
H.MODELS['BEST1J'] = (['A', 'B', 'C', 'A2', 'B2', 'C2', 'D2'], [9, -.5, -2.5, 9, 0, -2, -3], [14.5, 2, 1.5, 16, 5, 3, 3],
                      [[11.7481, .8178, -.7119, 12.06, 2.63, 2.5, 0.], [11.7481, .8178, -.7119, 12.06, 2.63, 2.5, 1.], [11.7, .8, -.8, 12.3, 2.3, 2.0, -1.]])
_gpred, _esd = H.gpred, H.esd
def mass(model, th, ms, mp, re, z):
    joint = NEW[model][2]; abc = th[:3] if joint else B.FROZEN['abc']; br = th[3:] if joint else th
    m = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]
    d2 = br[3] if NEW[model][1] else 0.
    return m + 10 ** br[0] * (ms / 1e11) ** br[1] * (re / 10) ** br[2] * (1 + z) ** d2
def gpred(model, th, a):
    if model not in NEW: return _gpred(model, th, a)
    m = mass(model, th, a['ms'], a['mp'], a['re'], a['z']); c = H.K.c_dm14(m, a['z'])
    if NEW[model][0]: c = c * a.get('fac', 1.)
    rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3); x = a['r'] * c / rv
    return a['gb'] + H.G * (1 - H.FB) * m * H.ff(x) / H.ff(c) / a['r'] ** 2 * H.ACC
def esd(model, th, st):
    if model not in NEW: return _esd(model, th, st)
    m = mass(model, th, st.ms, st.mp, st.re, st.z); c = H.K.c_dm14(m, H.K.Z0)      # isolated lenses: no subhalo factor
    rc = H.K.rho_crit(H.K.Z0); rv = (3 * m / (4 * np.pi * 200 * rc)) ** (1 / 3); rs = rv / c
    rhos = 200 / 3 * c ** 3 / H.ff(c) * rc
    e = (1 - H.FB) * rs[:, None, None] * rhos[:, None, None] * H.K.g_wb(st.R / rs[:, None, None]) / 1e12
    return st.w @ (np.sum(st.aw * e, axis=2) + st.bary)
H.gpred, H.esd = gpred, esd
H.dwarf_rows = dwarf_rows

if __name__ == '__main__':
    Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
