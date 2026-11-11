"""All-data what-if (not held out): the same IMF factor for the 16 SPARC early types (Lelli et al. 2017, Chabrier-based
M/L 0.8 at 3.6 micron), with sigma_e from Cappellari et al. 2013a (ATLAS3D XV, Table 1). Their mass models have no gas
term (Vgas = 0), so g_bar is all stars and scales by alpha at both radii. Caveat: the inner points are ATLAS3D JAM circular
velocities, and Posacki et al. 2015 fitted their IMF trend partly to ATLAS3D JAM models, so the inner check is not
independent; the outer (HI) points are. MOGF at the all-data coefficients."""
import io, zipfile, numpy as np, harness_T21 as F
H = F.H; B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
sig = {l.split()[0].replace('UGC06176', 'UGC6176'): 10 ** float(l.split()[1]) for l in open(H.S / 'ext/atlas3d/Cappellari2013a_Table1_16ETG_rows.txt') if not l.startswith('#')}
with zipfile.ZipFile(H.AST / 'inputs/data/ETG/Rotmod_ETG.zip') as zf:
    vg = {n.split('/')[-1].replace('_rotmod.dat', ''): np.loadtxt(io.BytesIO(zf.read(n)))[:, 3] for n in zf.namelist() if n.endswith('_rotmod.dat')}
print('gas terms all zero:', all(np.all(v == 0) for v in vg.values()), '| names matched:', sorted(vg) == sorted(sig))
a = H.arr(H.etg_rows()); th = np.array([11.182185, .173541]); s = np.array([sig[n] for n in a['name']])
def sc(b, kind=None):
    k = np.ones(len(b['r']), bool) if kind is None else b['kind'] == kind; c = H.sub(b, k); l = H.losses('MOGF', th, c)
    return np.sqrt(np.mean(list(l.values()))), np.mean(np.log10(c['go'] / H.gpred('MOGF', th, c)))
print('sigma_e range %.0f-%.0f km/s' % (s.min(), s.max()))
for lab, delta in (('as published', None), ('Chabrier 0.215', .215), ('Kroupa 0.180', .180)):
    b = dict(a)
    if delta is not None: b['gb'] = a['gb'] * 10 ** (.38 * np.log10(s / 200) - .06 + delta)
    print('%-15s all %.4f (%+.3f) | inner %.4f (%+.3f) | outer %.4f (%+.3f)' % ((lab,) + sc(b) + sc(b, 'ETG_inner') + sc(b, 'ETG_outer')))
