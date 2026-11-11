"""T26: T25's formula (MOGF0: gas switch for the Moster mass + halo floor, no contraction; 4 coefficients) with the IMF
factor of T24 (Posacki et al. 2015 trend, Chabrier scale, delta = 0.215) applied to every spheroid sample with a measured
dispersion: the 50 MaNGA BCGs (<sigma_los>, Tian et al. 2024) as T24, and the 16 SPARC early types (sigma_e from
Cappellari et al. 2013a, ATLAS3D XV, Table 1; their Lelli et al. 2017 mass models are Chabrier-based, M/L 0.8 at 3.6
micron, with no gas term, so g_bar scales by alpha at both radii). Caveat: the early types' inner points are ATLAS3D JAM
circular velocities, which Posacki et al. used, so only their outer (HI) points test the trend independently.
Usage: python3 harness_T26.py <mask> <out> MOGF0IE@0.1   (or through chunk_run.py)"""
import sys
import numpy as np
import harness_T24 as I
I.apply('C')
F = I.F; H, Y = I.H, I.Y
F.FLOOR['MOGF0IE'] = F.FLOOR['MOGF0']; H.MODELS['MOGF0IE'] = H.MODELS['MOGF0']
SIGE = {l.split()[0].replace('UGC06176', 'UGC6176'): 10 ** float(l.split()[1])
        for l in open(H.S / 'ext/atlas3d/Cappellari2013a_Table1_16ETG_rows.txt') if not l.startswith('#')}
assert len(SIGE) == 16
_etg = H.etg_rows
def etg_rows():
    out = _etg()
    for r in out: r['gb'] = r['gb'] * I.alpha(SIGE[r['name']], I.DELTA['C'])
    return out
H.etg_rows = etg_rows
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
