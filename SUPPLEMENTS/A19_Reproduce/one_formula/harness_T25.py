"""T25: the simpler formula with T24's IMF account: T21's MOGF0 (gas switch for the Moster mass + halo floor, no
contraction; 4 coefficients: A, B, C frozen on disks + log10 M_t) with the MaNGA stellar g_bar raised by the Posacki
et al. (2015) IMF factor on a Chabrier scale (harness_T24.py, delta = 0.215).
Usage: python3 harness_T25.py <mask> <out> MOGF0I@0.1   (or through chunk_run.py)"""
import sys
import harness_T24 as I
I.apply('C')
F = I.F; H, Y = I.H, I.Y
F.FLOOR['MOGF0I'] = F.FLOOR['MOGF0']; H.MODELS['MOGF0I'] = H.MODELS['MOGF0']
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
