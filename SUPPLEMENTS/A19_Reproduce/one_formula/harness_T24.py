"""T24: is the MaNGA offset the stellar initial mass function? T21's MOGF with the MaNGA BCGs' stellar g_bar raised by the
IMF factor that massive early-type galaxies show: log alpha = 0.38 log10(sigma_e / 200 km/s) - 0.06 relative to Salpeter
(Posacki et al. 2015, ATLAS3D + SLACS, rms 0.12 dex), converted to the IMF of the published masses with
log10(M_Salpeter / M_IMF) = 0.215 for Chabrier or 0.180 for Kroupa (Madau and Dickinson 2014: divide Chabrier masses by
0.61, Kroupa by 0.66). Each BCG's <sigma_los> from Tian et al. (2024), Table 1, stands in for sigma_e. Only g_bar of the
50 MaNGA points changes (gas is under 3.2 percent of their baryons); stellar masses that enter the halo formula stay on
the published scale, like the disks the rule was set on. No coefficient is added; nothing else changes.
Use through harness_T24C.py (Chabrier) or harness_T24K.py (Kroupa), which set the factor before the data load."""
import re
import numpy as np
import harness_T21 as F
H = F.H; Z = F.Z; Y = F.Y
DELTA = {'C': .215, 'K': .180}
for k_ in ('MOGFI', 'MOGFK'):
    F.FLOOR[k_] = F.FLOOR['MOGF']; H.MODELS[k_] = H.MODELS['MOGF']
def sigmas():
    sig = {}
    for l in open(H.PKG / 'data/Clusters/manga_bcg_table1_extracted.txt'):
        if not l.strip(): continue
        ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
        assert len(v) == 14; sig['MaNGA_' + ident.replace('∗', '').replace('*', '')] = v[8]
    assert len(sig) == 50; return sig
def alpha(s, delta): return 10 ** (.38 * np.log10(s / 200.) - .06 + delta)
_rows = H.cluster_rows
def apply(which):
    delta = DELTA[which]; sig = sigmas()
    def rows():
        out = _rows(); n = 0
        for r in out:
            if r['kind'] == 'MaNGA_BCG': r['gb'] = r['gb'] * alpha(sig[r['name']], delta); n += 1
        assert n == 50; return out
    H.cluster_rows = rows
