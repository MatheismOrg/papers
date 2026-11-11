"""T35 (PROTOCOL_T35.txt): the Milky Way with the contracted T18 formulas against the uncontracted frozen formula."""
import sys, io, contextlib
import numpy as np
with contextlib.redirect_stdout(io.StringIO()): import diag_milkyway as MW      # reuse baryons, data and helpers (its prints suppressed)
H, K = MW.H, MW.K; G = MW.G; R, vobs = MW.R, MW.vobs; v2b = MW.v2b; mstar, mp = MW.mstar, MW.mp
k = R < 21
def F(eta): return .45 + .38 * (eta + 1.16) ** .53
for lab, (A, B_, C, s) in {'frozen, no contraction': (11.7481, .8178, -.7119, 0.), 'CAUD, s = 0.25 (disks)': (11.6762, .7774, -.6769, .2468),
                           'CAU5, s = 0.5': (11.6039, .7599, -.6638, .5), 'CAU1, s = 1 (published)': (11.4187, .7019, -.6724, 1.)}.items():
    m = mp / H.FB + 10 ** A * (mstar / 1e10) ** B_ * (5.75 / 3) ** C
    c = K.c_dm14(m, 0.); rv = (H.G * m / (100 * H.hz(0.) ** 2)) ** (1 / 3); mnfw = m * H.ff(R * c / rv) / H.ff(c)
    mbar = v2b * R / G; eta = mbar / (H.FB * mnfw)
    v2d = G * (1 - H.FB) * mnfw * np.maximum(F(eta), 1.) ** s / R if s > 0 else G * (1 - H.FB) * mnfw / R
    vm = np.sqrt(v2b + v2d); rat = vm[k] / vobs[k]
    print('%-26s M_halo %.2e  r_s %.1f kpc  mean ratio %.3f  rms %.3f dex  v(8.19) %.0f  v(20.27) %.0f' % (lab, m, rv / c, np.mean(rat), np.sqrt(np.mean(np.log10(rat) ** 2)), vm[R == 8.19][0], vm[R == 20.27][0]))
