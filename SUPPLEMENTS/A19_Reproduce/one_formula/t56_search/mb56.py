"""SA73 T56 look: the final formula (T48) against all 14 Mandelbaum et al. (2016) lensing halo masses of SDSS centrals
(7 red used to set M_t, 7 blue not used), with z-scores, and the three terms of the formula. Usage: python3 mb56.py <dmfit dir>"""
import sys, csv, json; sys.path.insert(0, sys.argv[1])
import numpy as np
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; S = M.S; K = H.K; G43 = M.P.G; L8 = M.L8 if hasattr(M, 'L8') else None
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
if L8 is None:
    sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
(p1, p0), (q1, q0) = H.pop_coeff(H.disks())
rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
def pred(lm, z, red, abc=None, lmt=None, extra=None):
    ms = 10 ** lm; re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
    m = float(G43.mass(np.array([th[0] if lmt is None else lmt]), np.array([ms]), np.array([mp]), np.array([re]), np.array([z]), 0.)[0])
    return L8.to_m200m(m, K.c_dm14(m, z), z), m, ms, mp, re
if __name__ == '__main__':
    c2 = {'red': 0., 'blue': 0.}
    print('color logM*   obs  +err -err  formula   dev    z   | log terms: Mp/fb  stellar  q*M_E | M200c')
    for r in rows:
        lm, z, obs, ep, em = float(r['logMs_eff_VAGC']), float(r['z_eff']), float(r['logM200m_hinv']), float(r['err_plus']), float(r['err_minus'])
        red = r['color'] == 'red'; p, m, ms, mp, re = pred(lm, z, red); dev = p - obs; zs = dev / (ep if dev > 0 else em); c2[r['color']] += zs ** 2
        abc = B.FROZEN['abc']; t1 = mp / H.FB; t2 = 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2]; t3 = max(m - t1 - t2, 1.)
        print('%-5s %5.2f %6.2f %.2f %.2f %8.2f %+6.2f %+5.1f  |  %6.2f %7.2f %7.2f | %6.2f' % (r['color'], lm, obs, ep, em, p, dev, zs, np.log10(t1), np.log10(t2), np.log10(t3), np.log10(m)))
    print('chi2 red (7, used for M_t):', round(c2['red'], 1), ' blue (7):', round(c2['blue'], 1))
