"""T45 (PROTOCOL_T45.txt, sha256 193809a2...): T44's formula (MOGF0IESR) with log10 M_t set, in each partition, on the
seven red-central bins of Mandelbaum et al. (2016) instead of the benchmark: chi2 with asymmetric errors on a grid from
10.0 to 12.5 (step 0.01, refined to 0.001), red centrals with M_p = M*, R_e = 0.5 x the SPARC size relation, z_eff and
no selection, M200c -> M200m [h^-1] with the Dutton and Maccio concentration, as in the T38 check. MOGF0IESM: 3
coefficients fitted to the benchmark (A, B, C on the training disks)."""
import sys, csv
import numpy as np
import harness_T44 as P
S = P.S; H = P.H; K = H.K; Y = P.Y; B = S.B
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
MODEL = 'MOGF0IESM'
H.MODELS[MODEL] = H.MODELS['MOGF0']
MB = [r for r in csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')) if r['color'] == 'red']
assert len(MB) == 7
(P1, P0), (Q1, Q0) = H.pop_coeff(H.disks())
def red_pred(lmt, row):
    lm, z = float(row['logMs_eff_VAGC']), float(row['z_eff']); ms = 10 ** lm; re = 10 ** (P0 + P1 * lm) * .5
    m = float(P.G.mass(np.array([lmt]), np.array([ms]), np.array([ms]), np.array([re]), np.array([z]), 0.)[0])
    return L8.to_m200m(m, K.c_dm14(m, z), z)
def red_dev(lmt): return [red_pred(lmt, r) - float(r['logM200m_hinv']) for r in MB]
def red_chi2(lmt):
    return float(sum((d / float(r['err_plus'] if d > 0 else r['err_minus'])) ** 2 for d, r in zip(red_dev(lmt), MB)))
_CAL = {}
def calibrate():
    key = tuple(np.round(np.asarray(B.FROZEN['abc'], float), 10))
    if key not in _CAL:
        grid = np.round(np.arange(10.0, 12.5 + 1e-9, .01), 3); c2 = np.array([red_chi2(g) for g in grid]); k = int(np.argmin(c2))
        fine = np.round(np.arange(grid[k] - .01, grid[k] + .01 + 1e-9, .001), 4); fine = fine[(fine >= 10.0) & (fine <= 12.5)]
        cf = np.array([red_chi2(g) for g in fine]); j = int(np.argmin(cf))
        _CAL[key] = (float(fine[j]), float(cf[j]), k in (0, len(grid) - 1))
    return _CAL[key]
_fit = H.fit
def fit(model, weight, dtr, etr, lens, groups):
    if model != MODEL: return _fit(model, weight, dtr, etr, lens, groups)
    lmt, c2, edge = calibrate()
    return dict(th=[lmt], loss=c2, attempts=[], bound_hits=['logMt'] if edge else [], calibration='Mandelbaum2016_red_chi2',
                abc=[float(x) for x in B.FROZEN['abc']])
H.fit = fit
_gpred, _esd = H.gpred, H.esd
def gpred(model, th, a): return _gpred('MOGF0IESR' if model == MODEL else model, th, a)
def esd(model, th, st): return _esd('MOGF0IESR' if model == MODEL else model, th, st)
H.gpred, H.esd = gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
