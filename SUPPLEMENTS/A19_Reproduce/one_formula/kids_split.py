"""Unfitted check: the KiDS-1000 early/late lensing split of isolated galaxies (Brouwer et al. 2021, Fig. 8), which no fit
in turns 15 to 19 used. Early and late stacks have the same (matched) stellar-mass distribution (Fig. 7 digitization);
early types are gas-poor with van der Wel sizes, late types carry cold gas, as in the paper's KiDS model (A19_KiDS_v9).
For a model and coefficients: the predicted split S = log10 mean(ESD_early / ESD_late) over the inner bins, chi-square of
the red, blue, bulge and disc profiles (inner bins, their covariance blocks), and of the measured early-minus-late
difference against the model's. Import after the harness that defines the model."""
import json
import numpy as np
import harness_T15_hashed as H
K = H.K
KD = H.PKG / 'data/Lensing/KiDS'
def setup(gs):
    pop = K.populations(json.load(open(KD / 'brouwer2021_fig7_digitized.json')), 'primary')
    coeff = H.pop_coeff(gs)
    st = {k: H.Stack(pop['fig8_' + k], coeff) for k in ('early', 'late')}
    rows = [r for _, rr in pop['fig8_late'] for r in rr]
    lmg = float(np.log10(sum(n * K.m_gal(lm) for lm, n, _ in rows) / sum(n for _, n, _ in rows)))
    k = np.asarray(K.inner_mask(lmg), bool)
    prof = lambda f: K.load_profile(KD / f)['esd']
    D = {'blue': prof('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'), 'red': prof('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'),
         'disc': prof('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'), 'bulge': prof('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt')}
    Cc, _ = K.load_cov(KD / 'Fig-8_Colorbins_covmatrix_compact.json'); Cs, _ = K.load_cov(KD / 'Fig-8_Sersicbins_covmatrix_compact.json')
    return dict(st=st, k=k, lmg=lmg, D=D, C={'colour': Cc, 'sersic': Cs})
def check(model, th, S):
    E = H.esd(model, th, S['st']['early']); L = H.esd(model, th, S['st']['late']); k = S['k']
    out = dict(S_model=round(float(np.log10(np.mean(E[k] / L[k]))), 4))
    for lab, lo, hi in (('colour', 'blue', 'red'), ('sersic', 'disc', 'bulge')):
        C = S['C'][lab]; dL, dE = S['D'][lo], S['D'][hi]
        cL, cE = C[:15, :15], C[15:, 15:]
        chi = lambda d, m, c: float((d - m)[k] @ np.linalg.inv(c[np.ix_(k, k)]) @ (d - m)[k])
        Cd = cE + cL - C[15:, :15] - C[:15, 15:]
        out[lab] = dict(S_data=round(float(np.log10(np.mean(dE[k] / dL[k]))), 4),
                        chi2_late=round(chi(dL, L, cL), 2), chi2_early=round(chi(dE, E, cE), 2),
                        chi2_difference=round(chi(dE - dL, E - L, Cd), 2), chi2_difference_null=round(chi(dE - dL, 0 * E, Cd), 2),
                        n_inner=int(k.sum()))
    return out
