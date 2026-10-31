A43 population-scatter continuation -- 29 October 2026

Run from a complete extracted A43 archive, Python with numpy/scipy/pandas:
  OPENBLAS_NUM_THREADS=1 python3 dm_population_scatter/audit_pre.py
  OPENBLAS_NUM_THREADS=1 python3 dm_population_scatter/run.py
  OPENBLAS_NUM_THREADS=1 python3 dm_population_scatter/audit_post.py
  OPENBLAS_NUM_THREADS=1 python3 dm_population_scatter/moment_sensitivity.py
The four workers in run.py operate locally. No network or new data is needed.

The frozen protocol precedes these fits. The subsequent nonlinear-moment
sensitivity was motivated by the audit, keeps all fitted coefficients frozen,
and is explicitly post-audit. No independent astronomical confirmation occurred.

INPUTS
149 galaxies, 2700 radii; byte-preserved previous input and original split records
are read from dm_feedback_profiles. Error and component restoration uses the
original SPARC table and Rotmod_LTG.zip in dm_fit/data. matched_observations.csv
records all restored values. It checks both observed and baryonic accelerations
against the archived selection to 1e-12 relative agreement. No Vflat is a predictor.

MODELS AND SCORES
run.py implements analytic source/nuisance derivatives and Gaussian covariance
C=diag(sigma_loggobs^2+s^2)+UU^T. Shared nuisance directions: distance, inclination,
and a common stellar multiplier. Fixed concentration scatter is 0.11 dex; the
optional excess-mass scatter has a fitted width. This is an empirical, first-order
conditional population model, not full nonlinear Bayesian inference. Hyperparameter
uncertainty is not integrated. The mean joint negative log density per radial
measurement is averaged equally across galaxies. This weighted composite score
preserves within-galaxy covariance. Ordinary equal-galaxy RMS is also retained.

The test galaxy's observed rotation errors enter the evaluation covariance, as
in a conditional observational likelihood; the prediction's central halo state
never uses its rotation curve. Disk/bulge stellar fluctuations are approximated
by one multiplier. Inclination only propagates the kinematic correction. Unmodeled
photometric covariance, gas errors and other systematics remain. Error estimates
are not an additional physical predictor. Bounds and starts are in the protocol.

RESULTS
summary.json contains every model's complete three-repeat, tail and dwarf scores,
all-data coefficients, conditional paired intervals and convergence accounting.
All 140 fitted predictor sets and all numerical attempts converged. The 440 stored
attempt records include 20 exact comparator endpoints, not 440 numerical runs.

NFW with mass/concentration scatter improves mean held-out log score over the
specified RAR-plus-offset model by 0.18981 natural-log units per radius, with
conditional whole-galaxy interval [-0.30374,-0.09123]. Central RMS instead worsens:
0.15883 versus 0.15461 dex. Absolute 95% coverage is 96.1%, but 68% coverage is
79.2%; do not claim uniform calibration. No full MOND/CDM likelihood ranking follows.

NONLINEAR CHECK
The full-sample excess-mass width is 0.548 dex, with split fits sometimes reaching
the 0.6 bound. The audit's 20/32-point Gaussian quadrature changes median log-field
means by 0.012 dex and halo-only SD by up to 17.5% versus linearization. Its full
row-level checks are saved. The additional moment_sensitivity.py uses 32x32 nodes
on each excluded galaxy, replaces the linear halo mean/covariance by nonlinear
moments, and keeps all coefficients frozen. It preserves the log-score difference
(-0.19285; interval [-0.30455,-0.09521]). This still treats the mixture as Gaussian,
and keeps nuisance covariance linearized; it is not an exact nonlinear likelihood.

AUDIT
Analytic derivatives match independent finite differences to 2.6e-10 absolute;
Woodbury scores match direct matrices to 2.3e-15; old NFW predictions match to
2.2e-15 dex. audit_post.py checks exclusions and scores independently. These are
numerical checks, not validation of physical assumptions or nuisance priors.

Published sources/credit: Dutton & Maccio 2014, https://arxiv.org/pdf/1402.7073;
Li et al. 2020, https://arxiv.org/html/2001.10538v2; SPARC/Lelli et al. 2016.
The source authors retain credit. No article PDFs or third-party software included.
