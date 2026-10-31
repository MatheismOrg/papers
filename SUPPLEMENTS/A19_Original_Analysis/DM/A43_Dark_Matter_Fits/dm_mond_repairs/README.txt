A43 targeted MOND comparison and finite-range repairs
28 October 2026

This directory contains a completed exploratory test of seven fixed candidate
families on 2,931 acceleration estimates in 286 systems (285 validation groups).
No claim of beating the best published gravity/dark-matter theories follows.
The same five population balance, source data, nuisance assumptions and
cross-validation grouping as the previous A43 stage are retained.

Main balanced log10-acceleration scores (dex, smaller is better):
                          group-CV  whole-population-exclusion
RAR                         .329913      .359343
simple MOND                 .329076      .358612
standard MOND               .336787      .369869
nu_n shape fit              .329715      .378094
previous mixed_regular      .221975      .289205
mixed_transient             .221968      .290574
mixed_potential_taper       .227937      .328599
two_log_bumps               .232828      .374428

The simple/standard/nu_n/RAR rows are common algebraic interpolations fitted
across all five populations; they do not include full AQUAL/QUMOND geometry,
external fields or additional unseen cluster matter. Thus they are not the
strongest fully specified versions of MOND. Dark-halo field models and
cosmological likelihoods are not evaluated by this module.

No numerical breakthrough over previous mixed_regular was found. The transient
repair's primary CV change is -0.000007228 dex (paired frozen-prediction 95%
interval [-0.000026556,+0.000009870]); population-transfer change is +0.001369072
dex (interval [-0.000451673,+0.002921965]). These intervals exclude refitting,
model selection, measurement-systematic and reconstruction uncertainty.

The physically useful result is a simpler correction of the previous asymptote:
 a_eff/a_ref = C0 + (C1 U^-alpha S W + C2 s^beta)/(1+r/R),
 g_pred = g_b/[1-exp(-sqrt(g_b/a_eff))],
 a_ref=1e-10 SI, U=(g_b/a_ref)(r/kpc)^2, s=(g_b/a_ref)(r/kpc),
 S=r^4/(r^4+ell^4), W=1/[1+(g_b/1e-8)^2].
Augmented all-data fit:
 C0=.8809886676623764,
 C1=.0007748296513830285,
 C2=.40173444675522413,
 alpha=.95, beta=.95, ell=.010 kpc, R=307.65765105499855 kpc.
Alpha and beta hit their upper bounds, ell its imposed 10 pc lower bound.
Moving the inverse-mass term inside the existing radial cutoff adds no parameter
and restores V_inf^4=GM*a_ref*C0 for every fixed isolated point mass. The previous
formula instead retained a mass-dependent far-field normalization. No negative
point-source inferred extra density occurs for any allowed parameters, as shown
analytically in ANALYTIC_CHECKS.txt. This is not a full field theory.

The asymptotic normalization is 8.80989e-11 SI, 26.6% below the common nominal
1.2e-10 SI. Restoring a mass-independent asymptotic coefficient does not validate
that normalization or the observed finite-radius baryonic Tully-Fisher relation.
For a 1e6-solar-mass point source, at 10/100/1000 kpc the predicted V^4/(GM*a_inf)
is about 1.445/1.343/1.107; convergence can be remote. The 1-solar-mass isolated
10,000 AU acceleration remains 1.78655 times Newtonian; external-field/binary
and Solar quadrupole likelihoods have not been fitted.

Matched disk-only recalibration of established interpolations (149 SPARC hosts,
2,700 rows, equal-host five-fold CV) gives:
 model     group RMS dex    median speed error
 RAR          .152898             8.8926%
 simple       .153016             8.9783%
 standard     .157322             9.9218%
 nu_n         .153743             8.8480%
The mixed_transient formula calibrated across all populations gives SPARC
.187290 dex and 15.2629% median speed error. Different calibration objectives
must be stated; this comparison does not claim a like-for-like training design.
These current matched-procedure disk benchmarks should not be confused with
historical A43 disk analyses using different weighting/folding conventions.

Failures retained:
- Potential taper does worse and its s0 parameter hits the lower bound 1e-8.
  Across a fixed-radius response grid of 144,981 points, 683/144,800 slope
  intervals are negative (minimum logarithmic slope -0.1268). This is an extreme
  extrapolation failure: no negative slopes occur at the 2,931 observed inputs.
- Two bounded log-potential bumps do worse in population transfer. Their same
  grid has 2,900 negative response intervals (minimum -2.3161), including
  negative derivatives at 37 of the 2,931 observed inputs; a point-source
  grid has 613/24,000 negative inferred-extra-mass slope intervals.
- Generalized nu_n selects n=.7836 in the all-data fit; n<1 implies a growing
  absolute high-g residual, despite a vanishing fractional correction. This is
  a useful negative benchmark, not a Solar-compatible modification.

There are 239 fit sets, 471 optimizer starts, zero failed starts, and 79
prediction CSV tables (63 five-population tables and 16 disk-only tables).
Every coefficient, start, score and split is in results/results.json. The
primary source input and pre-fit protocol/run hashes are retained. Families
were fixed before this stage's fits; the overall study is adaptive because
previous data/results informed their design.

Reproduction requires Python with NumPy and SciPy:
 python run.py
 python diagnostics.py
run.py reuses previous mixed_regular prediction files for paired comparisons
from sibling ../dm_unified_regular/results/augmented; all new fits are made
from this directory's input.csv. It verifies previous input and fold identities.
Diagnostics scripts use no external source or network. Zero underflowed Solar
excesses retain their finite log10 excesses in solar_stress.csv.

Files:
 PROTOCOL.txt, pre_fit_manifest.json: design frozen before numerical fitting.
 input.csv: unchanged earlier augmented input.
 run.py, run.log, results/: complete fit calculation and predictions.
 diagnostics.py, diagnostics.json: physical grids and asymptotic coefficients.
 solar_stress.csv, smooth_core.csv, finite_radius_btfr.csv: diagnostic curves.
 ANALYTIC_CHECKS.txt: limiting behavior and positivity proof with qualifications.
 primary_sources.txt: prior published interpolation definitions.
 manifest.json: hashes of deliverable research files (not itself).

No measured temporal history enters these seven functions. This module neither
finds gravitational memory nor uses failure of a smooth fit as evidence that
spacetime is discrete.
