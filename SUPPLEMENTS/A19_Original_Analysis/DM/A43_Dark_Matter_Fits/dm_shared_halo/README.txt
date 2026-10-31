SHARED PHOTOMETRY-TO-HALO RELATIONS: A43 RESEARCH TEST
29 October 2026

Result
------
This finite exploratory test did not improve on the shared refitted exponential
radial-acceleration relation (RAR) for entirely withheld SPARC galaxies. Familiar
NFW and Burkert halo families can be compressed into four or six common
photometric coefficients with roughly comparable predictions. Their original
individual-galaxy fitting success does not automatically become a better common
predictive law. This is useful negative evidence and a reproducible bridge from
object-specific halo templates to common source-state relations, not a new
all-scale gravity theory or evidence that analytic formulas do not exist.

Input and leakage safeguards
----------------------------
The unchanged sample has 149 galaxies and 2700 radial measurements. Independent
SPARC Table 1 luminosity, HI mass, and effective light radius are joined by name.
No Vflat, observed acceleration, observed rotation speed, dynamical mass, or
previous fitted halo parameter enters the new prediction features. The existing
radial baryonic acceleration remains an input to every competing formula.

M_p = 10^9 (0.5 L_3.6 + 1.33 M_HI) solar masses,
where the catalogue luminosity and HI columns are in 10^9 solar units.
m = log10(M_p / 10^10 solar masses), h = log10(R_eff / 3 kpc).
Sigma_p = M_p / (2 pi R_eff^2).

This is a photometric/gas mass proxy, not exact total baryonic mass. In particular
it uses stellar M/L=0.5 for all light, omits molecular gas, and does not implement
a separate bulge correction. The existing radial gbar calculation retains its
disk 0.5 / bulge 0.7 convention. R_eff is a stellar-light radius, not the baryonic
half-mass radius. Distances and other systematics are shared with rotation-curve
inference; target-independent feature construction is not statistical
independence of all measurement errors.

Calibration range: M_p = 4.906e7 to 2.6886829e11 solar masses;
R_eff = 0.26 to 12.36 kpc; Sigma_p = 5.81137e6 to 1.88854e9 solar masses/kpc^2.

Formulas
--------
g_pred = gbar + A f(z) / [z^2 f(1)], z = r/r_s.
A is the extra halo acceleration at r_s; A > 0 and r_s > 0.
NFW: f(z) = ln(1+z) - z/(1+z).
Burkert: f(z) = 0.5 ln(1+z) + 0.25 ln(1+z^2) - 0.5 atan(z).

Mass-only relation:
log10(A / 1e-10 m/s^2) = a0 + aM m;
log10(r_s / kpc) = b0 + bM m.
Mass+size relation adds aR h and bR h, respectively.

All-data coefficients below describe calibration, not coefficients used to score
withheld data. Each validation fit has its own exclusively training coefficients
in results/results.json and diagnostics/selected_fit_bounds.csv.

                      a0          aM          aR          b0          bM          bR
NFW_mass          -0.646347     0.149806     --         1.158997    0.382132      --
Burkert_mass      -0.358201     0.129593     --         0.632975    0.420772      --
NFW_size          -0.595287     0.397430   -0.536789     1.123854    0.070706    0.515723
Burkert_size      -0.320607     0.288899   -0.418629     0.609475    0.307773    0.226284

For example the NFW mass+size relation is approximately
A = 2.539e-11 (M_p/1e10 Msun)^0.39743 (R_eff/3kpc)^(-0.53679) m/s^2,
r_s = 13.300 (M_p/1e10 Msun)^0.07071 (R_eff/3kpc)^0.51572 kpc.
Use the machine-readable full precision coefficients for reproduction.

RAR is gbar / [1 - exp(-sqrt(gbar/a0))]. Fixed RAR uses a0=1.2e-10 m/s^2;
refitted RAR has one shared coefficient. Its all-data a0 is about 1.0368e-10.
The existing seven-coefficient A43 mixed_regular formula is also refitted on the
same disk data and weights; run.py preserves its exact expression and bounds.

Validation
----------
The objective is root mean squared log10 acceleration residual, with each galaxy
carrying equal total weight and radii equally weighted within each galaxy. The
three fivefold seeds are 20261027, 7331, and 48105. Entire galaxies are excluded.
The four tail holdouts exclude the lowest or highest 37 galaxies by M_p or Sigma_p.
Strata are defined from predictor values only. Every tail galaxy lies outside
the training convex hull in (log M_p, log R_eff), as expected for these tails.
These are covariate extrapolation checks, not independent surveys.

Primary metric: equal-galaxy RMS log10 acceleration (dex), lower is better.

Model             CV seed 1  CV seed 2  CV seed 3   Low mass High mass Low Sigma High Sigma
RAR_fixed         0.154556   0.154556   0.154556    0.200280  0.122613  0.200228  0.121202
RAR               0.152898   0.153013   0.153367    0.197281  0.122522  0.196530  0.124710
mixed_regular     0.154169   0.154865   0.156944    0.251312  0.126207  0.198595  0.129456
NFW_mass          0.156294   0.157426   0.160311    0.210548  0.129576  0.202988  0.132567
Burkert_mass      0.160531   0.161005   0.165277    0.220238  0.130968  0.211596  0.133269
NFW_size          0.155601   0.154317   0.158489    0.199438  0.140581  0.187161  0.146192
Burkert_size      0.159399   0.160618   0.164852    0.208953  0.138456  0.199116  0.159388

The primary NFW_size minus RAR difference is +0.002703 dex; the descriptive
galaxy-bootstrap 95% interval is [-0.009930,+0.014848]. Its low-Sigma tail has a
favorable point estimate (-0.009370 dex), but interval [-0.039363,+0.021357].
It loses on the high-mass and high-Sigma tails. No reliable improvement claim
follows from the favorable low-Sigma point estimate. These intervals condition
on frozen predictions: they omit refitting, model selection, measurement errors
and shared systematics, and are not simultaneous confidence statements.

Primary median pointwise speed errors are 8.89% for refitted RAR, 9.06% for
NFW_size, and 8.80% for mixed_regular. These are row medians, unlike the primary
equal-galaxy RMS objective; reporting the favorable median alone would reverse
the conclusion without justification. In-sample NFW_size RMS is 0.148933 dex
versus RAR 0.152488, but that improvement disappears for withheld galaxies.

All 320 optimization attempts converged, spanning 120 fit sets, with all attempts
retained. Convergence does not prove global optimality. No selected halo or RAR
fit is within 1e-4 of its coefficient bound range. All 20 selected mixed_regular
fits touch or approach at least one bound. Bounds were not widened after scoring.

Physical scope
--------------
These formulas specify positive spherical dark-halo templates added to baryonic
gravity. They are not modified-gravity field equations. NFW has a density cusp;
Burkert has a core. Both masses diverge logarithmically at infinite radius without
an outer cutoff; a physically global model must specify truncation and boundary
conditions. No cosmological concentration prior or abundance-matching constraint
is imposed here. This does not compare against a full cosmological halo model.

At fixed R_eff and nonzero r, the fitted all-data mass exponents formally make
the extra acceleration vanish as M_p tends to zero. This narrowly defined check
does not establish path-independent empty-source behavior, a relativistic theory,
or a Solar-System limit. It extrapolates the mass relation below its calibrated
domain. Stellar and Solar systems lack the calibrated galaxy halo population
and light-radius meaning; no ephemeris likelihood or cosmological observable is
tested here. Current light and gas cannot identify a causal gravitational memory
kernel or a formation history. No such kernel was fitted or inferred.

What can be concluded
---------------------
Known halo families can reproduce much of galaxy-to-galaxy diversity through
common relations to measured baryonic mass and size. Such relations already have
substantial literature; no priority claim is made. Within this sample, allowing
these four/six shared coefficients did not beat one shared RAR coefficient on
whole-galaxy prediction. The more promising next test is a physical halo-state
model with a published mass-concentration relation, baryon efficiency, intrinsic
scatter and consistent nuisance treatment, scored with the same withheld-galaxy
design. It must be evaluated as a further model, not silently substituted here.
New confirmation requires an independent dataset or fully nested model selection.

Files and reproduction
----------------------
PROTOCOL.txt: frozen before the first numerical fit; finite model set and bounds.
PARSER_CORRECTION.txt / parse_failure.log: initial input parser failed before any
fit; correction changed parsing only. Both pre-fit manifests are preserved.
run.py: fitting, exact folds, target exclusions, all prediction output and bootstrap.
input.csv / galaxy_features.csv: exact analyzed rows and independent inputs.
results/results.json: all scores, coefficients, starts, train rows, physical checks.
results/*.csv: 56 complete prediction tables, retaining every model and holdout.
diagnostics.py: read-only post-fit summaries and plot; does not perform fits.
diagnostics/: scores, per-galaxy errors, halo states, bounds, extrapolation and plot.
SOURCES.txt: primary source provenance and prior work.
manifest.json: hashes for all included files except itself.

From the original workspace, run `python tmp/dm_shared_halo/run.py` and then
`python tmp/dm_shared_halo/diagnostics.py`. Dependencies: Python, NumPy, SciPy,
Matplotlib. run.py needs sibling dm_unified_regular/augmented_input.csv and
dm_fit/data/SPARC_Lelli2016c.mrt. Their hashes are in results/results.json; input.csv
and galaxy_features.csv retain the analyzed data for independent reconstruction.
The raw original Table 1 is not copied here; it remains at the recorded sibling
path. The parser uses 19 whitespace fields because displayed byte ranges did not
match actual row spacing. Name, luminosity, Reff and HI are tokens 0,7,9,13.
