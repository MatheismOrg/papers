Exploratory formula search for SPARC galaxies and seven X-COP R500 estimates
2026-10-27

This analysis does not establish a new dark-matter law. It is a small reproducible
comparison of fitting families, including established RAR phenomenology, with
whole-system internal validation and a diagnostic for a mass-derived predictor.

RUN
  python search_formulas.py --input ../dm_fit/results
  python analyze_search.py
Required packages: numpy, scipy, matplotlib.
Inputs: held_out_predictions.csv and cluster_test.json from the earlier fit.
The first script retains their sample selection, units, distances, inclinations,
mass-to-light assumptions, cluster catalog masses, and outer fold assignments.
It records input SHA256 hashes. All generated files are in results/.

SAMPLE
149 galaxies, 2700 radii; seven clusters, one estimate at R500 each. No new data
are introduced. Galaxies span 0.19-108.31 kpc, while cluster radii span
1153-1424 kpc. There are no training observations in the intervening radius gap.
The cluster set is selected by availability of stellar profiles, not a complete
or randomly selected cluster population.

FORMULAS
x=log10(gbar/[1e-10 m/s^2]), z=log10(r/kpc), y=log10(gobs/[1e-10 m/s^2]).
Seven candidate families were specified before the initial search execution:
standard RAR; quadratic in x; full quadratic in x,z; RAR with a radius-dependent
acceleration scale; RAR with a free low-acceleration exponent; RAR with both;
and an explicitly piecemeal RAR with separate galaxy/cluster acceleration scales.
See search_results.json for exact formulas, coefficients, bounds and fold fits.
The RAR variants guarantee positive gpred>=gbar and approach gbar at high gbar
for fixed radius; this is a limited asymptotic property, not a relativistic theory.
Polynomial candidates do not impose those properties. A separate class label in
the two-population model is deliberate and must not be described as continuous
physical unification.

VALIDATION
Five outer folds hold out entire systems. Three inner folds within each outer
training set select among candidates. Inner score is the mean of (1) equal-galaxy
log-acceleration MSE and (2) equal-cluster log-acceleration MSE. Training is repeated
with each system weighted equally and with each population weighted equally.
Outer folds retain the previous analysis assignments. Inner folds are stratified
by population, with random seed 49173+outer_fold. No per-test-system parameter is
fitted. The samples and radius dependence had already been inspected, so this is
exploratory internal validation, not independent discovery confirmation.

Galaxy-only experiments additionally retain the original point-weighted objective
for a direct comparison with earlier results. Equal-system results are supplied
separately; their numerical differences must not be attributed to a changed formula.
All fitting fixes x values and omits full observational covariance/nuisance inference.

MAIN RESULTS
Original-objective galaxy-only standard RAR: 0.13339018 dex CV RMS and 9.21674%
median absolute speed error. Quadratic: 0.13516582 dex and 9.31051%.
Generalized RAR exponent: 0.13318907 dex and 9.26643%. Thus its tiny RMS reduction
comes with slightly worse median speed error. The descriptive paired-galaxy
bootstrap difference in RMS is -0.00020111 dex, interval [-0.00147659,0.00114745].
This bootstrap conditions on existing CV fits; it omits retraining, nuisance and
selection uncertainty and is not a confirmatory test.

The two-population RAR gives a0_gal=1.03681522e-10 m/s^2 and
a0_cluster=3.17876670e-10 m/s^2 with equal-system joint fitting. CV median galaxy
speed error is 9.00224%; cluster acceleration error is 6.35045%. The all-candidate
nested selection selects this model in four outer folds and the radius quadratic
in one. These numbers describe seven R500 points, not full cluster profiles.

CRITICAL DIAGNOSTIC
R500 is derived with mass through an overdensity definition:
M500=(4*pi/3)*500*rho_crit(z_cosmological)*R500^3.
Consequently gobs(R500)=(4*pi*G/3)*500*rho_crit(z_cosmological)*R500.
Thus radius carries information about the target by construction. Holding out
systems does not eliminate this within-system dependency. A radius-only null
g=C r, evaluated with the same five outer folds, has 1.19859% median absolute
cluster acceleration error and 0.00660460 dex RMS. A constant-acceleration null
has 8.56115% median error. These are target-construction diagnostics, not new DM
models. They were added after an independent audit identified this issue.
The radius-quadratic cluster error of 5.63240% cannot support physical unification.
An independent audit also checks the exact relation with cosmological redshift.

No history data were supplied, so no gravitational-memory model was fitted.
No Solar-System, lensing, CMB, structure-growth, or independent survey validation
was performed. There is no demonstrated improvement that warrants calling these
candidate modifications a new universal dark-matter formula.

FILES
search_formulas.py: all candidate fits, nested validation and diagnostic nulls.
analyze_search.py: conditional paired-system bootstrap and two figures.
results/search_results.json: full machine-readable report.
results/*predictions.csv: all outer-fold predictions and residuals.
results/paired_galaxy_bootstrap.json: paired comparisons with standard RAR.
results/galaxy_candidate_errors.png: matched-objective galaxy model comparison.
results/cluster_definition_null.png: cluster definitional-null comparison.

All fourteen joint candidate CSVs were checked independently against the JSON
metrics, row counts, and system counts. The figures were visually inspected.
