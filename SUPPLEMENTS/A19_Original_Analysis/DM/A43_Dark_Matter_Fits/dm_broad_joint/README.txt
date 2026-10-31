BROAD ACCELERATION FITS - 27 October 2026

PURPOSE
Compare a finite list of explicit formulas across 2,910 measurements or inferred
acceleration estimates in 279 objects, grouped into278 conservative validation
groups:149 SPARC rotating galaxies,44 dwarf spheroidals,16 early-type galaxies,
50 MaNGA brightest cluster galaxies, and20 CLASH lensing clusters. The CLASH data
contain84 inner/outer measurements; early-type galaxies have32 paired measurements.
These are heterogeneous derived quantities, not279 identical raw experiments.
The smallest selected radius is27.1pc; the largest is600kpc. Baryonic acceleration
spans5.18e-15 to8.13e-9m/s^2. This is not a fit to all DM observables.

INPUTS
build_input.py reads standardized sibling analysis products. Each input hash and
original path is retained in input_metadata.json. combined_input.csv is sufficient
to rerun the joint fitting; acquisition scripts remain in the source subpackages.
Source annotations/error columns are retained but not used as a full likelihood.
SPARC: fixed published distance/inclination and stellar M/L.5/.7.
dSph: catalog-based Wolf half-light estimates, nominal stellar M/L2; spherical
equilibrium, dispersion and binary/tidal assumptions remain. Source M/L1/3 and
dispersion-precision sensitivities are provided by the dwarf subpackage.
ETG: published inner JAM/outerHI estimates with photometric M/L.8.
MaNGA BCG: published Jeans/Abel isotropic estimates. Two nearby BCGs (8989-12703
and8991-6102) are conservatively folded/weighted together; their common cluster
identity is not established. The50 objects are not asserted to be50 independent
cluster environments.
CLASH: NFW spherical mass reconstructions from strong/weak lensing+magnification,
with published cosmology/GR assumptions; each cluster's BCG and outer points share
a validation group. MaNGA and CLASH redshifts do not overlap.

RESERVED OBSERVABLES
The earlier X-COP fixed-radius hydrostatic profiles remain a separate stress test.
The R500 points are excluded because their radii encode the target mass.
KiDS/GAMA stacked lensing has no unique radius per RAR bin and is tested separately.
Group velocity dispersions test a global virial relation, not an invented local
acceleration. Globular central-dispersion controls, Solar-System and cosmological
constraints are likewise not mixed into this scalar regression.

REPRODUCE
Install Python packages numpy scipy matplotlib.
python fit_broad.py
python analyze_broad.py
python taper_sensitivity.py
The first command uses the supplied combined_input.csv and input_metadata.json;
it does not require network or the data-acquisition libraries.
Optional input rebuilding: python build_input.py (requires sibling datasets).
Optional post-search higher-precision dwarf sensitivity:
python prepare_sensitivity.py
python fit_broad.py --quick --input combined_high_precision.csv --metadata high_precision_metadata.json --output results_high_precision

METHOD
See PROTOCOL.txt for the candidate families and bounds. Training gives equal total
weight to five populations, equal weight to validation groups within each, and
equal weight to radial measurements within a group. Classical and ultrafaint dSphs
share one training population. Five outer folds hold entire groups out. Three
inner folds select among shared forms. Whole-population and MW/M31 dwarf-host
holdouts provide harder extrapolation checks. Observational covariance, shared
calibration and nuisance parameters are not marginalized. These are exploratory
internal comparisons of previously published/inspected data, not independent
confirmation of a newly chosen physical model.

MAIN RESULTS
The original SPARC RAR's frozen external transfer has median acceleration errors
70.8%dSph,34.1%ETG,48.4%MaNGA BCG and68.8%CLASH. The SPARC data retain a tight
relation; the added populations expose substantial systematic differences/scatter.
A common RAR re-calibrated equally across populations sacrifices the SPARC fit.
The six-coefficient logarithmic quadratic improves balanced group-CV RMS from
0.3360 to0.2355dex, but worsens to0.4883dex when an entire population is withheld.
It is a descriptive interpolation, not a physically constrained force law.

A compact three-coefficient potential-proxy correction gives group-CV RMS0.24946
and population-holdout RMS0.31581dex:
  x=log10(gbar/[1e-10m/s^2]); u=log10(radius/kpc); q=x+u
  a_eff=1e-10 *10^(0.06672835+0.35008148*q+0.17209066*q*q) m/s^2
  gpred=gbar/[1-exp(-sqrt(gbar/a_eff))]
q measures log(gbar*r/[3.08567758149e9m^2/s^2]); it is a baryonic potential-scale
proxy, not the full gravitational potential and not a new independent input.
Group-CV median acceleration errors:29.10%SPARC,39.62%ETG,53.13%dSph,21.14%MaNGA
BCG,34.43%CLASH. SPARC median speederror14.29% is worse than the original9.2%.
The proxy correction's coefficients are common, but it does not fit all phenomena
to9%or establish a new theory. Its raw high-g asymptotics can be unphysical.

POST-SEARCH TAPER
taper_sensitivity.py divides log10(a_eff/1e-10) by1+(gbar/1e-8)^2, using the same
coefficients and no refit. The gcut1e-8m/s^2 is chosen, not inferred. Group-CV and
population-holdout RMS remain0.24951/0.31589dex. This restores a central high-g
Newtonian limit; it is not a Cassini/external-field/relativistic/CMB calculation.
It was added after the main search and is not an independently selected model.

POST-SEARCH DWARF PRECISION SENSITIVITY
Keeping25dSph with quoted dispersion errors at most25% and no listed caution,
while retaining all other populations, gives potential-proxy coefficients
[0.14102344,0.31211354,0.16034434]. Its group-CV/population-holdout RMS are
0.21322/0.28721dex, versus0.29968/0.33583 for common RAR. Dwarf median acceleration
error becomes31.59%, while SPARC is31.08%. Positive proxy curvature remains under
this quality cut, but this is an exploratory sensitivity, not validation on new
observations. Full data and results are in combined_high_precision.csv,
high_precision_metadata.json and results_high_precision/.

FILES
results/broad_results.json: coefficients,folds,allmetrics,nested selections,
whole-population and dwarf-host holdouts.
results/*predictions.csv: reproducible cross-validated rows.
results/compact_summary.json: concise main-model comparison.
results/taper_sensitivity.json: post-search regulator calculation.
results/verification.json: independent CSV metric/count/fold verification.
Three PNG figures show cross-population patterns, error comparisons and harder
population holdout performance. Every percentage in the model-error figure is
an acceleration percentage. Frozen SPARC calibration is in-sample on SPARC itself,
external on the other populations. Later model rows use internal group validation.

NO MEMORY CLAIM
No formation histories were supplied or fitted. Radius/potential/density proxies
are not temporal memory measurements. Data/model systematics and missing physical
inputs remain viable sources of residual patterns. The fit is an empirical target
for future theory and independent tests, not a discovery of all dark matter.
