A43 targeted MOND cluster repair: a published gas-associated source component
28 October 2026. Exploratory reuse of existing cluster mass reconstructions.

PRIOR ART AND EQUATIONS
B. Famaey, L. Pizzuti & I. D. Saltas, Physical Review D 111,123042 (2025),
https://doi.org/10.1103/dccw-srks, Eq17; https://arxiv.org/html/2410.02612v2 .
The authors' parameter archive is https://doi.org/10.5281/zenodo.15299349 .
Gas profiles are credited to Vincenzo Salzano; no downloaded source code runs.

rho_extra(r)=eta rho_hotgas(r) exp(-r/Rc)
M_extra(<r)=eta integral_0^r exp(-s/Rc) dM_hotgas(s)
g_source=g_b+G M_extra/r^2; g_pred=nu(g_source/a0) g_source.
a0 is fixed to1.2e-10 m/s2. Simple nu=0.5+sqrt(0.25+1/y) is primary, as in
Famaey et al.; exponential RAR nu=1/(1-exp(-sqrt(y))) is a sensitivity.
Only eta and Rc are refitted, shared among all training clusters. This is a
phenomenological additional source mass, not a theory explaining its origin.
Do not call this density ansatz new or treat it as baryons-only modified gravity.

SCOPE
16 CLASH hosts have new public gas-profile parameters and65existing A43 points.
Four other CLASH hosts are excluded solely for unavailable gas parameters.
Seven X-COP hosts supply21existing500/750/1000kpc points. Total23hosts86points.
The old g_b and g_obs are unchanged. The new CLASH gas model is an extra source
predictor: its enclosed gas mass / old enclosed total-baryon mass ranges below1
for every point, median0.6811 and maximum0.9478. This does not guarantee equality
of the two studies' baryonic inference. Cluster geometry and gas-profile priors
remain those of their source models. Target values still depend on lensing-NFW
or hydrostatic modeling; no direct shear, X-ray or velocity likelihood is fitted.
The limited point set samples14.3--1000kpc and does not retest every central
radius and outer turning point examined by the source paper.

X-COP MGAS is directly read, never FGAS multiplied by inferred total mass. The
stored gas radius R/R500 is converted to kpc using its archive header. This is a
coordinate recovery; evaluated radii remain fixed500/750/1000kpc, not target-
defined overdensity radii. Gas and HSE inference nevertheless share observations,
cosmology and systematic uncertainties. Inner gas below firstgridpoint is given a
constant-density completion; rigorously bounding its entire weightedmass between
M_min*exp(-r_min/Rc) and M_min changes log predictions by at most0.000558dex over
all testedcoefficients. Numerical quadrature doubling changes≤7.8e-15dex here.

RESULTS (simple MOND primary; acceleration residual RMS in dex)
                                    Combined host-out     CLASH→XCOP  XCOP→CLASH
Fixed galaxy a0, no extra gas         0.440585               0.358108    0.509891
Cluster-trained shared a0            0.140111               0.214358    0.214735
Published-scale eta10,Rc400kpc        0.114564               0.079125    0.141382*
Two shared gas coefficients           0.113307               0.061197    0.139847

*The published-scale parameters were informed by CLASH. Their fixed predictions
on CLASH are descriptive, not independent source exclusion. CSV split labels for
fixed models only indicate which rows are displayed; they are not new fits.

The source-held-out refitted model trains on CLASH alone for X-COP predictions:
eta=11.774..., Rc=288.85...kpc and medianaccelerationerror7.63%; no X-COP target
enters that coefficient fit. Reverse source transfer medianerror16.59%. In the
combined all-data fit eta=12.01466..., Rc=268.799...kpc. Exact logged coefficients
are authoritative. Combined host-out medianaccelerationerror15.10% versus25.13%
for cluster-trained a0. Exponential sensitivity gives essentially same results.
CLASH-only host-out favors one-a0 (0.112979) over gas repair (0.147941), while
X-COP-only favors gas repair (0.068404 vs0.096740). Do not hide this reversal.

Across nine X-COP target/stellar sensitivities, repair gives0.03548--0.06878dex.
Smaller NFW/Einasto-target residuals are agreement with smooth reconstructions,
not new accuracy on raw observations. Those targets are sensitivity cases only.

Paired host-resampling of fixed predictions gives combined repair-minus-one-a0
95% interval[-0.06199,+0.00682]dex, repair-minus-publishedfrozen
[-0.01304,+0.01164]dex. Neither establishes a clear combined improvement. The
CLASH→XCOP difference versus source-trained one-a0 is[-0.17641,-0.12036]dex;
versus publishedfrozen[-0.03129,+0.00115]dex. These descriptive intervals exclude
refitting, modelselection, measurementuncertainty and source-systematic errors.

This is a useful successful transfer test of an existing model family. It does
not show superiority to the best published cluster model, and has not been fitted
to all other populations. Missing hot-gas profiles in disks/dwarfs are never set
to zero to manufacture a cross-scale claim. The source law is shared but requires
full hot-gas information and remains additional mass phenomenology.

PHYSICAL LIMITS
For eta>=0 and positive gas density the added mass/density is nonnegative.
The exponential cutoff gives finite total extra mass for the supplied power-law
outer gas profiles. The central fitted gas slopes alpha<=1 have finite enclosed
mass and no divergent added Newtonian force (it is constant if alpha=1). This
statement concerns the extra component, not an assumed complete central star
profile. Scaling all sources to zero makes g_pred vanish. Finite total isolated
source mass still has the usual MOND1/r farfield, normalized by baryons PLUS
extra mass. This is not universal baryonic BTFR preservation for gas-richclusters.
No lensing metric completion, external-field effect or cosmological evolution is
calculated. Hot-gas association in equilibrium does not explain collisionless
merger offsets.

REPRODUCTION
Place this directory beside original dm_fit, dm_fixed_radius and
 dm_unified_regular directories in the A43 supplement. Python requires numpy,
scipy and astropy (local dm_fit/deps is used when available).
  python run.py
  python diagnostics.py
No network is required. clash_gas_parameters.csv is a frozen numerical extraction
from the author archive. source_provenance.json records the originalsourcehash;
extract_parameters.py can optionally re-fetch and re-extract to audit it. Article
PDFs and third-party source .py files are intentionally not redistributed.

Outputs: inputs.csv, gas_profiles.csv, quadrature.json, results.json (allfits),
predictions.csv (192 logical tables,6128rows), provenance.json,
numerical_diagnostics.csv, diagnostics.json. PROTOCOL.txt was written beforefits.
Results include460fitsets/1610optimizerstarts,zero failures,no activeboundhits.
