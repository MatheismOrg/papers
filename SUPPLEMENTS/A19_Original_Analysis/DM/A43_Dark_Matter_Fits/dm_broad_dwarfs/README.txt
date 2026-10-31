External dwarf-spheroidal benchmark
==================================

This is a reproducible population-transfer calculation, not a fit of a new
dark-matter theory. The frozen SPARC RAR is applied to a separately acquired
nearby-dwarf catalog.

INPUTS
McConnachie, A.W. 2012, AJ 144, 4; DOI10.1088/0004-6256/144/1/4,
with the author's January2021 FITS update and October2019 morphology table:
https://www.cadc-ccda.hia.nrc-cnrc.gc.ca/en/community/nearby/
The complete FITS, morphology table text, update log, comments and numbered original
references are retained in data/. The catalog is a heterogeneous compilation,
not a new uniform stellar-velocity analysis. The code does not substitute
rounded or guessed measurements for its entries.

Wolf etal2010, MNRAS406,1220; DOI10.1111/j.1365-2966.2010.16753.x
provides the approximate spherical equilibrium mass estimator. We calculate
new estimates from the updated catalog; these are NOT published Wolf2010
table masses and should not be described as such.

REPRODUCE
Install numpy and astropy. The frozen morphology table extraction is supplied:
  python build_dwarfs.py
No PDF conversion or network download is required for this packaged run.
The script reads the earlier frozen RAR calibration from the sibling path
../dm_fit/results/fit_results.json. Its acceleration scale is
1.1603274962269242e-10 m/s^2. No dwarf data determine that value.

OUTPUTS
results/dwarfs_benchmark.csv: 44 selected systems with standardized columns
population,system,host,radius_kpc,gbar_m_s2,gobs_m_s2 and original observables.
results/dwarfs_predictions.csv: same sample with frozen RAR predictions.
results/excluded_dwarfs.csv: every rejected morphology-matched row and reasons.
results/dwarfs_summary.json: source provenance, hashes, methods and metrics.
results/morphology_lookup.json: exact source lines used for classification.
results/fits_field_descriptions.xml: embedded catalog column documentation.

DEFINITIONS AND ASSUMPTIONS
The published rh is a projected major-axis angular half-light radius. We
circularize it using sqrt(1-ellipticity), convert it using the distance modulus,
and adopt a spherical 3D half-light radius rhalf=4Re/3. This is an approximation
for flattened systems, not an exact deprojection of each individual galaxy.
The quoted stellar dispersion approximates the global luminosity-weighted
dispersion required by the estimator. McConnachie warns that single numbers
can conceal radial trends or multiple components.

gobs=3*sigma_los^2/rhalf and gbar=G*(Mstar/2)/rhalf^2.
Mstar=2*L_V, in solar units. Apparent V magnitudes are generally already
foreground-extinction corrected (McConnachie2012 section4); no second reddening
correction is applied. The adopted solar absolute V magnitude is4.83.
M/L=1and3 alternatives are sensitivity choices, not uncertainty intervals.
Gas is ignored for these dSph systems; entries with catalog HI mass above
10% of nominal stellar mass are excluded when HI is measured. Missing HI
values are allowed by morphology and do not independently establish zero gas.

The operational label ultrafaint means L_V<100000Lsun, independently of fit
residuals; the 25remaining systems are labeled classical. This is not a claim
that every source paper used identical taxonomy. The sample has19ultrafaints.
Systems with ambiguous dSph? classification, missing ellipticity, missing or
unresolved dispersion, strong catalog rotation, explicit selected tidal
disruption/ambiguous-cluster flags, or the Segue2 upper-limit comment are
excluded. Other possible tidal systems are retained with a caution flag.
Higher-precision sensitivity selects max(sigmaerror)/sigma<=0.25 and no
tidal/classification caution flag; it contains25systems.

STATISTICS
At nominal M/L=2, frozen RAR external acceleration RMS is0.655dex for44systems,
0.318dex for25classical systems and0.928dex for19ultrafaints. Median absolute
acceleration fractions are70.8%,50.0%,86.1%, respectively. These do not have
the same meaning as median rotation-speed errors in SPARC. The matching
dispersion error is sqrt(gpred/gobs)-1 under the same estimator. No formal
significance, intrinsic scatter, or new-law exclusion is inferred here.

First-order propagated measurement errors, when available, use the mean of
the quoted asymmetric errors and assume independent inputs. They omit all
structural/Jeans estimator uncertainty, binaries, contamination, tides,
external-field effects and selection systematics. Distance cancels from
stellar gbar, and is not counted twice in its error. Stellar-population M/L
uncertainty is only bracketed by the explicit sensitivity cases.

No gravitational-memory variable has been measured or fitted. No free
per-system accelerations were used. The broad discrepancies are a diagnostic
of transfer limits for this simple estimator and assumed stellar masses;
interpreting them physically requires a more careful population analysis.
