External early-type, weak-lensing, and galaxy-group checks
2026-10-27

SCOPE
This directory extends the frozen SPARC RAR to other actual public observations.
No parameters are fitted to these external samples. The acceleration scale is
1.1603274962269242e-10 m/s^2 from the prior 149-galaxy/2700-point fit. Numerical
comparison to different kinds of observations requires different measurement
models. We therefore do not collapse every system into a local disk acceleration
pair or claim validation of all dark-matter phenomena.

REPRODUCTION
Python packages: numpy, scipy, matplotlib. The published group table extraction
is supplied as a table-only excerpt. Input URLs are listed in provenance.json.
  python validate_external.py
The script reads the frozen fit from ../dm_fit/results/fit_results.json.
The six KiDS/GAMA profile and covariance files used by this calculation are
already unpacked in data/kids/. Their release README is included unchanged.
No PDF conversion or archive download is required for this packaged run.
Input SHA256 hashes are recorded in the result JSON and curation manifest.

EARLY-TYPE GALAXIES
Official source ETG_Lelli2017.mrt has 16 galaxies, each with two acceleration
pairs: an inner stellar Jeans-anisotropic-model (JAM) estimate and an outer HI
measurement. These are 32 points, not 25 galaxies or all early-type systems in
the paper. Radii were obtained from the matching Rotmod_ETG.zip two-point files.
The catalog acceleration pairs are used directly (rounded logs); reconstructing
them from the rounded radii and velocities agrees within 0.03 dex.

Baryonic acceleration assumes the source population prescription M/L=0.8 at
3.6micron for both disk and bulge. This is photometric, not the JAM dynamical
mass-to-light estimate. The inner total acceleration nevertheless depends on
JAM assumptions: axisymmetry, nearly oblate velocity ellipsoid, total mass
following light, constant r-band total M/L. It is not a direct circular-velocity
measurement. The outer points are HI velocity measurements. Source gbar/gobs
uncertainties in dex are retained. No uncertainty-weighted or nuisance refit.

Frozen RAR: all32 median acceleration error34.053%, speed-equivalent18.793%,
RMS0.208134dex. OuterHI16 alone: acceleration23.341%, speed11.782%, RMS0.147941.
InnerJAM16: acceleration44.721%, speed-equivalent20.939%, RMS0.254466.
Radius range0.17-28.61kpc. The inner and outer results must be separated because
their measurement models differ. Descriptive proxy correlations are reported,
but these share gbar with the residual and are not causal evidence.

WEAK LENSING
Official KiDS release from Brouwer2021 provides stacked Excess Surface Density
profiles and full covariance. Each profile has15 acceleration bins, not15
independent galaxies. gbar is the published baryonic-acceleration bin coordinate.
gobs is derived using their SIS approximation, Eq7: gobs=4G DeltaSigma.
Multiplicative shear bias corrections are applied to signal, errors, and full
covariance exactly as in the release README. h70=1. No unique physical radius is
assigned to a stack's gbar bin; radius_kpc is intentionally blank.

The primary external comparison uses the isolated GAMA-lens profile because
spectroscopic redshifts improve isolation selection. All15 bins give median
acceleration error25.864%, RMS0.181078dex, statistical-covariance chi2=18.544
for15 bins and no refitted parameters. Broad errors make median fractions a
poor stand-alone goodness-of-fit measure.

The isolated KiDS-lens profile gives median25.367%, but belowgbar~1e-13 its
isolation selection becomes unreliable. Restricting to the source's7 bins
above that threshold gives22.085% median error, statistical chi2=46.527/7.
Adding an independent0.1dex SIS-conversion uncertainty as a sensitivity check
reduces that chi2 to5.324/7. This is not a fully validated systematic covariance.
The release's nominal hot-gas variant and +/-0.2dex baryonic mass shifts are
included as sensitivity cases. They are not independent observations or proof
that a discrepancy has been physically resolved. KiDS and GAMA lensing analyses
overlap in sky/source data and must not be summed as independent experiments.

This comparison retains Brouwer2021's original assumptions. Later work, e.g.
Mistele et al.2024 arXiv:2310.15248, revises deprojection, baryonic masses and
isolation treatment; no claim is made that the original conversion is unique.

GALAXY GROUPS
Milgrom2019 TableI reproduces observed member counts, velocity dispersions,
harmonic radii and K-band luminosities from Makarov & Karachentsev2011. We test
the separate deep-MOND global virial relation
sigma_pred^4=(4/81) G Mbar a0,
using the frozen SPARC a0 and fixed Mbar/L_K=1 in solar units. This is a global
pressure-supported-system prediction, not the local disk RAR with an arbitrary
sigma^2/r conversion. Harmonic radii span62-738kpc in the source-selected rows.
Their values are retained as metadata, not fitted inputs.

The paper/caption state67total/56selected, but v2TableI contains66printed rows
and55unstarred rows. Both PDF table pages were visually inspected. We include
only these actual rows, preserving every star exclusion marker; no missing row
is fabricated. The source excludes10richgroups/clusters and NGC4216 for bad
distance. Our55unstarred rows give18.845% median dispersion error and0.117275dex
RMS. All66rows give20.508% and0.137549dex. These errors are dispersion errors,
not acceleration or rotation-speed errors. No individual error bars for sigma
are supplied in the table. Stellar-only mass ignores unmeasured gas; source
sample selection, distances, anisotropy, virialization, and isolation matter.
Published rounded sigma predictions at the source a0=1.2e-10 were independently
reconstructed towithin1.6km/s, guarding the table extraction and unit conversion.

FILES
results/etg_external_pairs.csv: 32 rows with standard radius/gbar/gobs columns.
results/lensing_external_pairs.csv:45 rows for three overlapping/sensitivity
profiles; use primaryGAMA15 alone or separately compareKiDS7 reliable bins.
results/*gobs_covariance_m2_s4.csv: full corrected15x15 covariance matrices.
results/group_velocity_dispersion_test.csv:66 groups,55source-selected flags.
results/external_validation_results.json: metrics, hashes and limitations.
provenance.json: primary URLs and references.

Citation/acknowledgement: cite Lelli2017, the SPARC master release where
applicable, Brouwer2021 and the KiDS data release, Milgrom2019, and the underlying
Makarov & Karachentsev2011 group catalog. The KiDS release page requests that
users follow the acknowledgement section of Brouwer2021.

Bilek2026's17 early types were reviewed as context. Its new table contains
three-parameter NFW-plus-light profile fits, some with unphysical coefficients,
not raw acceleration pairs. No such profiles are added here as fresh independent
data. That paper motivates independently measured environment for future tests.
