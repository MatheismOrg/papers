EXTERNAL LENSING AND BRIGHTEST-CLUSTER-GALAXY DATA

Data obtained 27 October 2026 Pacific time (server timestamps 28 October UTC).

1. CLASH: Tian, Umetsu, Ko, Donahue, Chiu (2020), ApJ896,70.
   DOI10.3847/1538-4357/ab8e3d; arXiv2001.08340.
   Official CDS/VizieR catalogue J/ApJ/896/70, DOI10.26093/cds/vizier.18960070.
   84 numerical acceleration pairs downloaded as TSV from the official VizieR
   South Africa mirror, not digitized from an image. The other mirrors returned
   server errors during retrieval. Original URLs are in external_results.json
   and the cached TSV headers.

   20 central BCG locations (14.3-30.3 kpc) and64 outer cluster locations at
   fixed physical100/200/400/600 kpc. These share20 host clusters; all locations
   from a host must remain together in train/test splits. Cluster targets are
   spherical NFW mass-profile estimates from combined strong/weak lensing and
   magnification posteriors. Baryonic mass uses X-ray gas and stellar estimates.
   They depend on cosmology and on the gravitational framework used for lensing.
   They are independent of X-COP HSE reconstruction, not theory-independent raw
   accelerations. Radial bins share mass-profile parameters. The published
   analysis explicitly ignored cross-bin covariance for simplicity; the supplied
   table contains marginal log errors only, not a covariance matrix.

2. MaNGA BCGs: Tian, Ko, Li, McGaugh, Poblete (2024), A&A684,A180.
   DOI10.1051/0004-6361/202347868; arXiv2402.12016.
   50 numerical pairs transcribed programmatically from the actual Table1 text
   on PDFpage8 using pdftotext -layout. The table was visually checked against
   the PDF. No plotted points were digitized. A cached50-line numerical table
   is supplied, with the original printed precision and errors.

   Physical radius = Re[kpc] times rlast/Re, both published in Table1; it is
   based on observed photometry/kinematic extent, not inferred dark mass. The
   observed acceleration is an isotropic(beta=0), spherical Jeans/Abel inversion
   estimate from stellar kinematics; baryons use photometric Sersic masses.

All loggbar/loggobs values and errors in the standardized files are base10
logarithms of accelerations in m/s^2; errors are symmetric1-sigma dex as printed.
CLASH paper Equation8 instead uses natural logarithms: its intercept must be
converted before applying it to base10 data. No covariance matrix is invented.

Standardized files:
  external_standardized.csv: all134 rows.
  clash_standardized.csv:84 rows in20 hosts.
  manga_bcg_standardized.csv:50 rows in50 observed BCGs.
Fields include population, system, parent_system, radius_kpc, gbar_m_s2,
gobs_m_s2, loggbar_err_dex, loggobs_err_dex, original log values, source/method,
redshift, and environment. CLASH rows include full cluster names and positions.
The MaNGA parent_system labels denote observed BCGs, not verified unique host
clusters; host-cluster IDs are absent from the table.

Crossmatching:
  MaNGA DR17 positions were retrieved by public SDSS SQL query to mangaDrpAll.
  No MaNGA object overlaps an X-COP cluster center: the nearest is71.27arcmin
  away with a large redshift difference. CLASH names are disjoint from the seven
  X-COP targets. MaNGA redshifts .027-.148 are separate from CLASH .187-.686.
  The closest MaNGA pair,8989-12703 and8991-6102, is2.694arcmin apart with
  delta-z .002943. They are distinct objects but may share an environment;
  conservatively keep them together in cross-validation or report sensitivity.

Frozen-model benchmarks in external_results.json use pre-existing SPARC and
X-COP fits. No external target is used to adjust those formula coefficients.
The published_CLASH_sqrt comparison is sqrt(2.02e-9*gbar); its CLASH score is
in-sample to Tian2020, and Tian2024 already compared it with MaNGA. It must not
be described as a new or independently discovered relation.

Reproduce standardized files and frozen benchmarks:
  python prepare_external.py
Requires numpy and astropy. The script expects sibling dm_fit and dm_fixed_radius
folders for previous fits and X-COP input headers. The numerical inputs needed
here are clash_idia.tsv, clash_table1_idia.tsv, manga_bcg_table1_extracted.txt,
and manga_positions_sdss.csv. The original PDFs, arXiv tar archive,75MB DRPall
FITS, temporary PDF render and full-paper text files are research intermediates
and should not be included in the user-facing reproduction package.
