A43 local-response audit, 28 October 2026

Run from the complete A43 supplement, with sibling dm_broad_joint:
  python dm_identifiability/audit.py
  python dm_identifiability/verify.py
Dependencies: Python, NumPy, SciPy. The input CSV is not changed.

This is an exploratory diagnostic on the repeatedly inspected 2,910-row,
279-system/278-validation-group benchmark. It is not a discovery sample, a
new gravitational theory, or a mathematical bound on achievable predictive error.
The finite local protocol was recorded before this diagnostic ran.

ONE FLEXIBLE SHARED RESPONSE

The prediction is
  log10(gpred/a_ref) = log10(g_RAR/a_ref) +
     sum_j W_j exp(-||z-z_j||^2/(2h^2)) d_j /
     sum_j W_j exp(-||z-z_j||^2/(2h^2)),
  z=(log10(gbar/a_ref), log10(r/kpc)),
  d_j=log10(gobs_j/g_RAR_j), a_ref=1e-10 m/s^2.
The RAR scale and all residuals are derived using the training sample only.
W_j assigns equal training weight to each population, then to each validation
group within that population, then divides that weight across its radial rows.
Three inner group folds select h from the fixed set 0.15, 0.30, 0.60, 1.20 dex.
No held-out target or target-population intercept is used for prediction.
The model's complexity comes from the training residuals, despite its compact
written expression. It is not a five-parameter competitor or a physical law.

Validation                           Kernel RMS      Common-RAR RMS
Five whole-system folds              0.234867 dex    0.335997 dex
Entire populations excluded          0.406016 dex    0.367130 dex
All MW-host dwarfs excluded (26)      0.466659 dex    0.623190 dex
All M31-host dwarfs excluded (18)     0.441280 dex    0.475004 dex

The group-held-out kernel median SPARC speed error is 17.306%.
Bandwidth 0.30 dex is selected in all five outer system folds. In population
transfer it selects 0.15 dex for CLASH and 0.30 dex for the others. These are
descriptive comparisons under the existing population/group metric, not
intrinsic-scatter estimates. The pattern supports useful interpolation, with
substantial degradation when a source population is unavailable for training.
It does not establish a hard ceiling for other two-input predictors.

INPUT-ONLY NEAREST-NEIGHBOR CHECK

Each row is matched to the nearest (x,u) point from a different validation group;
a second version also requires a different population. Matching never uses y.
At 0.05 dex tolerance in EACH coordinate, the cross-population version covers
331 rows, 106 validation groups, and 25.49% of the population/group-balanced
sample weight. Within that eligible fraction, 36.91% of the balanced weight has
target accelerations differing by a factor of two or more. The corresponding
RMS pair difference is 0.33480 dex. Coverage and factors are not independent
counts or probabilities: many rows can share one counterpart.

One especially close input pair is:
                     CLASH_MACS1206       UGC09037 (SPARC)
gbar (m/s^2)         2.870780582e-11      2.948341089e-11
radius (kpc)         14.30                14.55
gobs (m/s^2)         4.698941086e-10      4.878764918e-11
Its input separation is 0.01158 dex in gbar and 0.00753 dex in radius; its inferred
target differs by 0.98369 dex, a factor of 9.6314. This does not refute a
deterministic law: the inputs differ, the mass geometries and environments differ,
and the reconstructed targets have different systematic assumptions. In
particular, the cluster target is inferred from lensing with a mass-profile model,
whereas the disk target comes from a rotation curve. The combined table does
not supply a homogeneous joint error likelihood; no significance is calculated.
Some close dwarf pairs with large target differences have very large quoted
dispersion-derived errors. Do not interpret any example as noise-free proof.

All nearest matches and the three prespecified tolerances (0.05/0.10/0.20 dex)
are recorded. The diagnosis motivates independently measured geometry/state and
observable-level inference rather than an unlimited search for sharp curves.

FIELD-EQUATION CONSISTENCY

For gN=-grad(PhiN), the naive vector rule g_alg=nu(|gN|/a0) gN satisfies
  curl(g_alg)=grad(nu) cross gN.
It is generally nonzero in nonspherical geometry; then g_alg cannot itself be
the gradient of one scalar potential. A scalar fit is therefore not automatically
a conservative field theory. This issue and the appropriate field equations
are prior knowledge, not a discovery of this audit.

Controlled local example (dimensionless):
  PhiN=(a*x^2+b*y^2)/2, gN=(-a*x,-b*y), nu=|gN|^(-1/2).
Then curl_z(g_alg)=a*b*(a-b)*x*y /
                  [2*(a*a*x*x+b*b*y*y)^(5/4)].
At a=1, b=2, x=y=1 it is -0.1337480609952844. A central finite difference with
step 1e-5 gives -0.1337480609986397 (absolute discrepancy below 4e-12).
This is an analytic test of local integrability, not a fit to a galaxy.

AQUAL instead specifies div[mu(|grad Phi|/a0) grad Phi]=4*pi*G*rho.
QUMOND specifies
  Laplacian(PhiN)=4*pi*G*rho,
  Laplacian(Phi)=div[nu(|grad PhiN|/a0) grad PhiN].
With appropriate boundary data, QUMOND projects the algebraic field onto its
conservative component. The PDE solutions can depend on complete baryon
geometry and external fields even where two systems share local (gbar,r).
Neither displayed nonrelativistic equation alone supplies all lensing and
cosmological predictions; a relativistic completion and forward observables
are separate requirements.

Primary sources (each source has a 200-word retrieval attribution limit):
1. Bekenstein, J., and Milgrom, M. (1984), Astrophysical Journal 286, 7-14,
   DOI 10.1086/162570. Equations 3-5 and their discussion give the field equation
   and curl correction. https://articles.adsabs.harvard.edu/pdf/1984ApJ...286....7B
   Browser source turn144view1, relevant original pages 8-9.
2. Milgrom, M. (2010), MNRAS 403, 886-895,
   DOI 10.1111/j.1365-2966.2009.16184.x.
   https://arxiv.org/abs/0911.5464 ; browser source turn146view0.
   Publisher search source turn140search1 explicitly describes projection onto
   gradient fields; direct publisher open failed, so cite arXiv or original PDF.

REPARAMETERIZATION CHECK

M_eq=gbar*r^2/G is a spherical-equivalent enclosed mass. It equals actual
enclosed baryonic mass only with the relevant spherical relation. In disks,
exterior material and flattening contribute to local acceleration. It must not
silently be called measured enclosed mass. With r retained, (M_eq,r) contains
exactly the same information as (gbar,r). Similarly gbar*r has units of potential
and equals the local baryonic circular-speed squared, but is not generally the
gravitational potential, which depends on global geometry and boundary data.
Changing coordinates can improve useful functional forms; it adds no new data.

VALIDATION AND CONTENTS

audit.py performs 48 common-RAR calibrations, including inner folds; all optimizer
calls succeeded. There are eight held-out prediction tables and two complete
nearest-match tables, plus results.json. verify.py independently recomputes the
prediction equations, metrics, nearest-match identity, hashes, group exclusions,
all 26/18 dwarf host exclusions and analytic curl check. All ten tables pass.

Implementation correction: the first host test mistakenly searched for 'MW' in
the source environment column. It was corrected to the previous pilot's explicit
mapping: MW=[Milky_Way,G/L], M31=[M31,A/L]. Final results exclude all 26/18
satellites. This routine code correction did not change the candidate set or the
other validation results. No old numerical files or manuscripts were modified.

ADDITIONAL INDEPENDENT AUDITS

verify_stress.py checks all 1,890 frozen X-COP predictions (nine sensitivity
cases, ten families, 21 fixed-radius entries), their metrics and coefficient/input
provenance. The fixed-radius input matches the previous delivered archive byte
for byte in this run; the portable verifier also records its frozen SHA-256.
stress_verification.json records PASS. Read the sibling dm_unified_stress
protocol for this separate, already-inspected transfer sample.

check_candidates.py separately checked the 70 first-stage shared-family tables
and generated illustrative point-Solar probes. Those probes are calculations
against the Newtonian point-source reference, not observational likelihood fits.
The first-stage mixed_power response gives g/gNewton=2.3018 at 30 AU for one
solar mass, showing why its abstract high-acceleration limit is not enough.
Its numerical table and candidate_verification.json retain all ten families.

regularity_audit.txt gives analytic checks of the final proposed stage3 family,
including a proof of nonnegative Newtonian-equivalent density outside an isolated
point source, and the limit of a regular baryonic core. The imposed l>=10 pc
screening bound is not inferred from stellar observations.

Historical correction: dm_one_formula/search_one_formula.py searched 'MW',
which was absent from the catalog environment field, so that historical MW test
was skipped. Its M31 sensitivity held out only 16 explicitly labeled objects.
The original dm_broad_joint tests used catalog labels Milky_Way=25 and M31=16,
excluding the G/L Leo1 and two A/L Andromeda assignments from those host subsets.
The new 26/18 full-host exclusions supersede those incomplete sensitivities.
This host-label issue does not alter the earlier system/population scores.

verify_regular.py independently checks all 16 final-stage prediction tables,
all population/source metrics, minimum-cost selection from the 144 successful
optimizer starts, all 24 descriptive paired-resampling comparisons, exact
host/source exclusions and frozen input/protocol/script hashes. Its direct
formula uses ordinary dimensionful products rather than the fit script's
log-sum-exp implementation. regular_verification.json records PASS.
regular_stellar_probes.csv records the six stated radii for each of the two fits.
The broad-grid monotonicity check spans log10(gbar/a_ref)=-15..15 (601 points)
and log10(r/kpc)=-12..6 (73 points), finding no negative slope beyond 1e-8.
This finite numerical grid is not a global monotonicity or field-ellipticity proof.

For the augmented all-data final-stage fit, isolated one-solar-mass predictions
match Newtonian acceleration at 0.4, 1, 10 and 30 AU to floating precision; the
ratios are 1.000274 at 1,000 AU and 1.78654 at 10,000 AU. The latter retains the
base MOND-like term, and is not a wide-binary observational success claim.
The fitted length l hits the imposed 10 pc lower bound in both all-data fits,
so it has not been independently identified by these data.
