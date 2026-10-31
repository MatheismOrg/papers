CURRENT CONTINUATION: 29 October 2026, population scatter and stellar mass

The most useful new point-prediction result is a three-coefficient stellar-mass
NFW mapping. Original-fold disk RMS changes 0.153234 -> 0.153292 dex, while
frozen 44-dwarf RMS changes 0.667569 -> 0.479136 dex (28.2% lower).
Classical dwarfs worsen, ultrafaints improve; three disk feature tails worsen.
The old 25-dwarf quality sensitivity gives 0.575016 -> 0.337799 dex, 41.3% lower.
These are adaptive, conditional benchmark results, not a new independent sample.

The corresponding RAR/stellar-NFW blend gives disk RMS 0.150506, median speed
error 8.78%, and dwarf RMS 0.505205. Its tiny disk change versus the prior blend
is not established; its dwarf RMS is 24.2% lower. See Sections 4.5-4.6 and the
dm_population_scatter and dm_stellar_halo_curve directories. All families,
negative results, protocols, code, measurements and validation records are saved.

A distribution of NFW halo properties also improves a conditional Gaussian log
score versus the specified RAR scatter comparator, with worse central disk RMS.
The approximation, sample selection and covariance limitations are explicit.

Earlier chronological summaries follow. Current manuscript and newest branch
records supersede them; no all-scale or best-in-literature claim is established.

CURRENT REVISION: 29 October 2026, feedback-profile continuation

A43 remains a standalone working benchmark paper. It is not a claim of a new
universal law or established superiority to the best published models.
Sections 4.3-4.4 add published DC14/coreNFW prescriptions, a finite halo response,
a shared RAR/NFW ensemble and frozen transfer to 44 dwarf spheroidals.

On the original 149-galaxy, 2700-radius partitions, mean held-out RMS is
0.153093 dex (refitted RAR) versus 0.150602 (blend): 1.63% lower, conditional
95% difference interval [-0.008801,+0.004033] dex, including zero.
DC14 and assumed-duration coreNFW do not improve the deterministic mass mapping.
No tested shared predictor gives an accurate dwarf transfer. No new cluster,
Solar-System or CMB likelihood was fitted in this extension.

Start with A43_Dark_Matter_Fits.pdf and dm_feedback_profiles/README.txt.
PUBLICATION_ASSESSMENT.txt explains what warrants a separate working paper
and what is still missing for a focused journal submission.

The earlier README entries below document the sequence of exploratory stages;
the current manuscript and branch-specific records supersede earlier summaries.

A43: DARK-MATTER PHENOMENOLOGY ACROSS SCALES
Explicit Fits and Transfer Tests
Exploratory working paper, 27 October 2026 (US Pacific date)

START HERE
A43_Dark_Matter_Fits.pdf / .md: complete paper with sources and two figures.
dm_broad_joint/results/compact_summary.json: main numerical comparison.
dm_broad_joint/results/broad_results.json: coefficients, folds and full metrics.
BROAD_INPUTS_CURATION.txt: included data, portability and reproduction audit.

MAIN FINDINGS
The original 39.46% to 9.22% galaxy speed improvement compares Newtonian baryons
alone with the already published RAR form; it is not an improvement over the
best prior dark-matter model. A logarithmic quadratic gives 9.31%.

The broader benchmark has 2,910 inferred acceleration estimates in 279 physical
systems / 278 conservative validation groups, spanning 27.1 pc to 600 kpc:
149 SPARC disks, 44 dwarf spheroidals, 16 early-type galaxies, 50 MaNGA brightest
cluster galaxies, and 20 CLASH lensing clusters. These are heterogeneous modeled
quantities, not identical raw experiments or independent radial points.

A shared three-coefficient potential/velocity-scale correction improves balanced
group-held-out RMS from 0.3360 to 0.2495 dex relative to a commonly recalibrated
RAR. Entire-population holdout gives 0.3671 versus 0.3158 dex. The unrestricted
six-coefficient polynomial is better internally (0.2355) but worse on population
holdout (0.4883). The compact correction sacrifices SPARC speed accuracy (14.29%
versus the original 9.22%) and leaves large dwarf and cluster residuals.
It is not a universal law, a 9% fit across all scales, or a detected memory effect.

Galaxy-group dispersions, stacked weak lensing, fixed-radius X-COP profiles and
globular controls are assessed separately. Solar-System asymptotics and
cosmological requirements are described, not jointly fitted to raw likelihoods.
The post-search high-g taper and higher-precision dwarf subset are sensitivities,
not independent confirmations. Selection and missing nuisance covariance limit
all regression interpretations. Source authors retain credit and data rights.

CORRECTION TO THE EARLIER CLUSTER CLAIM
The earlier 5.63% R500-radius polynomial result is withdrawn as unification
evidence: the radius is mass-defined and encodes the target acceleration through
g500=250 H(z)^2 R500. All original calculations are retained for reproducibility.
The defining identity, radius-only negative controls, and fixed-physical-radius
stress test expose the dependence. The older Dark_Matter_Curve_Fits report is not
included because its interpretation is superseded by this paper.

REPRODUCE
Python 3.12 recommended. Install requirements.txt (numpy scipy matplotlib astropy).
Keep all sibling analysis folders together. Numerical inputs are frozen locally;
no network is needed for the listed calculations. Runtime depends on hardware.

Original galaxy benchmark and corrected cluster investigation:
  python dm_fit/fit_curves.py
  python dm_fit/test_clusters.py
  python dm_fit/fit_two_scales.py
  python dm_paper_search/search_formulas.py
  python dm_paper_search/analyze_search.py
  python dm_audit_r500.py --fit-root dm_fit --output dm_audit_r500.json
  python dm_fixed_radius/stress_test.py

Rebuild the broad inputs and external fixed-formula tests:
  python dm_broad_dwarfs/build_dwarfs.py
  python dm_broad_globular/build_control.py
  python dm_broad_groups/validate_external.py
  python dm_broad_groups/plot_external.py
  python dm_broad_lensing/prepare_external.py
  python dm_broad_joint/build_input.py

Main broad fits, nested validation and population-transfer tests:
  python dm_broad_joint/fit_broad.py
  python dm_broad_joint/analyze_broad.py
  python dm_broad_joint/taper_sensitivity.py

Post-search dwarf-quality sensitivity:
  python dm_broad_joint/prepare_sensitivity.py
  python dm_broad_joint/fit_broad.py --quick --input dm_broad_joint/combined_high_precision.csv --metadata dm_broad_joint/high_precision_metadata.json --output dm_broad_joint/results_high_precision
The --quick option omits nested selection, retaining group CV and population/host
holdouts as in the supplied sensitivity results. The main run includes selection.

The paper source references local figures in this directory. The complete
portfolio also supplies its PDF/HTML rendering tools in TOOLS and source in
ACADEMIC/SOURCE_MD; the standalone analysis ZIP focuses on numerical reproduction.
No research-paper PDFs, full article texts or installed dependencies are bundled.
Input URLs, author release notices, numerical-table provenance, SHA-256 hashes,
selection exclusions and uncertainties are documented in each subpackage.
The curation audit reran six input builders and checked unchanged numerical CSVs.
The joint prediction verification independently reconstructed all nine main CSV
metrics. These checks confirm implementation consistency, not physical truth.


INTEGRATED FOLLOW-UP: ONE SMOOTH POWER-SUM FORMULA
The paper now includes Section 5.4 and dm_one_formula. Five global coefficients
apply to all systems; no population label enters prediction. Original group-CV
B=0.249630dex vs0.249461 for the previous potential quadratic; population-holdout
B=0.277628 vs0.315811. SPARC speed error worsens to16.06%. There is no uniformly
best fit. Alpha reaches its upper bound1 and leaves a nonzero low-baryon floor;
this must not be used as a vacuum force law. The follow-up reused previously
inspected data and is exploratory; choosing after population holdouts introduces
selection bias. No new memory, cosmological or Solar-System likelihood fit.
Reproduce: python dm_one_formula/search_one_formula.py
All14families/497fits, all42groupCVfiles and full coefficients are included.

PRIOR ART AND CREDIT
Cite the sources listed in A43 and the subpackage READMEs, including McGaugh,
Lelli & Schombert (RAR/SPARC), McConnachie and Wolf (dwarfs), Lelli et al. (ETG),
Tian et al. (CLASH/MaNGA), Ettori/Ghizzardi/Eckert (X-COP), Milgrom and the
Makarov-Karachentsev group observations, Brouwer et al. and KiDS/GAMA weak lensing,
and Baumgardt & Hilker's updated globular catalog. Follow the KiDS release and
Brouwer2021 acknowledgements for publication using their data. No new license
is assigned to third-party numerical catalogs. The paper makes no priority claim
for RAR, potential-dependent gravity, polynomial fitting, or nonlocal memory.

FOLLOW-UP 28 OCTOBER 2026
Section 5.5 and dm_physical_repair add three physical-limit variants, 63 fits,
12 held-out prediction tables, and paired whole-group resampling intervals.
No uniform improvement is claimed. See that package README and PROTOCOL.

ENVIRONMENT/HISTORY PILOT
Section 8.1 and dm_environment add 91-galaxy external-field and 15-dwarf age
tests. No robust new gain; all 2,314 fits are reproducible, with 49 validation
prediction tables and predictor-error sensitivity scores. Read that README.

UNIFIED REGULARITY FOLLOW-UP
Sections 5.6 and 8.2 add nine shared families across three adaptive stages,
a 286-system expanded benchmark, a final seven-coefficient regular response,
geometry/field-law diagnostics, and complete numerical verification. Start with
README.txt in dm_unified_regular and HOST_HOLDOUT_CORRECTION.txt. All unsuccessful
families remain supplied in dm_unified_limits and dm_unified_repair. Data fits
and mathematical limit checks do not establish a universal gravitational theory.
The independent audits are in dm_unified_audit and dm_identifiability.

ESTABLISHED-MODEL BENCHMARK EXTENSION
Sections 1.2, 4.1, 5.7, 6.1 and 8.3 compare scale-appropriate published models,
test contiguous radial prediction in104SPARC galaxies, compare7MOND/repair forms,
and transfer an existing gas-associated density correction across23clusterhosts.
No result establishes superiority over the best published theories. Transient
repair changes mathematical asymptotes without a meaningful empirical gain.
The cluster improvement implements Famaey,Pizzuti,Saltas2025,not a new densitylaw.
Read dm_cluster_repair/README.txt for the CLASH-informed frozen-reference caveat.
The new291logicalpredictiontables,1953fitsets/5843starts are independently checked
in dm_best_audit. All original inputs and earlier results remain preserved.

SHARED-SOURCE BRIDGE EXTENSION, 29 OCTOBER 2026
Sections4.2,6.2,8.4 start from established halos, gas-source corrections and
positive nonlocal weak-field kernels. Shared mass/size halo relations and a
published concentration-mass prior predict unseen SPARC galaxies. The latter
reaches RAR-like accuracy with3fitted coefficients but no robust RAR advantage.
Uniform cold/hot gas coupling fails the stated constant-weight disk extension;
this is not the uncomputed finite-cutoff disk convolution. Fixed positive kernels
have a necessary dwarf mass-budget limit under isolation/half-mass assumptions;
exact Plummer-source integrations and shared refits do not beat dwarf RAR.
No complete MOND/NLG/MOG/CDM theory is ruled out; no all-scale victory is claimed.
The new284logicaltables,426fitsets/1362starts and158920predictionrows are checked
in dm_bridge_audit. Full-disk gas diagnostics add32400rows. Older results and
original inputs remain preserved. See each analysis README and PROTOCOL.
Run dm_bridge_audit/verify.py for independent numerical verification.
