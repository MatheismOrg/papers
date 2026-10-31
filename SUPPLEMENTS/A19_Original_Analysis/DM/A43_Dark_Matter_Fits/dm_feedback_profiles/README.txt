A43 feedback-profile continuation, 29 October 2026

REPRODUCE
From the extracted A43 directory, with Python 3, numpy, scipy and pandas:
  OPENBLAS_NUM_THREADS=1 python3 dm_feedback_profiles/run.py
  OPENBLAS_NUM_THREADS=1 python3 dm_feedback_profiles/audit.py
The complete archive is required: fold assignments are read from the frozen
sibling dm_cosmological_halo/results/NFW_Planck_size_cv_*.csv files.
run.py uses four local worker processes. No external data download is required.
PROTOCOL.txt records the local plan fixed before this branch's fits. The overall
project and dataset were already explored; this is not a preregistered discovery.

INPUTS AND METRICS
149 SPARC galaxies, 2700 radii, the previous cuts and fixed baryonic mass-to-light
ratios. input.csv and galaxy_features.csv are byte-identical to their source in
dm_cosmological_halo; audit.json records hashes. The 44-dwarf transfer uses the
archived mass-estimator benchmark from dm_field_bridge. Original citations,
selection and systematic limitations remain in those directories and A43.
Loss: each galaxy has total weight 1/N, split equally among its selected radii;
residual is log10(gpred/gobs). RMS is sqrt(mean galaxy mean residual squared).
The speed statistic is a pooled radius median, averaged across complete repeats.
It is a separate metric and does not have equal galaxy weights.

METHOD
Mass mapping and published c(M) follow A43 Section 4.2. DC14 uses the original
Di Cintio et al. 2014 equations 2, 3 and 6, with an explicitly translated reference
virial boundary and capped out-of-domain behavior. Shape and concentration depend
on X=log10(Mstar/Mdark_vir). The DC14 concentration coefficient is 0.00003, as in
the original PDF, not a garbled HTML rendering of a later catalog.

Density is proportional to x^-gamma (1+x^alpha)^[-(beta-gamma)/alpha].
I(x)=x^(3-gamma)/(3-gamma) * hyp2f1((3-gamma)/alpha,
 (beta-gamma)/alpha,1+(3-gamma)/alpha,-x^alpha).
M(<r)=Mvir*I(min(r,Rvir)/rs)/I(Rvir/rs).
The NFW/DC14 mixture shares this boundary and total mass. Both densities are
zero outside Rvir; 0<=lambda<=1 ensures nonnegative density and finite mass.
This reference boundary is not asserted to be the actual final hydrodynamical
overdensity boundary. Planck reference parameters transfer a WMAP3 calibration.

coreNFW modifies the NFW enclosed mass by
 [tanh(min(r,Rvir)/(1.75Re))/tanh(Rvir/(1.75Re))]^n,
 n=tanh(0.04*14Gyr/tdyn), tdyn=2pi sqrt(rs^3/[G Mdark_NFW(rs)]).
The assumed 14 Gyr is the catalog prescription, not a measured history.
The blend separately trains its RAR/NFW bases then a common convex force weight;
no inner cross-fitting. No coefficient is fitted to test-galaxy kinematics.

Primary references:
Di Cintio et al. 2014: https://arxiv.org/pdf/1404.5959
  https://doi.org/10.1093/mnras/stu729
Read, Agertz & Collins 2016: https://arxiv.org/abs/1508.04143
  https://doi.org/10.1093/mnras/stw713
Li et al. 2020: https://arxiv.org/abs/2001.10538
  https://doi.org/10.3847/1538-4365/ab700e
Dutton & Maccio 2014: see dm_cosmological_halo and A43 references.
These original authors retain credit; no article PDFs are redistributed.

RESULTS AND LIMITS
summary.json contains original-fold repeats, four tails, frozen dwarf transfer,
conditional paired bootstrap intervals, all-data coefficients and domain flags.
Each *_record.json preserves split membership, fits, starts, convergence and
predictions. 20 cases, 120 fitted predictor sets, 500 recorded base-fit attempts;
two attempts exhausted max evaluations, neither selected. No selected fit failed.
All six predictors are retained. DC14 and assumed-duration coreNFW do not improve
the restrictive shared relation; this does not refute their published individual
fits or physical feedback. Calibration extrapolation is material: the all-sample
DC14 fit has one upper-X cap and 71 halos outside its simulated mass range.

Blend point gain versus fitted RAR: 1.63% RMS, only 0.046 percentage points of
median speed error. Conditional 95% difference interval includes zero. The earlier
JavaScript partitions gave 2.04% on the same observed galaxies; the ensemble was
introduced adaptively there. The current run freezes that ensemble algorithm but
is still exploratory. The blend does not transfer well to the 44 dwarfs and has
no finite-DM-mass asymptote or complete dynamics/lensing/cosmology formulation.
Dwarf figures compare frozen disk rules under the isolation and stellar-mass
assumptions; they are not full Jeans, external-field or cosmological likelihoods.

AUDIT
math_checks.json: 169 hypergeometric-vs-quadrature checks and exact nested NFW.
audit.json: independently recomputed metrics; original fold identity; no test
system in training; all fitted DC14 shapes have positive increasing enclosed
mass; 596 fitted mass checks versus direct quadrature; exact outer mass; baseline
NFW predictions agree with archived results to 4.7e-8 dex. Conditional bootstrap
resamples whole galaxies with all repeats together; it does not refit models or
correct for adaptive formula search. Requirements versions are in environment.json.
