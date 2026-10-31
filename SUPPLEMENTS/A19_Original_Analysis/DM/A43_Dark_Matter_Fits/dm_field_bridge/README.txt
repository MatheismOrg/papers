A43 physical-source bridge: published NLG and constant-parameter MOG
29 October 2026

RESULT

This test does not improve on the exponential RAR benchmark. It identifies a
specific limitation of transferring two published positive linear kernels with
one universal coefficient pair into the existing 44-dwarf sample. It is not a
complete test of either covariant theory, a new gravity theory, or an all-scale
fit. The source rule and both parameter pairs are from the literature.

At nominal stellar M/L_V=2, log-acceleration RMS errors in dex are:

entry                 all-data   fivefold*   withheld MW   withheld M31
NLG published          1.861705    1.861705     2.104039       1.441451
NLG shared refit       0.700081    0.728845     0.933839       0.763417
MOG published          1.905164    1.905164     2.139397       1.503737
MOG shared refit       0.700081    0.728845     0.933839       0.763417
RAR fixed a0           0.650849    0.650849     0.773815       0.412987
RAR shared refit       0.528790    0.548437     0.696312       0.600618

*Fivefold seed 20261027. Seeds 7331 and 48105 are also retained: fitted NLG
0.729518 / 0.707438, fitted MOG 0.729518 / 0.707438, fitted RAR
0.544917 / 0.535199. Published/fixed entries are untrained predictions;
their fold tables are descriptive. Withheld MW means train on all 18 M31
objects and test on all 26 MW objects; withheld M31 reverses this split.
No sample-error or covariance likelihood was fitted. These central-estimate
descriptive scores are not p values or evidence ratios. RAR here is the
specific isolated algebraic relation, not a full MOND field/host/Jeans fit.
The primary-fold median absolute acceleration errors are 98.12% / 85.62%
for published/refitted NLG, 98.27% / 85.62% for MOG, and 70.27% / 78.68%
for fixed/refitted RAR. The fitted RAR improves the chosen RMS loss while
worsening this median percentage metric; no improvement on every metric is
claimed. The conventional fractional error is asymmetric for under/overprediction.

The universal positive-kernel mass budget gives D=gobs/gbar <= 1+2 alpha at
the baryonic half-mass radius, without assuming a particular spherical radial
profile. Violations at nominal M/L=2 are 31/44 for published NLG alpha=10.94
and 33/44 for published MOG alpha=8.89. At M/L=3, 29/44 remain for both.
These are conditional central-estimate inconsistencies, not statistical theory
exclusions. Changing the kernel range mu cannot evade this total-mass bound.
Plummer imposes the stronger bound D<=1+alpha; all counts are archived.

All 108 selected shared-kernel fit sets (2 kernels x 3 M/L x 18 splits) reach
the predeclared upper bound mu=10^5/kpc. For nominal M/L=2 the all-data alpha
is 58.413217. This approaches a local rescaling g=(1+alpha)gbar; the exact
constant-rescaling optimum is alpha=58.413216. Thus there is no measured,
finite preferred memory length here. Extending the search bound is unnecessary
to establish that this limiting fit underperforms RAR on these scores. Raising
alpha to fit dwarfs does not establish a fit to disks, clusters, Solar tests,
lensing or cosmology. The RAR fitted a0 is 7.07944e-10 m/s^2 in this sample;
its improvement likewise does not demonstrate universal coefficient transfer.

METHOD AND INPUT

The source equation is Laplacian(Phi)=4 pi G [rho_b + q*rho_b]. Both q are
positive, radial and integrate to alpha. Exact formulas, proofs and limits are
in ANALYTIC_CHECKS.txt. Frozen values are NLG (10.94, .059/kpc), MOG
(8.89, .042/kpc). Each refit has two shared coefficients, never one pair per
object. No population flag or measured discrepancy enters prediction.

Input comprises the exact 44 dSph rows in the previously curated broad table,
with one inferred acceleration at half-light per dwarf. The source metadata
and selection are in ../dm_broad_dwarfs/build_dwarfs.py and its results. This
test uses the same circularized projected radius and deprojection rhalf=4Re/3,
gobs=3 sigma_los^2/rhalf and gbar=G Mstar/(2 rhalf^2). Nominal Mstar=2 L_V.
For this conditional test only, assume spherical Plummer baryons with
a=rhalf/sqrt(1/(2^(2/3)-1)), so the adopted radius encloses exactly half the
source. Plummer a is not silently equated to the three-dimensional half-mass
radius. No new mass estimate or full observed source profile is invented.

The convolution integrates the entire Plummer distribution, including source
mass outside rhalf. A point mass made from enclosed baryons is not used.
Only the observed tracer/Wolf acceleration proxy is compared, not a fit to
individual stellar velocities or the full line-of-sight dispersion profile.
All galaxies are treated in isolation. Host mass, tides, anisotropy and binary
contamination can affect those proxies; the actual objects are satellites.
The cleanest interpretation is an internal consistency stress test under the
declared source and equilibrium assumptions. M/L=1,2,3 is a sensitivity grid,
not a posterior or an exhaustive range of possible stellar populations.

WHY THIS IS USEFUL, AND ITS LIMITS

Starting from a field/source law makes the failure informative: a universal
positive linear kernel has a finite mass budget. The same law also scales the
field linearly with total baryonic mass at fixed source shape and size, while
its dimensionless enhancement at a homologous radius depends only on size.
To repair this class, one must change an assumption (source masses, host and
dynamical modeling, universal coefficients, positivity/linearity, or source
response), then test that change in independent systems. An arbitrary new
coefficient is not a demonstrated cure. In particular, a stationary spatial
convolution does not measure time-dependent gravitational memory, and this
test offers no evidence that spacetime discreteness is needed.

These dwarf tensions and source models have precedents. Haghi & Amiri 2016
reported fixed-MOG parameter problems with eight dwarf Jeans profiles.
De Martino et al. 2025 used the same NLG q0 and Plummer source for eight
classical dwarfs, including an 18-parameter joint analysis (two shared kernel
parameters plus two intrinsic parameters per galaxy). They reported broad
compatibility and some parameter tension. Our fixed-M/L half-mass sample,
loss and holdouts differ and do not supersede their full-profile analysis.
The useful addition is a reproducible broader conditional transfer audit and
an explicit necessary bound, not a claim of first discovery of dwarf tension.

The astronomical q0 limit omits NLG's short-distance scale. Solar-system
analyses retain that additional parameter; our q0 run is not a Solar test.
The MOG comparison is specifically the constant-parameter weak-field model,
not every source-mass-dependent or scale-dependent version. No disk geometry,
complete cluster gas+stellar sources, gravitational lensing or cosmological
perturbations were solved here. Missing source data were never replaced by
zero gas. The present numerical work is 44 dwarfs; the cross-scale aspect is
transfer of previously published galaxy/cluster coefficients.

REPRODUCTION AND AUDIT

Run run.py with Python, NumPy and SciPy from the directory layout of this
portfolio (it reads ../dm_unified_regular/augmented_input.csv). Then run
diagnostics.py (also reads the earlier dwarf benchmark for mass/subset checks).
PROTOCOL.txt was written before the fit; run.py retains all starts, coefficient
bounds and optimizer outcomes. results.json records 162 fitted model/split
sets, 540 optimizer attempts, zero optimizer failures and 108 row-level tables.
Hashes connect protocol, input, script, source and kernel table.

Independent analytic shifted-Plummer integration agrees with the production
64-point quadrature to 2.0e-15 relative over 361 displacements. Independent
log-distance integration at 133 points per kernel gives maximum interpolation
relative errors 6.7e-6 (NLG), 1.1e-5 (MOG). Kernel probability normalizations
agree with unity to floating precision, H is nonnegative and <=1/2, and source
mass reconstruction agrees with prior curation to 2.3e-16 relative. Production
quadrature emitted convergence warnings at extreme grid points; these were
retained in run.log rather than suppressed. Independent checks verify the
reported precision. The independent team audit is provided separately.

Files: input.csv, basis.csv, mass_budget_bounds.csv, results.json,
results/*.csv, run.py, PROTOCOL.txt, diagnostics.py, diagnostics.json,
quadrature_checks.csv, ANALYTIC_CHECKS.txt, primary_sources.txt.
No third-party article PDFs are distributed by this directory.
