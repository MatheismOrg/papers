A43 gas-coupling bridge: an exploratory compatibility stress test
29 October 2026

MAIN RESULT
The simple constant-weight extension of a cluster-calibrated gas-associated
source term does not transfer to HI-inferred gas-force templates in disks. This result
does not reject the published hot-gas-only model, a full finite-cutoff disk
convolution, or every possible common gravity/matter law.

MODEL AND PRIOR ART
Famaey, Pizzuti & Saltas (2025), Phys. Rev. D 111, 123042, Eq. 17:
https://doi.org/10.1103/dccw-srks ; https://arxiv.org/html/2410.02612v2
rho_extra = eta rho_hotgas exp(-r/Rc). The prior archive and gas-profile
provenance are documented in sibling dm_cluster_repair. The density ansatz
is established prior art and represents additional phenomenological source
mass; this test supplies no microscopic identification.

We keep a0=1.2e-10 m/s2 and compute gpred=nu(gsrc/a0) gsrc. For spherical
clusters, gsrc=gbar+eta G/r^2 integral_0^r exp(-s/Rc) dM_hotgas(s), using
the frozen source quadrature. The declared disk stress-test extension is
gsrc=gstar+(1+eta) ggas, where ggas is the inferred signed gas force. A
constant gas-density multiplier gives exactly this field in the original
gas geometry. It is only a leading uniform-weight approximation to any
finite-Rc extension. The r<=30 kpc selection cannot bound the error of a
full disk convolution, since material outside the evaluation radius also
contributes. This caveat is essential to the interpretation.

SAMPLE AND VALIDATION
The primary disk subset uses baryonic predictors only: positive gas force,
gas force / total baryonic force >=0.5, and radius<=30 kpc. It includes
402 points in 61 galaxies, selected from 2700 frozen SPARC points in 149
galaxies. Gas-force fraction is not gas-mass fraction. The original SPARC
quality selection is inherited. The cluster sample has 86 points in 23 hosts
(65 points/16 CLASH; 21 points/7 X-COP). Joint totals: 488 points, 84 systems.
Gas/stellar accelerations are inferred source components, not direct force
measurements. Fixed distances, inclinations and stellar mass-to-light inputs
are inherited. Cluster targets retain their lensing/hydrostatic reconstruction
and shared source assumptions; no raw-data or nuisance-parameter likelihood
is newly fitted. The targets and candidate family were inspected in previous work, so this
is exploratory compatibility evidence, not untouched discovery validation.

Fivefold CV excludes whole systems. Seeds are 20261027 (primary), 7331 and
48105. Disk and cluster domains have equal total weight; within clusters
CLASH and X-COP have equal weights; within each source hosts are equal and
within each host radii are equal. Domain/source names determine geometry,
data bookkeeping and validation weights, not a fitted switching law.

PRIMARY RESULTS: whole-system CV, RMS log10 acceleration residual (dex)
Model                              Disk       Cluster    Joint balanced
---------------------------------------------------------------------
No extra gas (eta=0)               0.151214   0.440585   0.329379
Frozen CLASH calibration           0.630917   0.107572   0.452564
Joint common eta,Rc                0.272364   0.262600   0.267527
Separate eta / eta,Rc calibrations 0.151214   0.111405   0.132810

Frozen CLASH coefficients were trained using CLASH targets in the earlier
stage; the frozen CLASH score on CLASH is descriptive, not held-out evidence.
The frozen model disk score is a coefficient-transfer test. The fitted rows
exclude every held-out host from coefficient training. The separate row
allows different disk/cluster coefficients and is not a unified prediction.

All-data simple-MOND coefficients: disk eta=1.31e-11 (effectively the
nonnegative boundary); clusters eta=12.01466, Rc=268.799 kpc; joint eta=1.57550,
Rc=10000 kpc (upper boundary). Rc is unidentifiable in the disk approximation
and is never fitted to disks alone. All 15 disk CV fits have eta<=1e-10.
The frozen CLASH eta=11.77438, Rc=288.844 kpc gives disk median acceleration
error 298.61%, or median circular-speed error 99.65%, versus 22.89% and
11.60% for eta=0. These are different statistics and must not be interchanged.

The common compromise lowers the weak joint eta=0 score, but makes both
domains worse than their separately calibrated descriptions: disk 0.27236
versus 0.15121 dex; cluster 0.26260 versus 0.11141 dex. The cluster Rc reaching
the declared upper bound indicates tension within this search, not an
estimated universal length. Exponential-RAR sensitivity yields the same
pattern (disk eta~0; joint disk/cluster RMS 0.27294/0.26592 dex).

Across three CV seeds, simple separate-cluster B spans 0.10993--0.11141 dex;
joint B spans 0.26733--0.26781 dex. These are split sensitivities, not confidence
intervals. Paired 5000-host resamplings of frozen CV predictions give 95%
percentile intervals for joint-minus-separate B: disk [0.09928,0.14033] dex;
clusters [0.11289,0.18647] dex. These exclude refitting, model selection,
measurement errors and shared source systematics. They do not infer a new
physical temperature dependence or gravitational memory.

ENTIRE-SOURCE AND DOMAIN TRANSFERS (simple MOND)
Train disks+X-COP -> CLASH: eta=1.19505, Rc=10000 kpc; B=0.37208 dex.
Train disks+CLASH -> X-COP: eta=2.14121, Rc=10000 kpc; B=0.10726 dex.
Train clusters -> all selected disks: eta=12.01466, Rc=268.799; B=0.63580 dex.
Train disks -> clusters: eta~0 with Rc externally fixed at400 kpc; B=0.44059.
The last length is an explicit external assumption, not measured from disks.

ALGEBRAIC CHECK OF THE DECLARED EXTENSION
On the primary subset, f=ggas/gbar>=0.5; hence the frozen eta=11.774 source
force grows by t=1+eta*f>=6.887. Simple MOND satisfies sqrt(t)<=g(tx)/g(x)<=t
for t>=1. Therefore this extension must increase every predicted acceleration
by at least2.624 relative to eta=0 (circular speed by at least1.620). This
input-only bound is verified in analytic_checks.json and explains the large
transfer penalty without a fitting artifact. It is a comparison to baseline
predictions, not a bound on observed residuals, and applies only to the
constant-weight extension, not the missing finite-cutoff disk convolution.

FULL-DISK DIAGNOSTIC AND MISSING-DATA CAUTION
All 2700 rows are also retained as diagnostics. At eta=10, 8 rows have
nonpositive net source forces; at frozen simple CLASH eta, 9 rows do. These
arise from the signed gas-force contribution in the supplied disk geometry.
Their predicted inward MOND acceleration is left undefined, never clipped
or silently discarded. Metrics labelled valid_conditional describe only
valid rows and are not whole-sample scores when any rows are invalid.

In inputs.csv, ggas_m_s2, gstar_m_s2 and gas_force_fraction are DISK-ONLY
fields. Cluster zeros in ggas_m_s2/gstar_m_s2 are unused computational
placeholders; they are not measurements and never mean absent gas. The
cluster kernel is read from cluster_quadrature.json and overwrites the
placeholder. gbar_m_s2 is the existing total baryonic source for each row.
The separate schema file marks the cluster disk-only fields as inapplicable.
No unmeasured hot-gas profile in a disk has been set to zero. The experiment
instead asks whether one coefficient could multiply two distinct measured
components: HI-based gas in disks and intracluster hot gas.

FILES AND REPRODUCTION
PROTOCOL.txt was frozen before these fits. pre_fit_manifest.json records
its hash and the script hash. run.py and the numerical results are unchanged
since that run; this README, schema, summary and final manifest are post-run
documentation only. run.py requires Python with NumPy and SciPy and sibling
directories dm_halo_benchmark and dm_cluster_repair. Execute python run.py.
No network is used. Script copies its source inputs and logs hashes.

results.json: full coefficients, start outcomes/bounds/train rows, CV
metrics, source/domain transfers, bootstrap intervals and full-disk failures.
predictions.csv: 32380 rows in104 logical (scope,nu,model,split) tables.
full_disk_predictions.csv: 32400 rows, including every invalid prediction.
inputs.csv: primary488-row calculation; full_disk_components.csv and
cluster_inputs.csv preserve original source tables; cluster_quadrature.json
preserves integration nodes and weights. summary.csv is a convenient index.
All104 fit sets and382 optimizer starts converged; failures remain a logged
field rather than an omitted category. An independent audit is separate.

LIMITS
No claim of beating the best published model is made. No full disk cutoff
convolution, external-field calculation, Solar-System likelihood, galaxy
cluster merger, lensing-metric completion, or cosmological evolution is fitted.
What fails is this particular common uniform gas-source extension. A law
that distinguishes physical gas states or matter histories would require
measured state/history inputs and a fresh held-out test. The observed
coefficient tension alone cannot identify such a mechanism.
