Independent audit of the A43 halo/source bridge round
29 October 2026

Run: python /path/to/dm_bridge_audit/verify.py
Dependencies: Python 3, NumPy, SciPy. Keep this directory adjacent to the
dm_shared_halo, dm_cosmological_halo, dm_gas_bridge, dm_field_bridge and
underlying dm_fit, dm_unified_regular, dm_halo_benchmark and dm_cluster_repair
directories. Relative paths are resolved from the scripts, not the shell cwd.

verify.py reruns four independent scripts, writes per-branch logs and reports,
and creates verification.json. It does not import or execute analysis fitting
code, refit models, alter original outputs or choose new model families.
Every script independently recomputes predictions and scores from frozen
inputs and recorded coefficients. Source and script hashes are checked.

The scope is 284 logical forecast tables, 158,920 forecast rows, 426 fit sets,
1,362 optimizer attempts, and 94 descriptive paired bootstrap comparisons.
There are separately 32,400 full-disk extrapolation diagnostic rows and 49
exact comparator-table copies; these are not double-counted as new forecasts.
Repeated rows reflect different models, splits or assumptions. They are not
158,920 independent astronomical observations.

What is independently checked

1. Shared photometric halo maps: raw 175-row SPARC photometric-catalog parse;
   149-galaxy joins; no Vflat predictor; 56 forecasts; equal-galaxy loss;
   integrated NFW/Burkert profiles; all starts, minimum-cost choice, bounds,
   whole-galaxy folds, predictor-only tails and 42 bootstrap comparisons.
   The original parser failed before fitting; a retained correction manifest
   documents its replacement with the actual 19-whitespace-field layout.

2. Cosmological structural prior: exact input identity with the audited
   photometry; mass convention and cosmic baryon fraction; reference r200
   independently constructed from critical density; published c(M) h factor;
   dark-mass factor 1-fb; integrated/truncated NFW; 16 forecasts, 40 fit sets;
   all structural diagnostics, 49 comparator hashes and 42 bootstrap contrasts.
   r200 is a DMO-equivalent reference radius. Once actual baryons are added,
   it need not be the assembled system's true overdensity-200 radius.

3. Common cold/hot gas coefficient: independent component selection and sums;
   516 adaptive cluster-density or mass-derivative integrals at six Rc values;
   all 104 forecast tables and 382 starts; source/domain/host exclusions;
   10 paired host bootstrap contrasts; nonpositive source forces preserved
   in all 32,400 full-disk diagnostic rows. Missing hot gas is not represented
   as measured zero: cluster gas enters via density quadrature, while unused
   disk-only component columns contain documented placeholders.

4. Positive spatial kernels: analytical normalization and finite mass budget;
   exact, stable translated-Plummer-ball mass independently checked by shell
   solid angles; 1,602 infinite-support scalar integrals for the 401-point
   basis and 400 interval midpoints per kernel; every forecast, source join,
   optimizer objective and host exclusion; general and stronger Plummer
   bounds; selected fit saturation and analytic infinite-mu optimum.

Interpretation limits

PASS means the stated computation and validation are reproducible. It does
not establish physical validity, novelty, generalization to unseen surveys or
an all-scale theory. This entire round reuses previously examined targets.
No source profile, assembly history, cosmological likelihood, relativistic
lensing solution, Solar external-field test or causal temporal memory law is
supplied by an arithmetic audit.

The gas extension is a constant multiplier of measured cold-gas disk force,
not an exact finite-Rc disk field. An evaluation-radius cutoff cannot bound
the error caused by unmodeled exterior disk sources. Its failure does not
exclude the original hot-gas-only cluster proposal.

The field bounds use isolated positive spherical sources and a true baryonic
half-mass radius. Observations use an equilibrium half-light estimator and
assumed stellar masses. Tides, binaries, anisotropy, membership and host mass
can change the inference. The stronger 1+alpha bound also requires radially
nonincreasing baryon density. No confidence-level theory exclusion follows.

The field basis has a maximum independent relative discrepancy 1.4 parts per
million at extremely small z for MOG; the all-data forecast error against
independent direct integrals is below 6.5e-8 dex. This is numerically negligible
compared with model residuals, and is reported rather than hidden.

Manuscript review: halo_section.txt, gas_section.txt and field_section.txt in
dm_bridge_integration were read for mismatch/overclaim. Their numeric values
and scope match this audit. The reference-radius clarification was sent to
the integrating author. No manuscript was changed by the auditing scripts.

ANALYTIC_PROOFS.txt gives the physical normalization, inequalities and the
independent closed translated-ball formula used by the audit.
