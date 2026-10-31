A43 single-family regularity repair, 28 October 2026

This is a third, explicitly adaptive exploratory stage. Read PROTOCOL.txt for
the one seven-coefficient mixed_regular definition, frozen bounds and starts.
It follows the dm_unified_limits and dm_unified_repair searches; their failures
are retained in those sibling directories. No gravitational law is established.

Reproduce from the complete supplement:
  python dm_unified_regular/run.py
  python dm_unified_regular/diagnostics.py
Dependencies: Python, NumPy, SciPy. The sibling dm_unified_repair inputs and
results are needed for exact comparator predictions. Inputs are copied unchanged;
matched comparator use requires identical input hashes, row/group identity and
fold assignments. Results are deterministic conditional on library versions.
No source-specific or population-specific coefficient occurs in the new model.

The fit performs 48 independent training-set fits with three fixed starts each:
144 optimization attempts, zero reported optimizer failures in this run. It
supplies 16 prediction tables, including all-data, three whole-group splits,
full-population, MW/M31-host, and CLASH/X-COP-source exclusions. Complete selected
and attempted coefficients, costs and bounds are in results/results.json.
Five thousand paired group-Poisson draws compare fixed predictions against the
six prior comparators on each dataset. Their intervals omit refitting, adaptive
model selection, reconstruction and measurement uncertainty.

The original dataset contains 2,910 rows from 279 systems / 278 groups. Its
new-model whole-group RMS values are 0.22644, 0.22438 and 0.22482 dex, and
population-held-out RMS is 0.28744 dex. The augmented dataset adds 21 X-COP
fixed-radius rows and gives equal cluster-host weight jointly with CLASH:
2,931 rows, 286 systems / 285 groups. Its new RMS values are 0.22198, 0.22168
and 0.22088 dex, with population-held-out RMS 0.28921 dex. Both input sets
and all test populations were previously inspected. These are not new-survey
confirmation scores, and different objective/sample compositions must not be
ranked directly against one another.

The augmented whole-group fit improves over the prior power-sum comparator,
but population transfer remains worse. It also improves over mixed_gate's
whole-group score while worsening its population-transfer score. It repairs
mixed_gate_range's point-source density failure while giving slightly worse
whole-group accuracy. There is no uniform winner. The augmented SPARC median
speed error is 15.26%, still worse than A43's galaxy-only RAR result near9.2%.

Physical checks and limits are explicit in ANALYTIC_CHECKS.txt. The radius^4
gate cures the smooth baryonic-core divergence and the power1 outer gate yields
nonnegative inferred extra density for point sources. These claims are analytic,
with numerical spot checks in diagnostics.json and smooth_core.csv. The point
source 1/r far field remains, with a mass-dependent Tully-Fisher normalization.
Positive density for an arbitrary extended source is not proved.

The screening length is restricted to at least10 pc by an imposed prior; both
all-data fits reach that boundary. That length is not a measured discovery. At
30 AU, the illustrated isolated-Solar correction is tiny, but at10,000 AU a
roughly79% acceleration enhancement remains from the ordinary low-acceleration
term. Zeros in solar_stress.csv mean the exponential correction underflowed
double precision, not exact equality of the force laws. These are isolated force
predictions, not fits to ephemerides, wide-binary populations or Cassini's Galactic
external-field quadrupole. No vector field, metric lensing or cosmology is solved.

pre_fit_manifest.json records protocol/script hashes before fitting; results
records them again. manifest.json seals the final numerical files. Independent
verification, when supplied separately, is distinct from these diagnostics.
