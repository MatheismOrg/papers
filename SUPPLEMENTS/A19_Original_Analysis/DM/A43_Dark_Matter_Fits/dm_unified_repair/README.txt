A43 stage-2 screened/range follow-up, 28 October 2026

Run with Python, NumPy and SciPy:
  python dm_unified_repair/run.py
  python dm_unified_repair/verify.py
  python dm_unified_repair/diagnostics.py
Requires sibling dm_broad_joint/combined_input.csv and
 dm_fixed_radius/fixed_radius_inputs.csv from the A43 supplement.

The protocol was fixed before these fits, AFTER inspection of stage-1 candidate
failures and the earlier X-COP radial stress tests. This is a selected exploratory
repair. It is not independent discovery or a universal gravity theory. Exactly
two new families were fitted; four earlier comparators are retained. No further
family was silently added. Input source/protocol hashes and every optimization
attempt are in results/results.json.

Definitions: a_ref=1e-10 SI; x=log10(gbar/a_ref), u=log10(r/kpc), m=x+2u,
q=x+u. m is the spherical-equivalent mass coordinate, not the actual enclosed
mass in disks. W=1/[1+(gbar/1e-8)^2] uses the fixed older A43 regulator scale.
Both responses use g=gbar/[1-exp(-sqrt(gbar/a_eff))].
  mixed_gate: a_eff/a_ref = 10^A + W*10^(B-alpha*m) + 10^(C+beta*q).
  mixed_gate_range: divide the final term by 1+(r/rc)^2.
The five/six shared coefficients never use population labels or per-object fits.
Alpha,beta are bounded in [0,0.95]; log10(rc/kpc) is bounded in [0,4]. All exact
starts and bounds are in run.py. Coefficients for every validation split are saved.

Original dataset: the exact frozen 2,910-row five-population benchmark.
Augmented dataset: those rows, with CLASH renamed to population clusters, plus
21 X-COP rows at independently specified 500,750,1000 kpc radii. The 27 cluster
hosts (20 CLASH,7 X-COP) have equal weight within one population; the other four
population classes are unchanged. Total 2,931 rows,286 systems,285 conservative
groups. source retains the original archive provenance; CLASH is canonicalized
from Tian2020_CDS_fig2 only for source-group reporting. Predictions preserve both
source and source_detail. X-COP uses the previous primary gbar (gas plus 0.75 of
projected stellar mass) and forward hydrostatic mass target. These heterogeneous,
correlated reconstructions remain conditional inputs, not new raw observations.

All fits use equal population, then equal group, then equal within-group row
weights. Three whole-group partitions, full-population exclusions, and all 26 MW /
all 18 M31 satellite exclusions are included. The augmented benchmark also holds
out all 21 X-COP or all 84 CLASH rows for source transfer. There are 288 fits and 624
optimizer calls; all converge under the declared optimizer criteria. There are 96
prediction tables. Independent direct-equation and metric checks pass all 96;
verification.json also checks group/source identities and exact original-row
preservation apart from the stated population-name change.

Within the ORIGINAL objective, group-held-out / population-held-out RMS (dex):
  old power_sum      0.249630 / 0.277628
  first mixed_power  0.228065 / 0.276962
  mixed_gate         0.228065 / 0.276808
  mixed_gate_range   0.224972 / 0.318877
The fixed gate repairs the finite planetary extrapolation without changing the
astronomical fit appreciably. The range term improves the internal score but
substantially harms prediction when an entire population is absent from training.

Within the AUGMENTED objective, group-held-out / population-held-out RMS (dex):
  RAR                0.329913 / 0.359343
  potential_quad     0.249254 / 0.379614
  old power_sum      0.247490 / 0.274927
  first mixed_power  0.227872 / 0.275862
  mixed_gate         0.227872 / 0.275626
  mixed_gate_range   0.219798 / 0.321211
Original and augmented scores have different sample composition/weights; do not
rank numbers across those two objectives as if measured on identical data.
The range model's alternate group scores are 0.219451 and0.218638 dex. Its group-
resampling difference from old power sum is -0.02769 dex with descriptive 95%
interval [-0.04329,-0.01134]. Population transfer instead worsens by 0.04628 dex,
interval [+0.02675,+0.06468]. These intervals exclude refitting, selection and
input uncertainty. The improved internal score does not establish dominance.

Source-held-out RMS for augmented old power_sum -> mixed_gate_range:
  X-COP: 0.247963 -> 0.061965 dex
  CLASH: 0.248138 -> 0.196943 dex
The range model had already been motivated by inspected radial behavior, so this
is conditional source transfer, not a genuinely blind discovery test. Its all-data
range scale is approximately 406 kpc, with alpha and beta at their 0.95 bounds.
Full per-source and per-population metrics are retained, including SPARC speed.

PHYSICAL FAILURES RETAINED
1. Both gates make the isolated solar point-mass correction negligible at 0.4, 1,
   10,30 AU in the recorded spot checks. This is not a full ephemeris/Cassini test.
2. They still give enormous isolated 1-solar-mass enhancements at 1000/10000 AU:
   roughly 50 and580 times the Newtonian acceleration as EXCESS, not total ratio.
   Galactic external field, vector binary dynamics and observational likelihoods
   were not modeled. This is a failure of the naive isolated extrapolation.
3. A smooth constant-density baryonic core has gbar proportional to r and
   M_eq proportional to r^3. The inverse-mass term with alpha=0.95 yields
   g proportional to r^(-0.925), divergent at a regular center. The fixed-high-g
   gate cannot repair a low-g center. diagnostics.json includes this separate
   limiting path, which a fixed-radius vanishing-source test misses.
4. The squared radial range introduces negative point-source inferred-density
   shells on the declared grid: 26 of 4,991 probes for the augmented fit, with
   min dlog(M_eff)/dlog(r)=-0.311 near 1e15 solar masses and 1277 kpc. A positive
   scalar acceleration is insufficient to ensure a nonnegative effective halo.
   The ungated and gate-only mixed responses have no negative values on that
   same finite grid, but still have failures 2 and 3. This is not an all-domain proof.
   Signed effective (phantom) density can occur in modified-gravity models; this
   finding excludes a nonnegative extra-mass interpretation at these probes, not
   every possible modified-gravity completion by itself.

The candidate improves empirical fits on the broader selected sample and exposes
specific next requirements. A regular center, physically controlled source/scale
response, conserved extended-source field equations, metric lensing and cosmology
remain unresolved. A separate later repair must retain this stage's failures and
be labeled post-selection; it must not retroactively overwrite these results.
