A43 constrained shared-response follow-up, 28 October 2026

Run from the complete A43 supplement with NumPy and SciPy installed:
  python dm_unified_limits/run.py
  python dm_unified_limits/verify.py
  python dm_unified_limits/diagnostics.py
The sibling dm_broad_joint/combined_input.csv is required. run.py fits all ten
listed families. verify.py independently evaluates their direct equations and
metrics. diagnostics.py calculates mathematical extrapolations, not observational
likelihoods. No article PDFs or external software are redistributed here.

The protocol was fixed before the new fits, but after earlier A43 residuals and
predictor ranges had been inspected. All six new families and four old comparators
are retained. The candidate comparison is exploratory and subject to selection
bias. The same benchmark has been examined repeatedly across A43 versions.

New families are mass_power, mass_bounded, mass_additive, mixed_power,
index_power and mass_screened. See PROTOCOL.txt and the frozen SPEC in run.py.
Mass means the spherical-equivalent predictor gbar*r^2/G. For disks it is not the
actual enclosed baryonic mass. It is an invertible rearrangement of the same
(gbar,r) inputs, not an independent measurement. No object labels enter a formula.

Results: 230 fits, 552 optimizer calls, 70 prediction tables (including all-data,
three group CV partitions, population exclusions, and two full-host exclusions).
One start for mass_screened without the dSph population hit the evaluation limit;
it was not selected. The selected lower-cost start converged. All diagnostics
including the failed start remain in results/results.json. The MW host is all
26 objects with environment Milky_Way or G/L; M31 is all 18 with M31 or A/L.
The earlier inherited 'MW' alias error was corrected before these final outputs.
An independent direct-formula check passes all 70 tables, numerical residuals,
metrics, input/protocol hashes and host identities; see verification.json.

Best primary internal score is mixed_power:
  s = gbar*r/(a_ref*kpc), M = gbar*r^2/(a_ref*kpc^2), a_ref = 1e-10 SI
  a_eff/a_ref = 0.826377808 + 0.000783719652*M^(-0.95)
               + 0.574763043*s^(0.796356624)
  g = gbar/[1-exp(-sqrt(gbar/a_eff))].
These are all-data illustrative coefficients. Every reported holdout uses its
training-fold coefficients. The inverse exponent reaches its upper bound 0.95;
it has not been identified freely away from the computational limit.

Primary balanced system-held-out RMS: old power sum 0.249630 dex; mixed_power
0.228065 dex. Two alternate seeds give 0.226197 and 0.226292. A descriptive
paired group-Poisson interval for mixed minus old power sum is [-0.03894,-0.00521]
dex. It excludes model selection, parameter refitting and reconstruction errors.
Population-held-out RMS is 0.276962 versus 0.277628; that tiny difference has
interval [-0.00512,+0.00346] dex. Thus the new fit does not improve all transfer
criteria. SPARC median speed error is 15.83%, still worse than the separate
9.2% galaxy-only RAR benchmark. The gain is predominantly in dwarfs:
  dwarf group RMS 0.41500 -> 0.34719 dex;
  MW host RMS 0.44055 -> 0.32334 dex;
  M31 host RMS 0.41108 -> 0.34865 dex.
Other populations are approximately tied and some metrics worsen.

All six new families have g->0 at fixed positive r as gbar->0; relative Newtonian
recovery at high gbar; and a fixed isolated point-mass 1/r far field. That is an
analytic improvement over the old constant-floor extrapolation, NOT sufficient
physical viability. Finite-Solar evaluation reveals a decisive problem with the
best fit: mixed_power predicts +22.085% extra acceleration at 10 AU and +130.18%
at 30 AU for an isolated solar point mass. These are plainly unsuitable as a
Solar-System law. Its asymptotic a_infinity is 2.03e-5 SI for one solar mass,
far above its massive-galaxy limit. The prespecified mass_screened family has
negligible values at these Solar radii, but worse CV 0.258514 and population
transfer 0.423130 dex. No additional screening was silently fitted afterward.

The mass_bounded family has negative dlog(g)/dlog(gbar) at three actual predictor
locations and in a larger test grid, despite bounded asymptotes. No new family's
point-mass inferred-density derivative is negative on the declared finite grid.
These spot checks do not prove all-radius positivity or extended-source physical
consistency. diagnostics.json declares the grids and derivative thresholds.

For fixed point-source mass M_phys the asymptote obeys
  v_flat^4 = G*M_phys*a_infinity(M_phys).
Therefore a mass-dependent response changes the baryonic Tully-Fisher slope.
The mixed model tends to slope log M/log v = 4 for massive sources, but tends
towards 4/(1-alpha)=80 as the inverse term dominates at very small masses.
This is a formal extrapolation, not a BTFR data fit. It further prevents treating
correct radial scaling alone as validation of a universal force law.

The useful result is improved empirical dwarf prediction and a concrete finite-
Solar failure exposed by broader testing. Neither a covariant action, a conserved
extended-source vector field, a lensing metric, nor cosmological perturbations
have been derived. All comparisons are conditional on the fixed heterogeneous
baryonic and dynamical inputs of the existing benchmark.

An additional audit probes a smooth constant-density baryonic center, where
gbar is proportional to r and M_eq is proportional to r^3. The inverse-mass
term then implies g proportional to r^[(1-3*alpha)/2]. With alpha=0.95 this
diverges as r^(-0.925), despite the fixed-radius vanishing-source limit. The
bounded-mass family avoids that divergence. High-acceleration screening alone
does not cure a low-acceleration central limit. This is another reason the best
response cannot be presented as a regular universal baryon-sourced field law.
