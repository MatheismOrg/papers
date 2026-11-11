A42 Pass36 — Four-dimensional quantum matter and self-consistent geometry
Research working package, 31 October 2026 UTC

RESULT AND SCOPE
A published nonsingular gravity model is coupled to the published four-dimensional
conformal-vacuum stress in closed homogeneous/isotropic FLRW. Quantum matter
changes the geometry self-consistently, including through its angular pressure.
This replaces the inadequate radial-source promotion in this restricted sector.
The gravity action and anomaly tensor are prior work. No novelty priority is claimed.

The combined background equation is
  Phi/sqrt(1-alpha^2 Phi^2) - eta Phi^2 = C_d/a^3 + C_r/a^4,
  Phi=H^2+1/a^2, eta=N_eff*hbar*G/(360*pi).

An analytic condition lambda=eta/alpha < 8/(3*sqrt(3)) ensures positive
algebraic response for every nonnegative density. A second analytic inequality
bounds both independent FLRW curvature combinations in the dust sector.
For the specified free SM-count proxy at alpha=hbar*G it gives K<=24/alpha^2.
The inequality is conditional on the stated model, state and coupling assumptions.

COMPUTED RESULTS
Passing checks: 16/16 (implementation and structural checks).
Evolution entries: 20; distinct parameter sets: 18.
Free-field effective count: 995.5.
Free SM-count lambda: 0.880215254711; response threshold: 1.53960071784.
Maximum sampled quantum density fraction in main proxy runs: 0.440107616.
Maximum normalized constraint error across all entries: 2.03631e-10.
Maximum bounce-radius discrepancy: 2.57572e-14.
Independent finite-difference conservation residual: 1.73834e-13.
Angular-equation residual: 1.16249e-15.
Covariant quantum-tensor component error: 1.14492e-16.

Fixed dust constant c=3.22370979547:
  Classical minimum radius: 1.02712279278; quantum: 1.01605252375 (sqrt(alpha) units).
  Classical maximum sampled K*alpha^2: 18.5300581532; quantum: 20.8949172953.
Both bounce; quantum feedback does not decrease curvature in this comparison.
These trajectories share the dust constant, not their complete initial data.

FAILURES RETAINED
At lambda=1.5, maximum sampled K*alpha^2 is 618.229647114.
Thus the classical ceiling 24 cannot be extended to every subcritical coupling.
At lambda=1.6, c=2.0, response vanishes at y=0.414831004049.
The scaled proper time from half that y is 0.290603806488.
The nonzero source divided by vanishing response makes acceleration diverge.
This is an effective-equation failure, not a rounding artifact.

ASSUMPTIONS AND UNRESOLVED WORK
The total renormalized Box R coefficient is set to zero. With dynamical gravity
this fixes an independent physical curvature coupling. Restoring it changes
the differential order and can change stability. That test remains open.
Matter is free, massless and conformal. The SM field-count proxy is not an
interacting Standard Model: masses, running, symmetry breaking and scattering
are absent, and conformal Higgs coupling is assumed. There is vacuum polarization,
not conformal particle production. alpha is specified, not derived.
Near-Planck-curvature examples are formal semiclassical solutions; control over
quantum-gravity corrections has not been demonstrated. Spatial perturbations,
inhomogeneous black-hole centers and evaporation are not solved here.
The implicit-midpoint step is a standard integrator, not a fundamental discrete
spacetime law, exact finite-step conservation proof, or shortest ToE program.

REPRODUCTION
Recorded environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0.
In a suitable Python environment:
  python3 -m pip install -r requirements.txt
  python3 run_all.py
run.py is the authoritative executable experiment specification. It derives the
SM proxy from field counts and writes inputs.json, results.json and trajectories.npz.
inputs.json records that specification; editing it alone does not change a run.
run_all.py captures run.log, regenerates the plot and regenerates this summary.
All observables are calculated, not supplied as expected answers. Assertions
check tolerances/identities; the physical failure control is expected to fail
regularity and is retained explicitly. No downloaded numerical data are required.
Exact floating-point results can vary slightly across environments.

FILES
DERIVATION.txt: equations, proofs, assumptions and interpretation.
PAPER_INSERT.txt: proposed A42 addition; the main manuscript was not overwritten.
PROTOCOL.txt and AUDIT_NOTES.txt: choices before/after runs and limitations.
SOURCES.txt: primary literature and attribution.
model.py: equations, root searches and discrete update.
run.py, inputs.json, results.json, run.log: experiment and numerical record.
trajectories.npz: sampled trajectories, named by case and variable.
plot.py, quantum_backreaction.png: scientific figure and generator.
report.py, run_all.py, requirements.txt: reproduction helpers.
SHA256SUMS: integrity hashes of the distributed files (before rerunning).

CHECK REGISTER
  response_threshold: PASS
  bounce_evolution: PASS
  constraint: PASS
  conservation: PASS
  anomaly: PASS
  sufficient_curvature_bound: PASS
  excessive_coupling_rejected: PASS
  low_curvature_limit: PASS
  discrete_bounce: PASS
  discrete_convergence: PASS
  de_sitter_anomaly: PASS
  ward_independent_finite_difference: PASS
  full_angular_equation: PASS
  analytic_SM_uniform_bound: PASS
  covariant_quantum_tensor: PASS
  conformal_casimir_vacuum: PASS
