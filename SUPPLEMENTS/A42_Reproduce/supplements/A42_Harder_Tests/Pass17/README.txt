A42 PASS 17: LOCAL SOURCE AND GRAVITATIONAL FEEDBACK
Completed 2026-10-29 UTC.

READ FIRST
A42_Shortest_Rules_Search.pdf preserves the previous 46 pages and adds pp47-50.
This pass is a new 3+1 homogeneous semiclassical benchmark. It is not the same
1+1 model from pass16 enlarged without additional assumptions. The standard
GR/scalar Hamiltonian is an input. Geometry is one classical global scale
variable. Matter is a free Gaussian quantum scalar field on a fixed cubic
spatial lattice. These tests neither identify a ToE nor resolve singularities.

RESULTS
- Direct local T/2-V-T/2 gate response and geometry impulse agree independently.
- 30 coupled runs: 3 initial/parameter cases, 2 clocks, 5 time steps.
- 10 deliberately omitted-feedback runs: 2 clocks, 5 time steps.
- Coupled last-halving error ratios: 3.996864 to 4.008415 across a,H,rho,p,
  and Friedmann constraint. No constraint projection after initialization.
- Finest relative scale-factor clock difference <=1.18363e-6.
- Finest maximum relative Friedmann residual <=2.51378e-5.
- All-run Gaussian covariance purity error <=1.24345e-14.
- Omitting feedback preserves purity but gives analytically incompatible
  clock limits, 7.587759% scale-factor difference at proper time2. The proper
  clock constraint residual stays near84.8% as the numerical step shrinks.
- Bare-vacuum cutoff study shows rho approximately epsilon^-4, p/rho->1/3;
  a constant Lambda canceling initial rho leaves p+rho. This is a regulator
  diagnostic, not a no-go result for renormalization.

MODEL / INITIAL STATES
hbar=c=M_P=1; 4^3 periodic sites, comoving spacing1, box length4, volume64,
mass.4. The 64 real normal modes are grouped into equal lattice eigenvalues,
with degeneracies retained. Cases: vacuum(Lambda0), squeezed(r.4,angle.6,
Lambda0), vacuum with vacuum energy density Lambda.05. Every initial a=1;
H0 is fixed by Friedmann constraint. The same initial Gaussian covariance is
held fixed as h changes. This is time refinement at fixed spatial cutoff,
not pass16's joint UV/time limit. h=.08,.04,.02,.01,.005. Clocks N=a^s,
s=0,1; coordinate increments are converted to proper time. Runs cover t=2;
stored trajectories can contain one endpoint beyond2. Error comparisons use
nodes up to t=2; reported final values/clock comparisons interpolate to common
proper times. The auxiliary clock uses trapezoidal endpoint half-gates.

FILES
code/backreaction.py: model, exact split steps, independent a,H ODE reference,
  40 runs, histories and metrics.
code/audit_backreaction.py and gravity_formula_audit.py: independently built
  full real-space 64-mode canonical reference; pressures, lapse derivatives,
  reversibility, analytic no-feedback and empty-field deSitter checks.
code/audit_gate.py: matrix derivatives, generating action, local support,
  shadow averages, impulse equality, finite-step lapse-response counterexample.
code/cutoff_test.py: independent fixed-volume bare-vacuum sums.
code/summarize.py: aggregate checks and scientific figure.
code/build_update.py: appended four PDF pages, source-page preservation check.
code/*.json: authoritative numerical results and audit reports.
code/*.npz: reference and test trajectories, including Gaussian covariances.
code/gate_derivation.txt and gravity_review.txt: detailed derivations/provenance.
source/A42_PrePass17.pdf: unchanged 46-page input.
previous_pass16.zip: prior cumulative tests/code, including earlier passes.
SHA256SUMS.txt: hashes of this package's other files.

REPRODUCE from this package root
Python with numpy, scipy, matplotlib, pymupdf, reportlab and fonttools.
  python3 code/backreaction.py
  python3 code/audit_gate.py
  python3 code/audit_backreaction.py
  python3 code/cutoff_test.py
  python3 code/summarize.py
  python3 code/build_update.py --source source/A42_PrePass17.pdf --output A42_Shortest_Rules_Search.pdf
PDF authoring needs the font paths listed in build_update.py. Adjust those to
installed Liberation Sans and DejaVu Sans fonts on other systems. The numerical
runs do not need PDF/font packages. No external network data or secrets needed.

INTERPRETATION
Convergence is to a finite-regulator semiclassical reference. Four-dimensional
inhomogeneous gravity, a regulator-independent renormalized source, quantum
geometry, interactions, finite-information state encoding and uniqueness of
an underlying Turing-machine rule are unestablished. The natural next test is
an action-consistent vacuum subtraction whose geometry variation preserves
the same clock and constraint checks as the spatial cutoff changes.
