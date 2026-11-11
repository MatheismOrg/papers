A42 pass 11 - interacting gauge benchmark and separate error budgets

Main paper: 31 pages. Section 11 (pages 30-31) records the new benchmark.
Earlier 29 pages are preserved. Earlier next-step statements are historical;
Section 11 records the completion and scope of the dynamical-gauge benchmark.

REPRODUCE
  OPENBLAS_NUM_THREADS=1 python3 code/gauge_feedback.py
  OPENBLAS_NUM_THREADS=1 python3 code/audit_and_plot.py
Dependencies: Python, NumPy, SciPy, Matplotlib. Versions in environment.json.
The first writes the complete results JSON beside itself. The second runs
independent energy and matrix-exponential checks and regenerates the figure.
The PDF construction script retains original workspace and font paths and
needs path adaptation elsewhere; it is not needed to run the physics tests.

MODEL
This is the established staggered U(1) lattice Hamiltonian, not a derivation
from the previous Dirac walk. Six sites, five electric links, w=1, m=0.4,
g=1.2, T=1.5, hbar=1, open boundary electric fields fixed at one.
The initial state is a bare staggered vacuum with occupations 010101 and
all link fields one. It is not the interacting ground state. A particle-hole
pair count is defined relative to that staggered reference configuration.
Hopping changes occupations and the intervening electric flux jointly.
The hard-truncated raising operator obeys [E,U]=U but is not unitary at the
cutoff. Hamiltonian evolution and each local evolution gate are unitary.
No wraparound is used in the intended model. No projection step is applied.

WHAT PASSED
Exact local Gauss-value preservation over the full 15552-state cutoff-one
Hilbert space. Full-space and Gauss-reduced circuit agreement. Independent
exact-Hamiltonian references. Time-step convergence, norm and energy checks.
A coupling sweep demonstrates changes in matter and field evolution.

IMPORTANT FAILURE AND CORRECTION
A cutoff of |E|<=1 removes physically allowed states and overestimates the
pair mean by 29.82%, despite exact Gauss conservation and tiny time-step error.
For this particular six-site, fixed-boundary, three-fermion problem, |E|<=3
contains every allowed configuration (20 states); raising the cutoff further
does not change this sector. The corrected circuit reaches state-norm error
5.716e-8 at 4096 time steps, with norm error about 1.1e-12 and energy error
1.55e-7. These are errors against the finite Hamiltonian, not against continuum
QED or observations. Increasing system size needs a new cutoff study.

NEGATIVE CONTROLS
Frozen links and cyclic electric wraparound preserve norm and net charge,
but violate the intended U(1) Gauss constraints. Their recorded sum <G_i^2>
is about 2.48 and 16.14. Small time-step errors cannot fix the wrong algebra.

SCOPE
Known-method benchmark, not publication novelty, an efficiency advantage,
a chiral Standard Model, a 3+1 photon sector, or quantum gravity. Neither a
continuum extrapolation nor coupling to the prior curved-background branch
was performed. That integration is the next requirement, followed eventually
by matter backreaction and gravitational constraint preservation.

The cumulative earlier archive is in prior/. SHA256SUMS.txt covers this package.
