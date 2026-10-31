A42: Curved interacting gauge tests, pass 12
Completed 2026-10-28. Section 12 is on pages 32-34 of the updated paper.

RESULT
One unchanged local splitting prescription was applied to the established 1+1
Schwinger model in supplied static and expanding conformal geometry. Its flat
limit reproduces pass 11. The electric Hamiltonian coefficient scales as Omega^2
in conformal time, while the mass coefficient scales as Omega. Two intentionally
wrong electric factors still pass unitarity, Gauss-law and global-clock checks
but yield different quantum states. These are diagnostics, not new laws.

LIMITS
This is a fixed six-site model with no spatial/volume extrapolation. Electric
cutoff 3 exhausts its 20-state physical sector for the specified boundaries.
The cutoff-1 full tensor-product test checks gate structure only. There are no
independent transverse photons in this open 1+1 gauge model. Geometry is supplied;
there is no gravitational backreaction, Einstein-constraint solver, black-hole
center treatment, new experimental prediction, or demonstration of a ToE.
Bare occupations are not identified with asymptotic interacting particles.

REPRODUCE THE SCIENTIFIC CALCULATIONS
Python 3 with NumPy, SciPy, Matplotlib. From the code directory run:
  OPENBLAS_NUM_THREADS=1 python curved_gauge.py
  OPENBLAS_NUM_THREADS=1 python audit_and_plot.py
The first script writes curved_gauge_results.json. The second writes the
independent audits and plot. Both raise an error if an asserted check fails.
gauge_feedback_pass11.py is the unchanged earlier reference for the flat-limit
cross-check. No helper, account credential, or external network is needed.

REFERENCE INDEPENDENCE
Local gates explicitly change fermion occupations and electric links. A separate
bit-word enumeration eliminates flux by Gauss's law and builds a dense Hamiltonian.
Adaptive DOP853 integrates this Hamiltonian; static cases use matrix exponentials.
The physical-time reference jointly integrates the clock, while the circuit uses
an independently inverted analytic time map. References are checked at tighter
precision and Hamiltonian work derivatives by finite differences. These checks
establish the finite-system implementation, not its continuum physical accuracy.

PDF REBUILD
PyMuPDF, ReportLab and FontTools are needed in addition to the above. The builder
uses installed Liberation Sans and DejaVu fonts (paths are near its top). From the
code directory, after adjusting those paths if needed:
  python build_update.py --source ../source/A42_PrePass12.pdf --output ../A42_Rebuilt.pdf
This appends three pages and asserts all 31 previous pages have unchanged text.
The supplied final PDF was rendered and visually checked.

CONTENTS
Main paper; source/pre-pass PDF for regeneration; code/results/plot/logs; runtime
versions; file hashes; and previous_pass11.zip containing earlier cumulative work.
Earlier sections retain the historical next-step statements; section 12 gives the
current boundary. The current work extends the interacting benchmark branch;
it is not a derivation from the earlier free Dirac walk.

PRIOR ART
Budd et al. (2026), JHEP 04, 183; arXiv:2601.02331v2, Eqs. 80, 89, 98.
Sellapillay, Arrighi and Di Molfetta (2022), Scientific Reports 12, 2198.
Momeni and Channuie (2021), European Physical Journal C 81, 534, section 2.
The metric Hamiltonian and gauge-preserving framework are established methods.

NEXT TARGET
Separate spatial-regulator and boundary errors; then choose a specified gravity
model before testing backreaction and gravitational constraint propagation.
In two spacetime dimensions the ordinary Einstein tensor vanishes identically;
a dilaton reduction or higher-dimensional model is needed for that next stage.
