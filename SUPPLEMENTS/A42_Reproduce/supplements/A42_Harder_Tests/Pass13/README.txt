A42: Spatial refinement and continuum comparison, pass 13
Completed October 28, 2026 (America/Los_Angeles). New section 13, pages 35-37.

WHAT CHANGED
The known interacting staggered Schwinger benchmark is tested at 6 through 20
sites (up to 184,756 neutral configurations). Complete matter sectors are used;
there is no electric truncation. The known mass improvement m_lat=m-g^2*a/8 is
compared to the original unshifted prescription. A wall-length convention is
fixed from the exactly solvable free-chain spectrum: ell=(N+1)*a. Separate runs
with L=N*a are also retained, including a case where boundary/mass errors cancel
and the unshifted gap is closer. No convention or correction is fitted to outputs.

STRONGER CONTINUUM TEST
Massless interacting fermions plus electric flux are compared with the known
bosonized free field of mass g/sqrt(pi) in the same supplied homogeneous expanding
metric. Continuum oscillator modes are solved independently, with 20/40/80/160
mode cutoffs. The 20-site vacuum probability is 0.932641928 versus 0.932804285,
but excitation energy still has a 2.13% discrepancy. Energy discrepancy was
31.8% at six sites. Probability error crosses zero and is not monotone in
magnitude: the smaller N=18 probability error does not certify continuum accuracy.

BOUNDARIES AND PARAMETERS
Static spectral test: target mass 0, g=1, original boundary flux B=0, ell=4 and 8.
Massless expansion: g=1, ell=8, tau=0.15, Omega=1+0.5*tanh(eta/tau), eta from
-10*tau to +10*tau. B=-1/4 makes the shifted field E_tilde=E-(-1)^i/4 zero at both
walls. The improved electric+mass energy equals J*sum(E_tilde^2) plus a scalar;
this is verified over every basis state. It matches the intended zero shifted
Dirichlet flux and is a different boundary problem from the earlier B=1 example.
Massive expansion: retains m=.4, g=sqrt(2.88), B=1, tau=.6; ell=3.5. N=6 exactly
reproduces pass 12. The old vacuum probability changes substantially with N.
Fixed-a=.5 volume sweeps are stored separately. Each run starts in its own
interacting ground state. The homogeneous improved mass coefficient is
m*Omega-g^2*a*Omega^2/8. This is not a general curved-space renormalization theorem.

INTERPRETATION
These reproduce published lower-dimensional physics and diagnose regulator
errors. They are not a new ToE, proof of fundamental discreteness, gravitational
backreaction, independent transverse photons, or a new observational prediction.
The original massive six-site value is correct for its finite Hamiltonian but
has not become a regulator- or volume-independent continuum result. Tight ODE
errors and spectral residuals must not be confused with those physical errors.
Sparse/adaptive evolution is used for the larger spatial sweeps. The local gate
prescription is separately checked at ten sites; it was not run at every N.
The known continuum bosonization and mass improvement are not claimed as novel.

REPRODUCE (from the code directory; Python + NumPy/SciPy/Matplotlib)
  OPENBLAS_NUM_THREADS=1 python lattice_scaling.py
  OPENBLAS_NUM_THREADS=1 python dynamics.py
  OPENBLAS_NUM_THREADS=1 python extend_dynamics.py
  OPENBLAS_NUM_THREADS=1 python audit_scaling.py
  OPENBLAS_NUM_THREADS=1 python bosonization_test.py
  OPENBLAS_NUM_THREADS=1 python audit_bosonization.py
  python make_figure.py
Dynamics initially runs to 16 sites; the extension resumes the saved file through
20 sites and does not repeat completed N=18/20 rows. Re-running dynamics.py starts
that sweep over, after which the extension should run again. All tests use a CPU;
no quantum hardware or network is needed. Expect several minutes at the largest
sizes and sufficient RAM for sparse matrices and vectors. Each assertion must pass.

VALIDATION
Dense/sparse eigenvalues, exact free-chain gaps, pass-12 Hamiltonian and evolution,
local gate convergence, tighter ODE tolerances, continuum mode cutoff/tolerance,
free no-production and sudden-change oscillator controls, complete-basis shifted
field identity, and boundary/mass-control variants are included. A nonsymmetric
seed prevents restricting the Lanczos iteration to a symmetry sector. No random
sampling of physical outcomes is involved. All final assertions passed.

PDF REBUILD
The builder needs ReportLab, PyMuPDF and FontTools. It uses Liberation Sans and
DejaVu font paths near its top; adjust for another installation. From code:
  python build_update.py --source ../source/A42_PrePass13.pdf --output ../A42_Rebuilt.pdf
The builder checks that all previous 34 pages retain their text. The delivered
three new pages have been rendered and visually inspected.

FILES
The main paper, source/pre-pass PDF, code, full result JSON, logs, figure, runtime
versions, checksums and previous_pass12.zip (earlier cumulative work) are included.
Scripts do not contain account tools, credentials or file-transfer helpers.
Earlier next-step statements in the paper are historical; section 13 is current.

REFERENCES
Dempsey, Klebanov, Pufu and Zan, Physical Review Research 4, 043133 (2022),
https://arxiv.org/abs/2206.05308 (v3, especially the open-boundary discussion).
Budd, Florio, Frenklakh and Mukherjee, JHEP 04 (2026) 183,
https://arxiv.org/abs/2601.02331v2 (curved Schwinger model and bosonization).

NEXT
Use the tractable bosonic reference to benchmark stress-energy and reciprocal
matter/geometry energy exchange in a specified reduced-gravity model. Ordinary
1+1 Einstein gravity needs additional structure. Massive calculations still need
spatial and volume control; a fixed fundamental-spacing proposal needs tests of
its physical deviations rather than identifying numerical convergence with proof.
