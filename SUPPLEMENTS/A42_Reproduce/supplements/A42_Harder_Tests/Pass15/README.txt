A42 - Pass 15: lattice stress and rule selection
===============================================

The updated A42 appends section 15, pages 41-43, preserving the previous
40 pages. This is a bosonic spatial-lattice branch, not yet the stress tensor
of the interacting fermion/gauge chain from passes 11-13.

Findings
--------
Three local positive harmonic kernels share the desired low-momentum wave
equation. Two approach the intended one-scalar stress reference. A squared
centered derivative introduces a second species (decoupled even/odd chains),
and approaches twice the intended source despite good conservation.

The exact identity is O_centered(epsilon)=2*O_nearest(2*epsilon) for matched
homogeneous stress, excitation and subtraction histories. The real-space
block decomposition and numerical histories independently verify this.
The centered rule is consistent as a TWO-field model; it is rejected only
for the intended ONE-field target. Scalar doubling is established prior art.

The lattice subtraction is integrated from second-order WKB terms over the
entire Brillouin zone. The continuum anomaly coefficient is not inserted
directly into that calculation. The subtraction prescription and quantum
oscillator structure are still assumed. The mass-compensated trace cancels
the evolved state algebraically: it is a spectrum/subtraction diagnostic.
Energy and pressure histories are separate dynamical tests.

Conventions
-----------
hbar=c=1; metric a(eta)^2*(-deta^2+dx^2); R=2*(a'/a)'/a^2.
Main scalar mass mu=1/sqrt(pi). Supplied geometries:
  a=1+0.5*tanh(eta/tau), or a=1+0.5*sech(eta/tau)^2; tau=0.6.
Initial instantaneous vacuum at eta=-12*tau; final eta=12*tau.
Spatial spacing epsilon is COMOVING, so physical spacing is a*epsilon.
Infinite lattice, not a finite ring; positive-k quadrature includes the
factor1/pi for both signs. All spectral zeros, including the BZ edge, are
resolved. The finite-ring checks are separate operator and time-step audits.

Kernels, with s=sin(k*epsilon/2):
  nearest:  4*s^2/epsilon^2
  improved: 4*s^2*(1+s^2/3)/epsilon^2
  centered: sin(k*epsilon)^2/epsilon^2

Pressure is obtained by metric variation of H_eta/a at fixed canonical
variables and comoving sites. Its spatial spectral factor is kappa^2,
not the continuum k^2 and not a momentum derivative of the dispersion.

Numerical reproduction
----------------------
Python 3 with NumPy, SciPy, Matplotlib, PyMuPDF, reportlab and fonttools.
Versions are listed in environment.txt. From code/:

    OPENBLAS_NUM_THREADS=1 python3 lattice_stress.py
    OPENBLAS_NUM_THREADS=1 python3 audit_precision.py
    OPENBLAS_NUM_THREADS=1 python3 audit_operators.py
    python3 make_figure.py
    python3 build_update.py --source ../source/A42_PrePass15.pdf --output ../A42_Shortest_Rules_Search.pdf

The PDF builder uses installed Liberation Sans and DejaVu Sans paths near
its beginning; adjust those font paths if necessary on another machine.
continuum_reference.py is a copy of pass14's stress_test.py, imported as a
module. Do not run it standalone to reproduce pass15. It is included so
the current scripts have no dependency on external folders or network.

lattice_stress.py runs3stencils*2geometries*6spacings=36cases, plus independent
continuum references. Adaptive continuous-time evolution isolates spatial
regulator effects; it is not a fixed fundamental space-time update.
Histories are named by geometry, stencil and spacing, replacing decimal
points with p. lattice_results.json lists their exact names and all metrics.

audit_precision.py varies lattice quadrature/evolution, extends continuum
cutoff128->256, independently refines work quadrature and computes the
massless anomaly coefficients using adaptive momentum integration.
audit_operators.py builds real-space matrices, checks pressure by metric
variation, audits WKB conservation and exact doubling, and implements
finite-range Gaussian kick-drift-kick steps on a16-site periodic ring.
The time-step audit compares covariances with an independent adaptive
real-space map. It does not establish a simultaneous space-time limit of
the renormalized source at fixed step-to-spacing ratio.

Error meanings
--------------
rho_sup_error and p_sup_error are maximum absolute differences against the
continuum history over the sampled time interval. Relative sup errors divide
by the maximum magnitude of the corresponding continuum observable, not
pointwise values that may cross zero. They are not observational errors.
The work residual has model units and is not a percentage. All three rules
pass conservation, so it cannot select the correct field count by itself.
The finest improved errors approach the numerical reference floor; use the
separate precision checks when interpreting them. Main table uses epsilon=.1.
The massless coefficient calculation first integrates with small positive
mass, resolving every infrared region, then examines the limit.

Physical and ToE scope
---------------------
This verifies known lower-dimensional quantum-field behavior under a
prescribed source-renormalization scheme. It does not derive GR, quantum
probabilities, Standard Model matter, a fixed Planck scale, or black-hole
singularity resolution. No observational data were fitted.
Two passing stencils show that these matches do not determine a unique
program. A finite program can describe unbounded structures; this run does
not supply a fundamental state encoding, precision rule or exact TM update.
The local oscillators are untruncated. A complete fixed-step rule must
specify its vacuum and stress, and then connect to interacting matter and
dynamical geometry before it can be assessed as a unified candidate.

Contents
--------
A42_Shortest_Rules_Search.pdf:43-page updated paper.
source/A42_PrePass15.pdf:previous40-page paper.
code/:source, JSON outputs, NPZ histories, figure and logs.
previous_pass14.zip:the previous cumulative archive. Earlier numerical
packages are retained there; the original paper records the first searches.
SHA256SUMS.txt:integrity hashes for all other included files.

Primary references
------------------
A. Buonanno, G. Cella, G. Curci, Phys.Rev.D51,4494(1995),
https://arxiv.org/abs/hep-lat/9506008.
S. Bellucci, C. Sochichiu(2004),
https://arxiv.org/abs/hep-th/0404143, section5.
Adiabatic-renormalization sources remain U17,U18 from pass14.
