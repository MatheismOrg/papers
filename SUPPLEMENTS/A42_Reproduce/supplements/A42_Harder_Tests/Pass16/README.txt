A42 - Pass 16: joint space-time and vacuum tests
================================================

The paper appends section16, pages44-46, retaining the first43 pages.
This continues the scalar branch of passes14-15. It does not yet couple
dynamical geometry or derive the source from the gauge/fermion branch.

What changed
------------
The main runs now use one local kick-drift-kick (KDK) step with h=nu*epsilon.
No independent adaptive continuous-time solver evolves those states.
The scalar geometry is sampled at each step's time midpoint. The full
Brillouin zone is retained in momentum quadrature for the infinite chain.

The invariant Gaussian vacuum of the stable frozen step differs from the
continuous-time Hamiltonian vacuum. Importing the latter produces apparent
excitations relative to that same continuous vacuum, including an exact
massless-nearest first-step energy64*nu^6/(35*pi*epsilon^2). This increases
as the lattice is refined at fixednu. It is a vacuum/update mismatch, not
physical creation from the step's invariant vacuum.

Matched invariant vacua plus an explicitly chosen shadow-generator source
approach the previous continuum energy and pressure under joint refinement.
The60cases cover2geometries*2spatialkernels*3nu*5spacings, with no fitting.

Physics conventions
-------------------
hbar=c=1, signature(-,+), metric a(eta)^2*(-deta^2+dx^2).
Mass mu=1/sqrt(pi). Geometries a=1+0.5*tanh(eta/tau) and
a=1+0.5*sech(eta/tau)^2, tau=.6. Initial eta=-12*tau, finaleta=12*tau.
The initial state is the finite-time invariant Gaussian vacuum for the
frozen step at that geometry, not an asserted exact past-infinity state.
epsilon is comoving; proper spacinga*epsilon changes during expansion.
nu=.25,.5,.75; epsilon=.4,.2,.1,.05,.025. Uniform stability is checked.

For omega^2=kappa^2+mu^2*a^2, define
  y=h*omega/2, d=1-y^2>0, theta=2*asin(y), Omega=theta/h,
  r=omega*sqrt(d).
The invariant covariance is diag(1/(2*r),r/2). The selected positive
quadratic generator is HF_eta=Omega*(r*q^2+p^2/r)/2. This is a specified
lift of the stable map. It is not an absolute many-body quasienergy ground
state. The logarithmic generator is generally spatially nonlocal although
the KDK update factors into local gates.

The pressure convention is -d(HF_eta/a)/da per comoving length, holding
conformalh, comovingepsilon and canonicalvariables fixed. The zero-point
terms and second-order adiabatic terms of the smooth shadow generator are
subtracted. Applying that subtraction to the frozen midpoint sequence is
an explicit convention tested in the limit. The finite-step source is not
proved exactly local, covariant or conserved. Work residuals are retained.
mass_response in the data is the generalized a*dHF/da=mu*dHF/dmu response;
it is not simply mu^2*<phi^2> atfiniteh. The compensated trace is a
spectrum/subtraction diagnostic in which the evolved state cancels.

Reproduce
---------
Python3 with NumPy,SciPy,Matplotlib,PyMuPDF,reportlab,fonttools.
See environment.txt for exact versions. From code/:

    OPENBLAS_NUM_THREADS=1 python3 step_scalar.py
    OPENBLAS_NUM_THREADS=1 python3 audit_shadow.py
    OPENBLAS_NUM_THREADS=1 python3 audit_steps.py
    OPENBLAS_NUM_THREADS=1 python3 audit_precision.py
    python3 make_figure.py
    python3 build_update.py --source ../source/A42_PrePass16.pdf --output ../A42_Shortest_Rules_Search.pdf

The PDFbuilder's installed LiberationSans and DejaVuSans paths may need
adjustment on another machine. No numerical input requires network access.
continuum_reference.py and spatial_reference.py are imported helpers copied
from earlier passes. Do not run them standalone for this pass. --pilot in
step_scalar.py overwrites results with only5cases; the paper requires the
full60-case run (complete=true,pilot=false).

step_results.json contains all60rows and continuum reference metadata.
Each row names its NPZhistory, with decimalpoints replaced by p.
The histories include energy,pressure,subtractions,normalization and work.
static_results.json contains the imported-vacuum controls and directly
computed matched-vacuum single-step residuals.

Independent audits
------------------
audit_shadow.py imports no main evolution implementation. It checks exact
KDKmatrix rotation, invariantvacuum, scipy matrixlog, particle formulas,
metricderivatives, adiabatic Wardidentity, continuous-time subtraction limit,
and the analytic first-step integrated massless energy/number densities.
audit_steps.py directly evolves q,p by local kicks, separately differentiates
the metric energy, measures static invariantvacuum stability over100steps,
and includes a symplectic but unstable negativecontrol.
audit_precision.py doubles momentumquadrature, extends/tightens the
continuumreference, and checks the last joint refinements. It also refreshes
the direct staticcontrols after the finalimplementation.

Interpretation
--------------
The numerical error tables are differences from a specified continuum scalar
reference in modelunits, not observationalaccuracy. Relative pressure errors
divide the maximum difference by the maximum continuum pressure magnitude
on the corresponding sampled times; they do not divide by values crossing0.
The first-step mismatch table changes elapsedtime withh; a separate column
compares at the same time1.2. The exact epsilon^-2 one-step law is massless;
the tabulated massive runs approach its ultraviolet behavior.
Spatial improvement need not win at fixed h/epsilon because temporal error
and cancellations also contribute. No parameter is fitted to rescue a rule.

A geometry-dependent scalar phase can leave all tested field covariances
and particleprobabilities unchanged while shifting the chosen gravity source.
The paper exhibits exp(-i*h*L*Lambda*a^2). Retaining this term in the lifted
generator and total source while holding the finite vacuum counterterm fixed
gives constant densityLambda and pressure-Lambda. Redefining that counterterm
can absorb the shift; field-only data cannot choose this finite vacuum term.
This illustrates missing phase/action information, not a prediction or
resolution of the cosmologicalconstant.
The next construction must specify local metric response, phase and lapse,
then evolve geometry and connect the interacting matter sector.

Contents
--------
A42_Shortest_Rules_Search.pdf:updated46-pagepaper.
source/A42_PrePass16.pdf:previous43-pagepaper.
code/:complete current scripts,JSON,NPZ,figure and logs.
previous_pass15.zip:prior cumulative archive.
SHA256SUMS.txt:hashes of every other included file.

Primary references
------------------
V.Moncrief,Phys.Rev.D28,2485(1983),
https://doi.org/10.1103/PhysRevD.28.2485.
L.Ni,Z.Hu,On the relation between the velocity- and position-Verlet integrators,
https://arxiv.org/abs/2407.18726.
These are established quantum-leapfrog and harmonic-generator constructions.
This pass tests our implementation and its joint limit; it claims neither a
new integration algorithm nor a completed discrete theory of everything.
