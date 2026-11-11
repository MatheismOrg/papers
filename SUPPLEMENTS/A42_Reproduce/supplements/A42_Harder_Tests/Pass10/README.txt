A42 pass 10 - expanding background and two time coordinates

The updated main paper has 29 pages. Section 10 is new (pages 28-29).
Section 9.7 now points to this completed test. The other 26 earlier pages are
unchanged. The cumulative pass-9 package is included, including prior tests.

Run from the extracted package directory:
  python3 code/expanding_dirac.py
  python3 code/audit_expansion.py
Requires Python, NumPy and SciPy. The versions used are in environment.json.
The first script recomputes all 40 circuit cases and writes JSON beside itself.
The audit independently evolves a full real-space chain and checks it against
the exact Bloch reduction, then checks reference consistency and limiting controls.
The PDF build script records the original workspace layout and installed fonts;
it needs path adaptation to rebuild elsewhere, and is not needed for the tests.

PHYSICS
This is a prescribed 1+1 FLRW background, not dynamical quantum gravity.
The scale factor is A(eta)=1.75+1.25*tanh(eta/tau), with mass 0.7 and momenta
1 and 2. Eta endpoints are +/-10*tau. The two clocks obey dt=A*deta.
Conformal Hamiltonian: k*sigma_z + m*A*sigma_x.
Proper-time Hamiltonian: (k/A)*sigma_z + m*sigma_x.
Both use chi=sqrt(A)*psi. Both lattice parity branches remain in the test.
The eight selected modes are lifted to the 70-state four-fermion sector;
Fock quantization is assumed. This is not a total cosmological particle density,
a full vacuum sum, a chiral Standard Model, an interacting field calculation,
or a spacetime backreaction model.

RESULT
The same gate prescription reproduces selected-mode production, approaches
agreement between clock descriptions, and passes independent normalization,
massless, static, sudden-limit and slow-limit checks. The limiting sudden and
slow tests use the continuum solver. The supplied-geometry normalization term
is necessary; omitting it gives physical norm growth approximately sixfold.
An incorrect conformal mass term remains unitary but does not converge to the
specified physics. These are implementation checks, not a new physical law.

INTERPRETATION OF FOURFOLD ERROR REDUCTION
For second-order error E(a) approximately C*a^2, halving spacing gives E/4.
This is not fourfold improvement over another code at equal work. No matched
comparison against a published implementation or novel superiority is established.
At finite spacing the two time-coordinate discretizations differ. Refinement
reduces that discrepancy, but a fundamental discrete theory still needs a
physical time/update prescription or a consistent relational alternative.

PRIOR WORK
Fulgado-Claudio, Sanchez Velazquez and Bermudez, Quantum 7, 1042 (2023),
https://arxiv.org/abs/2212.01355, studies fermionic production in expanding
spacetimes and lattice regularizations. Our calculation uses a different supplied
smooth expansion and is a validation study; no claim of priority is made.

NEXT TARGET
Dynamical gauge fields, Gauss-law preservation and interacting evolution provide
stronger further tests than repeating free-field matches. Geometry backreaction
and gravitational constraint preservation remain later requirements.
