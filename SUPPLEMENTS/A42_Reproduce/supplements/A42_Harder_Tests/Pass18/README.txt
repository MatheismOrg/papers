A42 PASS 18: BLACK-HOLE ENDPOINT TESTS
29 October 2026

READ FIRST
The outcome is a bounded candidate screen, not a discovered Theory of Everything,
full Standard Model simulation, new evaporation solver, or resolved singularity.
The existing A42 paper is preserved through page 50. Pages 51-53 report this pass.
No changes were made to A41 or the separate DM/A43 project.

WHAT WAS DONE
1. Reproduced the exact RST shock/evaporation benchmark in reduced 1+1-dimensional
   semiclassical dilaton gravity. It includes a matter anomaly; geometry is not
   fully quantized. The strong-coupling endpoint is outside its controlled regime.
2. Updated the energy benchmark to the published 2026 complete-action Bondi
   charge, while retaining the traditional matter flux as a distinct quantity.
   Each charge has its own balance flux. Mixing them creates a false violation.
3. Tested 72 finite-step ledgers: six masses, four steps, three pairings. Of these,
   48 correctly paired cases conserve their budget to 1.43e-14 in energy units A.
   The 24 deliberately mixed cases retain their analytically predicted error.
4. Screened old-mass-zero deletion, inverse-branch switching, gap flooring and
   reflection. No complete endpoint law follows from these prescriptions alone.
   Floors/reflections could be new physics but need changed dynamics and sources.
5. Derived a conditional obstruction for keeping the unspliced RST bulk solution
   and changing only the smooth local reconstruction rho=phi=g(w). Generic late
   constant-w observers retain logarithmically growing curvature. A concrete
   softplus reconstruction was checked with 75-digit Decimal mixed derivatives.
   This fails a global bounded-curvature/late-vacuum goal, not a claimed theorem
   of finite-time geodesic incompleteness or a ban on other quantum gravity laws.

WHAT IS OUR CONTRIBUTION HERE?
The charge correction and RST solution are published results, not our discoveries.
The tests, their organized negative controls, explicit candidate calculations and
conditional reconstruction screen are this pass's work. No literature-priority
claim has been established for the generalized screen. Small numerical errors
verify calculations; they do not measure progress toward a ToE or superiority
over established physics codes. This work narrows candidate selection.

ARCHIVE LAYOUT
tmp/toe_pass18/endpoint_ledger.py       Main energy ledgers and clipping screen
tmp/toe_pass18/audit_rst.py             Independent charge/geometry audit
tmp/toe_pass18/audit_endpoint.py        Independent endpoint/reconstruction audit
tmp/toe_pass18/rst_derivation.txt       Detailed conventions and source map
tmp/toe_pass18/endpoint_review.txt      Chain-rule proof and candidate limitations
tmp/toe_pass18/*results.json            Main ledger results
tmp/toe_pass18/audit_*.json             Independent numerical results
tmp/toe_pass18/ledger_*.npz             All 72 sampled ledgers
tmp/toe_pass18/build_appendix.py        Three-page PDF append and preservation check
tmp/toe_pass18/source/                 Unaltered 50-page pre-pass A42
tmp/toe_pass18/appendix_text.txt        Extracted new-page text for review
tmp/toe_pass18/pdf_validation.json     Page count and prior-page preservation
tmp/toe_pass18/environment.json        Recorded dependency versions
tmp/toe_pass18/CONTINUATION.txt         Compact scientific handoff
output/toe_pass18/                     Updated 53-page A42 and standalone appendix
SHA256SUMS.txt                         Integrity manifest of package contents

REPRODUCE THE SCIENTIFIC TESTS
Unzip into an empty directory and run from the extracted archive root:

  python3 tmp/toe_pass18/endpoint_ledger.py
  python3 tmp/toe_pass18/audit_rst.py
  python3 tmp/toe_pass18/audit_endpoint.py

Dependencies: Python 3, NumPy, SciPy. The high-precision calculation uses Python's
standard-library decimal module. No external data or network is needed. Each
program writes its JSON result beside its source; the first also writes NPZ files.
Check all_checks_passed in ledger_results.json and audit_endpoint.json, and
all_pass in audit_rst.json. All three were true in the recorded run.

To rebuild the paper:

  python3 tmp/toe_pass18/build_appendix.py

Additional dependencies: ReportLab, PyMuPDF, and DejaVu Sans fonts in the usual
Linux directory (the script has a runtime-specific Liberation Sans fallback).
This creates three new pages and requires the original fifty pages to have
identical extracted text and identical half-scale raster renders after appending.
The three new pages were also visually inspected. A different font or library
version can alter layout; the script asserts exactly three added pages.

KEY CONVENTION WARNINGS
- A=kappa*lambda/4; our mu=M/A, unlike the mu symbol in the 2026 source.
- U=log(q) in endpoint_ledger.py is shifted affine retarded time for lambda=1.
- u in audit_endpoint.py is instead the Kruskal incoming coordinate x+.
  The paper uses s for that coordinate to avoid this clash.
- W=z-log(z) is normalized; audit_endpoint.py uses W for the unnormalized field
  called w in the new paper pages. See derivation files for exact conversions.
- The old-zero clipping postshock horizon claim requires mu>=1.536591355283385.
- The historical negative junction matter impulse is not interchangeable with
  the positive corrected Bondi endpoint charge.
- Numerical resolution is not a chosen fundamental physical lattice spacing.
  A physical modification scale, if proposed, must remain fixed during solver
  refinement. This does not require taking a physical discrete theory's scale
  to zero.

PRIMARY SOURCES
Russo, Susskind and Thorlacius (1992), End point of Hawking radiation.
  https://doi.org/10.1103/PhysRevD.46.3444
  https://arxiv.org/abs/hep-th/9206070
Eune, Gim and Kim (2014), Note on explicit form of entanglement entropy in the
RST model, Phys. Rev. D 89, 127501.
  https://arxiv.org/abs/1401.0590
Barenboim, Frolov and Kunstatter (2026), Quantum-corrected Bondi mass for 2D
Hawking radiation, Phys. Rev. D 114, 024067.
  https://doi.org/10.1103/q7hb-xjtf
  https://arxiv.org/abs/2510.18298v2
The third source changes the benchmark accounting, not the singular geometry.

PREVIOUS ARTIFACTS
Canonical A42 (same identity, updated by this pass):
  https://chatgpt.com/library?file_id=libfile_5c6a9820e2148191a10e6ac400f53138&account_id=personal
Pass 17 gravitational feedback archive (not duplicated in this archive):
  https://chatgpt.com/library?file_id=libfile_af479bee01b8819199e22d8abf97fb9a&account_id=personal
