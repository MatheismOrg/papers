# A42 engines in progress (November 11, 2026)

Working code behind results that are not yet in the paper. Nothing here is cited by a published version of A42. Each piece
lists what it was checked against.

## regge_chip: the gravity lattice on a graphics chip

The 4D Regge lattice of the paper (A42_SA74_Regge.c and its multi-core version regge_r2_par.c, included unchanged) with the
single-edge update written once for three builds: CPU double precision, CPU single precision and Apple Metal single precision.

- `rk.h`, `regge_gpu.c`: the update with the original engine's lookup tables.
- `rk2.h`, `regge_tab.c`: every neighbor found by formula. The state is the squared edge lengths plus three numbers per
  triangle; the only tables are two small ones (27 KB and 19 KB) built at start.

Build: `gcc -O2 -DWITH_ORIG -o regge_tab_chk regge_tab.c -lm -lpthread` (checks), `gcc -O2 -o regge_tab regge_tab.c -lm`
(CPU), and on a Mac `clang -x objective-c -fobjc-arc -O2 -DUSE_METAL -o regge_tab regge_tab.c -framework Metal -framework Foundation`.

Checks (`regge_tab_chk check L k lambda a eps seed cfg|-`):

- formulas against the original engine's tables for every edge, simplex and triangle: no mismatch at 5^4, 6^4, 8^4, 9^4, 12^4
- recomputation against the original engine: differences at most 4e-15
- formula update against table update in double precision with the same random numbers: the same decision on every edge in
  every sweep tested
- single against double precision: no differing decision in the first sweeps; later a rare knife-edge decision differs

Measured on an Apple M4 (seconds per sweep, update plus recomputation): 1.8 at 24^4, 6.1 at 32^4, 34 at 48^4; memory on the
chip 0.22, 0.69 and 3.5 GB. The eight-thread CPU engine takes 14.3 seconds per sweep at 32^4.

Not yet in the formula engine: step-size tuning, global dilation moves, the original output columns, checkpoints.

## cdt: causal dynamical triangulations in 3+1 dimensions

- `cdt4c.cpp`: the engine (SHA-256 d62bbfc3bc8aaeb4...), written from the published description of the model; point
  (kappa0, Delta) = (2.2, 0.6). Option KAVG freezes kappa4 at its mean over the second half of thermalization.
- `an5.py`, `planC.py`: the width of the universe from the volume profile (fit to cos^3). `planC.py` is the analysis fixed
  before the plan C runs.
- `sd.py`: spectral dimension by diffusion on the dual graph.
- `cdt_gate.py`: the paper's three dimension tests on the causal order of a CDT configuration (exploratory).

Status: a reproduction in progress. The width scales as the fourth root of the volume; its coefficient is below the
published 0.59 at the sizes run so far.

## rule_classes

`run_classes.py` and `summarize.py` run complete rule classes with the TOE@home client (see ../A42_TOE_at_Home), without a
coordinator.
