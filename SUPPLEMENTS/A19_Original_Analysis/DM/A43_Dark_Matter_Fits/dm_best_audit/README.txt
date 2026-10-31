Independent A43 benchmark and targeted-repair audit, 28 October 2026

No fitting-module prediction routines are imported by these audit programs.
All main analyses passed independent equation, score, split, parameter-bound,
optimization-record and source-input checks. Running verify.py executes all
checks and writes verification.json and the four detailed reports.

1. verify_mond.py: 79 prediction tables; 146,828 predictions; 239 fit sets and
471 optimizer attempts. Checks all five-fold seeds, full populations, complete
26/18 dwarf hosts and full CLASH/X-COP sources. Selected coefficients must match
the best declared start. Metrics and all recorded attempt costs are independently
computed. verify_mond_bootstrap.py checks 26 paired descriptive intervals.

2. verify_halo.py: 20 prediction tables; 15,760 predictions; 1,254 fit sets and
3,762 attempts. Independently matches all 2,700 SPARC rows to raw rotation files,
checks signed gas/star terms, both contiguous radial holdouts and all training
indices. NFW and Burkert enclosed masses are checked against adaptive density
integrals; all predictions, scores, fit costs and 54 bootstrap comparisons pass.

3. verify_cluster.py: 192 prediction tables; 6,128 predictions; 460 fit sets and
1,610 attempts. Checks independent direct MGAS input, all frozen A43 targets,
literal gas parameters when original source is present, all case/source/host
splits, all costs and coefficients. 344 adaptive density-integral checks span
three taper lengths plus unweighted mass. All 3,064 gas predictions have refined
quadrature and unresolved-inner-mass bound checks. verify_cluster_bootstrap.py
checks six descriptive intervals. Optional source code is not redistributed;
normal numerical reproduction uses the supplied parameter CSV. If the original
source is absent, the report honestly records primary literal verification false.

4. verify_physics.py and ANALYTIC_PRECHECK.txt check limiting behavior and the
reported physical diagnostic grids. The transient correction restores a common
formal point-source asymptotic BTFR normalization and retains nonnegative
point-source inferred extra density; it does not establish observed finite-radius
BTFR agreement, arbitrary-profile positivity, a conservative action, or a metric
and cosmological theory. Potential taper fails monotonic response only in the
extended grid, not at current measured coordinates. The two-bump family also
fails at 37 measured coordinates.

Total: 291 tables; 168,716 predictions; 1,953 fit sets; 5,843 attempts; no failures.
All 86 bootstrap comparisons condition on frozen predictions; none include
refitting, adaptive model-selection or observational/systematic uncertainty.

These audits validate computation under declared assumptions. They do not turn
previously inspected data into independent confirmation or establish best-in-
literature superiority. Frozen cluster eta=10,Rc=400 kpc is CLASH-informed and
its repeated CLASH-source table is not an independent source-exclusion test.
