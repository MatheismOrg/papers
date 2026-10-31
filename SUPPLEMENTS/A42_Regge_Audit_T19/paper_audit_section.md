### 8.9 Independent execution and equilibration audit

A new audit reran the supplied Regge C engine without changing its action, measure or update rules. The protocol was fixed before the main runs. At each of L = 3 and 4 and k = 0.05 and 0.065, three seeded chains began flat and three began from independently prepared rough configurations. Every rough preparation used the same prescribed 3,000-sweep evolution at k = 0.1; no configuration was selected for its observed outcome. Each of the 24 main chains used 3,000 equilibration sweeps followed by 12,000 measurement sweeps. All per-sweep data, starting preparations and ending configurations are retained. These are new executions and independently written diagnostics, not a second implementation of the gravitational dynamics.

The geometry self-tests reproduced, including the flat-start deficit angles and incremental/full recomputation check. The statistical audit used ordinary split R-hat across all six chains at each parameter pair, with a prespecified screening threshold of 1.05, and batch-mean, drift and start-dependence diagnostics. Passing this screen would not prove equilibrium; failing it prevents treating the pooled chains as a settled ensemble. Curvature means below are averages over each start class and are descriptive when the chains have not mixed.

| L | k | Flat-start mean R | Rough-start mean R | Split R-hat | Screen |
|---|---|---|---|---|---|
| 3 | 0.05 | -1.0470 | -0.6731 | 1.573 | fails |
| 3 | 0.065 | -0.6022 | 44.1443 | 13.032 | fails |
| 4 | 0.05 | -1.0655 | 31.2166 | 5.498 | fails |
| 4 | 0.065 | -0.6793 | 45.4855 | 17.716 | fails |

All four parameter groups fail the prespecified R-hat screen. The table and trajectories distinguish initial-state dependence from an established equilibrium phase. Negative mean curvature after a flat start is insufficient to demonstrate that the ensemble has equilibrated. These finite-run tests neither locate a critical point nor establish the nonexistence of an equilibrium limit; they restrict what can be claimed from the sampled histories.

All 24 chains satisfy the scaling identity in their measurement averages to within 0.09 percent, despite the failed mixing screen. This separates a successful overall-scale consistency check from equilibration of the geometry. It provides a concrete example of an implementation check that cannot establish the sampled phase.

![Regge audit trajectories. Each curve averages 100 sweeps for display; all diagnostics use the unaveraged measurement histories. Gray shading marks equilibration sweeps.](regge_audit.png)

The next tests are normalizability of the specified ensemble and convergence from both start classes at longer durations, followed by finite-size scaling and curvature correlations only where convergence is supported. A Newtonian potential and a Lorentzian matter-coupled theory require separate tests. The original histories remain unavailable; this audit supplies its own complete histories rather than claiming to reconstruct them.
