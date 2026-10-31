A43 physical-limit and uncertainty follow-up, 28 October 2026

Reproduce from this supplement root:
  python dm_physical_repair/test_repair.py
  python dm_physical_repair/verify.py
Requires Python, NumPy and SciPy; the sibling dm_one_formula and dm_broad_joint
packages supplied in the complete supplement must remain in their directories.

The three candidates were fixed before this follow-up's runs, after extensive
prior inspection of the benchmark. This is exploratory, not a preregistered
independent confirmation. PROTOCOL.txt records the local choices. results.json
contains starts, parameter bounds, all-data coefficients, every fold coefficient,
all population metrics, convergence warnings and 5000-resample paired intervals.
No new astronomical observations or full uncertainty likelihood are introduced.
The bootstrap uses identical group multipliers across populations and normalizes
within each population. It conditions on fitted out-of-fold predictions and
excludes refitting, selection and measurement uncertainty. It is not a posterior.

New families: softened inverse power (6 coefficients), bounded tanh quadratic
(5), and original power sum with inverse exponent fixed at one-half (4).
The bounded form has better low/high-acceleration and isolated finite-mass limits.
No candidate uniformly improves fit. The half-power form's lower population-
holdout score is not resolved by the descriptive interval and it retains rising
asymptotic isolated circular speed. Boundary hits are reported in the paper.

63 fitted models; 12 prediction files; no convergence warnings. verify.py checks
prediction-file residuals, metrics, group assignment and separately implements
the all-data formulas. verification.json records the result.

Source review: Li et al. 2018 (doi:10.1051/0004-6361/201732547), Desmond 2023
(doi:10.1093/mnras/stad2762), Zhao and Famaey 2012 (arXiv:1207.6232), Bilek et al.
2026 (arXiv:2603.23591v1), and Desmond et al. 2024 (arXiv:2401.04796). Their metrics
are not interchangeable with this cross-population score. Bilek et al. supports
an environment-focused next test, not a measured correction in this analysis.
