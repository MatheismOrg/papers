ONE SHARED FORMULA: EXPLORATORY FOLLOW-UP TO A43
27 October 2026, US Pacific date

All shared candidates use one coefficient vector for every system and no object
or population label in prediction. The labeled piecewise_RAR is a comparator.
14 model families / 497 fit calls; no convergence warnings. Same 2910 rows,
279 systems,278 conservative groups as A43. Three whole-group fold partitions,
whole-population holdouts and nested family selection. All 42 CV files checked.

A smooth power-sum acceleration scale essentially matches the previous three-
coefficient candidate internally and improves the whole-population-holdout score:
                         original-fold B   population-holdout B
Previous potential quad      0.249461           0.315811
New power sum                0.249630           0.277628
Potential plus radius        0.244682           0.308098
Piecewise five-pop RAR       0.281737           not defined for unseen class
Units are log10 acceleration RMS (dex), population/group-balanced, not percent.
Power sum's SPARC median speed error16.06% vs previous14.29% and galaxy-only RAR9.22%.
All-data coefficients and per-fold coefficients differ as required for validation.
The new family has5fittedcoefficients; alpha is at its imposed upper bound1.
At fixed r, gbar->0 leaves an acceleration floor, so this is not a vacuum theory.

These are exploratory results after repeated inspection of the same data.
Selecting a promising model after seeing holdouts means its score is not an
unbiased estimate of future discovery performance. No covariance/nuisance
likelihood; inherited dwarf, stellar, geometry and mass-inference assumptions.
No new CMB, Solar-System, shear-map or memory fit. No superiority-to-all claim.
The nested-selected internal B is0.242290; it selects other flexible models.

REPRODUCE
Install numpy scipy. From this archive's root:
  python dm_one_formula/search_one_formula.py
Keep sibling dm_broad_joint/combined_input.csv. The main A43 code/data package
contains original catalogs and selection provenance. This file repeats its exact
combined input; SHA-256 is recorded in results.json.

ATTRIBUTION
Lelli etal2017 Eq14 RAR plus acceleration floor is published prior art:
https://arxiv.org/html/1610.08981v2
Santos &Felizardo2026 reciprocal-log form is independently refitted here:
https://doi.org/10.3390/particles9030070
No original-paper score is equated with our heterogeneous benchmark. The other
families are empirical extensions, not newly established physical theories.

SEARCH PROMPT FOR PERPLEXITY / CONSENSUS
Find primary research comparing a single shared acceleration law across dwarf
spheroidals, SPARC disks, early-type galaxies, brightest cluster galaxies, cluster
lensing and weak galaxy lensing. Prioritize held-out prediction, shared parameters,
uncertainty modeling and numerical comparisons. Include EMOND, nonlocal gravity,
symbolic regression, two-component virial models and LambdaCDM predictions.
Distinguish shared-parameter predictions from per-object fits. Supply equations,
code/data, parameter counts, validation, limitations and DOI/arXiv links. Do any
published models depend on gbar*r through a sum of positive and negative powers?
