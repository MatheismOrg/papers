A43 environment and history pilot, 28 October 2026

Run from the complete A43 supplement:
  python dm_environment/prepare.py
  python dm_environment/run_tests.py
  python dm_environment/verify.py
Dependencies: Python, NumPy, SciPy. Keep the sibling dm_fit, dm_broad_joint,
dm_broad_dwarfs and dm_one_formula directories supplied in the package.

Data: Chae et al. (2021) ApJ 921:104, doi:10.3847/1538-4357/ac1bba, Table 3 only.
Author archive downloaded from:
https://prof.sejong.ac.kr/chae/resource-links/notice.do?articleNo=5205&attachNo=1496&mode=download
The Table 2 external-field values inferred FROM rotation curves are not predictors
and are deliberately omitted from this working data folder. The retained Table 3
fields are baryonic-environment reconstructions under max/no-clustering models.

Weisz et al. (2015) ApJ 804:136, doi:10.1088/0004-637X/804/2/136; arXiv:1503.05195v1,
Table 1, page 3. The 15 matched log10(tau90/yr) values, second (total) asymmetric
errors and host distances were transcribed and visually checked from that table.
These photometric age summaries are not complete star-formation histories or
gravity-memory measurements. Attempts to retrieve the earlier 2014 full electronic
histories from CDS/mirrors failed; no missing time series were reconstructed.
Source article PDFs and screenshots are not redistributed.

Matched subsets: 91 SPARC galaxies / 1,347 radial rows; 15 dwarfs in two hosts
(11 MW, 4 M31). No new target observations, no homogeneous all-population history
catalog. The original 279-system benchmark is preserved. Missing predictors are
not imputed. Galaxy folds and dwarf host splits are explicit in predictions.
The SPARC distance bands are a spatial stress test, NOT a complete host catalog.

There are 2,314 fits, 49 validation-prediction tables, 100 field-uncertainty
realizations per scenario, and 500 age-uncertainty realizations. Results JSON
contains main-fit coefficients and metrics plus sensitivity-score arrays;
intermediate sensitivity-fit coefficients are reproducible from the script.
All numerical optimizer calls report success; no global-optimum claim.

History sensitivity uses independent standard normals, multiplies negative or
positive variates by the respective total log-age error, and truncates ages at
13.8 Gyr. These are descriptive sensitivity draws, not full correlated posterior
samples. Dynamical errors, distance errors and correlated stellar-model
systematics are not propagated. The primary 5 Gyr history gain survives only 37.2%
of draws; it is not robust. The 2 Gyr and 10 Gyr variants are sensitivities, not
selected timescales. No inference of gravitational memory is established.

verify.py checks all 49 prediction tables, equations independently implemented,
input hashes, unsplit groups, equal-galaxy/equal-host RMS and the 44 full-benchmark
satellites excluded by the two broad-model host fits. provenance.json records
source interpretation and numerical hashes. PROTOCOL.txt was written before
fits, after prior inspection of the old benchmark; it is not preregistration.
