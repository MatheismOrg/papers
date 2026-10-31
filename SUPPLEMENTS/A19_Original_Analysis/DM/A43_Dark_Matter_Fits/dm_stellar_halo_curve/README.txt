A43 stellar-to-halo mapping continuation -- 29 October 2026

REPRODUCE
Requires the complete A43 archive and numpy/scipy/pandas:
  OPENBLAS_NUM_THREADS=1 python3 dm_stellar_halo_curve/run.py
  OPENBLAS_NUM_THREADS=1 python3 dm_stellar_halo_curve/ensemble_audit.py
  python3 dm_stellar_halo_curve/dwarf_sensitivity.py
The code reuses the measured-input loader and NFW constants from the sibling
population-scatter module. No new astronomical download is needed.

DESIGN
This branch was proposed after inspecting the population-scatter transfer results.
PROTOCOL.txt preserves that chronology and the four locally fixed families. The
stellar blend was added after the four-family outcomes, using the previously
established weight-fitting algorithm. These are adaptive analyses of reused data.

Keep the old finite NFW profile, c(M), baryon floor and exact original splits.
Use stellar mass in the fitted power-law mapping instead of total baryonic proxy
mass. Gas remains in gbar and the mass floor. Published and curved stellar-to-halo
relations are controls. Inverting Moster's median relation is not a calculation
of p(Mhalo|Mstar): selection, scatter and the halo mass function matter. Its quoted
z=0 coefficients are transferred to the shared Planck reference convention here,
not obtained from a rerun of abundance matching or from a satellite model.

KEY RESULTS
Model                       Disk CV RMS   Frozen 44-dwarf RMS
Old baryonic NFW             0.153234      0.667569
Published inverse Moster     0.171114      0.555907
Scaled inverse Moster        0.171280      0.566329
Stellar power NFW (3 coeffs)  0.153292      0.479136
Curved stellar NFW (5)       0.153805      0.478741
Old RAR/NFW blend (5)        0.150602      0.666126
RAR/stellar-NFW blend (5)     0.150506      0.505205
RMS is log10 acceleration, dex. Disk score averages three complete held-out
partitions of 149 galaxies/2700 rows. Dwarf coefficients are frozen from full
SPARC, not refitted to dwarfs. No stellar/halo state is fitted to a test galaxy.

The three-coefficient stellar relation improves frozen dwarf RMS 28.2% while
leaving average disk RMS essentially unchanged. It improves 30/44 individual
log errors. Conditional dwarf resampling gives a difference of -0.18843 dex,
interval [-0.28081,-0.08311], vs old NFW. This interval omits correlated host
systematics, measurement reconstruction and adaptive search uncertainty.

FAILURES AND SENSITIVITIES
Classical dwarfs worsen, ultrafaints improve. Milky Way improves; Andromeda worsens
for the pure stellar rule. Conservative host matching includes G/L and A/L and
accounts for all 44 dwarfs (26 Milky Way,18 Andromeda). Three of four photometric
disk tail scores worsen versus old NFW. The old 25-object high-precision selection
(max dispersion error/sigma<=0.25, no source caution flag) retains the aggregate
gain, 0.57502 -> 0.33780 dex. It is an earlier sensitivity, not new validation.

The stellar blend improves dwarf RMS 24.2% over the old blend, with disk accuracy
approximately unchanged. Its median disk speed error is 8.78%. The tiny disk RMS
difference versus the old blend has an interval spanning zero. It improves both
host subsets and narrowly improves classical-dwarf RMS vs old pure NFW, but does
not beat disk-trained RAR on the classical subset. It lacks a finite dark-halo
mass asymptote whenever its RAR weight is nonzero and has no full field law.

AUDIT AND REPRODUCIBILITY
All 60 numerical base fits,20 ensemble weights and20 fixed Moster evaluations
are retained; 200 numerical multistart runs and40 exact endpoints (240 records).
No numerical optimizer run failed. The curved transition hits its upper bound in
12/20 cases; additional curvature does not help disk prediction. Independent
scripts recompute scores, solve the published mass identity, check exact nesting,
positive finite mass, baryon floor and galaxy exclusions. The tables preserve all
fitted coefficients, starts, predictions and conditional intervals.

These are conditional mass-estimator transfer improvements, not a solution of
all dark matter, complete Jeans inference, or proof of new physics/memory.
Primary sources: Moster, Naab & White 2013, https://doi.org/10.1093/mnras/sts261;
Li et al. 2020, https://doi.org/10.3847/1538-4365/ab700e; Dutton & Maccio 2014;
SPARC and the dwarf catalog/Wolf estimator credited in the parent A43 paper.
