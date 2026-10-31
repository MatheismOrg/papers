# Published concentration–mass prior: a restricted NFW prediction test

29 October 2026. Exploratory analysis of 2,700 retained SPARC rotation-curve
points in 149 galaxies. The two model families and all validation splits were
fixed in `PROTOCOL.txt` before fitting this branch. Earlier analyses of this
sample already existed, so this is not a prospectively blinded confirmation.

**Result:** imposing a published cosmological concentration–mass relation makes
a compact, three-coefficient NFW predictor approximately as accurate as the
one-coefficient fitted RAR on unseen galaxies. It does **not** establish an
improvement over RAR/MOND, a literature-best model, or an all-scale law.

## What was tested

The shared photometric proxy is

\[
M_p/M_\odot=10^9(0.5L_{3.6}+1.33M_{\rm HI}),
\]

with the catalog luminosity and HI mass in their stated billion-unit columns.
It is an approximate cold-baryonic mass proxy; fixed total stellar M/L and
omitted molecular gas matter. Reff is an independently measured effective
radius. Neither feature uses Vflat or the target rotation velocities.

The two predeclared reference-mass mappings are

\[
M_{\rm ref}=M_p/f_b+10^A M_\odot
       (M_p/10^{10}M_\odot)^B
       (R_{\rm eff}/3\,\mathrm{kpc})^C.
\]

The mass-only model fixes C=0; the size model fits C. Bounds were A∈[8,14],
B∈[0.1,2], C∈[−1.5,1.5]. All selected optima were interior. The first term
imposes a proxy baryon-budget restriction. These mappings are empirical and
are not a published abundance-matching prescription.

For halo structure we use the NFW z=0 relation of Dutton & Macciò (2014),
their Eq.8 and Planck cosmology:

\[
\log_{10}c_{200}=0.905-0.101\log_{10}(M_{\rm ref}h/10^{12}M_\odot),
\quad h=0.671,
\quad f_b=0.0490/0.3175.
\]

The units are physical solar masses. Mref is a **dark-matter-only-equivalent
reference total-matter mass**, with r200=[GMref/(100H0²)]^(1/3), rs=r200/c200.
We assign dark mass (1−fb)Mref and add the supplied observed baryonic field.
This is an explicitly stated no-response approximation: the resulting total
system need not have its actual 200-critical radius at reference r200.

Writing f(z)=ln(1+z)−z/(1+z), the prediction is

\[
g_{\rm pred}(r)=g_b(r)+
 {G(1-f_b)M_{\rm ref}\over r^2}
 {f[\min(r/r_s,c_{200})]\over f(c_{200})}.
\]

The spherical halo is truncated at reference r200, guaranteeing finite dark
mass. No observed row in any fitted or validation split reached that edge;
truncation did not create the reported accuracy. The halo contribution is
positive and vanishes as Mp→0 at fixed nonzero radius and fixed Reff. It retains
the conventional idealized central NFW density cusp.

## Held-out results

Entries below are equal-galaxy RMS log10 acceleration errors in dex; lower is
better. No parameter was fitted separately to a held-out galaxy. Each galaxy's
radial points stay together. All four feature-tail tests hold out 37 galaxies.

| Model | Shared fitted coefficients | CV 20261027 | CV 7331 | CV 48105 |
|---|---:|---:|---:|---:|
| Fitted RAR | 1 | 0.152898 | 0.153013 | 0.153367 |
| Shared empirical NFW, mass | 4 | 0.156294 | 0.157426 | 0.160311 |
| Shared empirical NFW, mass + size | 6 | 0.155601 | 0.154317 | 0.158489 |
| NFW + published c(M), mass | 2 | 0.154580 | 0.155359 | 0.157385 |
| NFW + published c(M), mass + size | 3 | 0.152357 | 0.152790 | 0.154556 |

| Model | Lowest mass | Highest mass | Lowest surface proxy | Highest surface proxy |
|---|---:|---:|---:|---:|
| Fitted RAR | 0.197281 | 0.122522 | 0.196530 | 0.124710 |
| Shared empirical NFW, mass + size | 0.199438 | 0.140581 | 0.187161 | 0.146192 |
| NFW + published c(M), mass | 0.207471 | 0.121022 | 0.205450 | 0.129349 |
| NFW + published c(M), mass + size | 0.200314 | 0.131606 | 0.191984 | 0.127752 |

The published-c(M) size model improves on the six-coefficient shared NFW in all
three ordinary CV splits while using three coefficients. This supports testing
physical restrictions, but the small differences are not a broad superiority
claim: its paired 95% galaxy-bootstrap intervals against shared NFW cross zero
in two of the three splits. Tail behavior is mixed.

Against RAR, the size-model primary difference is −0.000541 dex, with a
descriptive 95% paired interval [−0.011441,+0.010199] dex. It improves by only
0.000223 dex in the second seed and loses by 0.001189 dex in the third. All three
intervals cross zero. Its primary median circular-speed-equivalent absolute
fractional error is **9.294%**, compared with **8.893% for RAR**. A tiny RMS
advantage in one split therefore should not be advertised as “beating MOND.”

The 5,000-replicate paired Poisson galaxy bootstrap is descriptive. It holds
predictions fixed and does not refit, account for model search, or marginalize
measurement and astrophysical systematics. The three CV splits reuse the same
149 galaxies and are not three independent experiments.

## Coefficients and physical checks

All-data fitted coefficients, for descriptive reproduction only:

| Model | A | B | C |
|---|---:|---:|---:|
| NFW + published c(M), mass | 11.426411792745512 | 0.8094059079901207 | 0 (fixed) |
| NFW + published c(M), mass + size | 11.4672833579147 | 0.9937487126154467 | −0.6355007557861937 |

The size relation assigns a larger reference halo mass to more compact systems
at fixed Mp. This is a fitted conditional trend, not evidence for gravitational
memory or a causal formation mechanism.

In the size all-data fit Mref spans 3.23×10^9–5.42×10^12 Msun, concentrations
7.05–14.93, and rs=2.10–52.86 kpc. The approximate proxy-retention fraction
Mp/(fb Mref) spans 0.043–0.356. Eight galaxies (26 rows) lie below the published
c(M) calibration's approximate 10^10 h^-1 Msun lower boundary. The three CV
splits extrapolate for 8/9/9 galaxies (26/36/36 rows). This is explicitly flagged
in every prediction table; masses were not clipped. The greatest observed
r/r200 among all size-model validation predictions is about 0.378.

Both candidates together used 40 fit sets and 120 deterministic optimizer
attempts, with zero reported failures or active parameter bounds. Results
include all starts, coefficients, training row IDs, and physical diagnostics.

## Scope and source conventions

The Dutton–Macciò relation is a published structural prior, calibrated on
relaxed collisionless simulations. Its intrinsic concentration scatter is
approximately 0.11 dex; this deterministic baseline omits that scatter,
halo-response physics, assembly and environment. Primary source:
[Dutton & Macciò 2014, arXiv:1402.7073](https://arxiv.org/pdf/1402.7073).

Moster, Naab & White (2013) were checked directly. Their paper calls its masses
“virial” but explicitly defines them at 200 times critical density; it uses
WMAP7, and satellites require infall mass/redshift. We did not silently splice
its frozen abundance-matching coefficients into the different Planck cosmology
or substitute present masses for infall masses. Primary source:
[Moster, Naab & White 2013, arXiv:1205.5807](https://arxiv.org/pdf/1205.5807).

Neither these candidates nor the empirical shared NFW comparators represent a
full ΛCDM prediction. No fitted individual halo scatter, feedback, observational
nuisance marginalization, satellite stripping, lensing, cluster, CMB, or growth
likelihood is included. The mass proxy contains cold gas; applying its fitted
mapping unchanged to hot-gas-dominated clusters would be unjustified. Failure
of these deterministic mappings would not reject ΛCDM. Success here would not
establish a unified gravity theory or novelty relative to all existing work.

## Files and reproduction

- `PROTOCOL.txt`: exact pre-fit families, starts, bounds, and validation rules.
- `run.py`: executable analysis; requires numpy and scipy, with the sibling
  `dm_shared_halo/results` comparator directory in the full research package.
- `input.csv`, `galaxy_features.csv`: frozen inputs copied exactly from that
  independently extracted photometry branch.
- `results/results.json`: complete coefficients, fit logs, metrics, tails,
  physical checks, comparator hashes and paired bootstrap intervals.
- `results/metrics_summary.csv`: all nine comparator/candidate models over
  all seven held-out designs (63 rows).
- `results/NFW_Planck_*.csv`: 16 complete new prediction tables, eight/model.
- `comparators/`: exact copies of the 49 source held-out comparator tables.
- `manifest.json`: file checksums; source/protocol/script hashes are also in
  the immutable numerical result JSON.

Reproduce from the full package with `python dm_cosmological_halo/run.py`.
Input rows, group labels, folds, radii, accelerations and features were checked
against every imported comparator table before comparisons were computed.
