29 October 2026. **Status: independent observational reconstruction completed.** The original Einasto and common-core candidates do not improve mean galaxy prediction. A subsequently tested shared RAR–NFW blend lowers mean held-out RMS by 2.04%, with uncertainty spanning no gain and slightly worse median speed error. Results below are exploratory and do not establish an all-scale gravity law. The main execution environment remains unavailable; the existing A43 PDF and portfolio ZIP were not changed.

# Completed observational comparison

**A small numerical gain was found only in an additional exploratory blend.** Its mean held-out logarithmic RMS is 2.04% lower than refitted RAR, but the uncertainty interval includes no improvement and median speed error is slightly worse. The two originally proposed profile changes do not improve the average prediction.

| Predictor | Shared coefficients | Mean held-out RMS, dex | Mean median speed error |
|---|---:|---:|---:|
| Fixed published RAR | 0 | 0.154556 | 8.996% |
| Refitted RAR | 1 | 0.153736 | 8.932% |
| Constrained NFW | 3 | 0.153947 | 9.467% |
| Cosmological Einasto | 3 | 0.160093 | 9.067% |
| Common core response | 4 | 0.154341 | 9.528% |
| Exploratory force blend | 5 | 0.150597 | 8.963% |

Each model predicts every retained radius of each excluded galaxy. RMS gives each galaxy equal total squared-log-error weight; the headline averages three complete fivefold repeats. Speed-error medians pool radial points within each repeat, then average the three medians. A 2.04% relative reduction in log RMS is not a two-percentage-point reduction in speed error.

## What changed in the execution

The main execution environment remained unavailable. A separate JavaScript calculation reconstructed the complete primary [SPARC catalog](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt) and [radial table](https://astroweb.case.edu/SPARC/MassModels_Lelli2016c.mrt) from research-tool text: 175 galaxies and 3,391 radii, with no missing lines or conflicting overlaps. Applying the paper’s stated selections recovered exactly 149 galaxies and 2,700 radii. Its all-data NFW coefficients reproduce the archived A, B and C to a maximum absolute difference of 1.31×10⁻⁷.

The archived Python fold assignments were inaccessible. These runs therefore use newly fixed Mulberry32/Fisher–Yates partitions, not the old assignments, even though their seed integers match. They are an independently implemented reconstruction on previously explored observations, not a new survey or independent discovery sample. Both baselines were refitted within every new training partition; no old fold predictions were substituted.

## One shared empirical blend

After the five initial models had been evaluated, one additional candidate was frozen and tested to examine whether RAR and NFW make complementary errors:

```math
g_{\rm blend}=(1-w)g_{\rm RAR}+w\left[g_b+\frac{GM_{\rm NFW}(<r)}{r^2}\right],\qquad 0\le w\le1.
```

RAR’s acceleration scale and NFW’s three mass-mapping coefficients are fitted separately inside each outer training set; one common weight is then fitted on those training predictions. The test galaxy supplies no rotation data to any fitted coefficient. This estimates five shared coefficients in total. Both endpoint weights are evaluated explicitly. No inner cross-fitting is used for the weight, and the candidate was chosen after examining the initial results, so its comparison remains exploratory.

The all-sample descriptive weight is w = 0.532088, with RAR a₀ = 1.036815×10⁻¹⁰ m s⁻² and the original all-sample NFW mapping. Cross-validation weights range from 0.3830 to 0.6885. The weight is not a dark-matter mass fraction or a probability that either theory is correct. All reported validation predictions use the corresponding training-only coefficients.

## Repeat and transfer results

| Predictor | Repeat 20261027 | Repeat 7331 | Repeat 48105 |
|---|---:|---:|---:|
| Fixed published RAR | 0.154556 | 0.154556 | 0.154556 |
| Refitted RAR | 0.155100 | 0.153339 | 0.152768 |
| Constrained NFW | 0.154485 | 0.154180 | 0.153178 |
| Cosmological Einasto | 0.161222 | 0.160486 | 0.158570 |
| Common core response | 0.155040 | 0.154595 | 0.153388 |
| Exploratory force blend | 0.151532 | 0.151055 | 0.149205 |

The blend reduces RMS relative to both fitted RAR and NFW in all three repeats. A paired 10,000-resample galaxy bootstrap gives blend minus RAR = −0.003139 dex, with a descriptive 95% interval [−0.009060, +0.002861]. Every resampled galaxy carries its predictions from all three repeats; the calculation does not treat 447 galaxy-repeat entries as independent. It excludes refitting and model-search uncertainty, so it does not establish superiority.

Four additional tests exclude 37 galaxies at each extreme of mass proxy or mass proxy/(2πR_e²). All entries below are held-out equal-galaxy RMS in dex:

| Predictor | Low mass | High mass | Low surface density | High surface density |
|---|---:|---:|---:|---:|
| Fixed published RAR | 0.200280 | 0.122613 | 0.200228 | 0.121202 |
| Refitted RAR | 0.197281 | 0.122522 | 0.196530 | 0.124710 |
| Constrained NFW | 0.200314 | 0.131606 | 0.191984 | 0.127752 |
| Cosmological Einasto | 0.205317 | 0.129811 | 0.194500 | 0.127393 |
| Common core response | 0.200428 | 0.131933 | 0.191984 | 0.127633 |
| Exploratory force blend | 0.194649 | 0.125834 | 0.189054 | 0.121555 |

The blend improves the RAR point estimate in three of these four tails but worsens the high-mass tail. Every blend-versus-RAR tail interval includes zero. Fixed published RAR also remains a competitive reference; no candidate dominates all metrics and subsets.

## What the failed profile repairs teach

The common core exponent fitted to all observations is n = 0.035683, close to the unmodified NFW limit. Across the 15 validation training folds, n ranges from 0 to 0.09761; three folds choose zero. Coring lowers all-data training RMS from 0.150692 to 0.150637 dex, yet increases held-out RMS in each complete repeat. That is evidence against promoting this particular common correction as an improvement.

Radial diagnostics locate part of Einasto’s shortfall. For r < 0.5R_e, covering 66 galaxies, its mean logarithmic residual is +0.13675 dex and RMS 0.28468, versus +0.04069 and 0.23986 for constrained NFW. The small core response reduces NFW’s central bias but does not improve the overall benchmark. These conditional residuals do not establish that all Einasto halos or all feedback models fail.

All 340 halo solution records, including 20 exact nested references and 320 optimizer runs, are preserved. Every optimizer run converged; feasible coordinate perturbations found no loss reduction exceeding 10⁻¹⁰. The core training loss never exceeded its nested NFW reference. Actual fitted halo boundary masses agree within 6.44×10⁻¹⁵ relative; 57 independent density-integral checks of fitted Einasto profiles agree with the gamma implementation within 3.16×10⁻¹⁰ relative.

## Physical and publication limits

The blend is one shared prediction formula, but it is an empirical ensemble rather than a derived gravitational field equation. For any fixed w < 1 and an isolated finite baryonic mass, its RAR term eventually falls as 1/r while a truncated NFW force falls as 1/r². It therefore does not preserve finite effective halo mass if reinterpreted as Newtonian dark matter. No new cluster, lensing, Solar-System or cosmological likelihood has been fitted in this reconstruction, and no gravitational-memory signal has been detected.

Publication value remains in a reproducible comparison of restricted shared predictors and their failures, not a demonstrated advance over MOND or a new universal theory. Existing work already investigates shared galaxy–halo models, broad formula searches and the difficulty of tightening the SPARC RAR: [Boreiko et al. 2026](https://arxiv.org/html/2601.07799v1), [Desmond, Bartlett and Ferreira 2023](https://arxiv.org/abs/2301.04368), and [Stiskalek and Desmond 2023](https://arxiv.org/abs/2305.19978). A stronger submission needs a clearly distinct question, uncertainty treatment and validation on genuinely new observations.

The next useful test is to freeze this blend and the original candidates, then evaluate independent observations and correlated distance, inclination, stellar-mass and concentration uncertainties. If pursuing a physical response model, independently measured assembly history or environment should supply new information; another freely adjusted curve on these same held-out results would weaken the evidence.

[Executable calculation, selected data and complete fit records](https://chatgpt.com/space/page_88a2001f480c8191ab9e086fd9fe872f). Saved code and records were read back exactly. PDF/ZIP integration and a rendered document preview remain pending.

# Starting point

The current [A43 paper](https://chatgpt.com/library?file_id=libfile_34fdf19e4f048191bbb268e4ecd18595&account_id=personal), version 6, Section 4.2, contains the verified baseline. Its 149 SPARC galaxies and 2,700 radial measurements use independently measured stellar luminosity, HI mass and effective radius. Three shared coefficients predict a reference halo mass; a published concentration relation fixes its shape. Across the three recorded galaxy-validation splits, the earlier mean RMS was 0.153234 dex for constrained NFW and 0.153093 dex for fitted RAR. These are previous results, not new runs.

Keep the existing mass proxy and reference-mass mapping, with masses in solar units and effective radius in kpc:

```math
M_p=0.5L_{3.6}+1.33M_{\rm HI},\qquad
\frac{M_{\rm ref}}{M_\odot}=\frac{M_p}{f_bM_\odot}+10^A\left(\frac{M_p}{10^{10}M_\odot}\right)^B\left(\frac{R_e}{3\,\mathrm{kpc}}\right)^C.
```

Here f_b = 0.0490/0.3175, h = 0.671, and H₀ = 67.1 km s⁻¹ Mpc⁻¹. Reference radius R = [G M_ref/(100 H₀²)]^(1/3), and dark mass M_h = (1−f_b)M_ref. This preserves the earlier dark-matter-only-equivalent reference convention; R need not equal the assembled baryon-plus-dark-matter system’s true overdensity radius. Stellar mass-to-light ratios, molecular-gas omission, distances and inclinations retain the earlier benchmark assumptions.

# Candidate 1 Cosmological Einasto

Replace the NFW density with Einasto, using its own concentration and shape relations. Dutton and Macciò’s redshift-evolution fit gives the following z = 0 prescription; this protocol chooses it before fitting, instead of switching between their evolution fit and direct z = 0 table. Their published relations are fixed, so only A, B and C are fitted. [Primary source, Sections 6 and Appendix A](https://arxiv.org/html/1402.7073v2).

```math
\rho_E(r)=\rho_{-2}\exp\!\left[-\frac{2}{a_E}\left(\left(\frac{r}{r_{-2}}\right)^{a_E}-1\right)\right],\qquad r_{-2}=R/c_E.
```

```math
m=\log_{10}(M_{\rm ref}h/10^{12}M_\odot),\qquad
\log_{10}c_E=0.977-0.130m,
\quad a_E=0.155+0.0095\nu^2,
\quad\log_{10}\nu=-0.11+0.146m+0.0138m^2+0.00123m^3.
```

Normalize the dark mass exactly at R. Writing P for the regularized lower incomplete gamma function, the resulting enclosed mass is:

```math
M_E(<r)=M_h\frac{P\!\left(3/a_E,\,(2/a_E)[\min(r,R)/r_{-2}]^{a_E}\right)}{P\!\left(3/a_E,\,(2/a_E)c_E^{a_E}\right)},\qquad
g_{\rm pred}=g_b+\frac{GM_E(<r)}{r^2}.
```

Density is positive and finite at the center. The adopted truncation retains constant mass beyond R. Concentration scatter and baryonic halo response remain unmodeled in this candidate; low-mass extrapolations must be flagged. In-memory arithmetic checked 126 combinations of shape, concentration and radius. Incomplete-gamma masses agreed with independent direct density integration to a maximum relative difference of 4.79×10⁻⁹. This verifies implementation mathematics, not fit quality.

# Candidate 2 One common core response

Read, Agertz and Collins provide a coreNFW family in which feedback history can flatten an NFW cusp. Their calibrated core size is proportional to projected stellar half-mass radius. Here R_e approximates that radius, and a single common exponent n is fitted for the whole training sample. This is an empirical profile-response test, not a measured-history fit. [Primary source, equations 16–21](https://arxiv.org/html/1508.04143v3).

Let M₀(<r) be the existing finite NFW dark halo, set r_c = 1.75 R_e, and define:

```math
M_n(<r)=\begin{cases}
M_0(<r)\left[\dfrac{\tanh(r/r_c)}{\tanh(R/r_c)}\right]^n,&0<r<R,\\
M_h,&r\ge R,
\end{cases}\qquad 0\le n\le1.
```

All four coefficients A, B, C and n are shared. The finite-boundary denominator is an explicit normalization choice added to the published family. With f = tanh(r/r_c), differentiation gives:

```math
\rho_n(r)=\frac{f^n\rho_0(r)+n f^{n-1}(1-f^2)M_0(<r)/(4\pi r^2r_c)}{\tanh^n(R/r_c)}.
```

This establishes positive density, exact recovery of NFW at n = 0, an inner slope ρ ∝ r^(n−1), and exact conservation of total dark mass. Exterior force is unchanged. For fixed halo mass and concentration, the interior force can only decrease. Consequently this correction can address overprediction, while underprediction requires another change. Density positivity alone does not prove an admissible equilibrium distribution function.

Two independently written in-memory checks covered 21,780 and 30,060 dimensionless test points. Both verified the normalization, positive density and absence of enclosed-mass excess. Analytic mass derivatives agreed with centered finite differences within 1.10×10⁻⁹ relative. For a complete core and a distant boundary, halo-force multipliers at r/r_c = 0.1, 1 and 2 are 0.099668, 0.761594 and 0.964028. These are halo-only mathematical changes, not observed error reductions.

# A diagnostic before optimizing

For the recorded equal-galaxy objective L = Σ wᵢ eᵢ², define eᵢ = log₁₀(gᵢ/g_obs,i) and hᵢ as the baseline halo acceleration. The derivative at n = 0 is:

```math
\left.\frac{dL}{dn}\right|_0=\frac{2}{\ln10}\sum_i w_i e_i\frac{h_i}{g_i}\ln\!\left[\frac{\tanh(r_i/r_{c,i})}{\tanh(R_i/r_{c,i})}\right].
```

A negative value identifies an initially improving direction on that training fold. A positive value means infinitesimal coring initially worsens it at the baseline coefficients. Neither sign determines global optimization or held-out performance. Calculate and preserve this diagnostic using training rows only, and report the predeclared core candidate regardless of its sign.

# Archived intended comparison

Use the archived inputs from [A43 Code and Data](https://chatgpt.com/library?file_id=libfile_ad3b757ae6ac81919116d5024496b395&account_id=personal). Preserve their hashes, selections, weights and original baryonic forces.

The six items below preserve the originally intended protocol. The completed reconstruction above uses separately defined JavaScript partitions because the archived assignments could not be loaded. It does not claim to have executed those exact old folds.

1. Compare fixed RAR, fitted RAR, the existing three-coefficient NFW predictor, three-coefficient Einasto, and four-coefficient core response. Reuse verified frozen baseline predictions only after row and split identity checks.

2. Retain the three fivefold whole-galaxy partitions with seeds 20261027, 7331 and 48105. Every radius of a test galaxy stays excluded from fitting. Retain the four original predictor-only mass/surface-density tail tests of 37 galaxies each.

3. Refit A, B and C inside every training split. Bounds remain A ∈ [8,14], B ∈ [0.1,2], C ∈ [−1.5,1.5]; the core response adds n ∈ [0,1]. Use the previous three deterministic mass-mapping starts. For the core family cross these starts with n = 0, 0.5 and 1, and additionally evaluate the fitted baseline at n = 0 as an exact nested reference. Keep all starts, failures, coefficients and boundary hits.

4. Report equal-galaxy logarithmic RMS, median speed error, every split, their mean, all four tail scores, and radial residual patterns. Preserve paired galaxy uncertainty comparisons; descriptive resampling of fixed predictions excludes refitting and selection uncertainty.

5. Do not select a favorable split or metric as evidence of superiority. Report both new candidates. Choosing a family after this comparison remains exploratory; a prospective new sample is required for stronger confirmation.

6. If either candidate improves the galaxy benchmark, freeze its rules before attempting dwarf, early-type and cluster transfer. Keep those results separate. Dwarf Jeans profiles and cluster lensing need the appropriate source geometry, nuisance treatment and observables; a galaxy score alone cannot establish an all-scale law.

# History and the broader search

The published coreNFW relation uses n = tanh(0.04 t_SF/t_dyn), where t_SF is a star-formation duration. A quenching lookback time t₉₀ does not determine that duration. A genuine history test therefore needs independently reconstructed histories; fitting one common n does not measure gravitational memory. Ordinary baryonic-potential evolution is already a physical way for history to alter collisionless halos.

DC14 remains a worthwhile subsequent comparison. Its original halo-mass, virial-overdensity and concentration conventions require conversion before attaching it to the present reference-M200 model. In particular, its profile scale radius is not generally the slope-minus-two radius. The published SPARC comparison by Li et al. already includes DC14, Einasto and coreNFW with galaxy-specific nuisance inference. Its better profile fits motivate these tests but do not establish better prediction of entire excluded galaxies. [DC14 primary paper](https://arxiv.org/abs/1404.5959); [SPARC halo comparison](https://arxiv.org/html/2001.10538v2).

# Reproducible arithmetic check

These synthetic checks use no astronomical targets. The root core grid used c = 3, 5, 8, 12, 20, 40; r_c/R = 0.0001, 0.001, 0.01, 0.03, 0.1, 0.3; n = 0, 0.25, 0.5, 0.75, 1; and 121 logarithmic radii from 10⁻⁷R to 10⁻⁰·⁰¹R. Centered derivatives used multiplicative steps exp(±10⁻⁵). The independent grid extended r_c/R to 1 and used different concentration/exponent samples.

The Einasto check used a_E = 0.1, 0.155, 0.2, 0.3, 0.5, 1, 2; c = 3, 10, 40; and r/R = 10⁻⁵, 0.001, 0.01, 0.1, 0.5, 1. One implementation evaluated regularized incomplete gamma by series/continued fraction. The other integrated density directly with 8,192 Simpson intervals in logarithmic radius, using s = x exp(v), v ∈ [−80,0]. All normalized masses were positive and at most the boundary mass. These were the preliminary mathematical checks. The completed observational fits and fitted-profile checks are reported at the top of this updated document.
