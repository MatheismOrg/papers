# Dark Matter Across Galaxies: Better Dwarf Predictions and Limits on Fixed-Kernel Gravity

Working research paper, revised 29 October 2026. Portfolio paper A43.

## Abstract

A shared stellar-mass halo relation improves transfer from rotating galaxies to dwarf spheroidals, while a mass-budget bound exposes failures of published fixed-kernel gravity prescriptions. Trained on 149 SPARC galaxies, the three-coefficient stellar relation lowers frozen 44-dwarf log-acceleration RMS from 0.66757 to 0.47914 dex, a 28.2% reduction, with essentially unchanged galaxy-held-out RMS. The gain concentrates in ultrafaint dwarfs; classical dwarfs worsen. Independently, positivity and finite kernel normalization imply g/g_bar <= 1+2 alpha at the baryonic half-mass radius of an isolated spherical source. Published nonlocal-gravity and constant-parameter MOG normalizations fall below central inferred discrepancies in 31 and 33 of 44 dwarfs; changing their range cannot repair this budget deficit. A wider 286-system benchmark identifies further failures of transfer between populations. These results establish a conditional predictive improvement and a necessary constraint on fixed linear kernels. Adaptive model selection and uncertain satellite mass estimates require independent validation and uncertainty-aware dynamical tests.

Keywords: dark matter; radial acceleration relation; dwarf galaxies; gravitational lensing; galaxy clusters; cross-validation; gravitational memory

## 1. What is being tested, and what already exists

Two results organize this paper. First, changing the mass variable in a compact halo mapping improves dwarf predictions without adding coefficients or materially degrading average disk predictions (Section 4.6). Second, a positive kernel with a finite integral imposes a mass ceiling that excludes the central dwarf estimates from the tested published parameterizations, regardless of kernel range (Section 8.4). The first result identifies a better predictor in this benchmark; the second identifies a structural limitation of a class of proposed explanations.

There is no absence of smooth formulas for dark-matter phenomenology. Many galaxy rotation curves, halo profiles and statistical relations have compact descriptions. The harder problem is predictive transfer: one parameter set must describe systems with different sizes, baryonic distributions, environments and histories, and ultimately predict different kinds of measurements. A flexible curve can interpolate a finite dataset without supplying that physical explanation.

McGaugh, Lelli and Schombert (2016) established a tight radial acceleration relation in rotating galaxies. We use their empirical form

$$g_{\mathrm{RAR}}=\frac{g_{\mathrm{bar}}}{1-\exp[-\sqrt{g_{\mathrm{bar}}/a_0}]}.$$

At high baryonic acceleration it approaches the Newtonian baryonic value; at low acceleration it approaches sqrt(a0 times g_bar). Its success is prior knowledge, not a discovery of this paper. Navarro, Frenk and White (1997) likewise provide a common halo-profile family with varying characteristic scales. Comparing a fit only against baryons alone is therefore a substantially weaker test than comparing it with existing dark-matter or modified-gravity models.

Formula searches and additional state variables also have precedent. Desmond, Bartlett and Ferreira (2023) applied exhaustive symbolic regression to the RAR. Zhao and Famaey (2012) investigated acceleration and potential dependence to connect galaxies and clusters. Tian et al. (2020, 2024) measured distinct cluster and brightest-cluster-galaxy relations. Rahvar and Mashhoon (2014) explored nonlocal gravity. A recent broad comparison by Dantas, Ribeiro and Capelato (2026) considers a two-component virial approach across early-type galaxies, brightest cluster galaxies, clusters and weak lensing. Inspection of its full paper shows that its constant-interaction approximation varies across systems and radial apertures, while its lensing treatment adds a halo-based interaction model. That broad comparison is not one fixed coefficient set across every population; its quantitative pipeline is not reproduced here.

The contribution here is a reproducible empirical comparison, including failed transfers and an audit of a misleading radius fit. The particular coefficients below are fitted in this analysis; no priority claim is made for the mathematical families, potential dependence, cross-scale unification, or gravitational memory.

### 1.1 What this adds, and what stronger published results mean

The value of this study is the common, population-balanced evaluation of compact shared formulas, with whole-system and whole-population exclusions, transparent failure cases and reproducible numerical records. This is a useful exploratory benchmark, not a claim to be the first cross-scale study or the best explanation of dark matter. A lower score against one recalibrated RAR baseline establishes a narrower result than superiority to every published model.

Li et al. (2018) report 0.057 dex residual scatter after fitting individual SPARC galaxies while marginalizing over stellar mass-to-light ratios, distances and inclinations. Desmond (2023) infers intrinsic scatter of about 0.034 dex with a joint nuisance-parameter analysis. Neither number can be ranked directly against our balanced held-out RMS: they concern different samples, statistical targets and fitted quantities. Conversely, the 16.1% median disk speed error of our shared power sum is worse than the 9.2% of the galaxy-only RAR in our own comparison. A broader domain has not produced a better disk-galaxy fit.

Potential-dependent unification is already present in EMOND (Zhao and Famaey 2012). Symbolic equation searches have precedent in Desmond, Bartlett and Ferreira (2023), and relativistic constructions with calculated linear CMB and matter spectra in Skordis and Zlosnik (2021). This paper supplies neither their physical derivations nor comparable cosmological tests. Its specific coefficients and test records are contributions of this analysis; priority and literature-wide superiority remain unestablished.

### 1.2 What must be beaten at each scale

There is no single meaningful ranking of rotation-speed residuals, cluster mass profiles, planetary tracking and CMB spectra. The observation, allowed source information, nuisance parameters and prediction task must match before comparing scores. A43 has not demonstrated superiority over the strongest published models.

| Domain | Established comparison | What a stronger result must predict |
| --- | --- | --- |
| Solar System and strong gravity | GR with ordinary sources | Orbital/relativistic observables and external-field effects |
| Disk galaxies | Shared MOND/RAR; NFW, Burkert, DC14 and Einasto halos | New galaxies separately from new radii in measured galaxies |
| Dwarf spheroidals | Jeans/phase-space halo models and MOND with host fields | Stellar dispersions with anisotropy, tides, binaries and membership |
| Clusters and lensing | NFW/Einasto plus baryons; MOND plus residual mass | Profiles and projected observables across independent sources |
| CMB and large-scale structure | GR with cold dark matter, including LambdaCDM | Expansion, perturbations, CMB and growth likelihoods |

Li et al. (2020) find that cored Burkert/DC14 profiles generally fit SPARC better than NFW after allowing halo and observational nuisance parameters. Eckert et al. (2022) find average NFW cluster-profile deviations below 10%, with remaining shape diversity. These are fitted profiles, not parameter-free predictions from local baryonic acceleration. ACT DR6 supplies a current quantitative cosmological benchmark (Louis et al., 2025). Relativistic MOND completions such as Skordis and Zlosnik (2021) already calculate CMB and linear matter spectra; the algebraic RAR alone does not. Sections 4.1, 5.7 and 6.1 now test concrete, restricted comparisons and repairs. They do not reproduce the full published likelihoods.

## 2. Coverage from stellar systems to cosmology

The main regression uses five populations with independently specified physical radii, covering 27.1 pc to 600 kpc and baryonic accelerations from 5.18 times 10^-15 to 8.13 times 10^-9 m/s^2. Its 2,910 rows represent 279 physical systems, grouped conservatively into 278 validation groups. They are not 2,910 independent experiments. Other observables are tested separately rather than converted into artificial local acceleration points.

| Main fitted population | Systems | Rows | Quantity and principal assumptions |
|---|---:|---:|---|
| SPARC rotating galaxies | 149 | 2,700 | Resolved circular speeds; fixed distance, inclination and stellar mass-to-light ratios |
| Dwarf spheroidals | 44 | 44 | Half-light mass estimates from stellar dispersions; spherical equilibrium |
| Early-type galaxies | 16 | 32 | Inner stellar Jeans models and outer neutral-hydrogen rotation |
| MaNGA brightest cluster galaxies | 50 | 50 | Published isotropic Jeans/Abel estimates |
| CLASH lensing clusters | 20 | 84 | Central and outer NFW mass reconstructions from lensing |

### 2.1 Rotating and pressure-supported galaxies

The public SPARC archive (Lelli, McGaugh and Schombert, 2016) contains 175 galaxies and 3,391 radial measurements. We select quality flag Q below 3, inclination at least 30 degrees, positive radius and observed speed, quoted relative speed error at most 0.1, and positive net inward baryonic acceleration. Four of the 153 catalog-selected galaxies have no surviving radial measurements, leaving 149 galaxies and 2,700 rows. This differs from the original RAR paper's selection.

Disk and bulge mass-to-light ratios remain fixed at 0.5 and 0.7 in solar units at 3.6 micrometers, with catalog distances and inclinations. Signed squared components preserve the archive convention:

$$V_{\mathrm{bar}}^2=V_{\mathrm{gas}}|V_{\mathrm{gas}}|+0.5V_{\mathrm{disk}}|V_{\mathrm{disk}}|+0.7V_{\mathrm{bulge}}|V_{\mathrm{bulge}}|.$$

$$g_{\mathrm{bar}}=V_{\mathrm{bar}}^2/r,\qquad g_{\mathrm{obs}}=V_{\mathrm{obs}}^2/r.$$

For dwarf spheroidals, we use McConnachie's author-maintained January 2021 nearby-galaxy catalog and morphology information. Selection leaves 25 classical and 19 ultrafaint systems; the operational division is V-band luminosity below 10^5 solar luminosities. We calculate approximate Wolf et al. (2010) half-light estimates, rather than treating them as published Wolf masses:

$$r_{1/2}=\frac{4}{3}R_e,\qquad g_{\mathrm{obs}}\simeq\frac{3\sigma_{\mathrm{los}}^2}{r_{1/2}},\qquad g_{\mathrm{bar}}=\frac{G(M_\star/2)}{r_{1/2}^2}.$$

The circularized projected radius is obtained from angular major-axis size, distance and ellipticity. Nominal stellar mass is twice the V-band luminosity in solar units; ratios of one and three are supplied as sensitivities. We exclude unresolved dispersions, selected tidal cautions and measured gas fractions above 10% of nominal stellar mass; missing gas measurements are allowed for the selected morphologies. The frozen selection and exclusion reasons are supplied. Binary stars, tides, flattening, dispersion sampling and heterogeneous catalog uncertainties remain important, especially for ultrafaint systems.

The 16 early-type galaxies use Lelli et al. (2017): one inner Jeans anisotropic-model estimate and one outer HI estimate per object. Baryonic acceleration uses the photometric mass-to-light ratio 0.8 at 3.6 micrometers, not a dynamical mass-to-light ratio fitted from the target. The inner target still carries the source's dynamical-model assumptions. Radii span 0.17 to 28.61 kpc. These galaxies do not overlap the selected SPARC systems.

### 2.2 Cluster centers and lensing clusters

The 50 MaNGA brightest cluster galaxies use the numerical table of Tian et al. (2024), including the observed outer-to-effective-radius ratio and photometric Sersic masses. Physical radii span 3.87 to 24.78 kpc. Their accelerations are inferred with isotropic spherical Jeans/Abel modeling, not measured force vectors. Two nearby objects, 8989-12703 and 8991-6102, are conservatively kept in the same validation group; common host identity is not established. These 50 objects are not asserted to represent 50 independent cluster environments.

CLASH provides 84 central and outer estimates from 20 hosts (Tian et al., 2020): 20 central points and 64 points at fixed 100, 200, 400 or 600 kpc. The total masses are spherical NFW reconstructions using strong and weak lensing plus magnification, with the source's cosmology and gravity assumptions. Baryonic masses retain the published X-ray and stellar modeling. All radii from one host remain in one validation group. CLASH and MaNGA redshift ranges do not overlap; the selected X-COP hosts are also distinct from CLASH.

### 2.3 Observables retained outside the joint regression

Seven X-COP clusters contribute 21 fixed-radius hydrostatic estimates at 500, 750 and 1,000 kpc as a separate stress test. The 55 selected galaxy groups contribute global velocity dispersions. KiDS/GAMA supplies stacked weak-lensing acceleration bins with covariance but no single physical radius per bin. A frozen globular-cluster catalog supplies 151 conservatively selected stellar-dynamical controls. Solar-System, wide-binary and cosmological studies supply additional constraints and requirements; their raw likelihoods are not fitted here.

This covers a much wider range than ordinary disk galaxies and cluster outskirts alone, but it is still a selected benchmark. Merger maps, all satellite populations, gas-poor intermediate systems and full cosmological likelihoods are not exhaustively analyzed. Results from different measurement pipelines can differ because of both astrophysics and inference assumptions.

## 3. Metrics, weighting and validation

Let a_ref = 10^-10 m/s^2 and define

$$x=\log_{10}(g_{\mathrm{bar}}/a_{\mathrm{ref}}),\quad u=\log_{10}(r/\mathrm{kpc}),\quad y=\log_{10}(g_{\mathrm{obs}}/a_{\mathrm{ref}}).$$

The broad fit gives each of the five populations equal total weight and each conservative system group equal weight within its population; radial samples divide that group's weight. Thus 2,700 SPARC measurements cannot overwhelm 44 dwarfs. The aggregate score B is the square root of the mean, over populations, of their group-averaged mean squared log10 acceleration residuals. It is measured in dex, not percent. Equal population weighting is a declared modeling choice, not a cosmic abundance weighting.

Five outer folds hold out entire groups, with random seed 20261027. Three inner group folds choose among seven shared formula families. The population-specific RAR is a separate descriptive comparator. An additional test leaves out each entire population and refits on the other four; no target-population calibration is allowed. Conservative Milky Way and M31 dwarf host holdouts are supplied in Sections 5.6 and 8.1 because holding out individual satellites does not hold out their common environment. Earlier catalog-label host sensitivities were incomplete; the correction is specified in Section 5.6.

The candidate list was written before these broad fits were run, but the datasets and earlier residuals had already been inspected. This is exploratory internal validation, not a preregistered discovery or a genuinely untouched survey confirmation. Comparing the candidates after seeing population holdouts also prevents treating the preferred candidate's holdout score as a final unbiased selection estimate.

We report acceleration error as median absolute predicted-to-observed fractional deviation. Rotation-speed error instead takes the square root of that ratio before subtracting one. These percentages are not interchangeable. Pooled radial medians do not imply equal accuracy in every object. The original SPARC-only benchmark uses equal point weights, so it is reported separately from the broad balanced objective.

No full errors-in-variables likelihood, covariance propagation, distance/inclination refitting, anisotropy marginalization or baryonic-mass uncertainty inference is performed in the main regression. Desmond (2023) demonstrates why nuisance-parameter treatment matters. Our residual scatter is not an estimate of intrinsic scatter; the scores are point estimates conditional on the fixed inputs.

## 4. What the original galaxy improvement means

Holding out entire galaxies, the original point-weighted SPARC calculation gives:

| Model | Coefficients | Log acceleration RMS | Median speed error |
|---|---:|---:|---:|
| Newtonian baryons alone | 0 | 0.5145 dex | 39.46% |
| Logarithmic quadratic | 3 | 0.1352 dex | 9.31% |
| Degree-nine polynomial | 10 | 0.1367 dex | 9.58% |
| Published RAR form | 1 | 0.1334 dex | 9.22% |

The all-data quadratic is y = 0.235654 + 0.716863x + 0.066349x^2. The all-data RAR calibration is a0 = 1.16033 times 10^-10 m/s^2. Fold predictions use coefficients refitted without the test galaxies. Higher polynomial powers do not improve prediction. About 16.1% of the RAR radial speed predictions still have errors above 20%.

A matched-objective extension search finds only a 0.000201 dex RMS reduction from freeing the RAR exponent, with a slightly worse median speed error. A descriptive paired whole-galaxy bootstrap interval for the RMS difference is [-0.001477, 0.001147] dex; it excludes refitting and selection uncertainty. It supplies no clear improvement over standard RAR. The reduction from roughly 39% to 9% is real for this benchmark, but it is relative to baryons alone, not a new advance over the best known formulas.

Freezing the SPARC RAR and applying it elsewhere gives median acceleration errors of 50.0% for classical dwarfs, 86.1% for ultrafaint dwarfs, 34.1% for early-type galaxies, 48.4% for MaNGA cluster centers and 68.8% for CLASH. Those are conditional transfers through each source's measurement model, not formal exclusion probabilities. The large diversity is precisely why a galaxy-only fit cannot settle the broad question.

### 4.1 Established halos and MOND on withheld radial intervals

We now test a different prediction task: omitted radii inside galaxies whose other radii are observed. Of the 149 selected SPARC galaxies, 104 have at least nine retained radii, giving 2,470 eligible measurements. Each of two declared tests withholds 788 points: either the outer third of each galaxy, or its middle third. All models refit without those targets; the remaining 45 galaxies can inform shared disk coefficients. This is not whole-galaxy exclusion.

The familiar NFW and Burkert density forms are

$$\rho_{\mathrm{NFW}}=\frac{\rho_s}{z(1+z)^2},\qquad \rho_{\mathrm{Burkert}}=\frac{\rho_s}{(1+z)(1+z^2)},\qquad z=r/r_s.$$

We integrate each halo, add its Newtonian acceleration to the fixed baryonic acceleration, and fit two positive halo coefficients per galaxy. A common RAR and the Section 5.6 seven-coefficient A43 response are separately trained on the same disk-only targets and equal-galaxy objective. The broad-trained A43 is a distinctly weighted comparison. A declared sensitivity also fits one common stellar multiplier per galaxy, bounded between 0.2 and 2, on the original disk and bulge contributions; gas, distance and inclination remain fixed. The fixed-a0 RAR uses 1.2 times 10^-10 m/s^2.

| Model | Fitted freedom | Outer-third RMS | Middle-third RMS |
| --- | --- | ---: | ---: |
| Fixed-a0 RAR | 0 | 0.11494 | 0.12070 |
| Disk-trained RAR | 1 shared | 0.11346 | 0.12440 |
| Disk-trained A43 | 7 shared | 0.11360 | 0.12515 |
| Broad-trained A43 | 7 shared | 0.16288 | 0.18356 |
| Separate-a0 RAR diagnostic | 1 per galaxy | 0.10482 | 0.05643 |
| NFW halo | 2 per galaxy | 0.09006 | 0.05271 |
| Burkert halo | 2 per galaxy | 0.10982 | 0.04818 |
| Fixed-a0 RAR + stellar scale | 1 per galaxy | 0.09517 | 0.07337 |
| NFW + stellar scale | 3 per galaxy | 0.09483 | 0.04937 |
| Burkert + stellar scale | 3 per galaxy | 0.11950 | 0.04493 |

Scores are equal-galaxy log10 acceleration RMS in dex. NFW/Burkert with two coefficients use 208 fitted quantities across 104 galaxies, versus one common RAR constant. Allowing a separate a0 is a diagnostic of additional freedom, not a claim of universal MOND. The A43 disk-only response does not improve on shared RAR. Halo state and stellar-mass freedom help several within-galaxy predictions, particularly interpolation through the middle block. More freedom does not always improve outer extrapolation.

These halos lack cosmological mass-concentration priors and full nuisance marginalization. In particular, NFW scale radii reach the imposed 1,000-kpc upper bound in 41 outer and 24 middle fits. Descriptive paired galaxy bootstraps give NFW minus shared-RAR intervals [-0.04865, 0.00282] dex outside and [-0.08955, -0.05446] dex in the middle: the former includes zero. The intervals condition on fitted predictions. No inference here ranks all CDM against all MOND implementations. The supplement preserves 20 prediction tables, 1,254 fit sets and all 3,762 converged optimizer attempts.

### 4.2 Predicting halo state without fitting the test galaxy

The radial successes in Section 4.1 motivate a stricter question: can common relations predict halo properties before seeing any rotation measurements from the test galaxy? We retain all 149 selected SPARC galaxies and 2,700 radii. Independent luminosity, neutral-hydrogen mass and photometric effective radius provide a deliberately approximate mass proxy M_p = 0.5 L_3.6 + 1.33 M_HI, in solar units, and R_e. This assumes one stellar mass-to-light ratio even for bulges and omits molecular gas. Neither observed flat speed nor a dynamically inferred halo mass is a predictor.

For each NFW or Burkert shape, define A_h as its halo acceleration at r_s. The first shared family makes log10(A_h/10^-10 m s^-2) and log10(r_s/kpc) affine functions of log10(M_p/10^10 solar masses), requiring four coefficients. A six-coefficient version adds log10(R_e/3 kpc) to each. The complete radial response remains the integrated halo plus measured baryonic acceleration. Coefficients are refitted after excluding every radius of each test galaxy. Three declared fivefold assignments and four predictor-only tail tests are retained; the latter withhold the lowest/highest mass or surface-density quartile, 37 galaxies each. Individual galaxy weights are equal.

We then impose a published structural relation instead of fitting both halo variables independently. Dutton and Maccio (2014) give the Planck, redshift-zero NFW relation

$$c_{200}=10^{0.905}\left(\frac{M_{\rm ref}}{10^{12}h^{-1}M_\odot}\right)^{-0.101},\qquad h=0.671.$$

Here M_ref is explicitly a dark-matter-only-equivalent reference mass within a radius enclosing 200 times critical density. With f_b = 0.0490/0.3175, our two predeclared empirical mass mappings are

$$\frac{M_{\rm ref}}{M_\odot}=\frac{M_p}{f_bM_\odot}+10^A\left(\frac{M_p}{10^{10}M_\odot}\right)^B\left(\frac{R_e}{3\ {\rm kpc}}\right)^C.$$

The mass-only version fixes C = 0; the size version fits A, B and C. The first term enforces M_p <= f_b M_ref. We assign dark mass (1-f_b)M_ref, use r_200 = [G M_ref/(100 H_0^2)]^(1/3), set r_s = r_200/c_200 and add the observed baryonic field separately. The spherical halo is truncated at r_200; none of the observed radii reaches that truncation. This is a conditional deterministic halo approximation, with no response to baryonic assembly, concentration scatter or full cosmological likelihood. The actual baryons-plus-dark-matter system need not have its 200-critical radius at this reference radius. Eight galaxies in the all-data size fit, and eight or nine across its validation splits, extrapolate the low-mass end of the published concentration calibration.

| Common galaxy predictor | Fitted coefficients | Galaxy-held-out RMS, dex | Median speed error |
| --- | ---: | ---: | ---: |
| Fixed published RAR | 0 | 0.15456 | 9.00% |
| Refitted RAR | 1 | 0.15290 | 8.89% |
| Disk-trained A43 | 7 | 0.15417 | 8.80% |
| NFW, mass relation | 4 | 0.15629 | 9.38% |
| NFW, mass and size | 6 | 0.15560 | 9.06% |
| Burkert, mass relation | 4 | 0.16053 | 9.24% |
| Burkert, mass and size | 6 | 0.15940 | 9.29% |
| Published concentration + mass mapping | 2 | 0.15458 | 9.09% |
| Published concentration + mass and size mapping | 3 | 0.15236 | 9.29% |

The three-coefficient all-data mapping has A = 11.46728, B = 0.99375 and C = -0.63550. Validation uses fold-specific coefficients instead. The published concentration relation improves the point estimate over our six-coefficient NFW relation in all three random splits while using fewer fitted constants. Its improvement over RAR is only 0.00054 dex in the primary split and 0.00022 in the second; it is worse by 0.00119 in the third. The primary descriptive paired bootstrap interval for its difference from RAR is [-0.01144, 0.01020] dex. Median speed errors are worse. These results establish comparable conditional accuracy, not a reliable gain over MOND.

Transfer to the four feature tails remains mixed. For the lowest surface-density quartile, three-coefficient NFW gives 0.19198 dex versus RAR's 0.19653; it is worse than RAR in the other three tails. Even the isolated tail gains have intervals spanning zero. The common halo relations therefore compress object-specific fits into a few measured-input coefficients, but do not yet explain the diversity or unify other populations. Baryon-halo scaling and its connection to RAR already have substantial prior art, including Navarro et al. (2017) and Li et al. (2019). This is a matched predictive comparison of such ideas, not their first proposal.


### 4.3 Starting from published feedback profiles

We next test DC14 (Di Cintio et al. 2014) and coreNFW (Read, Agertz and Collins 2016), both of which encode a baryonic response in a halo profile. Li et al. (2020) compare these and other established profiles while fitting individual galaxies and marginalizing over nuisance parameters. Our question is narrower: can those profile prescriptions improve the common measured-input predictor of Section 4.2, with no rotation measurement from the test galaxy entering its fit? This is not a reproduction of their full inference or a ranking of complete theories.

The mass mapping, three starting points, bounds, 149 galaxies and 2,700 radii are retained. We reuse the original three NumPy fivefold assignments, plus the four 37-galaxy feature tails. The table reports means across the three complete held-out prediction sets, rather than the primary split alone. A local protocol was frozen before these fits on data already explored by the project.

For DC14 we translate the reference NFW halo to the original virial convention, with overdensity 92.8, and use its published stellar-to-dark-halo-mass dependence of the density slopes and concentration. Both shapes are normalized to the same dark mass at that reference radius. Values below the recipe's stellar-to-halo mass domain revert to NFW; above it the recipe is capped. This is a declared Planck-reference transfer of a WMAP3 calibration, not an exact simulation match. In the full-sample DC14 fit, one galaxy reaches the upper-domain cap and 71 lie outside the original simulated halo-mass interval. The protocol records those extrapolations and the full conversion.

We also fit one common response strength, with nonnegative densities and the same finite total dark mass:

$$\rho_{\rm response}(r)=(1-\lambda)\rho_{\rm NFW}(r)+\lambda\rho_{\rm DC14}(r),\qquad 0\leq\lambda\leq1.$$

Exact endpoint fits are retained. For coreNFW we use the assumed-duration prescription described by Li et al. (2020): core radius 1.75 R_e and core strength tanh(0.04 times 14 Gyr divided by the halo dynamical time). This is an assumed common duration, not measured history or a detection of gravitational memory. The supplement supplies the precise enclosed-mass equations, normalization, quadrature checks and every optimizer attempt.

| Shared predictor | Global coefficients | Mean held-out RMS, dex | Median speed error |
| --- | ---: | ---: | ---: |
| Refitted RAR | 1 | 0.15309 | 8.93% |
| Constrained NFW | 3 | 0.15323 | 9.30% |
| DC14-informed | 3 | 0.18878 | 10.88% |
| NFW/DC14 response | 4 | 0.15376 | 9.29% |
| coreNFW, assumed duration | 3 | 0.16253 | 10.02% |
| RAR/NFW force blend | 5 | 0.15060 | 8.88% |

The feedback profiles do not improve this deterministic shared mapping. The response selects lambda = 0.04894 on all galaxies, close to NFW, and is slightly worse than NFW on held-out galaxies. This constrains this particular coupling and mass mapping, not baryonic feedback in general. Allowing halo scatter, object-specific histories and nuisance parameters could change the comparison. A separate JavaScript reconstruction used different galaxy partitions and found no improvement from a cosmologically constrained Einasto profile or a single fitted common core strength; its complete negative results are retained in dm_fallback_validation.

### 4.4 A small ensemble gain and a failed dwarf transfer

The last row combines two existing predictors with a common force weight:

$$g_{\rm blend}=(1-w)g_{\rm RAR}+w(g_{\rm bar}+g_{\rm NFW}),\qquad 0\leq w\leq1.$$

Within each training split we fit the RAR constant and three halo-mapping coefficients separately, then fit w on their training predictions. All five quantities are fixed before predicting the excluded galaxies. There is no inner cross-fitting of the weight. The full-sample value is w = 0.53209; held-out predictions use their training-only values. This is an empirical ensemble of established formulas, without a priority claim or a field-theoretic derivation.

The blend lowers mean held-out RMS by 0.00249 dex, or 1.63%, relative to refitted RAR, and lowers it in all three repeats. Its average median speed error changes only from 8.93% to 8.88%. A paired bootstrap resamples whole galaxies with all repeats together; the 95% interval for blend minus RAR is [-0.00880, +0.00403] dex. It includes zero and conditions on the fitted predictions, omitting refit and model-search uncertainty. The blend improves RMS in three feature tails but worsens the highest-mass tail, 0.12583 versus 0.12252 dex. Thus the gain is small and not established as a reproducible population advantage.

The earlier JavaScript partitions gave a 2.04% RMS gain, 0.15060 versus 0.15374 dex, with an interval likewise spanning zero and a slightly worse median speed error. The blend was introduced after inspecting those initial model results. Applying its frozen algorithm to the original NumPy assignments checks implementation and partition sensitivity on the same already inspected galaxies; it is not independent astronomical confirmation.

As a further transfer test, freeze the full-SPARC coefficients and apply all six predictors to the archived 44 dwarf-spheroidal half-mass estimates. Use stellar mass and photometric size, with no dispersion or inferred halo mass as a predictor. RMS values are 0.66937 dex for RAR, 0.66757 for NFW, 1.26536 for DC14-informed, 0.65052 for the NFW/DC14 response, 0.65977 for assumed-duration coreNFW, and 0.66613 for the blend. None gives an accurate unified disk-and-dwarf relation. These are isolated mass-estimator stresses; they do not substitute for a Jeans likelihood, orbital/environment modeling, or a full MOND external-field calculation.

The blend also lacks a finite-mass dark-halo interpretation at arbitrarily large radius: for an isolated finite baryonic source and w < 1, its RAR contribution eventually falls as 1/r while the truncated halo falls as 1/r squared. This asymptote is not by itself an objection to modified gravity, but it prevents interpreting this blend as a complete finite dark halo. No lensing or cosmological field law follows from the weighted force formula.

### 4.5 Predicting a distribution of halo states

An exact halo assignment ignores both population variation and observational uncertainty. We restore the selected SPARC velocity errors, catalog distance and inclination errors, and stellar/gas components, then compare conditional predictive distributions. Dutton and Maccio (2014) supply a 0.11-dex concentration scatter; the stellar-amplitude prior has 0.1-dex width, following Li et al. (2020). We use one common disk/bulge multiplier and linearize the source and observational transformations. Distance changes radii and photometric masses together; inclination changes the inferred speed. These approximations omit other correlated systematics and are not a reproduction of full published MCMC analyses.

For each galaxy, the covariance of its log-acceleration residual vector is

$$C_j=\operatorname{diag}(s_{v,jk}^{2}+s^{2})+U_jU_j^{T},\qquad s_{v,jk}=\frac{2\sigma_{v,jk}}{v_{jk}\ln 10}.$$

Columns of U carry the galaxy-shared nuisance responses and, where included, halo concentration and excess-mass variation. We fit an independent radial residual scale s. The mass-scatter model additionally fits a Gaussian width for the logarithm of the positive excess-mass term in Section 4.2; this is an empirical conditional width, not the inverse of a published stellar-mass scatter. RAR comparisons include the same nuisance treatment, with or without an additional marginalized galaxy-wide residual offset. No held-out galaxy is assigned a latent state inferred from its own kinematics.

The objective is the mean per-galaxy joint Gaussian negative log density divided by that galaxy's radius count. This weighted marginal log score penalizes excessive uncertainty through the covariance determinant. It is distinct from RMS and remains a first-order, plug-in population approximation. Three original fivefold assignments give:

| Conditional predictor | Mean negative log score | Central RMS, dex | 95% marginal coverage |
| --- | ---: | ---: | ---: |
| RAR plus nuisance errors | -0.9803 | 0.15592 | 87.8% |
| RAR plus shared residual offset | -1.0414 | 0.15461 | 96.4% |
| NFW plus concentration scatter | -1.1056 | 0.16400 | 92.2% |
| NFW plus mass/concentration scatter | -1.2313 | 0.15883 | 96.1% |
| NFW plus mean response and scatter | -1.2281 | 0.15990 | 96.1% |
| Previous force blend plus offset | -1.0420 | 0.15060 | 96.6% |

Lower log scores are better; density units are log10 acceleration and comparisons use identical observations. NFW with both halo scatters improves the score relative to offset-RAR by 0.1898 natural-log units per radius, averaged equally over galaxies; its conditional paired interval is [-0.3037, -0.0912]. This is an improvement in the tested predictive distribution, not in the central galaxy RMS. Its full-sample mass width is large, 0.548 dex, and some split fits reach the 0.6-dex bound. Nominal 68% coverage is 79.2%, so calibration is not uniformly correct. The result does not establish superiority to full MOND or cosmological inference.

Nonlinear quadrature of the halo states reveals a median shift of 0.0120 dex from the nominal central prediction and up to 17.5% discrepancy in the halo-only standard deviation. Replacing the linear halo moments by nonlinear moments, with all coefficients frozen, retains the log-score improvement: -0.1929 relative to offset-RAR, interval [-0.3045, -0.0952]. This sensitivity still approximates the resulting mixture by a Gaussian and keeps the other nuisance responses linearized. It is not an exact nonlinear marginal likelihood. A separate four-coefficient concentration-offset repair, optimized for the original RMS, gives 0.15353 dex and provides no point-prediction gain.

### 4.6 A better disk-to-dwarf transfer from stellar mass

The population fits motivate an explicitly adaptive follow-up. Their disk-trained NFW parameters give dwarf RMS around 0.50-0.51 dex, better than the previous 0.66757, but worsen disk RMS. We therefore test the stellar-to-halo relation directly. Such relations and their curvature have substantial prior art, including Moster, Naab and White (2013); they are not new physical laws. The follow-up is fixed locally after seeing the preceding results on already inspected data.

Four finite-NFW prescriptions are compared: the inverse published median Moster relation; a version with a common normalization and size dependence; a three-coefficient stellar power law; and a smooth, five-coefficient curved stellar law. All retain the published concentration relation and baryon-budget floor; the fitted extensions include photometric size. The quoted Moster coefficients are transferred to our common Planck reference convention, without rerunning abundance matching. The simple stellar rule differs from Section 4.2 only in which mass enters its fitted power:

$$\frac{M_{\rm ref}}{M_\odot}=\frac{M_p}{f_bM_\odot}+10^{11.74805}\left(\frac{M_\star}{10^{10}M_\odot}\right)^{0.81783}\left(\frac{R_e}{3\ {\rm kpc}}\right)^{-0.71187}.$$

Here M_star = 0.5 L_3.6 and M_p = M_star + 1.33 M_HI as before. Gas still contributes to the measured baryonic field and the mass floor. The three coefficients shown are fitted to all SPARC galaxies; validation uses training-only coefficients. Both the original and stellar rules use three shared coefficients. The same finite NFW prescription then applies to disks and dwarf spheroidals without a population-specific parameter set.

| Shared point predictor | Coefficients | Disk held-out RMS | Frozen dwarf RMS |
| --- | ---: | ---: | ---: |
| RAR, disk-trained | 1 | 0.15309 | 0.66937 |
| Previous baryonic-mass NFW | 3 | 0.15323 | 0.66757 |
| Published Moster relation + NFW | 0 | 0.17111 | 0.55591 |
| Stellar power + NFW | 3 | 0.15329 | 0.47914 |
| Curved stellar relation + NFW | 5 | 0.15380 | 0.47874 |
| Previous RAR/NFW blend | 5 | 0.15060 | 0.66613 |
| RAR/stellar-NFW blend | 5 | 0.15051 | 0.50521 |

All RMS values are dex. Disk columns average the three complete galaxy-held-out prediction sets. Dwarf coefficients are frozen from full-SPARC fits before transfer to the same 44 half-mass estimates. The scaled Moster alternative gives 0.17128/0.56633 dex and does not improve either column over the fixed version. The curved relation does not improve disk prediction over the stellar power; its transition reaches the upper bound in 12 of 20 evaluation cases. Inverting a population median is also not equivalent to predicting a conditional halo-mass distribution with selection and scatter.

The stellar power lowers dwarf RMS by 28.2% versus the previous NFW mapping while changing disk RMS by only +0.04%. Median disk speed error improves from 9.30% to 9.09%. A conditional paired dwarf bootstrap gives a difference of -0.18843 dex with interval [-0.28081, -0.08311]; it does not include formula-search, host-correlated or mass-reconstruction uncertainty. Thirty of 44 individual dwarf log errors improve. The prior 25-object higher-precision sensitivity gives 0.33780 versus 0.57502 dex, a 41.3% RMS reduction, but this is another reused subset.

The gains are not uniform. Classical-dwarf RMS worsens from 0.34763 to 0.39482; ultrafaint RMS improves from 0.93435 to 0.57142. Conservative host matching retains all catalog variants, including G/L and A/L: the Milky Way subset improves, while Andromeda worsens. Three of four disk predictor-tail scores also worsen relative to the previous NFW rule. Subclasses and known failure cases must therefore remain visible beside the aggregate gain.

After this comparison we apply the unchanged training-only blend algorithm with the stellar NFW base. Its 0.15051-dex disk RMS is only 0.00010 below the previous blend, with an interval spanning zero. Its dwarf RMS improves by 24.2%, and its full-sample weight is w = 0.58209. This blend improves both host subsets and narrowly improves classical-dwarf RMS versus the previous pure NFW rule, but remains worse for classical dwarfs than disk-trained RAR. The extra blend test is post-result exploration, not independent confirmation.

This is a useful improvement in a specific shared disk-to-dwarf prediction task. The remaining dwarf errors are large, satellite evolution is not modeled, and no new cluster, Solar-System or cosmological likelihood has been fitted. The result strengthens the case for independently measured source composition in halo-state prediction; it neither detects gravitational memory nor supplies an all-scale field theory.

## 5. Shared formulas across the wider sample

### 5.1 Candidate families and aggregate performance

We test standard RAR, RAR with a free exponent, linear or quadratic log-radius dependence of its effective acceleration scale, quadratic potential-proxy or density-proxy dependence, and a full quadratic in x and u. The potential proxy is q = x + u; the spherical density proxy is x - u. Both are rearrangements of the same two inputs, not independent new measurements. A piecewise comparator gives one RAR scale to each of the five populations.

The free exponent replaces sqrt(g_bar/a_eff) inside the exponential by (g_bar/a_eff)^p. For RAR variants we fit log10(a_eff/a_ref), with intercept in [-3,3], linear coefficient in [-2,2], quadratic coefficient in [-1,1], and free exponent where used in [0.1,0.9]. The two-variable polynomial has six unconstrained coefficients. Bounds are computational choices; they do not make these relations complete physical laws.

| Candidate | Free coefficients | Group-held-out B | Population-held-out B |
|---|---:|---:|---:|
| Common recalibrated RAR | 1 | 0.3360 | 0.3671 |
| RAR with free exponent | 2 | 0.3400 | See supplement |
| RAR with linear log-radius scale | 2 | 0.3395 | See supplement |
| RAR with quadratic log-radius scale | 3 | 0.2820 | See supplement |
| RAR with potential-proxy quadratic | 3 | 0.2495 | 0.3158 |
| RAR with density-proxy quadratic | 3 | 0.3382 | See supplement |
| Full two-variable quadratic | 6 | 0.2355 | 0.4883 |
| Population-specific RAR | 5 | 0.2817 | Undefined for unseen class |

Scores are logarithmic RMS in dex with the same population/group balancing. The common recalibrated RAR has all-data a0 approximately 5.34 times 10^-10 m/s^2, quite different from the galaxy-only value. Its larger scale sacrifices disk-galaxy accuracy. The frozen galaxy-calibrated RAR has broad B = 0.4172 dex, but that mixed score includes its SPARC training sample and is not a fully held-out baseline.

The flexible polynomial wins the internal comparison, yet deteriorates when an entire population is absent from training. Nested inner selection chooses that polynomial in four outer folds and the potential-proxy form in one, with combined B = 0.2451 dex. The less flexible potential-proxy form has steadier population transfer. This favors further testing of that structure, not a claim that the lowest internal residual identifies the correct law.

![Broad acceleration comparison](broad_acceleration_planes.png)

Figure 1. The broad benchmark separates populations, compares them with the frozen galaxy RAR, and shows its residuals versus physical radius. Inferred accelerations retain each source's dynamical or lensing assumptions. Curves and residual patterns compress heterogeneous data; they do not establish identical measurement likelihoods.

### 5.2 An explicit common candidate

A compact fitted expression is

$$g_{\mathrm{pred}}=\frac{g_{\mathrm{bar}}}{1-\exp[-\sqrt{g_{\mathrm{bar}}/a_{\mathrm{eff}}}]}.$$

$$q=\log_{10}\left[\frac{g_{\mathrm{bar}}}{10^{-10}\,\mathrm{m\,s^{-2}}}\frac{r}{\mathrm{kpc}}\right].$$

$$\log_{10}(a_{\mathrm{eff}}/a_{\mathrm{ref}})=0.066728+0.350081q+0.172091q^2.$$

These are the all-data coefficients; validation uses fold-specific estimates. No object-specific adjustment or population label enters this expression. The product g_bar times r has potential units, but does not generally equal the full baryonic gravitational potential. The proxy is especially incomplete for disks, extended distributions and external fields. The measured pattern can also reflect selection, stellar mass calibration and differing dynamical reconstructions.

| Population | Held-out median acceleration error |
|---|---:|
| SPARC rotating galaxies | 29.10% |
| Early-type galaxies | 39.62% |
| Classical and ultrafaint dwarfs combined | 53.13% |
| MaNGA brightest cluster galaxies | 21.14% |
| CLASH lensing clusters | 34.43% |

For SPARC this corresponds to 14.29% median speed error, compared with the earlier galaxy-optimized 9.22%. Classical and ultrafaint dwarf acceleration errors are 36.96% and 61.82%, respectively. There is therefore no common 9% solution here. The broad correction improves the balanced comparison while leaving substantial and population-dependent residuals.

The linear radius coefficient is nearly zero in its own fit, whereas curvature in q helps. Excess inferred acceleration occurs in both faint, small dwarfs and massive cluster environments; it is not a monotonic function of size. Wolf et al. (2010) already documented a broad U-shaped mass-to-light trend across systems. The present curvature is an empirical way to represent part of that diversity, not a discovery of the trend or proof of a potential-dependent force.

As a post-search sensitivity, retaining only 25 dwarfs with relative dispersion uncertainty at most 25% and no selected caution flag gives potential-proxy coefficients (0.1410, 0.3121, 0.1603). Positive curvature remains. Its balanced group score is 0.2132 dex versus 0.2997 for common RAR; population-holdout scores are 0.2872 versus 0.3358. Dwarf median error falls to about 31.6%, while SPARC acceleration error remains about 31.1%. This sensitivity reduces the influence of the least precise dwarfs without eliminating modeling or selection uncertainty.

![Transfer comparison](population_holdout_comparison.png)

Figure 2. The polynomial's advantage within represented populations does not persist when a whole class is withheld. All fold coefficients and predictions are supplied.

### 5.3 Extrapolation and a limited high-acceleration repair

An unconstrained polynomial is not a universal force law. Evaluated far outside its galaxy domain along the central-Sun acceleration-radius relation, the original galaxy quadratic gives acceleration ratios to Newton of 110.57 at 1 AU, 6.77 at 9.58 AU and 2.58 at 30.07 AU. These are asymptotic checks, not fits to ephemerides. The positive q-squared term in the new candidate also makes unrestricted extreme extrapolation unsafe.

A post-search mathematical regulator replaces the candidate's right-hand side by

$$\log_{10}(a_{\mathrm{eff}}/a_{\mathrm{ref}})=\frac{a+bq+cq^2}{1+(g_{\mathrm{bar}}/g_{\mathrm{cut}})^2},\qquad g_{\mathrm{cut}}=10^{-8}\,\mathrm{m\,s^{-2}}.$$

The cutoff is chosen, not inferred. Applying it with the existing fold coefficients and no refitting changes B only from 0.249461 to 0.249511 dex, and population-transfer B from 0.315811 to 0.315886. It retains a negligible isolated central-force correction at the three illustrative Solar-System radii while repairing the unrestricted high-acceleration asymptote. This preserves the empirical fit while repairing one limiting behavior. It does not calculate an external-field quadrupole, planetary ephemerides, relativistic effects or a valid low-acceleration limit at arbitrary radius. The regulator is not independently validated or part of the originally selected model list.

### 5.4 One smooth power-sum formula

The earlier three-coefficient candidate is already one shared formula: the same constants apply to every system. The piecewise RAR is a different comparator that uses a separate scale for each population. To test whether another single expression could match or improve the common fit, we specify a follow-up list of thirteen shared families and that one piecewise comparator. The data are unchanged. Besides the earlier baselines, this list includes an acceleration floor, polynomial extensions, a learned acceleration-radius combination, a positive power sum, and a reciprocal-logarithm family. The complete definitions, bounds and start points are supplied.

Two competitors have direct prior art. Lelli et al. (2017, equation 14) add an exponentially suppressed acceleration floor to the RAR. Santos and Felizardo (2026) report a reciprocal-logarithm family for galaxy data. We independently refit those functional forms under the present selection and weighting; this is not a reproduction of their original pipelines or a comparison with their quoted paper statistics. The proposed power sum is a conventional empirical combination of powers, with no priority claim for that construction.

The follow-up retains the original five whole-group folds and adds two fixed alternative partitions, with seeds 7331 and 48105. Each candidate also faces whole-population exclusion. Three inner group folds select among shared families within the original outer folds. There are 497 fitting calls with no reported convergence failures; all 42 group-validation prediction files were independently checked against the recorded balanced scores. These checks establish numerical consistency, not independence from prior inspection of the dataset.

One smooth candidate uses the dimensionless input

$$s=\frac{g_{\mathrm{bar}}}{a_{\mathrm{ref}}}\frac{r}{\mathrm{kpc}},\qquad a_{\mathrm{ref}}=10^{-10}\,\mathrm{m\,s^{-2}},$$

with an effective acceleration scale

$$a_{\mathrm{eff}}=a_{\mathrm{ref}}\left(C_0+C_1s^{-\alpha}+C_2s^\beta\right),$$

and the same RAR response

$$g_{\mathrm{pred}}=\frac{g_{\mathrm{bar}}}{1-\exp[-\sqrt{g_{\mathrm{bar}}/a_{\mathrm{eff}}}]}.$$

The all-data coefficients are

$$C_0=0.72239028,\quad C_1=0.00248832,\quad C_2=0.63776531,$$

$$\alpha=1.00000000,\qquad\beta=0.77129797.$$

All five coefficients are shared. There is no object-specific fit or galaxy-versus-cluster switch. Population labels are used only to balance the objective and define validation groups, not as predictors. Writing intermediate quantities does not make the expression piecemeal: their substitution gives one smooth function of positive baryonic acceleration and radius. The three terms contribute continuously. Parameter bounds are log10(C0) in [-3,3], log10(C1) and log10(C2) in [-5,3], alpha in [0,1] and beta in [0,0.95]. Alpha reaches its imposed upper bound, so its fitted value is not an independently established exponent. Validation always refits coefficients without the held-out systems.

| Selected comparison | Coefficients | System-held-out B | Population-held-out B |
|---|---:|---:|---:|
| Common RAR | 1 | 0.3360 | 0.3671 |
| Previous potential quadratic | 3 | 0.2495 | 0.3158 |
| Potential quadratic plus radius | 4 | 0.2447 | 0.3081 |
| Smooth power sum | 5 | 0.2496 | 0.2776 |
| Full two-variable response quadratic | 6 | 0.2355 | 0.4883 |
| Population-specific RAR | 5 | 0.2817 | Undefined for unseen class |

All scores use the same population/group-balanced log10 acceleration RMS in dex. The power sum essentially matches, rather than improves, the previous candidate on the original system folds. Its two alternative-partition scores are 0.2457 and 0.2455 dex, compared with 0.2475 and 0.2468 for the previous potential quadratic. Population transfer improves in this experiment. However, choosing a candidate after examining those transfers makes that score a selection result, not an unbiased confirmation on a new survey.

Nested selection gives B = 0.2423 dex and chooses other flexible forms: the quadratic effective RAR scale in three folds, the response quadratic in one and the potential cubic in one. The cubic illustrates the extrapolation danger: its internal B of 0.2455 rises to 1.3488 on population holdout. All fourteen candidates, including weaker results, remain in the supplement. The power sum is highlighted for its compact positive form and steadier transfer, not because it wins every selection criterion.

![Shared-formula comparison](one_formula_comparison.png)

Figure 3. The power sum trades disk accuracy for steadier transfer across the represented populations. The right panel reports SPARC speed error, a different quantity from balanced logarithmic acceleration RMS. The galaxy-only RAR benchmark is shown as a separate reference.

| Population | Previous candidate: acceleration error | Power sum: acceleration error |
|---|---:|---:|
| SPARC rotating galaxies | 29.10% | 33.55% |
| Early-type galaxies | 39.62% | 40.56% |
| Dwarf spheroidals | 53.13% | 56.16% |
| MaNGA brightest cluster galaxies | 21.14% | 19.19% |
| CLASH lensing clusters | 34.43% | 29.34% |

These are median absolute acceleration errors in the original whole-system folds. The corresponding SPARC speed error rises from 14.29% to 16.06%; the separately optimized galaxy RAR still gives 9.22%. This prevents a claim of uniform improvement, even though the power sum is a common expression and improves the balanced population-transfer score.

The expression also has a substantive physical limitation. At fixed radius and arbitrarily small positive baryonic acceleration, its fitted inverse-power term yields

$$g_{\mathrm{pred}}\longrightarrow a_{\mathrm{ref}}\sqrt{C_1\frac{\mathrm{kpc}}{r}}\simeq4.99\times10^{-12}\sqrt{\frac{\mathrm{kpc}}{r}}\,\mathrm{m\,s^{-2}}.$$

This nonzero floor is an empirical extrapolation, not a derivation of gravity in empty space. At high baryonic acceleration and fixed radius, beta below one gives the Newtonian limit, but neither limit supplies a complete vector, relativistic or cosmological theory. No new Solar-System, CMB, shear-map or temporal-memory fit was performed in the follow-up. Existing RAR, potential-dependent, nonlocal and halo models have not all been compared with matched observable-level likelihoods. The defensible claim is a useful shared fitting family with an explicit domain and known tradeoffs.

### 5.5 Testing physical repairs and uncertainty

A follow-up on 28 October 2026 fixed three additional candidates before running them on the same previously inspected data. This local protocol is supplied for transparency; it is not an independent preregistration or a new test sample. The same population weights, three five-fold group assignments and whole-population exclusions were used. All 63 fitted models converged under the stated optimizer criteria; this does not establish global optimality. All attempted families are reported.

The first repair softens the inverse power using a positive transition scale s0:

$$a_{\mathrm{eff}}=a_{\mathrm{ref}}[C_0+C_1(s+s_0)^{-\alpha}+C_2s^{\beta}].$$

This six-parameter form is finite as s approaches zero. At fixed positive radius, its predicted acceleration therefore vanishes as sqrt(g_bar), instead of leaving the previous finite acceleration floor. Its all-data fit drives log10(s0) to the lower permitted bound -5, so the data do not identify a resolved turnover scale. Population transfer deteriorates; the extra parameter is not supported as a predictive improvement.

The second repair bounds the entire effective scale through a five-parameter smooth function, with q = log10(s):

$$\log_{10}(a_{\mathrm{eff}}/a_{\mathrm{ref}})=A+BT+CT^2,$$

$$T=\tanh[(q-q_c)/w],\qquad w>0.$$

Here the effective scale is finite at both ends. The RAR response tends to zero at low g_bar and to the Newtonian value at high g_bar, at fixed radius. The fitted tuple (A, B, C, q_c, w) is (0.14622, -3.00000, 6.00000, -2.61914, 4.70749). B and C reach bounds; these coefficients are conditional empirical estimates, not measured constants. The limiting effective scales are about 0.14 m/s^2 at q approaching minus infinity and 1.4 times 10^-7 m/s^2 at plus infinity. In particular, the low-end value is far above ordinary galactic acceleration scales. Correct qualitative limits therefore do not establish a physically acceptable amplitude. The observed whole-population transfer is also worse.

The third variant fixes alpha at one-half in the original power sum, leaving four fitted coefficients. It is still one shared expression, with no population label or object-specific coefficient. The all-data tuple (log10 C0, log10 C1, log10 C2, beta) is (-3.00000, -0.78140, 0.03653, 0.64401). C0 reaches its lower bound. Its low-g_bar acceleration is proportional to g_bar raised to the one-quarter power at fixed radius, which vanishes but does not reproduce the usual square-root limit.

| Shared model | Coeff. | System-held-out B (dex) | Population-held-out B (dex) | SPARC speed error |
| --- | ---: | ---: | ---: | ---: |
| Earlier power sum | 5 | 0.2496 | 0.2776 | 16.06% |
| Softened inverse power | 6 | 0.2516 | 0.3498 | 16.09% |
| Bounded tanh quadratic | 5 | 0.2465 | 0.3399 | 18.33% |
| Fixed half-power | 4 | 0.2609 | 0.2640 | 16.48% |

The two alternate group-fold seeds give B = 0.2481/0.2476 for the softened form, 0.2429/0.2431 for the bounded form and 0.2594/0.2580 for the half-power form. These repetitions do not remove the tradeoff between internal fit and population transfer. There is no justified replacement that is better on all stated criteria.

A separate analytic check considers an isolated spherical finite baryonic mass, with g_bar = GM/r^2 at large r. For a dominating inverse-power term s raised to -alpha, the RAR response gives g proportional to r raised to (-1 + alpha/2), and circular speed proportional to r raised to alpha/4. The earlier alpha = 1 form therefore extrapolates to speed proportional to r raised to 1/4, and the half-power form to r raised to 1/8. In contrast, a finite positive limiting effective scale gives g proportional to 1/r and asymptotically constant circular speed. This is a useful structural constraint, not a fit to new outer-halo data. None of these radial scalar formulas supplies a vector field equation, momentum conservation, lensing metric or cosmology by itself.

To quantify conditional sampling uncertainty, we performed 5,000 paired whole-group Poisson resamples of the fixed out-of-fold residuals. A group's multiplier is shared wherever that group occurs, and each population is renormalized separately. These descriptive 95% percentile intervals exclude refitting, model selection, observational uncertainty and reconstruction systematics. Population-exclusion intervals resample systems within the same five populations; they do not quantify uncertainty over all possible future populations.

For system holdouts, power sum minus common RAR gives a B difference of -0.0864 dex, with interval [-0.1154, -0.0576]. Power sum minus potential quadratic gives +0.00017 dex, with interval [-0.00266, +0.00287]: no resolved internal advantage. For population exclusions, power sum minus potential quadratic gives -0.0382 dex, with interval [-0.0577, -0.0148], conditional on the already inspected population tests. The new half-power variant minus power sum gives -0.0137 dex, with interval [-0.0350, +0.0082]. Its apparent improvement is not resolved by this conditional resampling. The bounded model's small internal gain likewise has an interval crossing zero, [-0.00874, +0.00224] dex.

The actionable finding is that a shared fit can outperform a single recalibrated RAR across this mixed benchmark while retaining substantial tradeoffs. A physically preferable limit and a better numerical transfer score are distinct requirements. Recording both, including unsuccessful repairs, improves the study more than promoting another selected minimum as a universal solution.

### 5.6 Mass dependence, broader cluster coverage and a regularized shared response

A further search treats physical regularity as a criterion alongside prediction error. Six new families use the dimensionless mass coordinate U = (g_bar/a_ref)(r/kpc)^2, alongside s = (g_bar/a_ref)(r/kpc). The corresponding mass a_ref kpc^2 U/G is a spherical-equivalent quantity. It equals enclosed baryonic mass under spherical symmetry, but not generally for disks. These coordinates reorganize the same two inputs; they do not add independent measurements. Shared mass-dependent responses have prior art, including Moffat and Toth (2009).

The first stage tests mass-based power sums, a bounded mass response, an additive response, a mixed mass/potential form, a learned index and an acceleration-screened variant. The mixed form gives system-held-out B = 0.2281 dex versus 0.2496 for the earlier power sum, predominantly improving dwarfs. Its population-held-out score, 0.2770 versus 0.2776, is essentially tied. The inverse mass exponent reaches its upper bound 0.95. All six candidates and four comparators are retained.

**Physical failures determine the repairs.** The mixed form predicts 22% excess isolated-Sun acceleration at 10 AU and 130% at 30 AU. Correct limits at arbitrarily high acceleration do not prevent that finite-radius failure. It also diverges in a smooth baryonic core: g_bar proportional to r implies U proportional to r^3, so the inverse leg yields g_pred proportional to r^((1-3 alpha)/2), or r^-0.925 at the fitted exponent. A second stage gates that leg at high acceleration and tests a squared radial cutoff for the positive potential term. Planetary spot checks improve, but the core singularity and large stellar low-acceleration boost remain. The squared cutoff can also produce negative inferred extra-density shells around a point source; this prevents a nonnegative dark-halo interpretation there, although signed effective density alone does not exclude all modified gravity.

**Expanded joint sample.** We add the 21 fixed-radius X-COP estimates described in Section 7. The 20 CLASH and seven X-COP hosts share one clusters class, with equal host weights; the other four classes retain equal total weights. The expanded sample has 286 systems, 285 conservative groups and 2,931 rows. Its weighting differs from the original benchmark, so the two objectives are reported separately. The added measurements were inspected previously; this is exploratory reuse, not new-survey confirmation. No source-specific coefficient is fitted. Complete X-COP and CLASH exclusions additionally test transfer between hydrostatic and lensing reconstructions.

A final, separately specified repair tests exactly one seven-coefficient family:

$$U=\frac{g_{\mathrm{bar}}}{a_{\mathrm{ref}}}\left(\frac{r}{\mathrm{kpc}}\right)^2,\qquad s=\frac{g_{\mathrm{bar}}}{a_{\mathrm{ref}}}\frac{r}{\mathrm{kpc}}.$$

$$S(r)=\frac{r^4}{r^4+\ell^4},\qquad W(g_{\mathrm{bar}})=\frac{1}{1+(g_{\mathrm{bar}}/g_c)^2}.$$

$$\frac{a_{\mathrm{eff}}}{a_{\mathrm{ref}}}=C_0+C_1U^{-\alpha}S(r)W(g_{\mathrm{bar}})+\frac{C_2s^\beta}{1+r/R}.$$

$$g_{\mathrm{pred}}=\frac{g_{\mathrm{bar}}}{1-\exp[-\sqrt{g_{\mathrm{bar}}/a_{\mathrm{eff}}}]},\quad a_{\mathrm{ref}}=10^{-10}\,\mathrm{m\,s^{-2}},\quad g_c=10^{-8}\,\mathrm{m\,s^{-2}}.$$

The powers four and one and g_c are fixed. We fit the three positive amplitudes, alpha, beta, ell and R, with the same shared values for every object. Parameter bounds, three starts and all fold coefficients are supplied. The imposed ell range is 10 pc to 1 kpc; its lower bound is chosen to suppress the inverse-mass excess at stellar separations. It is not a measurement from binary observations. R ranges from 1 to 10,000 kpc; alpha and beta lie in [0,0.95]. For the expanded all-data fit,

$$C_0=0.880950439,\quad C_1=0.00077405352,\quad C_2=0.401566357,$$

$$\alpha=0.95,\quad\beta=0.95,\quad\ell=10\,\mathrm{pc},\quad R=308.221\,\mathrm{kpc}.$$

Both exponents and ell reach their respective bounds; their displayed values therefore are constrained fit results, not independently established constants. All reported validation refits the coefficients without the held-out objects. The family is smooth for positive inputs and has no galaxy/cluster branch or population label. Its two length scales describe one expression, not a separate fit for each system.

| Shared model | Original system B | Original population B | Expanded system B | Expanded population B |
|---|---:|---:|---:|---:|
| Common RAR | 0.3360 | 0.3671 | 0.3299 | 0.3593 |
| Earlier power sum | 0.2496 | 0.2776 | 0.2475 | 0.2749 |
| Mass/potential with high-acceleration gate | 0.2281 | 0.2768 | 0.2279 | 0.2756 |
| Gate plus squared radial cutoff | 0.2250 | 0.3189 | 0.2198 | 0.3212 |
| Regularized seven-coefficient response | 0.2264 | 0.2874 | 0.2220 | 0.2892 |

These are equal-population/equal-group logarithmic RMS scores in dex. The final response improves expanded whole-system RMS by 10.3% relative to the earlier power sum, with scores 0.2217 and 0.2209 in the two alternative partitions. A 5,000-draw paired group-resampling interval for its primary difference is [-0.0411,-0.0093] dex. This descriptive interval conditions on the fitted predictions and omits refitting, candidate selection and measurement-systematic uncertainty. Whole-population transfer instead deteriorates by 0.0143 dex, with interval [-0.0053,+0.0326]. No uniform predictive gain is established. The final family trades a little internal accuracy for repaired limits relative to the squared-cutoff fit.

![Regularized shared-response comparison](unified_regular_comparison.png)

Figure 4. The expanded objective weights five classes, with CLASH and X-COP together in the clusters class. Physical regularity and whole-system accuracy improve, but entire-population transfer remains a separate test. Dwarfs retain the largest residuals.

On the expanded sample, the final model's median acceleration errors are 31.4% for disks, 51.7% for dwarfs, 40.0% for early-type galaxies, 21.2% for brightest cluster galaxies and 28.2% for the combined clusters. Disk speed error is 15.26%, still worse than the separate 9.22% galaxy-only calibration. Excluding every X-COP host gives 0.1432 dex RMS and 29.3% median acceleration error there; excluding every CLASH host gives 0.2158 dex and 32.6%. The earlier power sum gives 0.2480 and 0.2481 dex on those same source exclusions. Common RAR still performs better on excluded X-COP, at 0.1096 dex and 17.5% median acceleration error. Thus improved transfer relative to the earlier power sum is not dominance over every comparator. These conditional transfers retain the different source reconstruction assumptions. Adding the hydrostatic sample does not make a common observable-level likelihood.

**The repaired limits are explicit.** At fixed positive radius, removing baryons gives g_pred proportional to g_bar^((1-alpha)/2), which vanishes for alpha below one. High baryonic acceleration returns both the relative and absolute Newtonian limits. In a smooth core the inverse leg now behaves as r^(4-3 alpha) and vanishes, leaving g_pred proportional to sqrt(r). The central force is finite and vanishes, although its inferred density can have an integrable cusp. For a fixed isolated point mass at large radius, a_eff tends to a positive constant and g_pred is proportional to 1/r. The supplement proves nonnegative point-source extra density for all allowed coefficients: each positive a_eff term has logarithmic radial slope above -2, so inferred extra enclosed mass increases with radius. This proof does not extend automatically to arbitrary baryonic geometry.

Illustrative isolated-Sun corrections at 0.4-30 AU are negligible for the final fit. At 1,000 AU its acceleration ratio is 1.00027, but at 10,000 AU it is 1.7865. The constant C0 contribution still produces a low-acceleration enhancement; stellar screening is partly imposed and no binary or ephemeris likelihood has been fitted. The asymptotic speed relation is V_infinity^4 = GM a_ref(C0+C1 U^-alpha); its low-mass normalization remains mass-dependent. Solar external fields, global baryonic Tully-Fisher data, lensing metrics and cosmology therefore remain necessary tests.

**Selection and host correction.** The three stages are sequential and adaptive: six initial families, two repairs after failures, and one final regularity repair. Their 566 fit sets use 1,320 optimizer starts; one unsuccessful start was not selected, and all selected fits converged. The 182 prediction tables are independently checked. Historical host sensitivities used exact catalog labels: one earlier Milky Way test was omitted, and other tests left one Milky Way and two Andromeda members in training. Those host interpretations are superseded by conservative exclusions of all 26 Milky Way or 18 Andromeda dwarfs. The original system/population scores are unaffected. The old files are retained as historical records with a correction notice.

### 5.7 Targeted MOND comparisons and an asymptotic repair

A further fixed local protocol compares seven families on the same expanded 286-system sample: exponential RAR; the published simple and standard MOND interpolations; their generalized shape family; and three finite-range/window repairs. All population/group weights, three system-fold seeds, complete host exclusions and source exclusions are retained. The known interpolation family is

$$\nu_n(y)=\left[\frac{1+\sqrt{1+4y^{-n}}}{2}\right]^{1/n},\qquad g=\nu_n(g_b/a_0)g_b.$$

Here n = 1 gives the simple form and n = 2 the standard form. Their strong-field limits differ: the simple function leaves an absolute excess approaching a0, and n below one makes that excess grow. The all-population fitted n is 0.784. Merely fitting this shape is therefore not a Solar-System repair.

| Shared model | Coeff. | System-held-out B | Population-held-out B |
| --- | ---: | ---: | ---: |
| Exponential RAR | 1 | 0.32991 | 0.35934 |
| Simple MOND | 1 | 0.32908 | 0.35861 |
| Standard MOND | 1 | 0.33679 | 0.36987 |
| Generalized nu_n | 2 | 0.32972 | 0.37809 |
| Section 5.6 regular response | 7 | 0.22197 | 0.28920 |
| Transient correction | 7 | 0.22197 | 0.29057 |
| Potential-taper correction | 7 | 0.22794 | 0.32860 |
| Two potential windows | 7 | 0.23283 | 0.37443 |

Scores are in dex. The transient correction makes the smallest useful mathematical change: both nonconstant terms in Section 5.6 now share its finite-range denominator,

$$\frac{a_{\mathrm{eff}}}{a_{\mathrm{ref}}}=C_0+\frac{C_1U^{-\alpha}S(r)W(g_b)+C_2s^{\beta}}{1+r/R},\qquad g=\frac{g_b}{1-\exp[-\sqrt{g_b/a_{\mathrm{eff}}}]}.$$

Using the earlier definitions of U, s, S and W, the all-data coefficients are C0 = 0.88098867, C1 = 0.00077482965, C2 = 0.40173445, alpha = beta = 0.95, ell = 10 pc and R = 307.65765 kpc. It remains one expression with seven shared coefficients and no population switches. Alpha, beta and the imposed ell minimum again contact bounds.

For an isolated point source its added effective-scale terms now vanish at infinity, giving a common asymptotic normalization V_inf^4 = GM a_ref C0. The previous mass-dependent normalization is removed. Empty-source, high-acceleration and vanishing smooth-core force limits are retained, and the point-source inferred extra mass remains nondecreasing. This establishes neither arbitrary-profile positivity nor a conservative field equation. The fitted asymptotic scale is 8.81 times 10^-11 m/s^2, 26.6% below the canonical 1.2 times 10^-10; at finite radii corrections persist. A one-million-solar-mass point-source probe has V^4/(GM a_ref C0) = 1.445 at 10 kpc. Thus an observed finite-radius baryonic Tully-Fisher fit has not been established.

The fit gain itself is negligible: system RMS changes by -0.00000723 dex, with descriptive paired interval [-0.00002656, 0.00000987]. Population transfer slightly worsens. This is a physical extrapolation repair, not a new accuracy result. The potential taper and two windows do not improve aggregate prediction; numerical probes find decreasing force response in the potential taper only in extreme extrapolation, while the two-window candidate fails this local derivative check at 37 observed input locations and can imply negative point-source extra density. All failed candidates remain supplied.

Disk-only recalibration under the same equal-galaxy weighting gives RAR speed error 8.89%, versus the transient model's cross-population disk error 15.26%. The earlier 9.22% disk reference used point weighting. Neither comparison supports superiority to MOND on disks. This stage supplies 79 prediction tables, 239 fit sets and 471 successful optimizer attempts, with independent prediction, split, metric and limit checks.

![Comparison of distinct prediction tasks](best_models_comparison.png)

Figure 5. The asymptotic repair leaves broad accuracy essentially unchanged. Established halo fits help predict radial structure within measured galaxies, a different task with more object-specific information. Lower bars are better; the panels must not be combined into one score.

## 6. Separate checks of groups, lensing and stellar controls

Galaxy groups test a global virial prediction rather than a local circular acceleration. Using the printed numerical table of Milgrom (2019), based on the Makarov-Karachentsev group catalog, we retain 55 source-selected groups from 66 printed entries. The article's prose counts differ by one; the frozen extraction documents the actual rows. Harmonic radii span 62 to 738 kpc. With the frozen SPARC a0 and stellar K-band mass-to-light ratio one, the deep-MOND expression

$$\sigma_{\mathrm{pred}}^4=\frac{4}{81}GM_{\mathrm{bar}}a_0$$

gives 18.85% median dispersion error and 0.1173 dex dispersion RMS. Hot gas is not added. The relation assumes the relevant deep-MOND, isolated virial regime; group membership, equilibrium and gas inventory remain limitations. This dispersion percentage cannot be pooled with acceleration or circular-speed errors.

For Brouwer et al. (2021), we use the official KiDS release with shear calibration and covariance. Its conversion g_obs = 4G DeltaSigma uses a singular-isothermal-sphere approximation. With no refit, the isolated GAMA profile's 15 bins give 25.86% median acceleration error and statistical chi-squared 18.54 for 15 bins. Seven KiDS bins with g_bar at least 10^-13 m/s^2 give 22.08% median error but statistical chi-squared 46.53. Adding an approximate independent diagonal 0.1 dex conversion uncertainty reduces the latter to 5.32; that is a sensitivity assumption, not a validated systematic covariance model. GAMA and KiDS samples overlap and are not independent confirmations. Baryonic-mass and hot-gas sensitivities are supplied.

The stacked bins lack a unique physical radius, so the new radius-conditioned formula is not evaluated by inventing one. A proper test requires returning to the lens galaxies and projected shear likelihood, with a specified mapping from the theory to lensing. The datasets are credited to KiDS and the authors; the original release README accompanies the frozen numerical inputs.

The Baumgardt-Hilker modeling catalog provides 151 selected globular-cluster controls from a current frozen table of 171 rows. Their projected half-light radii span 1.05 to 29.04 pc; median dynamical V-band mass-to-light ratio is 1.81, with 133 of 151 between one and three. These are model-based stellar-system controls. Central dispersions cannot be substituted for the global dispersion in the dwarf estimator, and dynamical mass cannot supply an independent baryonic predictor. We therefore report no invented RAR error for this catalog.

Together these checks show why radius alone is inadequate. Some galaxy groups comparable in size to cluster cores approximately follow a galaxy-calibrated virial relation, whereas small but massive cluster centers favor a larger acceleration normalization. Tiny dwarf spheroidals and small globular clusters have different inferred mass budgets. Geometry, baryonic inventory, environment, history and inference assumptions need testing alongside scale.

### 6.1 Repairing cluster MOND with an observed gas profile

Famaey, Pizzuti and Saltas (2025) propose an extra density that follows hot gas with an exponential cutoff. We implement that published density-level family rather than inventing an enclosed-mass multiplier:

$$\rho_{\mathrm{extra}}(r)=\eta\rho_{\mathrm{hot}}(r)\exp(-r/R_c),\qquad M_{\mathrm{extra}}(<r)=\eta\int_0^r\exp(-s/R_c)\,dM_{\mathrm{hot}}(s).$$

The spherical source acceleration is g_src = g_bar + G M_extra/r^2; the predicted acceleration is nu(g_src/a0)g_src. We hold a0 at 1.2 times 10^-10 m/s^2 and fit only two shared coefficients, eta and Rc. Both simple MOND and exponential RAR are tested. This explicitly adds phenomenological missing mass; it does not explain its physical identity.

The public gas-parameter release matches 16 CLASH hosts with 65 existing A43 points; seven X-COP hosts add 21 fixed-radius points. Four CLASH hosts lack the required gas parameters and are excluded by availability, not fit quality. X-COP uses direct enclosed MGAS, not gas fraction multiplied by the target mass. These gas and gravitational reconstructions still share observational assumptions and uncertainties. The combined objective balances the two sources, then hosts and radii. It is a separate 23-host test, not the 286-system objective.

| Simple-MOND comparison | Combined host-held-out RMS | CLASH to X-COP RMS | X-COP to CLASH RMS |
| --- | ---: | ---: | ---: |
| Fixed galaxy a0 | 0.44059 | 0.35811 | 0.50989 |
| One a0 fitted on training clusters | 0.14011 | 0.21436 | 0.21473 |
| Published illustrative eta=10, Rc=400 kpc | 0.11456 | 0.07913 | 0.14138 |
| Two shared gas coefficients refitted | 0.11331 | 0.06120 | 0.13985 |

Scores are in dex; every fitted coefficient excludes the test host or entire test source. CLASH-trained eta = 11.77 and Rc = 288.85 kpc give X-COP median acceleration error 7.63%, versus 53.05% for the CLASH-trained one-a0 comparator. The reverse transfer gives 16.59% versus 37.85%. Exponential interpolation produces nearly the same combined RMS, 0.11308 dex.

This is a promising transfer of an existing published family. The illustrative literature constants already give combined RMS 0.11456 dex, close to the refitted 0.11331. They were CLASH-informed, so their CLASH score is a descriptive reference, not independent source-held-out validation. Descriptive combined-host bootstrap intervals for the repair minus the one-a0 comparator [-0.06199, 0.00682] dex and minus the frozen literature reference [-0.01304, 0.01164] both include zero. Moreover, within CLASH alone the gas repair's host-held-out RMS is 0.14794, worse than the one-a0 comparator's 0.11298. Thus there is no uniform victory even inside this selected cluster benchmark. Nine X-COP target/stellar sensitivities, quadrature convergence and unresolved inner-gas bounds are supplied. The extra density is positive and has finite total mass for these gas profiles; its distant MOND normalization involves baryons plus that extra mass. No raw shear/temperature likelihood, merger evolution or causal gas-to-dark-component mechanism has been fitted. Missing hot-gas measurements in other populations were not set to zero; no all-scale result follows from this subset.

### 6.2 Can the successful cluster source correction also act on galactic gas?

We test a concrete attempted unification of Section 6.1. If the same extra component follows ordinary galactic gas with a spatially constant multiplier, linear Newtonian source integration gives

$$g_{\rm src}=g_\star+(1+\eta)g_{\rm gas},\qquad g_{\rm pred}=\nu(g_{\rm src}/a_0)g_{\rm src}.$$

This is exact for a constant density multiplier. It is only a proposed leading approximation to extending the cluster exponential profile into disks. A disk's force depends on gas beyond the evaluation radius, so small r/R_c alone does not bound the error of the finite-cutoff convolution. The existing rotation-component tables do not supply all required gas-density fields. We therefore test the constant-weight extension itself, not an exact finite-range disk field law.

The primary disk subset is selected entirely from predictors: positive gas acceleration, gas fraction g_gas/g_bar at least 0.5, and radius at most 30 kpc. It contains 61 galaxies and 402 points. Stars and measured cold gas retain their original SPARC components. The cluster calculation uses the same 23 hosts and exact hot-gas mass integrals as Section 6.1. The joint objective gives disks and clusters equal weight, with equal galaxies or source-balanced cluster hosts inside each domain. a0 is fixed at 1.2 times 10^-10 m/s^2.

| Simple-MOND gas comparison | Disk RMS, dex | Cluster RMS, dex |
| --- | ---: | ---: |
| Separate domain calibration, held out by system | 0.15121 | 0.11141 |
| Common eta and R_c, held out by system | 0.27236 | 0.26260 |
| Frozen CLASH eta transferred to disks | 0.63092 | Not this test |

The disk calibration prefers eta approximately zero, with 11.60% median speed error. Transferring CLASH eta = 11.77 raises that speed error to 99.65%. The joint all-data compromise has eta = 1.575 and R_c at the imposed 10,000-kpc upper bound; both domains deteriorate relative to separate calibration. Bootstrap intervals for those domain-specific RMS deteriorations exclude zero, conditional on these fitted predictions. Exponential RAR gives the same qualitative failure. Source exclusions and alternate galaxy splits are supplied.

The full 2,700-point disk diagnostic retains invalid cases: the frozen simple-MOND CLASH multiplier makes nine source accelerations nonpositive where the gas contribution points outward. Those rows are not clipped or silently dropped to report a whole-sample score. This further limits the algebraic extension. The result rejects naive uniform coupling to measured cold and hot gas within this test; it does not reject Famaey et al.'s original hot-gas-specific proposal, or a correctly solved finite-cutoff disk field. A distinction based on gas temperature, environment or formation could be investigated, but no such measured predictor or mechanism is established here.


## 7. Correcting the earlier cluster result

The original comparison used seven X-COP estimates at R500 and appeared to obtain 5.63% median cluster acceleration error with a common acceleration-and-radius quadratic. That radius is defined by the same reconstructed mass that supplies the target:

$$M_{500}=\frac{4\pi}{3}500\rho_{\mathrm{crit}}(z_c)R_{500}^3,\qquad g_{500}=250H(z_c)^2R_{500}.$$

Here z_c is cosmological redshift. A radius-only leave-one-cluster-out fit gives 0.8375% error; including the known redshift factor reduces this to 0.000144%. The defining identity with the catalog cosmology gives 0.000203%, essentially rounding. Holding out whole clusters does not break this target dependence. We withdraw the earlier 5.63% result as evidence of physical unification, while preserving its reproducible numerical calculation.

A fixed-physical-radius stress test uses A1795, A2029, A2142, A2319, A644, A85 and ZW1215 at 500, 750 and 1,000 kpc, within their tabulated domains. The primary total mass is the archived forward hydrostatic reconstruction, with NFW and Einasto alternatives. Gas plus 0.75 times projected stellar mass is an approximate middle baryon estimate. That factor is justified by the source at R500 only; zero-star and full-projected-star cases test sensitivity outside it, not confidence limits.

| Fixed-radius model | Test | Median acceleration error |
|---|---|---:|
| Galaxy-calibrated RAR | Frozen transfer | 56.16% |
| Earlier R500 radius quadratic | Frozen radial transfer | 39.02% |
| Cluster-calibrated RAR | Leave one cluster out | 15.98% |
| Joint galaxy-cluster quadratic | Leave one cluster out | 21.02% |
| Radius-only power law | Leave one cluster out | 17.96% |

The frozen R500 formula previously saw these same clusters, so its test concerns radial transfer, not new hosts. The joint quadratic trains on all galaxies and six clusters and is evaluated only on the omitted cluster. The 21 profile points are correlated mass-model outputs. The cluster RAR calibration on all 21 is 7.15914 times 10^-10 m/s^2. Its median error is 15.57-16.99% across stellar sensitivities and 13.61% or 17.29% with NFW or Einasto targets. Even in the primary case the radius-only power law has lower logarithmic RMS than cluster RAR, 0.09002 versus 0.09678 dex. No single summary statistic establishes a uniformly superior model.

## 8. Memory, geometry and the remaining universal test

The biological analogy is useful when stated precisely: common growth rules can generate different trees from different starting conditions, histories and environments. Likewise, a common dynamical theory can generate a family of galactic structures. The analogy does not require one arbitrary formula per object. A useful family must predict its variation from independently measured conditions or a generative model, and must earn its complexity on held-out data.

History already matters in ordinary dynamics. Pontzen and Governato (2012) model how changing baryonic potentials redistribute collisionless matter; Read, Walker and Steger (2019) relate star-formation histories to inferred dwarf central densities. A new temporal dependence in the gravitational law is a separate hypothesis. Present-day residuals cannot be used to define each object's memory and then counted as successful predictions.

One explicit research design would augment a baseline with a dimensionless state m, relaxation time tau > 0, and independently specified driving history s:

$$\tau\frac{dm}{dt}+m=s(t),\qquad m(t)=e^{-T/\tau}m(t-T)+\frac{1}{\tau}\int_{t-T}^{t}e^{-(t-t')/\tau}s(t')\,dt'.$$

A bounded coefficient multiplying m could alter log10(a_eff/a_ref), with high-acceleration suppression specified in advance. The driver might be a normalized star-formation, baryonic-potential or interaction history, but its meaning, timescale, initial state and measurement uncertainties must be fixed before fitting. Several relaxation modes can approximate more general causal kernels. Complete time-resolved driving histories are not in the joint table. Section 8.1 tests a photometric quenching-age proxy instead; it does not infer this dynamical memory kernel.

A discrete update rule can produce this kind of effective state: m_(n+1) = (1-epsilon)m_n + epsilon s_n, for 0 < epsilon < 1, retains a weighted past. This is a numerical illustration, not a spacetime derivation. Collective long timescales can emerge from many local updates, while a smooth coarse-grained equation can still describe the result. Failure to find a simple curve would not prove discreteness; the present fits actually demonstrate substantial smooth regularity.

A18's geometric/jet and harmonic proposals and A41-A42's discrete-rule program can use these datasets as targets. A credible connection must derive the observable acceleration, effective states and parameter relations from its specified geometry or update rule. Additional harmonic modes or dimensions cannot be assigned after observing every residual without a predictive penalty. This paper supplies empirical targets, not that derivation.

| Remaining domain | Required prediction |
|---|---|
| Solar-System tracking and wide binaries | Vector dynamics, Galactic external field, orbital populations and observational likelihoods |
| Lensing and merging clusters | Both metric potentials, projected maps, baryonic geometry and time evolution |
| CMB anisotropies | Expansion, photon-baryon dynamics, perturbations, recombination and lensing transfer |
| Large-scale structure and BAO | Distances, sound horizon, scale-dependent growth and matter power with covariance |

Cassini constraints studied by Desmond, Hees and Famaey (2024) show why an isolated Newtonian high-acceleration limit alone is insufficient. Planck's inferred cold-dark-matter density is a parameter of a cosmological model, not another local acceleration point. DESI distances and clustering likewise require dynamical predictions. Skordis and Zlosnik (2021) illustrate a proposed relativistic MOND completion with calculated CMB and matter spectra; an arbitrary empirical RAR extension does not inherit those results. The supplied source assessment describes these tests and the unsettled wide-binary literature. Their likelihoods remain unfit here.

Independently measured environment is a motivated extension to test before another flexible radius correction. Bilek, Renaud and Samurovic (2026) report distinct acceleration behavior for central versus noncentral and isolated early-type galaxies, with the transition associated with host-group mass. This motivates testing external baryonic potential, host mass, central/satellite status and baryonic geometry under identical observable-level likelihoods. It does not establish that any one correction works here. Their analysis revisits earlier galaxy profiles, so it must not automatically be treated as an independent validation set. Host properties inferred from the same dynamical target also require an explicit leakage audit.

### 8.1 Independent environment and stellar-history pilots

We obtained predictors independently of the acceleration residuals and fixed a local comparison protocol before running the new fits. These are subset tests on previously inspected acceleration data, not an independent discovery sample or a history-corrected law for all five populations. Missing predictors were not set to zero. The analysis performed 2,314 fits and supplies 49 held-out prediction tables plus sensitivity-score arrays; independent code checked the equations, metrics and group separation.

**External fields in 91 rotating galaxies.** Chae et al. (2021), Table 3, supplies environmental Newtonian-field estimates from surrounding baryonic matter, with maximum-clustering (MC) and no-clustering (NC) assumptions. The matched sample contains 1,347 radial points in 91 of our selected SPARC galaxies. We use only those environmental estimates, not the external fields fitted to rotation curves in their Table 2. These are model-dependent reconstructed predictors, not direct force measurements.

For positive external field ge, their idealized one-dimensional equation can be written

$$g=F(g_{\mathrm{bar}}+g_e;a_0)-F(g_e;a_0),$$

$$F(z;a_0)=\frac{z}{2}+\sqrt{\frac{z^2}{4}+a_0z}.$$

We compare it with the same simple interpolation at ge = 0 and with the standard RAR. Each fits one common a0; a separate two-coefficient variant also fits one global multiplier of ge. The published dimensionless environmental fields are converted with their fixed normalization 1.2 times 10^-10 m/s^2, independently of our fitted a0. The model is prior art, not a new environmental force law.

| Model on the same 91 galaxies | Galaxy-held-out RMS | Distance-band-held-out RMS |
| --- | ---: | ---: |
| Standard RAR | 0.14057 dex | 0.14151 dex |
| Simple interpolation, no external field | 0.14098 dex | 0.14196 dex |
| Published MC environmental field | 0.14512 dex | 0.14772 dex |
| Published NC environmental field | 0.14048 dex | 0.14195 dex |
| MC field with fitted global multiplier | 0.14158 dex | 0.14184 dex |

All radii of a galaxy remain together. Equal galaxy weights prevent long curves from dominating. Two additional random fold assignments give the same qualitative result. A stricter distance-transfer check withholds whole bands separated at 10, 20, 40 and 80 Mpc. Those bands contain 27, 42, 7, 6 and 9 galaxies; they are not verified host memberships, and spatially correlated errors can remain. The minute NC gain in the original galaxy folds does not persist in the distance-band score. Allowing a field multiplier does not produce a consistent improvement.

One hundred per-galaxy field-error perturbations, shared across all radii and followed by refitting, give 5th-95th percentile distance-band RMS ranges of 0.1466-0.1579 dex for MC and 0.1404-0.1456 dex for NC. Spatial error covariance and other nuisance parameters are not included. This fixed-baryon benchmark does not reproduce or refute Chae et al.'s full inference, which treats galaxy nuisance parameters differently.

**Stellar-history test in 15 dwarfs.** Weisz et al. (2015), Table 1, provides the lookback epoch t90 when 90% of stellar mass had formed, reconstructed from Hubble color-magnitude diagrams. Fifteen objects match our dwarf selection: eleven Milky Way satellites and four Andromeda satellites. We transcribed t90, its total asymmetric errors and published host distances, and checked the table visually. The ages use stellar photometry, not the velocity dispersions that supply the acceleration targets. A quenching epoch is only a summary of history and does not reconstruct earlier gravitational potentials.

We test the deliberately simple proxy

$$m_\tau=\exp(-t_{90}/\tau),\qquad \tau=5\ {\rm Gyr},$$

$$\log_{10}(a_{\mathrm{eff}}/10^{-10}{\rm\ m\ s^{-2}})=a+b\,m_\tau,$$

with the standard RAR response. The five-Gyr scale is fixed, not inferred; two and ten Gyr are separately declared sensitivities. This is a history-dependent empirical predictor, not a fitted causal memory integral. The all-data coefficients are a = 1.22335 and b = -1.85925. Validation refits them on one host and tests all matched satellites of the other, then reverses the split. The history coefficients are negative in both host fits, but their magnitudes differ substantially.

| Dwarf model | Coefficients | Host-held-out RMS |
| --- | ---: | ---: |
| Constant-scale RAR | 1 | 0.23598 dex |
| RAR plus five-Gyr history proxy | 2 | 0.20259 dex |
| Two-Gyr history sensitivity | 2 | 0.19633 dex |
| Ten-Gyr history sensitivity | 2 | 0.20851 dex |
| RAR plus log host distance | 2 | 0.28445 dex |
| History plus host distance | 3 | 0.21139 dex |
| RAR plus log stellar mass | 2 | 0.38745 dex |
| History plus stellar mass | 3 | 0.38958 dex |

This score weights the two hosts equally, then their satellites equally. The five-Gyr point-estimate improvement is 14.1% in RMS, with pooled median acceleration error falling from 44.1% to 35.5%. It is not shared equally: the Milky Way RMS falls from 0.2649 to 0.2011 dex, whereas Andromeda changes from 0.2029 to 0.2040 dex. The two-Gyr fit reaches a coefficient bound. The stellar-mass control transfers poorly despite giving the best leave-one-dwarf-out score, 0.2001 dex. This illustrates why individual-object validation can miss failure across hosts; four training satellites cannot reliably separate several correlated effects.

For context, refitting the original broad power-sum family after excluding every satellite of the target host removes 26 Milky Way or 18 Andromeda dwarfs from its 2,910-row input. On the same 15 history matches it gives 0.2775 dex RMS but a slightly smaller pooled median acceleration error, 35.4%, than the history model. The training populations and objectives differ; neither statistic establishes universal superiority.

Crucially, 500 independent asymmetric-error perturbations of the published log ages, capped at 13.8 Gyr and followed by refitting both host folds, give a median history-model RMS of 0.2617 dex and a 5th-95th percentile range of 0.1770-0.4400 dex. Only 37.2% of those trials improve on the constant-scale baseline. These are conditional sensitivity draws, not posterior probabilities; shared stellar-model systematics, distance errors and dynamical-estimator uncertainty are unmodeled. The central-estimate gain is therefore not robust to the available age uncertainties.

The result identifies a testable history association but supplies no detected modification of gravity. A dependence on stellar history could arise from ordinary baryonic feedback, changing stellar mass-to-light ratios, environmental processing or modeling errors. Better resolved histories in additional independent hosts, combined with observable-level uncertainty inference and halo/feedback comparisons, are needed before promoting the proxy into the common five-population formula.


### 8.2 From one radial response to one field law

A shared law in the sense of general relativity is a rule for fields and matter with common constants, completed by source, initial and boundary data. Its different solutions need not collapse onto one pointwise curve of g_bar and r. The biological analogy in Section 8 can therefore be made precise: variation should follow from the same dynamics acting on measured state, rather than from a new coefficient for each object.

Our additional diagnostic uses one flexible Gaussian-kernel correction in the two existing inputs. Three inner folds select its bandwidth from four fixed choices. Its whole-system score is 0.2349 dex, but its whole-population score deteriorates to 0.4060 dex, versus 0.3671 for common RAR. The correction has training-sample complexity despite its compact written expression; it is not a few-parameter theory. This test does not set a mathematical ceiling on other formulas, but shows that added interpolation flexibility alone has not solved transfer.

An input-only nearest-neighbor check makes the missing information concrete. CLASH MACS1206 and the disk UGC09037 have g_bar values 2.871 and 2.948 times 10^-11 m/s^2, at radii 14.30 and 14.55 kpc. Their inferred accelerations differ by a factor 9.63. This is not a noise-free contradiction: the inputs differ, geometry and environment differ, and lensing-model and rotation-curve targets carry different assumptions. It motivates full-source geometry and uncertainty treatment. A one-law theory can allow different outcomes here without assigning arbitrary labels to the two systems.

A basic field consistency condition also matters. If the Newtonian vector field is g_N = -grad(Phi_N), a naive rescaling has

$$\nabla\times[\nu\mathbf{g}_N]=\nabla\nu\times\mathbf{g}_N.$$

It is generally nonzero, so the rescaled field need not be the gradient of a scalar potential. An independently checked quadratic-potential example in the supplement verifies this directly. QUMOND (Milgrom, 2010) instead solves

$$\nabla^2\Phi_N=4\pi G\rho_b,\qquad \nabla^2\Phi=\nabla\cdot[\nu(|\nabla\Phi_N|/a_0)\nabla\Phi_N].$$

These common field equations incorporate source geometry and boundary conditions. They do not guarantee that a selected nu fits every dataset, or give a radial fit an automatic relativistic completion.

**A concrete derivative/jet route.** Milgrom's (2023) generalized QUMOND supplies prior art for an action depending on higher Newtonian-potential derivatives. Write H_ij = partial_i partial_j Phi_N and choose dimensionless invariants

$$X=\frac{|\nabla\Phi_N|^2}{a_0^2},\qquad Y=\frac{\ell_*^2H_{ij}H_{ij}}{a_0^2}.$$

A specialization of that established action class is

$$\mathcal{L}_g=-\frac{1}{8\pi G}[2\nabla\Phi\cdot\nabla\Phi_N-a_0^2Q(X,Y)].$$

Including the usual matter coupling and varying both potentials gives the Newtonian source equation above and

$$\nabla^2\Phi=\partial_i(Q_X\partial_i\Phi_N)-\ell_*^2\partial_i\partial_j(Q_YH_{ij}).$$

The Hessian term is essential whenever Q depends on Y. Our independent variation check verifies its sign and normalization. Local derivative proxies ell_J = sqrt(6)|grad(Phi_N)|/sqrt(H_ij H_ij) and M_J = 6|grad(Phi_N)|^3/(G H_ij H_ij) reduce to r and M for an isolated point mass. For general geometry they are proxies, and degenerate fields require regularity analysis. These constructions connect derivative-based geometry to explicit testable equations; neither the class nor this point-mass normalization establishes new physics. A43's radial coefficients do not fit Q, and A18's geometric ideas have not derived it.

Park et al. (2026) report the updated Cassini external-field quadrupole constraint Q2 = (1.6 +/- 1.8) times 10^-27 s^-2. Their tested classical modified-gravity families face galaxy/Solar-System tension. An isolated high-acceleration limit does not calculate that quadrupole. A decisive next step therefore requires full baryonic density fields, a specified action function, converged field solutions and predictions of observed velocities, dispersions and both lensing potentials. The present table does not contain those complete source fields. The new formulas are empirical constraints for such a construction, not a completed substitute for it.

### 8.3 One law can require source state and history

Different fitted halo parameters do not imply different fundamental laws. Already in ordinary collisionless dynamics, the shared evolution of a mass phase-space density f obeys

$$\partial_t f+\mathbf{v}\cdot\nabla_{\mathbf{x}}f-\nabla\Phi\cdot\nabla_{\mathbf{v}}f=0,\qquad \nabla^2\Phi=4\pi G\left(\rho_b+\int f\,d^3v\right).$$

This nonrelativistic example needs an initial distribution and evolving baryonic sources. Formation and interactions can leave different present-day structures under the same dynamics. It is not an A43 discovery or a complete cosmological calculation. The radial test illustrates the practical distinction: information inferred from one part of a galaxy can predict another part, whereas a local g_bar,r pair discards much of that state. A useful next unification must predict or independently constrain that state, not assign each residual an arbitrary memory coefficient.

History-dependent alternatives also have specific prior art. Milgrom (2023, modified inertia) discusses trajectory/frequency-dependent responses with external-field predictions different from AQUAL/QUMOND, including a weaker inner-Solar Galactic response in the models considered. Such a route requires orbital histories and an explicit dynamics, not an age-only correction.

Berkimbayev (2026) proposes a stationary projected-density response as a possible endpoint of deformation memory. Its reported galaxy-held-out comparison does not significantly outperform RAR. The paper explicitly lacks a universal multi-source normalization rule and finds its unscreened Solar extrapolation unacceptable. This is relevant prior art, not evidence that memory has solved the cross-scale problem. A43 has not reproduced that full pipeline or measured a temporal kernel. Alongside history, independently measured source geometry and hot-gas structure remain concrete inputs to test.

### 8.4 Fixed NLG and MOG kernels cannot meet the inferred dwarf mass budgets

We also begin with published field-level alternatives. In the stationary astronomical weak-field approximation, Rahvar and Mashhoon's nonlocal gravity (NLG) and constant-parameter MOG can be represented by positive spatial kernels:

$$\nabla^2\Phi=4\pi G(\rho_b+\rho_{\rm eff}),\qquad \rho_{\rm eff}(\mathbf{x})=\int q(|\mathbf{x}-\mathbf{y}|)\rho_b(\mathbf{y})\,d^3y,$$

$$q_{\rm NLG}(d)=\frac{\alpha\mu(1+\mu d)e^{-\mu d}}{8\pi d^2},\qquad q_{\rm MOG}(d)=\frac{\alpha\mu^2e^{-\mu d}}{4\pi d}.$$

Both integrate to alpha over all space, so a finite isolated source generates total effective mass alpha M_b. For any positive spherical baryonic source at its true baryonic half-mass radius, this immediately implies the necessary bound

$$\frac{g}{g_b}=1+\frac{M_{\rm eff}(<r_{1/2})}{M_b/2}\ \leq\ 1+2\alpha.$$

This elementary mass budget does not depend on the range mu^-1 or the detailed radial profile. It applies to the stated isolated, constant-coefficient kernels; external sources, nonlinear response or scale-dependent coefficients need separate treatment. For a radially decreasing Plummer source, a translated ball cannot enclose more mass than the centered half-mass ball. The stronger bound is then g/g_b <= 1+alpha. Neither bound is claimed as a new fundamental result.

The published frozen values (alpha, mu in kpc^-1) are (10.94, 0.059) for NLG and (8.89, 0.042) for MOG. The general half-mass ceilings are 22.88 and 18.78. At nominal stellar M/L_V = 2, the existing 44-dwarf inferred acceleration ratios exceed these ceilings in 31 and 33 objects, respectively. At M/L_V = 3, 29 objects still exceed each ceiling. These counts concern model-based central estimates, not significance levels: the stellar mass, half-mass conversion, membership, binary stars, anisotropy, tidal disturbance and equilibrium assumptions remain material.

For a more specific prediction, we integrate each kernel against a spherical Plummer baryonic source with total stellar mass and scale fixed by the photometric inputs. Its three-dimensional half-mass radius is the radius used in the earlier dwarf estimator. The convolution includes exterior baryonic matter; it is not a point-mass substitution. We compare frozen coefficients with two common coefficients refitted after excluding each test dwarf, and with RAR on exactly the same 44 objects. M/L values one and three, two alternate fold assignments and complete Milky Way/Andromeda host exclusions are retained.

| Model, nominal M/L_V = 2 | Galaxy-held-out acceleration RMS, dex |
| --- | ---: |
| Published fixed-a0 RAR | 0.65085 |
| One common a0 refitted on dwarfs | 0.54844 |
| Published fixed NLG kernel | 1.86170 |
| Published fixed MOG kernel | 1.90516 |
| Two common NLG coefficients refitted | 0.72885 |
| Two common MOG coefficients refitted | 0.72885 |

All 108 selected kernel fit sets reach the upper mu bound, 10^5 kpc^-1. This collapses the response toward a short-range kernel and an approximately constant enhancement, rather than identifying a useful nonlocal length. The all-data nominal fit has alpha approximately 58.4. Even this freedom does not beat the matched dwarf RAR comparator; complete host transfers remain poor. RAR itself has large dwarf errors, and the refitted-a0 RMS improvement does not improve its pooled median relative error. No model in this table solves the dwarf problem.

Related tensions are already published. Roshan and Mashhoon (2022) explicitly discuss the relatively small effective dark-matter content predicted for dwarf galaxies in NLG. Thus the broad dwarf-deficit mechanism is prior art; the contribution here is its explicit half-mass application and reproducible 44-object test.  Haghi and Amiri (2016) find no common MOG parameter pair for eight classical dwarfs with plausible stellar masses. De Martino et al. (2025) use the same NLG kernel and Plummer sources in full Jeans-profile analyses, including a joint shared-kernel fit with galaxy-specific stellar and anisotropy parameters. That richer inference can fit profiles while leaving parameter tensions. Our fixed-input half-mass test does not supersede it. Nor does it test running-parameter MOG, a full covariant theory, or a causal time-memory kernel. Solar applications of NLG retain an additional short-distance kernel parameter (Chicone and Mashhoon, 2016), absent from this astronomical approximation. Its specific contribution is an auditable transfer stress and necessary mass-budget condition for proposed universal linear kernels.

![Tests of shared halo relations, gas coupling and nonlocal kernels](bridge_comparison.png)

Figure 6. Three distinct attempts to connect established successful descriptions. Cosmologically constrained halo relations approach RAR accuracy with common measured inputs. Uniform cold/hot gas coupling and the tested fixed linear kernels fail their declared transfers. Panel A shows differences from RAR, with descriptive 95% paired galaxy-bootstrap intervals conditional on fitted validation predictions. Panels B-D show absolute logarithmic RMS on their own samples; C and D separately evaluate cold disks and hot clusters under simple MOND. Each panel uses its own sample and validation task; scores must not be pooled into an all-scale ranking.


## 9. Conclusions

**Stellar mass improves disk-to-dwarf prediction at fixed model complexity.** Replacing the total baryonic-mass power in a three-coefficient NFW mapping by stellar mass lowers frozen dwarf RMS by 28.2%, from 0.66757 to 0.47914 dex. Average galaxy-held-out RMS changes from 0.15323 to 0.15329 dex. The RAR/stellar-halo blend lowers dwarf RMS by 24.2% relative to the previous blend while preserving average disk accuracy. These are improvements over explicitly matched shared predictors, with the same coefficient counts. Their gains are concentrated in ultrafaints; classical dwarfs and several disk tails expose remaining failures. The adaptive search makes an uninspected external sample the decisive next test.

**A finite positive kernel has a finite missing-mass budget.** For the stated isolated spherical source, kernel normalization alpha gives total effective mass alpha M_b and the necessary half-mass bound g/g_bar <= 1+2 alpha. At the published NLG and constant-parameter MOG normalizations, 31 and 33 of 44 central dwarf estimates exceed their respective ceilings. Raising stellar M/L_V from two to three leaves 29 exceedances for each. No adjustment of kernel range can remove these violations while keeping the other assumptions fixed. The result specifies what a repair must change: source mass, environment/dynamical reconstruction, normalization, or the linear response itself. It tests the quoted prescriptions; the counts are conditional inconsistencies, not statistical exclusions of complete theories. Earlier dwarf tensions and effective-mass restrictions are prior art.

**Successful fits do not automatically transfer between populations.** A hot-gas-associated cluster correction fails when applied uniformly to cold galactic gas. Flexible cross-population formulas improve some aggregate scores while sacrificing disk accuracy or excluded-population prediction. Population-scatter modeling improves the tested predictive log score even when its central disk RMS worsens. These comparisons separate accurate central prediction, calibrated uncertainty and physical consistency, which require different tests.

The resulting research target is concrete: predict halo or response state from independently measured stellar content, geometry, gas phase and formation history, then test that state on uninspected velocity and dispersion profiles. The shared stellar relation supplies a better starting predictor; the kernel bound supplies a necessary condition for alternatives. Full dynamical, lensing and cosmological likelihoods are the further tests needed for an all-scale theory.

## Data and computational availability

The supplement contains frozen numerical inputs, selection/exclusion tables, provenance and hashes, all coefficients and fold predictions, numerical scripts, sensitivity results and reproduction instructions. It credits the originating SPARC, nearby-galaxy, X-COP, CLASH, MaNGA, KiDS/GAMA and globular-cluster studies, plus the Chae environmental and Weisz stellar-history catalogs. This work collects no new astronomical observations. Article PDFs and installed software environments are not redistributed. Regression scores are conditional estimates. Section 5.5 adds descriptive group-resampling intervals; reconstruction, refitting and model-selection uncertainty remain unquantified. The dm_physical_repair supplement supplies its protocol, all three attempted candidates, coefficients, predictions, bounds and bootstrap results. The dm_environment supplement adds source provenance, matched predictors, the environment/history protocol, main-fit coefficients, held-out predictions, sensitivity scores and independent numerical checks. The dm_unified_limits, dm_unified_repair and dm_unified_regular supplements add all three adaptive stages, frozen expanded inputs, 182 prediction tables, mathematical diagnostics and bootstrap comparisons. The dm_identifiability, dm_unified_stress and dm_unified_audit directories contain the flexible-response, nearest-input, source-transfer and independent verification records. Historical host-label omissions are identified in HOST_HOLDOUT_CORRECTION.txt; corrected tests supersede those earlier sensitivities.

The current benchmark extension adds dm_mond_repairs (79 prediction tables), dm_halo_benchmark (20 tables), dm_cluster_repair (192 model/case/split tables in one CSV) and dm_best_audit. The three analyses together retain 1,953 fit sets and 5,843 optimizer starts. Protocols were fixed before each local run on repeatedly inspected data, not before the overall exploratory project. The cluster supplement freezes extracted gas parameters, source URLs/hashes and direct gas profiles; downloaded article PDFs and third-party software are not redistributed. Independent numerical audits and their reproducible scripts are supplied. The literature comparison notes are in dm_best_literature.

The 29 October bridge extension adds dm_shared_halo, dm_cosmological_halo, dm_gas_bridge, dm_field_bridge and dm_bridge_audit. Together these retain 284 new logical prediction tables, 158,920 predictions, 426 fit sets and 1,362 optimizer attempts, plus 32,400 full-disk gas diagnostics. Comparator copies reused by the cosmological branch are not counted again. Independent scripts reproduce formulas, source integrals, weights, exclusions and numerical summaries. New protocols were frozen before their local runs on already inspected datasets; these remain exploratory comparisons. The supplied manifests preserve exact inputs, coefficients and all failed alternatives.

The feedback extension adds dm_feedback_profiles: 20 evaluation cases, 120 fitted predictor sets, 500 recorded optimization/endpoint candidates, 18 complete cross-validation tables and a six-model frozen dwarf-transfer table. Two attempts exhausted their evaluation budget; neither supplied a selected fit. All selected base fits converged. Independent checks reproduce baseline predictions to 4.7 times 10^-8 dex, recompute every complete-set score, verify galaxy exclusions, and compare 596 fitted-profile masses with direct quadrature to 1.4 times 10^-14 relative error. The preceding dm_fallback_validation preserves the independent JavaScript calculation, its distinct partitions, recorded outputs, selection data and provenance. Both branches use already inspected SPARC data and are exploratory.

The subsequent dm_population_scatter and dm_stellar_halo_curve directories add restored observational errors, linearized marginal distributions, a nonlinear-moment sensitivity, concentration and stellar-mass repairs, and 39 complete cross-validation prediction tables. Their records retain 620 numerical optimizer runs, 60 exact comparator endpoints and 20 additional scalar blend-weight fits. No selected or unselected numerical optimizer run failed. Independent scripts check derivatives, direct covariance calculations, fold exclusions, score arithmetic, source matching and the inverse published mass relation. Dwarf subset sensitivities preserve the original quality cut and conservative host variants. The stellar-law and blend stages were selected adaptively after inspecting prior outcomes; their conditional bootstrap intervals do not remove that selection uncertainty.

## References

Baumgardt, H., and Hilker, M. (2018). A catalogue of masses, structural parameters and velocity dispersion profiles of 112 Milky Way globular clusters. Monthly Notices of the Royal Astronomical Society 478, 1520. https://doi.org/10.1093/mnras/sty1057 . Updated author catalog: https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/parameter.html

Berkimbayev, D. (2026). A deformation memory model for galactic rotation curves. Physics Letters B 881, 140925. https://doi.org/10.1016/j.physletb.2026.140925

Bilek, M., Renaud, F., and Samurovic, S. (2026). Deviations from the radial acceleration relation in the central galaxies of clusters, subclusters, and groups. Astronomy & Astrophysics 708, A258. https://doi.org/10.1051/0004-6361/202558664

Brouwer, M. M., et al. (2021). The Weak Lensing Radial Acceleration Relation: Constraining Modified Gravity and Cold Dark Matter theories with KiDS-1000. Astronomy & Astrophysics 650, A113. https://doi.org/10.1051/0004-6361/202040108 . Numerical release: https://kids.strw.leidenuniv.nl/sciencedata.php

Burkert, A. (1995). The structure of dark matter halos in dwarf galaxies. Astrophysical Journal Letters 447, L25. https://arxiv.org/abs/astro-ph/9504041

Chae, K.-H., Desmond, H., Lelli, F., McGaugh, S. S., and Schombert, J. M. (2021). Testing the Strong Equivalence Principle. II. Relating the External Field Effect in Galaxy Rotation Curves to the Large-Scale Structure of the Universe. Astrophysical Journal 921, 104. https://doi.org/10.3847/1538-4357/ac1bba

Chicone, C., and Mashhoon, B. (2016). Nonlocal Gravity in the Solar System. Classical and Quantum Gravity 33, 075005. https://doi.org/10.1088/0264-9381/33/7/075005

Dantas, C. C., Ribeiro, A. L. B., and Capelato, H. V. (2026). The Radial Acceleration Relation in Galaxies and Clusters from a Two-Component, Virial-Motivated Framework. https://arxiv.org/abs/2607.29314 . Full paper inspected; quantitative pipeline not reproduced here.

De Martino, I., Della Monica, R., and Roshan, M. (2025). Dwarf galaxies in non-local gravity. Physics of the Dark Universe 49, 102044. https://doi.org/10.1016/j.dark.2025.102044

DESI Collaboration (2025). DESI 2024 VII: Cosmological Constraints from the Full-Shape Modeling of Clustering Measurements. Journal of Cosmology and Astroparticle Physics 07, 028. https://arxiv.org/abs/2411.12022

DESI Collaboration (2025). DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints. https://arxiv.org/abs/2503.14738

Desmond, H. (2023). The underlying radial acceleration relation. Monthly Notices of the Royal Astronomical Society 526, 3342-3351. https://doi.org/10.1093/mnras/stad2762

Desmond, H., Bartlett, D. J., and Ferreira, P. G. (2023). On the functional form of the radial acceleration relation. Monthly Notices of the Royal Astronomical Society 521, 1817-1831. https://doi.org/10.1093/mnras/stad597

Desmond, H., Hees, A., and Famaey, B. (2024). On the tension between the Radial Acceleration Relation and Solar System quadrupole in modified gravity MOND. Monthly Notices of the Royal Astronomical Society 530, 1781. https://arxiv.org/abs/2401.04796

Di Cintio, A., Brook, C. B., Dutton, A. A., Maccio, A. V., Stinson, G. S., and Knebe, A. (2014). A mass-dependent density profile for dark matter haloes including the influence of galaxy formation. Monthly Notices of the Royal Astronomical Society. https://doi.org/10.1093/mnras/stu729 . Recipe: https://arxiv.org/abs/1404.5959

Dutton, A. A., and Maccio, A. V. (2014). Cold dark matter haloes in the Planck era: evolution of structural parameters for Einasto and NFW profiles. MNRAS 441, 3359. https://doi.org/10.1093/mnras/stu742

Eckert, D., Ettori, S., Pointecouteau, E., van der Burg, R. F. J., and Loubser, S. I. (2022). The gravitational field of X-COP galaxy clusters. Astronomy & Astrophysics 662, A123. https://doi.org/10.1051/0004-6361/202142507

Ettori, S., et al. (2019). Hydrostatic mass profiles in X-COP galaxy clusters. Astronomy & Astrophysics 621, A39. https://doi.org/10.1051/0004-6361/201833323

Famaey, B., Pizzuti, L., and Saltas, I. D. (2025). Nature of the missing mass of galaxy clusters in MOND: The view from gravitational lensing. Physical Review D 111, 123042. https://doi.org/10.1103/dccw-srks . Gas parameters: https://zenodo.org/records/15299349

Ghizzardi, S., et al. (2021). Iron in X-COP: tracing enrichment in cluster outskirts with high accuracy abundance profiles. Astronomy & Astrophysics 646, A92. https://doi.org/10.1051/0004-6361/202038501

Haghi, H., and Amiri, V. (2016). Testing modified gravity with dwarf spheroidal galaxies. MNRAS 463, 1944. https://doi.org/10.1093/mnras/stw2140

Lelli, F., McGaugh, S. S., and Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. Astronomical Journal 152, 157. https://doi.org/10.3847/0004-6256/152/6/157

Lelli, F., McGaugh, S. S., Schombert, J. M., and Pawlowski, M. S. (2017). One Law To Rule Them All: The Radial Acceleration Relation of Galaxies. Astrophysical Journal 836, 152. https://doi.org/10.3847/1538-4357/836/2/152

Li, P., et al. (2019). A constant characteristic volume density of dark matter haloes from SPARC rotation curve fits. MNRAS 482, 5106. https://doi.org/10.1093/mnras/sty2968

Li, P., Lelli, F., McGaugh, S. S., and Schombert, J. M. (2018). Fitting the radial acceleration relation to individual SPARC galaxies. Astronomy & Astrophysics 615, A3. https://doi.org/10.1051/0004-6361/201732547

Li, P., Lelli, F., McGaugh, S. S., and Schombert, J. M. (2020). A comprehensive catalog of dark matter halo models for SPARC galaxies. Astrophysical Journal Supplement Series 247, 31. https://doi.org/10.3847/1538-4365/ab700e

Louis, T., et al. (ACT Collaboration, 2025). The Atacama Cosmology Telescope: DR6 Power Spectra, Likelihoods and LambdaCDM Parameters. Journal of Cosmology and Astroparticle Physics 11, 062. https://doi.org/10.1088/1475-7516/2025/11/062

McConnachie, A. W. (2012). The Observed Properties of Dwarf Galaxies in and around the Local Group. Astronomical Journal 144, 4. https://doi.org/10.1088/0004-6256/144/1/4 . January 2021 author update: https://www.cadc-ccda.hia.nrc-cnrc.gc.ca/en/community/nearby/

McGaugh, S. S., Lelli, F., and Schombert, J. M. (2016). The Radial Acceleration Relation in Rotationally Supported Galaxies. Physical Review Letters 117, 201101. https://doi.org/10.1103/PhysRevLett.117.201101

Milgrom, M. (2010). Quasi-linear formulation of MOND. Monthly Notices of the Royal Astronomical Society 403, 886-895. https://doi.org/10.1111/j.1365-2966.2009.16184.x

Milgrom, M. (2019). MOND in galaxy groups: a superior sample. Physical Review D 99, 044041. https://doi.org/10.1103/PhysRevD.99.044041 . Table I supplies the observed group quantities used here.

Milgrom, M. (2023). Generalizations of quasilinear MOND (QUMOND). Physical Review D 108, 084005. https://doi.org/10.1103/PhysRevD.108.084005

Milgrom, M. (2023). MOND as manifestation of modified inertia. https://arxiv.org/abs/2310.14334

Moffat, J. W., and Rahvar, S. (2013). The MOG weak field approximation and observational test of galaxy rotation curves. MNRAS 436, 1439. https://doi.org/10.1093/mnras/stt1670

Moffat, J. W., and Rahvar, S. (2014). The MOG weak field approximation II. Observational test of Chandra X-ray Clusters. MNRAS 441, 3724. https://doi.org/10.1093/mnras/stu855

Moffat, J. W., and Toth, V. T. (2009). Fundamental parameter-free solutions in Modified Gravity. Classical and Quantum Gravity 26, 085002. https://doi.org/10.1088/0264-9381/26/8/085002

Moster, B. P., Naab, T., and White, S. D. M. (2013). Galactic star formation and accretion histories from matching galaxies to dark matter haloes. Monthly Notices of the Royal Astronomical Society. https://doi.org/10.1093/mnras/sts261 . https://arxiv.org/abs/1205.5807

Navarro, J. F., et al. (2017). The origin of the mass discrepancy-acceleration relation in LambdaCDM. MNRAS 471, 1841. https://doi.org/10.1093/mnras/stx1705

Navarro, J. F., Frenk, C. S., and White, S. D. M. (1997). A Universal Density Profile from Hierarchical Clustering. Astrophysical Journal 490, 493-508. https://doi.org/10.1086/304888

Park, R. S., Hees, A., Famaey, B., Desmond, H., and Durakovic, A. (2026). Improved constraints on modified Newtonian gravity from Cassini radio tracking data. Physical Review D 114, 024066. https://doi.org/10.1103/R7N8-KW38

Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. Astronomy & Astrophysics 641, A6. https://doi.org/10.1051/0004-6361/201833910

Pontzen, A., and Governato, F. (2012). How supernova feedback turns dark matter cusps into cores. Monthly Notices of the Royal Astronomical Society 421, 3464-3471. https://arxiv.org/abs/1106.0499

Rahvar, S., and Mashhoon, B. (2014). Observational tests of nonlocal gravity: Galaxy rotation curves and clusters of galaxies. Physical Review D 89, 104011. https://doi.org/10.1103/PhysRevD.89.104011

Read, J. I., Agertz, O., and Collins, M. L. M. (2016). Dark matter cores all the way down. Monthly Notices of the Royal Astronomical Society. https://doi.org/10.1093/mnras/stw713

Read, J. I., Walker, M. G., and Steger, P. (2019). Dark matter heats up in dwarf galaxies. Monthly Notices of the Royal Astronomical Society 484, 1401. https://doi.org/10.1093/mnras/sty3404

Santos, R., and Felizardo, M. (2026). Symbolic Discovery of a Non-Linear Acceleration Scaling Relation in Galaxy Rotation Data. Particles 9(3), 70. https://doi.org/10.3390/particles9030070

Skordis, C., and Zlosnik, T. (2021). A new relativistic theory for Modified Newtonian Dynamics. Physical Review Letters 127, 161302. https://doi.org/10.1103/PhysRevLett.127.161302

Tian, Y., Ko, C.-M., Li, P., McGaugh, S., and Poblete, S. L. (2024). A distinct radial acceleration relation across the brightest cluster galaxies and galaxy clusters. Astronomy & Astrophysics 684, A180. https://doi.org/10.1051/0004-6361/202347868

Tian, Y., Umetsu, K., Ko, C.-M., Donahue, M., and Chiu, I.-N. (2020). The Radial Acceleration Relation in CLASH Galaxy Clusters. Astrophysical Journal 896, 70. https://doi.org/10.3847/1538-4357/ab8e3d

Weisz, D. R., Dolphin, A. E., Skillman, E. D., Holtzman, J., Gilbert, K. M., Dalcanton, J. J., and Williams, B. F. (2015). The Star Formation Histories of Local Group Dwarf Galaxies. III. Characterizing Quenching in Low-Mass Galaxies. Astrophysical Journal 804, 136. https://doi.org/10.1088/0004-637X/804/2/136

Wolf, J., et al. (2010). Accurate masses for dispersion-supported galaxies. Monthly Notices of the Royal Astronomical Society 406, 1220-1237. https://doi.org/10.1111/j.1365-2966.2010.16753.x

Zhao, H., and Famaey, B. (2012). Unifying all mass discrepancies with one effective gravity law? Physical Review D 86, 067301. https://doi.org/10.1103/PhysRevD.86.067301

Roshan, M., and Mashhoon, B. (2022). Characteristics of Effective Dark Matter in Nonlocal Gravity. Astrophysical Journal 934, 9. https://doi.org/10.3847/1538-4357/ac75d6 . https://arxiv.org/abs/2201.12852
