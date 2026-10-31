A43 halo / MOND radial-holdout benchmark
=======================================
This is a restricted numerical baseline exercise, not a claim of best-in-literature performance or a global test of LCDM. All calculations use the prior A43 fixed selected SPARC points and initially fixed baryonic components. It answers whether models can predict a contiguous block of radii inside a galaxy whose other radii are observed. It cannot be compared directly to entirely unseen-system or population holdouts.

104 galaxies with at least 9 retained radii; 2470 eligible points. Each outer-third or middle-third test contains 788 points; 1682 eligible points remain for per-galaxy fits. Shared disk fits use 1912 SPARC training points (including the 45 ineligible-for-test galaxies); broad A43 fits use 2143 training rows. All 788 corresponding gobs targets are excluded from every training objective. The two tests are not independent data sets. Radius ordering is deterministic; no model selected which radii to predict.

Held-out equal-galaxy RMS log10 acceleration [dex]:
Model                                     Outer third   Middle third   Fitted parameters
Published exponential RAR, a0=1.2e-10       0.114936      0.120702        0
Common RAR refit to disk train             0.113463      0.124402        1 shared
Prior A43 regular refit to disk train      0.113600      0.125148        7 shared
Prior A43 regular refit to broad train     0.162879      0.183561        7 shared
RAR with independent a0 per galaxy         0.104819      0.056428        104 (not universal a0)
NFW + fixed baryons                        0.090056      0.052708        208
Burkert + fixed baryons                    0.109817      0.048176        208
Published-a0 RAR + stellar multiplier      0.095172      0.073374        104
NFW + stellar multiplier                  0.094826      0.049370        312
Burkert + stellar multiplier              0.119500      0.044934        312

The extra seven-parameter A43 shape does not beat simple common RAR in either matched disk-only objective. The global five-population compromise is substantially worse on these disk predictions. Conventional fitted halo shapes and the one-nuisance MOND/RAR sensitivity perform better here on many comparisons, especially for the omitted middle block. This is evidence against presenting A43 as already superior to established models; it is not an adjudication of MOND vs dark matter.

Outer NFW minus shared-RAR RMS is -0.023407 dex, descriptive paired galaxy bootstrap 95% interval [-0.048650, +0.002819]; middle difference -0.071693, interval [-0.089555,-0.054457]. These intervals condition on fitted predictions, do not refit, and omit measurement/systematic/model-selection uncertainty. They should not be treated as final significance tests. Middle Burkert minus NFW interval also spans zero [-0.011541,+0.002530]. Do not declare a unique best halo based on these two probes.

Important limitations:
- No cosmological halo mass/concentration prior, no baryonic-feedback prescription, no SIDM evolution. Thus these are NFW/Burkert shape fits, not full LCDM/SIDM predictions.
- NFW scale radius hits its imposed 1000-kpc upper bound in 41/104 outer fits and 24/104 middle fits; Burkert in 10 and 2. Extrapolated virial masses are not inferred. Halo parameter degeneracy and prior sensitivity are material.
- One common stellar multiplier per galaxy rescales baseline disk M/L=0.5 and bulge=0.7 together by [0.2,2]; gas is unchanged. No positivity adjustment was actually needed. Distance, inclination, gas corrections, separate bulge M/L and uncertainties remain fixed. This does not reproduce the strongest published nuisance-marginalized MOND pipeline.
- Minimize mean squared log acceleration, not chi-square of velocities; no measurement covariance. Shared disk fits equal-weight galaxies. The separate broad fit equal-weights populations then galaxies; its objective intentionally differs and is labeled.
- Per-galaxy a0 is a diagnostic, not a MOND model with universal a0. Per-galaxy halo state is additional inferred information and should not be described as 208 different fundamental laws.

Files:
PROTOCOL.txt was written before fitting; its hash is in results.json.
run.py generates all 20 prediction CSVs and results/results.json, with 1254 selected fit sets and all 3762 optimizer attempts. All converged. It includes each exact training-row list, coefficient vector, objective, start, bounds, cost, convergence status, and selected ID. It requires numpy/scipy and sibling directories dm_unified_regular and dm_fit already included in the A43 supplement.
matched_components.csv independently matches the raw SPARC components to all 2700 selected rows. Gas uses signed vgas*abs(vgas), stellar components .5*vdisk*abs(vdisk)+.7*vbul*abs(vbul).
Each prediction CSV records the original row and selected-fit ID. Shared mixed_regular is refitted, not frozen from a data-leaking previous full fit.
source_notes.txt describes primary references. summary.json provides compact results and bound counts. External independent verification is supplied separately in the new A43 audit directory.
