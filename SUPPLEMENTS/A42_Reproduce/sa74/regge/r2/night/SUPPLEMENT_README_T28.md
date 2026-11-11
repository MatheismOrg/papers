# SA74 night supplement: the curvature-squared Regge ensemble on larger grids (6-7 November 2026)

## What was asked
Does the supplied-lattice gravity of A42 Section 8.10 (Regge action with Hamber's curvature-squared term, a = 0.005, DeWitt measure dl², lambda = 1, the 15-edge, 50-triangle, 24-simplex hypercubic triangulation) show the finite-size signs of a critical point near Hamber's k_c = 0.2443(11)? How small can a be while the ensemble still mixes, and does a = 0.02 behave the same way?

## Answer, in short
- **No critical growth from 3⁴ to 10⁴.** Reweighted across couplings (multiple histograms, Ferrenberg and Swendsen 1989), the largest curvature fluctuation between k = 0.2415 and 0.258 is 5.6(6), 6.3(5), 5.3(5), 6.4(8), 4.7(5), 5.9(4) on 3⁴, 4⁴, 5⁴, 6⁴, 8⁴, 10⁴: fitted power -0.04(11) in L, where a critical point with delta = 0.62 predicts +0.94 (results/multihist.json, results/chi_reweighted.png).
- **Smooth zero crossing** of the average curvature at k = 0.2517 to 0.2533 on 4⁴ to 10⁴, a little above Hamber's fitted k_c; its slope d<R>/dk = cov(R, Σ δA) (2.2 to 2.9) matches finite differences between neighboring couplings and does not steepen with size.
- **No slowing with size:** median autocorrelation times near the crossing 130 to 290 sweeps from 3⁴ to 10⁴.
- **Hamber's own fit is not excluded:** his 1993 fit of the fluctuation, (0.16 + 1.11/L) L^0.64(11) at k = 0.244, rises only 1.27 times from L = 4 to 16 (1.11 times from 4 to 10).
- **Rough starts on 8⁴:** at k = 0.250 and 0.255 they reach the flat-start values once the long edges of the rough preparation shrink; at k = 0.245 and 0.2475 about 60 of 61,440 edges stay two to three times the typical length, in clusters on about 30 of 4,096 sites, and R stays higher while they last (results/defects.json). Plan Q6: 6,000 more sweeps of the k = 0.245 rough chain; its long edges shrank to flat-start lengths after about 7,000 sweeps in all, and R fell with them to the flat-start value (last third -0.019; the three flat starts average -0.017 to -0.027). A k = 0.245 start from the second rough preparation reached -0.013 after 4,500 sweeps, with 21 long edges on 10 sites left at its end. Flat starts form such clusters too: within flat chains near the crossing, 100-sweep means of the longest edge and of R correlate (+0.52 on 4⁴, +0.35 on 6⁴, +0.29 on 8⁴), and on 6⁴ at k = 0.2475 four flat chains spread 3.6 times their internal errors because of them.
- **Small a:** a = 0.001 is start-dependent near its crossing (split R-hat 1.21, flat and rough means 4.6 standard errors apart); a = 0.002 is slow but unbiased (1.136, 0.68 standard errors); a = 0.003 passes the screen (1.037, 0.33 standard errors); at a = 0.005 T25's stage B passed on 3⁴ at k = 0.24 (1.047) and, continued, on 4⁴ (stage B2, 1.029) (results/as5_results.json, as3_results.json).
- **a = 0.02:** the fluctuation is a tenth as large (0.44 to 0.69) and does not grow from 3⁴ to 6⁴ (fitted powers -0.11 and -0.15 at k = 0.42 and 0.44); at k = 0.44, past Hamber's k_c = 0.4244(36), R is still negative on 4⁴ and 6⁴ (results/au_results.json). AU2 (the a = 0.02 crossing on 4⁴ to 8⁴) was still running when this supplement was built; its chains follow in a later update.
- **12⁴ and 14⁴** flat starts still drift in mean edge length and R after 2,000 sweeps; they are left out of the fits. The Mac's runner 2 (plan A1) continues 12⁴ from tonight's configurations.

## Stages (each plan fixed and hashed before its runs; times PST; NIGHT_LOG.md has the order of events)
| Stage | Where | What |
|---|---|---|
| FSS8 | cloud | 8⁴ flat starts at k = 0.24, 0.26 |
| AS, AS2 | cloud | a = 0.001, 0.002: start independence at k = 0.05, 0.065 on 4⁴; k-scans on 4⁴ |
| AS3 | cloud | a = 0.001 near its crossing: start independence at k = 0.14 on 4⁴; 6⁴ k-scan |
| AS4 | cloud | a = 0.002 start independence at k = 0.175 on 4⁴ |
| AS5 | cloud | a = 0.003 start independence at k = 0.20 on 4⁴ |
| AU1, AU2 | cloud | a = 0.02: 3⁴, 4⁴ at k = 0.40 to 0.46, 6⁴ at 0.42, 0.44 (AU1); 4⁴ to 8⁴ at 0.46 to 0.50 (AU2) |
| M1b, N2b | Mac (VM chunks, then native) | 8⁴, 10⁴, 12⁴ flat starts at k = 0.245 to 0.255 |
| S1 | server | 12⁴ to 16⁴ (stopped at 00:27: it halved SA73's share of the server) |
| Q1 | Mac queue | 14⁴ and 16⁴ (16⁴ stopped at 01:31 as too slow to equilibrate before morning), 10⁴, 12⁴, 8⁴ points, and 8⁴ rough starts from an a = 0, k = 0.1 preparation |
| Q2, Q3 | Mac queue | continuations of 8⁴ chains at k = 0.255 and 0.26; rough starts at 0.245, 0.2475; continuation of the 0.250 rough chain; a second rough preparation |
| Q4, Q5 | Mac queue | replica flat starts on 8⁴ (k = 0.245, 0.2475, 0.2525) and 10⁴ (0.2475) |
| Q6 | Mac queue | 6,000 more sweeps of the k = 0.245 rough chain; a k = 0.245 rough start from the second preparation |
| A1 | Mac runner 2 (sequential edge order, Wilson lines; no end time) | 8⁴ and 12⁴ at k = 0.244 to 0.254 in steps of 0.002; 16⁴ at 0.248, 0.252 |
| HY1 | cloud | flat against rough starts at k = 0.250, 0.255 on 4⁴ (3 + 3) and 6⁴ (1 + 1) |
| FSS2, FSS3, FSS4 | cloud | 3⁴ to 6⁴ at k = 0.2425 to 0.265, up to five chains per point |
| V1, V2 | Mac VM | exploratory: per-slice curvature sums on 8⁴ from thermalized configurations |

## Observables
- R: the engine's average curvature 2 Σ δ_h A_h / Σ V (column 5); Hamber's R(k) = <l²><2 Σ δA>/<Σ V> (Hamber 2015, eq. 14) is computed in the analysis.
- χ = 4 var(Σ δ_h A_h) / <Σ V>, Hamber's curvature fluctuation; at a critical point χ grows as L^(2/ν - 4) (Hamber 2015, eq. 93).
- Reweighting: the engine samples exp(k Σ δA) times a k-independent weight, so each chain reweights to nearby k with weights exp((k' - k) Σ δA); multihist.py combines all chains on a grid, each weighted by n/ESS.
- Split R-hat across flat and rough starts (screen 1.05), effective sample size, autocorrelation time n/ESS.
- Longest squared edge (column 10) and, from final configurations, the edges longer than four times the mean squared length.

## Rebuild
analysis/night_final.py reruns every analysis from the chain files and writes results/night_summary.json and the figures.
