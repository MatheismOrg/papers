# SA74 night supplement: the curvature-squared Regge ensemble on larger grids (6-7 November 2026; T29 update, morning of 7 November)

## What was asked
Does the supplied-lattice gravity of A42 Section 8.10 (Regge action with Hamber's curvature-squared term, a = 0.005, DeWitt measure dl², lambda = 1, the 15-edge, 50-triangle, 24-simplex hypercubic triangulation) show the finite-size signs of a critical point near Hamber's k_c = 0.2443(11)? How small can a be while the ensemble still mixes, and does a = 0.02 behave the same way?

## Answer, in short
- **No critical growth from 4⁴ to 12⁴.** The curvature fluctuation at fixed coupling, combined across all chains on each grid by multiple-histogram reweighting (Ferrenberg and Swendsen 1989), is <<f248_list>> on 4⁴, 5⁴, 6⁴, 8⁴, 10⁴, 12⁴ at k = 0.248 (fitted power <<p248>> in L) and <<f252_list>> at k = 0.252 (<<p252>>); at Hamber's k_c, k = 0.244, <<f244_list>> (<<p244>>). A critical point with delta = 0.62 predicts +0.94 (a factor 2.8 from 4⁴ to 12⁴). The largest value near the crossing (k = 0.2415 to 0.258, and within 0.002 of a coupling run) is <<peak_list>>: power <<pw>>, or <<pw2>> from the second halves of the chains (results/multihist_T29.json, multihist_skip50_T29.json, a1_results.json, chi_fixedk.png). Grids below 4⁴ are left out, as Hamber (1993) left them out (the term couples next-nearest neighbors; the 3⁴ value also moved between the halves of its chains).
- **Hamber's own fit is not excluded:** his 1993 fit of the fluctuation, (0.16 + 1.11/L) L^0.64(11) at k = 0.244, rises 1.17 times from L = 4 to 12 (an effective power of 0.14) and 1.27 times to 16; these grids cannot tell that from no growth.
- **12⁴ (plan A1, Mac runner 2, sequential edge order):** <<n12>> chains at k = 0.244 to 0.254, each 3,000 sweeps at its coupling and then 4,000 measured (the r1 chains continue the night's 12⁴ configurations, which drifted after 2,000 sweeps; the k = 0.244 r2 chain is a fresh flat start and agrees with r1: <<r12_244>>). Between the halves of the measurement the mean squared edge moves by at most <<dq12>> and R by at most <<dR12>>.
- **14⁴:** two single chains of 800 + 1,600 sweeps at k = 0.2475 and 0.2525 show no drift by sixths in R or the mean squared edge; they are too short for the fits; their fluctuations, <<c14a>> and <<c14b>>, show no growth either. (T28's paper text grouped them with the drifting 12⁴ chains; T29 corrects that.)
- **Smooth zero crossing** of the average curvature at k = <<cr_lo>> to <<cr_hi>> on 4⁴ to 12⁴ (12⁴: <<cr12>>), a little above Hamber's fitted k_c; its slope d<R>/dk = cov(R, Σ δA) matches finite differences between neighboring couplings and does not steepen with size.
- **No slowing with size:** median autocorrelation times (n/ESS of R) of chains near the crossing are <<tau_list>> sweeps on 4⁴, 5⁴, 6⁴, 8⁴, 10⁴, 12⁴.
- **The slowest motion is local:** clusters of long edges form and anneal over thousands of sweeps and lift R while they last; within flat chains near the crossing, 100-sweep means of the longest edge and of R correlate (+0.52 on 4⁴, +0.35 on 6⁴, +0.29 on 8⁴, +0.12 on 10⁴), and the shift in R per unit of longest edge falls with volume (results/defects.json). On 6⁴ at k = 0.2475 four flat chains spread 3.6 times their internal errors because of them.
- **Rough starts on 8⁴:** at k = 0.250 and 0.255 they reach the flat-start values within 3,000 sweeps, from either of two rough preparations at k = 0.250 (the second: <<rough2_250>>). At k = 0.245 the rough chain needed about 7,000 sweeps, until its long edges shrank (last third -0.019; the three flat starts -0.017 to -0.027); a start from the second preparation reached -0.013 after 4,500 sweeps with 21 long edges on 10 sites left at its end.
- **Small a:** a = 0.001 is start-dependent near its crossing (split R-hat 1.21, flat and rough means 4.6 standard errors apart); a = 0.002 is slow but unbiased (1.136, 0.68 standard errors); a = 0.003 passes the screen (1.037, 0.33 standard errors); at a = 0.005 T25's stage B passed on 3⁴ at k = 0.24 (1.047) and, continued, on 4⁴ (stage B2, 1.029) (results/as5_results.json, as3_results.json).
- **a = 0.02 (AU1, AU2):** the fluctuation is a tenth as large and grows by at most a fifth from 4⁴ to 8⁴ (at k = 0.46: <<au46>> on 4⁴, 6⁴, 8⁴; power <<au_pow46>>), where critical growth with Hamber's delta = 0.628 would nearly double it; the curvature crosses zero at k = <<au_cr4>> (4⁴), <<au_cr6>> (6⁴) and about <<au_cr8>> (8⁴, from R = <<au_R8>> at k = 0.46 and the fluctuation slope), drifting down with size but above Hamber's k_c = 0.4244(36) (results/au_results.json).
- **16⁴ (plan A1):** two chains at k = 0.248 and 0.252 (3,000 + 3,000 sweeps) were still running at T29.

## Stages (each plan fixed and hashed before its runs; times PST; NIGHT_LOG.md has the order of events)
| Stage | Where | What |
|---|---|---|
| FSS8 | cloud | 8⁴ flat starts at k = 0.24, 0.26 |
| AS, AS2 | cloud | a = 0.001, 0.002: start independence at k = 0.05, 0.065 on 4⁴; k-scans on 4⁴ |
| AS3 | cloud | a = 0.001 near its crossing: start independence at k = 0.14 on 4⁴; 6⁴ k-scan |
| AS4 | cloud | a = 0.002 start independence at k = 0.175 on 4⁴ |
| AS5 | cloud | a = 0.003 start independence at k = 0.20 on 4⁴ |
| AU1, AU2 | cloud | a = 0.02: 3⁴, 4⁴ at k = 0.40 to 0.46, 6⁴ at 0.42, 0.44 (AU1); 4⁴ at 0.48, 0.50, 6⁴ at 0.46, 0.48, 8⁴ at 0.46 (AU2, done 05:28) |
| M1b, N2b | Mac (VM chunks, then native) | 8⁴, 10⁴, 12⁴ flat starts at k = 0.245 to 0.255 |
| S1 | server | 12⁴ to 16⁴ (stopped at 00:27: it halved SA73's share of the server) |
| Q1 | Mac queue | 14⁴ and 16⁴ (16⁴ stopped at 01:31 as too slow to equilibrate before morning), 10⁴, 12⁴, 8⁴ points, and 8⁴ rough starts from an a = 0, k = 0.1 preparation |
| Q2, Q3 | Mac queue | continuations of 8⁴ chains at k = 0.255 and 0.26; rough starts at 0.245, 0.2475; continuation of the 0.250 rough chain; a second rough preparation |
| Q4, Q5 | Mac queue | replica flat starts on 8⁴ (k = 0.245, 0.2475, 0.2525) and 10⁴ (0.2475) |
| Q6 | Mac queue | 6,000 more sweeps of the k = 0.245 rough chain; a k = 0.245 rough start from the second preparation |
| A1 | Mac runner 2 (sequential edge order, Wilson lines; no end time) | 12⁴ at k = 0.244 to 0.254 in steps of 0.002 (r1 from the night's configurations, r2 flat), 16⁴ at 0.248, 0.252, then 8⁴; <<a1_status>> |
| HY1 | cloud | flat against rough starts at k = 0.250, 0.255 on 4⁴ (3 + 3) and 6⁴ (1 + 1) |
| FSS2, FSS3, FSS4 | cloud | 3⁴ to 6⁴ at k = 0.2425 to 0.265, up to five chains per point. The cloud chain stopped at its two-hour limit at 06:16: FSS4 was done (05:50); of FSS3's last runs, 6⁴ at k = 0.2525 r2 and r3 finished, two 4⁴ replicas at k = 0.2475 were cut off halfway (no done line; every analysis leaves them out), and the 4⁴ k = 0.2475 r3, 4⁴ k = 0.2525 and all 3⁴ replicas did not run |
| V1, V2 | Mac VM | exploratory: per-slice curvature sums on 8⁴ from thermalized configurations |

## Observables
- R: the engine's average curvature 2 Σ δ_h A_h / Σ V (column 5); Hamber's R(k) = <l²><2 Σ δA>/<Σ V> (Hamber 2015, eq. 14) is computed in the analysis.
- χ = 4 var(Σ δ_h A_h) / <Σ V>, Hamber's curvature fluctuation; at a critical point χ grows as L^(2/ν - 4) (Hamber 2015, eq. 93).
- Reweighting: the engine samples exp(k Σ δA) times a k-independent weight, so each chain reweights to nearby k with weights exp((k' - k) Σ δA); multihist.py combines all chains on a grid, each weighted by n/ESS, with bootstrap errors (chains within each coupling, then five blocks of each chain). T29 settings: A1INC=1 (plan A1's finished chains replace the earlier 12⁴ chains), LMAX=12, KPAD=0.002 (the largest χ only within 0.002 of a coupling run, so that no value is extrapolated far from the chains on 12⁴; the values on 4⁴ to 10⁴ are unchanged by it), χ at fixed couplings (KFIX) with bootstrap errors. The T28 runs of the same script are kept as multihist_T28.json and multihist_skip50_T28.json (script: multihist_T28.py).
- Split R-hat across flat and rough starts (screen 1.05), effective sample size, autocorrelation time n/ESS.
- Longest squared edge (column 10) and, from final configurations, the edges longer than four times the mean squared length.

## Rebuild
analysis/night_final.py reruns every analysis from the chain files and writes results/night_summary.json; fig_fixedk.py draws results/chi_fixedk.png. The T28 versions of the final pass and of this builder are kept as night_final_T28.py and night_supplement_T28.py.
