LEAN DISTRIBUTION NOTE: the per-sweep histories (histories.npz) and every final configuration are in SA74_T25_A42_REGGE_R2_SUPPLEMENT_2146_261106.zip; they regenerate byte for byte from the seeds in the plans with the included code.

# A42 T25: the Regge audit rerun with a curvature-squared term

**What this is.** The T19 audit (A42 Section 8.9) found that the supplied Regge engine does not forget its start: chains begun rough stayed rough, so all four parameter groups failed the split R-hat screen. T25 adds one term to the action, Hamber's lattice curvature-squared term, and reruns the same audit design.

**The term.** I = λ Σ_s V_s − k Σ_h δ_h A_h + a Σ_h δ_h² A_h² / V_h, with V_h a tenth of the summed volume of the 4-simplices that share triangle h (Hamber's volume per hinge), a = 0.005, λ = 1. The measure (uniform in squared edge lengths), the moves (single-edge Metropolis plus ten global rescalings per sweep) and everything else are unchanged. The term is invariant under the global rescaling, so the scaling identity ⟨2λΣV − kΣδA⟩ = N_edges still holds.

**Results.**

| Grid | k | Flat-start mean R | Rough-start mean R | Split R-hat | Split R-hat at a = 0 (T19) |
|---|---|---|---|---|---|
| 3⁴ | 0.05 | -0.2762 | -0.2717 | 1.004 | 1.573 |
| 3⁴ | 0.065 | -0.2649 | -0.2665 | 1.006 | 13.032 |
| 4⁴ | 0.05 | -0.2751 | -0.2737 | 1.003 | 5.498 |
| 4⁴ | 0.065 | -0.2616 | -0.2621 | 1.004 | 17.716 |

Stage A (the T19 design, 3,000 + 12,000 sweeps): all four groups pass; rough starts enter the flat-start band within 500 sweeps; effective sample sizes 140 to 334 per chain (T19: 3 to 29); scaling identity within 0.08 percent; incremental totals within 1e-13 of full recomputation after every sweep.

Stage B (k = 0.20 and 0.24, near Hamber's k_c = 0.2443 for a = 0.005; 2,000 + 6,000 sweeps): 3⁴ at k = 0.20: split R-hat 1.013 (passes), flat -0.1131, rough -0.1167; 3⁴ at k = 0.24: split R-hat 1.047 (passes), flat -0.0032, rough -0.0350; 4⁴ at k = 0.20: split R-hat 1.023 (passes), flat -0.0988, rough -0.1036; 4⁴ at k = 0.24: split R-hat 1.077 (fails), flat -0.0350, rough -0.0246. Effective sample sizes 17 to 84 per chain: mixing slows near the transition.

Stage B2 (the six 4⁴, k = 0.24 chains continued for 500 + 12,000 sweeps with new seeds): split R-hat 1.029 (passes), flat -0.0282, rough -0.0313, start difference -0.37 standard errors, effective sample sizes 27 to 75 per chain.

Scans (exploratory, from flat and rough starts): the average curvature rises from about -0.24 at k = 0.10 and crosses zero between k = 0.24 and 0.26 on 3⁴, and in the same interval on a shorter 4⁴ scan. At k = 0.30 on 3⁴ the two starts begin to differ.

**Lower bound.** On each triangle, λV_h − kδA + aδ²A²/V_h ≥ (λ − k²/4a) V_h (complete the square in δA). Summed, I ≥ (λ − k²/4a) Σ_s V_s, positive for k < 2√(aλ) = 0.141 at a = 0.005, λ = 1. Every recorded sweep of stage A satisfies it with margin. This bounds the action; it does not by itself show that the integral over degenerate configurations converges, and it does not reach k_c.

**Bits.** Adding the term to the supplied-lattice specification ("plus 0.005 (delta A)^2/V per triangle (V: the volume of its simplices / 10)") raises the LZMA count from 16,272 to 16,624 bits (bits_r2.json; three other wordings cost 376 to 400 bits).

**Files.**
- `make_r2.py` builds `regge_r2.c` from the unchanged audited engine `regge_audited_T19.c` (sha256 31dfe6fe...) by 19 exact, unique text replacements. The coupling a is read from the environment variable `R2A`; unset or 0 gives the audited engine.
- `run_audit_r2.py ENGINE.c [workers] [plan]` builds, self-tests and runs a plan: `plan.json` (stage A: the T19 plan with a = 0.005 in the 24 main chains; preparations stay at a = 0), `plan_B.json`, `plan_B2.json`. Each plan was written before its runs.
- `analyze_r2.py [plan] [runsdir] [out]` computes the T19 diagnostics with the T19 statistics functions copied verbatim (ordinary split R-hat, batch means, ESS, drift) plus the per-sweep check of the lower bound. Outputs: `results_A.json`, `results_B.json`, `results_B2.json`.
- `histories.npz`: per-sweep R, sumV, sumDA and sumR2 (float32) for every chain, preparation, scan and pilot run, with each run's number of equilibration sweeps. `runs_final_cfg/`: every final configuration. `logs/`: per-run logs and execution records. `COMMANDS.txt`: exact commands of the exploratory runs. `fig/`: figures and `make_figs.py`.

**Checks.** With R2A unset or 0 the new engine reproduces the audited engine's self-test output (3⁴ and 4⁴), per-sweep histories and saved configurations byte for byte, including from a rough loaded start. With a = 0.005 the incremental curvature-squared total and every V_h agree with full recomputation (self-test and every sweep), and a flat start has zero curvature-squared energy (6e-26).

**Reproduce.** `python3 run_audit_r2.py regge_audited_T19.c 2 plan.json`, then `plan_B.json` and `plan_B2.json`, then `analyze_r2.py` on each. About 90 minutes on two cores. The full 28-column histories regenerate byte for byte with gcc -O2 on x86-64.

**Limits.** Lattices of 3⁴ and 4⁴ sites only. R-hat below 1.05 is a screen, not proof of equilibrium. Near k_c mixing slows. The Weyl-squared part of the term carries a massive spin-2 ghost in perturbation theory (Stelle, 1977); Hamber (1992) found the critical exponents independent of a. No claim here of a continuum limit, a Newtonian potential or Lorentzian dynamics.

**Sources.** Hamber, H. W. (1993), Nuclear Physics B 400, 347-389 (a = 0.005, k_c = 0.2443(11)). Hamber, H. W. (1992), Physical Review D 45, 507-512. Hamber, H. W. (1993), arXiv:hep-th/9310152 (the action with the a-term). Utiyama, R., & DeWitt, B. S. (1962), Journal of Mathematical Physics 3, 608-618. Stelle, K. S. (1977), Physical Review D 16, 953-969.
