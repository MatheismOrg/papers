# A42_Regge_R2_Night_T37: the supplement Section 8.10 cites

Section 8.10 of the paper cites this supplement for its plans (fixed before their runs), per-sweep histories and analysis code. This package holds inputs and code only, so the supplement's parts sit here:

- **Plans.** Cloud stages: `sa74/regge/r2/audit_a005/plan.json` and the other `plan_*.json` files in that folder. Mac plans: `ext/night/jobs_native.txt` (N2b), `sa74/NIGHT/mac_pull/runner/queue.txt` (Q1 to Q6), `sa74/NIGHT/mac_pull/queue2.txt` (A1 to A4) and `sa74/NIGHT/mac_pull/queue3.txt` (A5, the 20⁴ chains). `sa74/NIGHT/NIGHT_LOG.md` gives the time each plan was fixed and its SHA-256.
- **Engines and runners.** `sa74/regge/r2/regge_r2.c` (cloud), `ext/night/regge_r2_big.c` and `sa74/regge/par/regge_r2_big_x.c` (Mac, one core), `sa74/regge/par/regge_r2_par.c` (Mac, several cores); the runner scripts are in `ext/night/` and `sa74/NIGHT/mac_pull/runner/`.
- **Analysis code.** `sa74/regge/r2/night/` (`multihist.py`, `analyze_a1.py`, `analyze_au.py` and the rest) and `sa74/regge/r2/fss/analyze_fss.py`.
- **Per-sweep histories and result files.** Not in the zip: `regenerate.py` reruns every chain from its plan, and `verify.py --regge` reruns the analyses and compares each number with the paper (`PAPER_NUMBERS.json`).
- **`README_NIGHT.md`** in this folder is the supplement's own account of its stages, plans and findings, as written for T31. The result files and figures it mentions are the outputs the commands above regenerate.

One addition since T31: `sa74/regge/r2/night/multihist.py` also writes the bootstrap spread of the zero crossing, which `verify.py` uses. Its other numbers are unchanged (a rerun on the original chains reproduced them exactly); `sa74/regge/r2/night/multihist_T31.py` is the earlier file. `sa74/regge/r2/night/reweight.py` now reads the cloud chains in the recorded order of the paper's run (`sa74/regge/r2/audit_a005/chain_order_T30.txt`); `sa74/regge/r2/night/reweight_T31.py` is the earlier file.


## T37 (November 11, 2026): plan A6, the 24⁴ and 32⁴ chains

Plan A6 was appended to `sa74/NIGHT/mac_pull/queue3.txt` on November 9 before its runs: A6_L20_k0.244_r2 (flat start, 2,500 + 3,000 sweeps), A6_L24_k0.244 (flat, 2,500 + 3,000), A6_L24_k0.248 and A6_L24_k0.252 (each from the previous chain, 1,500 + 3,000) and A6_L32_k0.244 (flat, 2,500 + 3,000), all with the multi-core engine. `regenerate.py run --plan queue3.txt` reruns them; the 32⁴ chain takes about 12 hours on eight cores and needs about 13 GB of memory. The analysis is `analyze_a1.py`, unchanged, with `MHFILE=multihist_T30.json OUTSUFFIX=_T37`, as `verify.py --regge` runs it. The original per-sweep histories of these five chains are in the papers repository (MatheismOrg/papers, SUPPLEMENTS/A42_SA74_Data/regge_night_A6) once its pull request is merged.
