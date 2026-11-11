# A42 verification package: inputs, code and paper (SA74 T32, November 10, 2026)

This zip holds the paper "Searching for a Short Discrete Theory of Everything: Rewriting Rules and Supplied-Lattice Tests" (A42), the code behind its claims, and the input data that code needs: run plans with their seeds and couplings, the specification texts, parameter files and rule lists. It holds no outputs. Every output (Markov-chain histories, result files, bit counts) is regenerated from the code and the inputs, and `verify.py` compares the regenerated numbers with the numbers the paper states.

## Quick start

    unzip SA74_T32_A42_VERIFY_*.zip && cd A42_VERIFY
    python3 verify.py                     # about 15 minutes

Requirements: Python 3.9 or later with NumPy and SciPy (SymPy and matplotlib for a few scripts), and a C compiler (cc, gcc or clang). The report goes to `VERIFY_REPORT.txt`. The checks work in `_work/`, a copy in which the absolute paths the scripts ran with (`/home/claude/sa74`, `/mnt/user-data/outputs/night`) are replaced by this folder's path.

## What verify.py checks

| Step | What it does | Expected result |
|---|---|---|
| 1 files | SHA-256 of every file against `MANIFEST.sha256` | all match |
| 2 engines | builds the four Regge engines; runs their self-tests (cell counts, flat deficit angles below 1e-12, the regular simplex's dihedral angle arccos(1/4), the Schlaefli identity, incremental against full recomputation); checks that the multi-core engine without threads gives the sequential engine's run and that 1 and 2 threads give the same run | all pass |
| 3 generators | reruns 47 fast scripts of the supplements from their code and inputs (`GENERATOR_JOBS.json`), and compares three headline numbers with the paper: 2,935 shortest rules, 2,488 bits for the Dirac step, 16,272 bits for the known laws on a supplied lattice; the other scripts' outputs are compared by hand with the sections the table below names | all run; the three numbers match exactly |
| 4 regge (`--regge`) | after the chains are regenerated (below): reruns the analyses of Section 8.10 and compares each quoted number with the paper's, within the combined statistical errors (next section): the fluctuation on every grid at three couplings and its fitted powers, no critical growth, the comparison with Hamber's curve, the largest fluctuation, the zero crossing on every grid, the 16⁴ start check, the median autocorrelation time on every grid, and at a = 0.02 the crossings and fluctuations | every check passes (on the original chains every difference is zero) |

## How the statistical checks decide

A regenerated chain is a new sample of the same Markov chain, so its numbers differ from the paper's by their statistical errors. Step 4 uses one rule and no hand-set ranges:

- `PAPER_NUMBERS.json` lists each number Section 8.10 quotes with its standard error, from the original chains. `regge_checks.py` computes the same numbers and errors from the regenerated chains; `make_paper_numbers.py` shows how the first file was written.
- Each comparison is a z value: z = (regenerated - paper) / sqrt(se_regenerated² + se_paper²). A single number passes at |z| <= 3. A check that compares m numbers at once passes when the largest |z| is within the bound that gives the whole check the same false-alarm chance, 0.27%: 3.21 for 2 numbers, 3.46 for 5, 3.58 for 8.
- The report prints every z value, so a pass shows how close the rerun came.

What that allows, in the paper's units. Where many chains stand behind a number the check is tight: the zero crossing must repeat within 0.0005 on 16⁴ and 0.0014 on 12⁴, the median autocorrelation time within a factor of 1.7 on 16⁴. Where few chains stand behind it the check is as wide as the number's own error: the 20⁴ median, from three chains, within a factor of 5; at a = 0.02 the fluctuation on 10⁴, 0.62(15), within 0.72, and the 4⁴ crossing, 0.479(19), within 0.09. The paper states the a = 0.02 errors beside its numbers. These widths assume the rerun's errors equal the paper's.

With fifteen such checks, each with a 0.27% false-alarm chance when the errors are exact, about one full rerun in twenty-five fails one check by chance; errors estimated from few independent sweeps make that somewhat more likely. A failed check marks a number to look at again (the report gives its z value), for example by rerunning that chain from another seed.

## Regenerating the lattice-gravity chains (Sections 8.9 and 8.10)

    python3 regenerate.py list               # every chain: plan file, grid, coupling, sweeps, seed, start
    python3 regenerate.py cloud all          # the cloud stages (T25 audit, FSS, AS, AU, HY1), each with its own plan and runner
    python3 regenerate.py scan               # the 3⁴ pilot and the 3⁴ and 4⁴ scans
    python3 regenerate.py run --all          # the Mac plans: N2b, Q1 to Q6, A1 to A4, A5 (20⁴, multi-core engine)
    python3 verify.py --regge

Eight Mac chains stopped before their planned length in the original night and FSS3's last 16 cloud chains did not run; regenerate.py repeats those lengths, so the regenerated chains are the inputs the paper's analyses had (`--full` runs the eight to their planned lengths). This is days of CPU time: 4⁴ to 8⁴ chains take minutes, 12⁴ chains a few hours, 16⁴ chains about ten hours each, the three 20⁴ chains about 1.5 hours each on eight cores. One plan or chain at a time: `regenerate.py run --plan queue3.txt`, `regenerate.py run A5_L20_k0.244`. `--scale 0.01` runs every chain at a hundredth of its length into `_smoke/`, a quick test of the machinery. The cloud chains repeat byte for byte with gcc -O2 on x86-64. The Mac chains ran in reseeded 600-second segments with clang on Apple silicon, so a rerun here is a new sample of the same chain: its averages agree within their errors, its digits differ.

## Where each claim comes from

| Claim (paper section) | Code and inputs | How to regenerate |
|---|---|---|
| Census of the 2,935 shortest rules and 3,500 sampled longer ones; no rule reads four on all three dimension tests; near misses 4.7 to 21 standard deviations below four (Sections 2 to 4) | `supplements/A42_Rule_Search.py`, `A42_Compile_Results.py`, `A42_Search_README.md` | the commands in `A42_Search_README.md` (several hours); step 3 runs the enumeration and the validation of the tests |
| Update order belongs in a rule's length; the link bound; transitive and originary percolation (Sections 4.6, 5, 5.1) | `supplements/A42_SA74/` (Extensions, tp_big.c, Numbers, Classify; `long27_rules.json`; Protocol) | commands in `A42_SA74_README.md` |
| One local unitary step reproduces the Dirac equation in curved spacetime, horizon crossing and pair creation in 2,488 bits (Sections 8.1 to 8.5) | `supplements/A42_Step_Rule/`, `supplements/A42_Supplied_Physics/` | step 3 |
| Harder tests: pair creation, gauge fields, the source of gravity (Section 8.6) | `supplements/A42_Harder_Tests/Pass10` to `Pass36_R2`, each with its README | step 3 |
| The known laws on a supplied lattice compress to 16,272 bits (Section 8.7, Table 12) | `supplements/A42_SA74/A42_SA74_Supplied_Lattice.py`, `supplements/A38_Symbolic_Benchmark_Spec.txt`, the A38 bit-count code | step 3 |
| Independent Regge audit: all four groups fail the mixing screen (Section 8.9) | `supplements/A42_Regge_Audit_T19/` (protocol with SHA-256, plan, run and analysis code) | `PROTOCOL.txt` and `run_audit.py` there |
| With Hamber's curvature-squared term all four groups pass (Section 8.10) | `supplements/A42_Regge_R2_T25/` and `sa74/regge/r2/audit_a005/` (plans `plan.json`, `plan_B.json`, `plan_B2.json`) | `regenerate.py cloud plan.json` (and B, B2), then `analyze_r2.py` (see the T25 README) |
| Larger grids, 4⁴ to 20⁴: fixed-coupling powers -0.08(10), -0.11(7), -0.05(7); Hamber's fitted curve 2.8 (2.3) standard errors above; crossings 0.2516 to 0.2541; median autocorrelation 120 to 240 sweeps; the 16⁴ start check; a = 0.02 crossings 0.479(19), 0.467(3), 0.468(1), 0.467(1) (Section 8.10) | plans: `sa74/regge/r2/audit_a005/plan_*.json`, `ext/night/jobs_native.txt`, `sa74/NIGHT/mac_pull/runner/queue.txt`, `queue2.txt`, `queue3.txt`; engines and runners; analysis code in `sa74/regge/r2/night/` and `fss/` | `regenerate.py`, then `verify.py --regge` |
| How each number in the text follows from the result files | `sa74/scripts/make_sa74.py` (the SA74 edits applied to `sa74/src/A42_SA72_T372_package.md`, with asserts on each sentence's wording) | runs once every result file it reads is regenerated |

## Layout

- `paper/`: the paper (PDF and Markdown).
- `sa74/`: the working tree of the SA74 additions at its original relative paths: `scripts/` (the paper builder), `src/` (the base text), `regge/r2/` (engines, plans, runners, analysis code), `regge/par/` (the multi-core engine), `NIGHT/` (the run log with the times and hashes of the plans, and the Mac runner files).
- `ext/night/`: the Mac engines, runner scripts and plans as first written.
- `supplements/`: the code, inputs and documentation of the published A42 supplements and the A38 bit-count files. Result files their READMEs name are produced by the commands given there.

## Notes

- Plans were fixed and hashed before their runs; `sa74/NIGHT/NIGHT_LOG.md` gives the order of events (its November 7 entries labeled 13:05 to 13:40 were written between about 12:56 and 13:14 by the file clocks; corrected in its November 9 entries).
- Step 3 adjusts two published scripts' runs (see `GENERATOR_JOBS.json`): `Pass13/audit_scaling.py` imports `curved_gauge_pass12.py`, a copy of `Pass12/curved_gauge.py` that the Pass13 folder lacks, so the working copy gets one; `Pass36_R2/r2_bounce.py` runs with the 14 values of b behind the paper's scan (its default list has 11). `A42_SA74_Numbers.py` is left out of step 3: it reads the rule-search runs, which take hours to regenerate.
- Some scripts name the folders they ran in (`/home/claude/sa74/...`); to run one by hand, replace that prefix with the path of this folder's `sa74/` (and `/mnt/user-data/outputs/night` with `ext/night/`).
- A single chain's variance runs low by about 2/ESS; the paper quotes the fitted power with and without that correction.
- The supplement Section 8.10 cites, A42_Regge_R2_Night_T31, is `supplements/A42_Regge_R2_Night_T31/`: a guide to where its plans, engines and analysis code sit in this package, and its own README.
- `sa74/NIGHT/mac_pull/queue3.txt` and the run log are as they stood on November 9 (T31), the state the paper's results rest on; later lines plan runs the paper does not use.
- Replaced before zipping: the address of the author's server, the name of its key file and the name of an unrelated project, in the run log and in the scripts that launched runs on that server. Those launch scripts therefore no longer match the hashes the log quotes for them; no plan of a chain the paper uses is affected, and neither script nor server is needed here.
- `sa74/regge/r2/night/multihist.py` has one addition since T31, the bootstrap spread of the zero crossing; a rerun on the original chains reproduced every other number exactly. `multihist_T31.py` is the earlier file.
- The bootstrap in `multihist.py` draws chains in the order they are read, so its errors, and the powers fitted with them, depend on that order (with 40 bootstrap samples, by a few hundredths in the powers). T31's `reweight.py` read the cloud chains in the file system's order (`reweight_T31.py`); `sa74/regge/r2/audit_a005/chain_order_T30.txt` now records the order of the paper's run and `reweight.py` follows it, so the same chains give the same numbers on any machine.
