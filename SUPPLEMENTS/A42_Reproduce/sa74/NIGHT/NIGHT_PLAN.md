# SA74 night of 6-7 November 2026: overnight ToE improvement search (Gary 22:45)

Gary: "i give u access to keep using my mac to keep prompting u to keep searching for improvements all night till 5am la my time"

## Standing rules
- Nothing applied to A42 (or any paper) without Gary's "apply"; results become proposals (props_N*.json) with full paragraphs, bold changes, Score/WC/WI/page.
- Mac (gs-mac-mini-local, 4 aarch64 cores, 3.9 GB, VM via device_bash): compute only in the VM scratch (cwd of device_bash). Connected folders ~/scripts/matheism-portfolio and ~/papers: read only, never modify. No background processes survive between device_bash calls (tested 22:46): run in chunks of at most ~165 s, save configurations, continue next call.
- Cloud container: 2 cores; background jobs up to 2 h re-invoke on exit.
- Every plan fixed (written, hashed) before its runs; exploratory scans labeled as such.
- Morning summary at about 05:00 LA: phone PDF (TOOLS/reply_pdf.py), short, bottom line first, bold key ideas, sources.

## Queue (most valuable first)
1. Finite-size growth of curvature fluctuations near k_c, a = 0.005: 8^4 (cloud, running), then 8^4/10^4/12^4 on the Mac. Observable: R(k), chi = 4 var(sum dA)/<sum V> at k = 0.24-0.255; crossing k0(L).
2. Smallest a that still mixes: a = 0.001, 0.002 at 4^4 from the six rough preparations at k = 0.05, 0.065 (T19 design, shorter), then k-scans near k_c(a) to fit delta(a). Question: does delta move toward 0.30 (a = 0) while mixing holds?
3. Matter: nphi = 4 free scalars (the Higgs doublet's real components) on the a = 0.005 ensemble: smooth phase kept? shift of the zero crossing?
4. Paper: proposals only; T26 proposals pending (8.10 larger grids, Hamber 1992 sentence fix, costs sentence).

## Files
- Engine: /home/claude/sa74/regge/r2/regge_r2.c (R2A env = a). Mac copy built in the VM scratch.
- Results so far: regge/r2/audit_a005 (stages A, B, B2, fss), regge/r2/fss (analysis, fits), t5/T25, t5/T26.
- Log: NIGHT/NIGHT_LOG.md (append one line per event).
