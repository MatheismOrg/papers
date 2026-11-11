Noise floors and pattern tests for Section 5.5, Table 9 (session SA73, turn T54, 10 November 2026 PST). No protocol:
nothing is fitted to the benchmark; the formula is the frozen T48 candidate (A, B, C = 11.7481, 0.8178, -0.7119; log M_t = 11.339).
Run from the dmfit folder (the harness modules must be importable); <dmfit dir> is that folder.
  noise54_disks.py <dmfit dir> 2000   disks: noise from SPARC's e_D, e_Inc, e_L[3.6], errV and Lelli et al. (2017)'s 25 percent
                                      in M/L (0.11 dex) and 10 percent in HI flux -> res54_noise_disks.json
  patterns54_disks.py <dmfit dir> 1000  within-galaxy and galaxy-to-galaxy statistics against the same noise -> res54_patterns_disks.json
  diag54_within.py <dmfit dir>        stacked within-galaxy residuals and per-galaxy slopes, formula, MOND, standard halos (prints)
  audit54_lens_etg.py <dmfit dir>     lensing chi-square by stack with N and p; early types' noise floor -> res54_lens_etg.json
  audit54_rest.py <dmfit dir>         dwarfs, MaNGA, CLASH and early types (with distance errors on the inputs): full noise
                                      distributions, the three accounts placed in them, mean offsets, dwarf slope -> res54_rest.json
  lensfloor54.py <dmfit dir>          best non-negative, non-rising profile per stack (floor) and the 24-point set -> res54_lensfloor.json
  dwarfs54_peri.py <dmfit dir>        the dwarfs' slope with size by pericenter and host (Pace, Erkal and Li 2022) -> res54_dwarfs_peri.json
Explorations, not adopted:
  explore54_profiles.py <dmfit dir>   DC14 and coreNFW on the final formula, A, B, C refitted on all disks -> res54_profiles.json
  explore54_sizec.py <dmfit dir>      concentration scaled by stellar surface density, fitted on disks -> res54_sizec.json
