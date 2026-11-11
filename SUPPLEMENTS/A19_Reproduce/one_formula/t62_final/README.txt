The final formula with its core (session SA73, turn T62, 10 November 2026 PST). Gary's decision is in ../DECISION_T62.txt.
The model is PROTOCOL_T61's (../harness_T61.py, MOGF0IESBC). ../final_T62.py freezes it at the values from all 149
benchmark disks (A, B, C = 11.7481, 0.8178, -0.7119; t = 0.1501; kappa = 0.6229; log M_t = 11.339) behind the interface
the earlier check scripts use, so each script here is the T54 or T60 script with one import changed. In their outputs the
key "T48" or "formula" now means the final formula with its core.

  python3 -I noise62_disks.py .. 2000     Table 9, disks: noise 0.144 (0.122 to 0.176); shares of draws at least as large as
                                          the formula's 0.1531, MOND's 0.1546, standard halos' 0.1655: 24, 22, 8 percent
  python3 -I patterns62_disks.py .. 2000  Section 5.5, first pattern: within disks 0.066 against 0.045 (0.043 to 0.048);
                                          galaxy means 0.134 against 0.136; no trend with any property
  python3 -I diag62_within.py ..          within disks: formula 0.0664, MOND 0.0716, standard halos 0.0804; formula 0.037 dex
                                          too high inside 0.3 R_e; scatter of per-galaxy slopes 0.22 (velocity errors 0.07)
  python3 -I lensfloor62.py ..            lensing on all data: 63.04 (29 points), 38.08 (25), 21.25 (24 clean points; MOND
                                          37.18, standard halos 68.67)
  python3 -I check62.py ..                KiDS early-late split 0.0507 (0.0509 without the core); the negative lensing point,
                                          -4.10 sigma, costs 16.8; held-out lensing shares 0.016 and 3.4 percent; w = g L is
                                          zero for the blue SDSS centrals above 10^10.5 suns and for the Milky Way
  PYTHONPATH=.. python3 noise62_dropped.py 2000   Section 5.6: noise on the points the benchmark drops, the cut applied to
                                          the redrawn velocities: 0.273 (0.235 to 0.352); the formula's 0.360 is reached by 1.9
                                          percent of draws, MOND's 0.381 by 1.1, standard halos' 0.426 by 0.45; in the draws the
                                          formula beats MOND's relation by 0.015 (0.001 to 0.028), observed 0.021
  python3 ../compare62_boot.py (run in ..)        Table 8 intervals: disks against MOND -0.011 to +0.008, against standard
                                          halos -0.023 to -0.002; the other scales as before -> ../res62_boot.json
Held-out scores of the final formula: ../res61_MOGF0IESBC_primary.json and _strict.json (PROTOCOL_T61's run);
dropped and all points, held out and outside: ../t61_core/ (Table 10).
