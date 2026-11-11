Search for a fit to the pattern inside small gas-rich disks (session SA73, turns T56 to T58, 10 November 2026 PST).
The final formula stays T48 (A, B, C = 11.7481, 0.8178, -0.7119; log M_t = 11.339). Nothing here is adopted.
The harness modules and the protocols are one folder up: harness_T56.py, PROTOCOL_T56.txt, T56_script.sha256,
res56_MOGF0IESBU_{primary,strict}.json, res56_analysis.json; harness_T57.py (all-data look only), harness_T57d.py,
PROTOCOL_T57.txt, T57_script.sha256, res57_{MOGF0IESBD,MOGF0IESBD5}_{primary,strict}.json, res57_analysis.json;
FORECAST_T58.txt, res58_excluded.json. Run with PYTHONPATH set to that folder.

All-data looks, no protocol (what each showed is in the script's docstring and in work/OVERNIGHT_PLAN.md):
  mb56.py        the 14 Mandelbaum et al. (2016) halo masses: red chi-square 61.2 for 7, blue 12.0 for 7 (p 0.10)
  szc56.py       size acting on concentration: disks 0.1478 at k = 0.4, dwarfs 0.358 at k = 0.2
  ctilt56.py     concentration tilted with stellar mass (0.1487), gas share (no gain) or surface brightness
  sig56.py       concentration from stellar surface density: disks 0.1458 but dwarfs 0.46 to 0.53 (B falls to 0.26 to 0.45)
  core56.py      cores where cold gas remains: 0.1492; weighted by surface density 0.1474
  gate56.py      looser halos for diffuse disks with and without a cold-gas gate
  gasmu56.py, mu56.py, describe56.py   cold gas weighted mu times in the baryon term (the look before PROTOCOL_T56)
  lensgas56.py   lensing chi-square against the cold gas assumed for the stacks
  lscat56.py     halo-mass scatter in the lensing stacks: 63.2 to 58.8 at best
  describe57.py  gas making every halo bigger and looser: disks 0.1476, dwarfs 0.352, lensing 79.4
  diffuse57.py   the same confined to diffuse gas-rich galaxies (g L): disks 0.1476, within disks 0.0733 to 0.0679, lensing 62.4

Hashed runs:
  analyze56.py   PROTOCOL_T56 (hash 08d5fcb1...), held out: disks 0.1533 to 0.1542, dwarfs 0.3702 to 0.3565 (-0.035 to +0.008),
                 early types 0.1852 to 0.1912, lensing 64.7 to 75.7 (fails its tolerance); mu 2.3 to 7.2. FAILED.
  analyze57.py   PROTOCOL_T57 (hash 9ccdb892...), held out: disks 0.1533 to 0.1513 (-1.26 percent; paired range -0.0054 to
                 +0.0017), lensing 64.67 to 63.84, nothing else changed; t 0.24 to 0.53, s 0.10 to 0.21. Mark MISSED (the range
                 includes zero). Five coefficients fitted together: disks -0.02 percent, dwarfs +5.6 percent.
  excluded58.py  FORECAST_T58 (hash 52ce22af...): the 17 SPARC galaxies the benchmark excludes, all frozen (t = 0.361,
                 s = 0.146): 0.3169 to 0.3095 (-2.3 percent; -0.018 to +0.003). Lower, not significant.
The checks on galaxies outside SPARC and on the points the benchmark's cut removes are in ../t59_outside/.
