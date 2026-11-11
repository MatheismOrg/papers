A core in the halos of diffuse gas-rich disks (session SA73, turn T61, 10 November 2026 PST; Gary 15:41 "can u improve to
beat mond?"). PROTOCOL_T61.txt (hash b370fbdc..., written 15:48 PST before the held-out run), harness_T61.py,
T61_script.sha256, res61_MOGF0IESBC_{primary,strict}.json and FORECAST_T61.txt (the eighth forecast) are one folder up.
Run from this folder with PYTHONPATH set to the folder above.

Model MOGF0IESBC: T48 plus, for light-selected systems, M_ref x 10^(t w) and M_DM(<r) = M_NFW(<r) tanh(r / 1.75 R_e)^n,
n = tanh(kappa w), w = g L (gas share times diffuseness, as T57). A, B, C as T48 on the training disks' benchmark points;
t, kappa on all points of the training disks, weighted by 1 / (sigma_V^2 + 0.131^2).

All-data looks before the protocol (not tests):
  explore61.py    looser halos against cores, fitted on the benchmark's points or on all points
  explore61b.py   the trade between kept and dropped points as kappa grows; error-weighted fits
  explore61c.py   the gate: g L against g alone
  describe61.py   the protocol's fit on all 149 disks: t 0.1501, kappa 0.6229 -> res61_describe.json
Hashed run: python3 ../chunk_split.py harness_T61 <primary|strict> MOGF0IESBC@0.1 <parts file> 6000 0 1, then merge.
  Benchmark, held out: disks 0.15306 (T48 0.15326), lensing 64.52 and 39.37 (T48 64.67 and 39.43), all else identical.
  heldout61.py <res61 primary json>   marks (b), (c) -> res61_heldout.json
     dropped points of held-out disks (114 galaxies): T48 0.4007, T61 0.3439 (-0.077 to -0.039), MOND 0.3639
     (T61 - MOND -0.031 to -0.009), standard halos 0.4074; all points (149): 0.2181, 0.1926 (-0.038 to -0.014), MOND 0.1996
     (-0.015 to +0.002), 0.2289; slower than 30 km/s: 0.548, 0.418, 0.459. t 0.11 to 0.18, kappa 0.54 to 0.71, no bound hits.
     ALL PRESET MARKS MET: T61 is a candidate to replace T48, for Gary to decide.
  outside61.py <t59 folder>   outside checks (d), (e), everything frozen -> res61_outside.json
     (d) Iorio et al.'s eight outside SPARC, every point: T48 0.4108, T61 0.3079, MOND 0.414 (T61 - MOND -0.181 to -0.032)
     (e) 26 SPARC galaxies outside the benchmark, every point: 0.4351, 0.392, 0.4146 (-0.036 to -0.007)
     with the benchmark's cut: Iorio's three 0.1897, 0.1959, MOND 0.169 (worse); SPARC's 17: 0.3169, 0.3163, MOND 0.3376
  describe61b.py  description with the frozen all-data values: all 175 galaxies, dropped 582 points 0.4156 -> 0.3596
     (MOND 0.381), all 3,389 points 0.2603 -> 0.2310 (MOND 0.2438), scatter within disks 0.0733 -> 0.0664.
After the runs, not in the protocol, in MOND's favor:
  heldout61m.py <res61 primary json>, outside61m.py <t59 folder>   MOND's relation with a0 refitted on the same training
     points and weights (a0 0.96 to 1.03e-10): held out, kept 0.1532, dropped 0.3522 (T61 0.3439; -0.018 to +0.001),
     all 0.1954 (T61 0.1926; -0.009 to +0.005); outside (d) 0.3925 (T61 0.3079; -0.154 to -0.017), (e) 0.3963 (T61 0.392;
     -0.016 to +0.007). Against the refitted constant T61 is lower everywhere but beats it only on (d).
../t59_outside/drift61.py: SPARC's speeds at the dropped points against Iorio et al.'s drift-corrected ones (ratio 1.00).
