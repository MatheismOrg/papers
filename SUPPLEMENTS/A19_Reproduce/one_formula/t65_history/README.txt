Does the core follow the gas or the history of star formation? (session SA73, turn T65, 10 November 2026 PST; the
twenty-first protocol, ../PROTOCOL_T65.txt, hashed 23:22 PST before any score, SHA-256 0123d59c...; script and input
hashes in ../T65_script.sha256; the tenth forecast, ../FORECAST_T65.txt, hashed 23:26 PST, c9c75c03...).

Verdict under the protocol's preset reading: not decided by these data.
 Test 3, disks (../harness_T65.py, MOGF0IESBH: gate "any cold gas", w = L, in place of w = g L; t and kappa refitted as in
   T61): presence is enough. Held out: benchmark points 0.1530 (gas gate 0.1531), dropped 0.3399 (0.3439; -0.010 to
   +0.001), all points 0.1888 (0.1926; -0.009 to +0.001), lensing 64.42 / 39.34 (64.52 / 39.37); outside 0.3074 and
   0.3921 (0.3079, 0.3920). t 0.084 to 0.117, kappa 0.36 to 0.45; all data t = 0.1007, kappa = 0.4003.
 Test 1, dark-matter density at 150 pc (Read, Walker and Steger 2019), 8 dwarf spheroidals: formula (cusp) 0.212 dex,
   published history recipe (Read, Agertz and Collins 2016, nothing fitted) 0.133, formula's core with a history gate
   0.138, self-interaction at 1 cm^2/g 0.246. History recipe minus formula -0.079 (-0.227 to +0.064): not decided.
   History beats self-interaction (+0.020 to +0.189). 8 dwarf irregulars (not independent): no core 0.806, formula 0.378
   (0.32 dex too high on average), history recipe 0.258, self-interaction 0.250.
 Test 2, the benchmark's points of 20 dwarf spheroidals with tau_90 (Weisz et al. 2015): formula 0.2611, history recipe
   0.2675, formula's core with a history gate 0.2604, self-interaction 0.3254 (+25 percent, +0.012 to +0.121). Amplitude of
   the recipe's shift in the residuals 0.33 (-0.74 to +1.36). Rank correlation of residual with t_SF -0.22 (p = 0.33).

Run from this folder with PYTHONPATH=..:<t59 folder> (model59.py and lt59_inputs.json):
  python3 -I parse65.py <DM_heating.tex> <ms.tex> sfh65_inputs.json    the two published tables (SOURCES.txt)
  python3 unit65.py ; python3 smoke65.py        limits and inputs only; nothing compared with a measurement
  python3 alldata65.py                          all-data t, kappa under the gate "any cold gas" -> res65_alldata_fit.json
  (in ..) python3 chunk_split.py harness_T65 primary MOGF0IESBH@0.1 parts65_H_primary.jsonl 3000 0 1, then merge; and strict
  python3 heldout65.py ../res65_MOGF0IESBH_primary.json ../res61_MOGF0IESBC_primary.json     -> res65_heldout.json
  python3 outside65.py <t59 folder> res65_alldata_fit.json                                   -> res65_outside.json
  python3 dens65.py <t59 folder> sfh65_inputs.json res65_alldata_fit.json                    Test 1 -> res65_density.json
  python3 points65.py sfh65_inputs.json res65_alldata_fit.json ../t64_particles/res64_sidm_points.json   Test 2 -> res65_points.json
  python3 look65.py                             after the runs, an all-data look (not a test): the published recipe with
                                                nothing fitted in every disk: dropped 0.323, benchmark points 0.164
  python3 check65.py .. <paper .src> <report .src>      every quoted number read back
Readings are defined in lib65.py (N no core; G the formula; H1 published history recipe; H2 the formula's core with a
history gate; P self-interaction).
