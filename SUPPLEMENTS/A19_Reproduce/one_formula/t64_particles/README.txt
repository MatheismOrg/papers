Three kinds of dark-matter particle in the core's place (session SA73, turn T64, 10 November 2026 PST; the twentieth
protocol, ../PROTOCOL_T64.txt, hashed 22:46 PST before any held-out run, SHA-256 40d5bbdd...; script hashes in
../T64_script.sha256; the ninth forecast, ../FORECAST_T64.txt, hashed 22:57 PST, 6cb4c1d8...).

Question: the final formula's core is set by gas and diffuseness. Fuzzy, warm and self-interacting dark matter each change
the inside of a halo by a rule that follows the halo, with one number that must be the same in every halo. Does any of
them, with its number set on the training disks and carried unchanged to every population, do what the gated core does?
Answer: no. Fuzzy is not chosen (cold end in 60 of 60 partitions). Warm is chosen at 2.0 to 2.4 keV, helps the disks a
little (all points 0.218 -> 0.210; MOND 0.200; core 0.193) and more than doubles the dwarf error (0.370 -> 0.769).
Self-interaction is chosen at 1 cm^2/g in 13 of 15 disk partitions and meets marks (b) and (c) (dropped points -14
percent, all points -9 percent), but fails (a): benchmark points +5.7 percent, dwarfs +25 percent, CLASH centers +26
percent.

Models: ../harness_T64.py (MOGF0IESBF fuzzy, MOGF0IESBW warm, on T48's formula); sidm64.py with ../sidm46.py (unchanged).
Run from this folder with PYTHONPATH=.. (and ../t59_outside for the outside dwarfs' inputs):

  python3 unit64.py            cold limits reproduce T48; soliton core masses match Schive et al. (2014)
  python3 explore64.py         all-data look before the protocol (not a test) -> res64_explore.json
  python3 alldata64.py         all-data fit of the two numbers -> res64_alldata_fit.json (fuzzy: cold end; warm 10^9.16)
  (in ..) python3 chunk_split.py harness_T64 primary MOGF0IESBF@0.1 parts64_F_primary.jsonl 3000 0 1, then merge; the
          same for strict and for MOGF0IESBW -> ../res64_MOGF0IESB{F,W}_{primary,strict}.json (benchmark scores, mark a)
  python3 heldout64.py ../res64_MOGF0IESBF_primary.json ../res64_MOGF0IESBW_primary.json ../res61_MOGF0IESBC_primary.json
                               marks (b), (c) for fuzzy and warm -> res64_heldout.json
  python3 sidm64.py 2 <t59 folder>   self-interaction on all 344 hosts and 8 cross-sections (4 minutes) -> res64_sidm_points.json
  python3 sidm64_heldout.py res64_sidm_points.json ../res48_MOGF0IESB_primary.json
                               cross-section chosen on the training disks; every population held out -> res64_sidm_heldout.json
  python3 sidm64_check.py res64_sidm_points.json ../res48_MOGF0IESB_primary.json
                               solver failures (9 early types, 3 massive disks at s0 = 1) and scores without them
  python3 outside64.py <t59 folder> res64_alldata_fit.json     outside checks (d), (e) -> res64_outside.json
  python3 describe64.py        soliton joined at 2, 3, 4 core radii; what each population alone would choose; the warm
                               number in keV; inner residuals in halos of like mass -> res64_describe.json
  python3 check64.py .. <paper .src> <report .src>             every quoted number recomputed

Caveats fixed in the protocol: retrospective; the gated core is the incumbent (gate chosen on these data, two numbers
against one); the recipes are the published dark-matter-only ones (no pull of baryons on the soliton, no core collapse);
for self-interaction A, B, C are at their all-disk values. One finding against a literal gate: in halos of 10^9.7 to
10^10.3 suns the formula without a core predicts too much inside gas-rich diffuse disks by 0.26 +- 0.10 dex and in the
eight dwarf spheroidals there by 0.14 +- 0.07; the 67 smaller dwarf spheroidals show none (-0.02 +- 0.05).
