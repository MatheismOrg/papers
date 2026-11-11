# A19, the one formula: code, protocols, forecasts and saved results (protocols 10 to 21)

This folder holds what Section 6 of "Dark Matter Behaves Like Matter" rests on (version of 11 November 2026): the harnesses that build the one formula in steps and score it on held-out systems (Table 7), the comparison with MOND's relation and standard halos and with the noise of the measurements (Tables 8 and 9), the test on every rotation-curve point (Table 10), the tests of three kinds of dark-matter particle and of star-formation history (Section 6.5), and the ten forecasts (Table 11).

Protocols 1 to 9, the public data they use and the independent rebuild are one level up, in `A19_Reproduce`. The original analysis branches are in `A19_Original_Analysis`.

## What is here

| Files | What they are |
|---|---|
| `harness*.py`, `chunk_*.py`, `final_T62.py` and the other `.py` files | The code. Two modules, `harness.py` and `harness_T15_hashed.py`, name the input folders; every later harness imports them. |
| `PROTOCOL_Txx.txt` with `PROTOCOL_Txx.sha256`, `Txx_script.sha256` | Each protocol as it was written before its run, with its SHA-256. The script hashes fix the code that ran. |
| `FORECAST*.txt` with `.sha256` | The forecasts, each hashed before any test of it. |
| `FROZEN_Txx_candidate.json`, `DECISION_T62.txt` | The coefficients frozen at each step, and the note that made the core part of the formula. |
| `res*.json`, `parts*.jsonl`, `log*.txt`, `merge*.txt` | Saved results, per partition and merged, with the run logs. |
| `t53_f5/` to `t65_history/` | The later rounds, each with its own `README.txt`. |
| `ext/`, `ext_hecs/` | Small external tables. `DATA_SOURCES.md` gives each one's source and terms. |
| `astra_T63/T56_fold_maps.json` | The fold assignments of the dwarfs, early types, MaNGA galaxies and clusters. |
| `tools/make_runtree.py` | Makes a folder from which the scripts run (see "Run it"). |
| `README_OVERNIGHT.txt` | The working list of rounds T15 to T45. |

Round numbers (T15 to T65) are the turns of the working session in which each step was made. The paper numbers its protocols 1 to 21; the map below joins the two.

## From the paper to the files

**Sections 6.1 and 6.2, Table 7.** Sixty partitions: three assignments of the disks to five folds, each crossed with four held-out lensing stacks; both lensing masks.

| Table 7 row | Paper protocol | Round and files |
|---|---|---|
| 1, the stellar rule | | `harness_T15_hashed.py`; earlier rounds T15 to T18 |
| 2 and 3, switch, mean Moster halo, floor; then the IMF trend | 10 | T19, T21, T24 to T26 (`harness_T19.py`, `harness_T21.py`, `harness_T24.py`, `harness_T26.py`); checks T27 to T30 |
| 4, the expected halo with each sample's selection | 11 | T37 (`harness_T37.py`); checks T38, T39; all-data looks T31 to T34 |
| The blind test on 216 HeCS-omnibus clusters | 12 | T40 (`blind_hecs_T40.py`, `ext_hecs/`, `res40_posthoc_notes.txt`) |
| 5, the switch on for halo-selected samples | 13 | T43 (`harness_T43.py`) |
| 6, all-halo concentrations in the lensing stacks | 14 | T44 (`harness_T44.py`); the look before it, `explore_cscan.py` |
| 7, the threshold set on SDSS red centrals | 15 | T45 (`harness_T45.py`, `res45_MOGF0IESM_*.json`, `FROZEN_T45_candidate.json`) |
| 8, hot gas among the clusters' baryons | 16 | T48 (`harness_T48.py`, `res48_MOGF0IESB_*.json`, `FROZEN_T48_candidate.json`) |
| 9, cores in diffuse gas-rich disks | 19 | T61 (`harness_T61.py`, `res61_MOGF0IESBC_*.json`), `DECISION_T62.txt`, `final_T62.py` |

**Section 6.3, Tables 8 and 9.** `compare46.py`, `compare48.py` and `compare62_boot.py` (`res62_boot.json`) score MOND's relation and standard halos on the same benchmark, with the paired intervals. The noise tests and the five patterns: `noise50_dwarfs.py`, `noise50_clusters.py`, `t54_noise/` and `t62_final/`.

**Section 6.4, Table 10.** Protocols 17 and 18, the two terms that failed: `harness_T56.py` and `harness_T57d.py`, with `t56_search/`. Protocol 19, the core: `harness_T61.py` and `t61_core/` (`heldout61.py` gives the held-out rows, `res61_outside.json` the galaxies outside the benchmark, the two `*_mond_refit.json` files the refitted MOND rows). The galaxies outside SPARC, the dropped points and their noise: `t59_outside/`.

**Section 6.5.** The three alternatives on the benchmark's points: `sidm46.py`, `sidm46_bench.py`, `sidm46_heldout.py`, `contract46.py`, `sizec46.py`. Protocol 20, fuzzy, warm and self-interacting dark matter in the core's place: `harness_T64.py`, `res64_*.json`, `t64_particles/`. Protocol 21, the history of star formation: `harness_T65.py`, `res65_*.json`, `t65_history/`.

**Section 6.6, Table 11.** Forecasts 1 to 4: `FORECASTS_T51.txt`. Forecast 5: `FORECASTS_T52.txt`; its test: `t53_f5/`. Forecast 6: `FORECASTS_T53.txt`. Forecast 7: `FORECAST_T60.txt`, Part B (hash in `T60_forecast.sha256`). Forecast 8: `FORECAST_T61.txt`. Forecast 9: `FORECAST_T64.txt`. Forecast 10: `FORECAST_T65.txt`. `FORECAST_T58.txt` and `FORECAST_T59.txt` are the checks of the looser-halo term, which was not adopted.

## Run it

Python 3 with numpy and scipy; some scripts also need colossus (`pip install colossus`).

The scripts are byte-identical to the ones that were hashed, so five of them still name the base folder of the session in which they ran. `tools/make_runtree.py` copies the code to a folder of your choice and replaces that one path in the copies. Nothing in the repository is changed.

    python3 tools/make_runtree.py /path/to/new/folder
    cd /path/to/new/folder/dmfit
    PYTHONPATH=. python3 chunk_split.py harness_T61 primary MOGF0IESBC@0.1 my_parts.jsonl 6000 0 1
    PYTHONPATH=. python3 chunk_split.py merge harness_T61 primary my_res.json my_parts.jsonl
    python3 ../compare_runs.py my_res.json res61_MOGF0IESBC_primary.json
    cd t61_core && PYTHONPATH=.. python3 heldout61.py ../my_res.json

The first three commands rerun the last row of Table 7, the final formula held out on the main lensing mask, in a minute or two. The last command gives the held-out rows of Table 10. Use `strict` in place of `primary` for the stricter mask, and another harness module and model name for another row (each protocol names its own).

Checked on 11 November 2026 from a fresh copy of this repository: the twelve benchmark numbers of the final formula came out identical to the saved ones (disks 0.153058, dwarfs 0.37022, early types 0.185164, MaNGA 0.115982, cluster centers 0.10198, clusters at 100 to 600 kpc 0.100181, lensing chi-square 64.5199), and so did every score of Table 10's held-out rows.

## What the files still say

The protocols, forecasts and scripts are the files that were hashed, unchanged. They therefore carry the working notes of the sessions in which they were written: session and turn numbers, the author's requests quoted as given, and the machines used. As the earlier protocols in `A19_Reproduce/blind_test` state, the protocols and code were written with an AI model, Claude, at the author's direction.

## What is left out

Scripts that queued jobs on a private server, the tools that build the manuscript, and raw third-party source files. For the last, `t59_outside/SOURCES_bundle59.txt` and `t65_history/SOURCES.txt` give each file's address, size and SHA-256, and the values read from them are in `t59_outside/lt59_inputs.json` and `t65_history/sfh65_inputs.json`.
