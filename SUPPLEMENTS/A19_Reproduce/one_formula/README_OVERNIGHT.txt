A19 dark-matter formula search, overnight 6 to 7 November 2026 (session SA73). Code, protocols, results.

Inputs not included (unchanged from the earlier zips): the A19 data package (pkg/A19_Data_Package) and Astra's
SO73 T63 benchmark folder (so73_t63). The code expects them at the fixed path in harness_T15_hashed.py.

FINAL FORMULA
FROZEN_T45_candidate.json (frozen 10:42 PST, the final formula): T44's formula with the switch's threshold set on the
SDSS red centrals of Mandelbaum et al. (2016) instead of the benchmark; only A, B, C (disks) are fitted to the benchmark.
FROZEN_T44_candidate.json: T43's formula with the KiDS lens stacks' concentrations scaled to all haloes (Neto et al. 2007).
FROZEN_T43_candidate.json (frozen 04:29 PST): T37's formula with the switch on (q = 1) for samples selected by their halos.
The formula that took the HeCS blind test (T40) is T37's:
FROZEN_T37_candidate.json (SHA-256 in FROZEN_T37_candidate.sha256; frozen 03:55 PST). Halo formula: the paper rule +
for gas-poor galaxies above M_t, the expected halo given the galaxy's stellar mass and its sample's selection (Moster 2013
with 0.18 dex scatter x Tinker 2008 halo mass function x the survey's halo cut: CLASH 5e14, MaNGA BCGs 1e14, others none)
+ halo floor 1e8 Msun; no contraction; 4 coefficients (harness_T37.py). Baryons of spheroids with a measured dispersion
carry the Posacki et al. 2015 IMF trend (harness_T24.py, harness_T26.py). The previous candidate, with the Moster mass
inverted at its mean, is FROZEN_T26_candidate.json (01:20 PST).

ROUNDS (each PROTOCOL_Txx.txt was hashed, PROTOCOL_Txx.sha256, and locked before its run)
T15 to T16  earlier turns (see the T15 to T18 zips)
T17, T17b   switches for the giant-galaxy term; MOSS (size switch)          harness_T17.py, harness_T17b.py
T18         Cautun et al. 2020 contraction, rule refitted on disks          harness_T18.py
T19         gas-poor switch MOG; MOGC with contraction                     harness_T19.py
T20         robustness of the switch                                        harness_T20.py
T21         halo floor (MOGF, MOGF85, MOGF0)                                harness_T21.py
T22         hot gas halo (MOGH5, MOGHf, MOGFH5)                             harness_T22.py
T23         fitted switch steepness: drafted, never hashed or run           harness_T23.py
T24         IMF trend on MaNGA BCG stars (MOGFI Chabrier, MOGFK Kroupa)     harness_T24.py, _T24C.py, _T24K.py
T25         T24 without contraction (MOGF0I)                                harness_T25.py
T26         IMF trend also on the 16 early types (MOGF0IE)                   harness_T26.py
T27         T26's formula, M_t fitted without MaNGA and CLASH (MOGF0IEX)      chunk_run_lopo.py, describe27.py
T28         T26's formula, M_t fitted without CLASH (MOGF0IEY)               chunk_run_loco.py, describe28.py
T29         T26's formula at population weight 1 (MOGF0IE@1)                 chunk_run.py
T30         T26's formula against Mandelbaum et al. 2016 central halos       diag_mandelbaum.py (res30_mandelbaum.json)
T31         T30 with the Moster term inverted with scatter + Tinker 2008 HMF diag_mandelbaum_scatter.py (res31_mandelbaum_scatter.json)
T32         CLASH, scatter-aware inversion with the survey selection         diag_clash_selection.py (res32_clash_selection.json)
T33         MaNGA, scatter-aware inversion with the program's stated aim    diag_manga_selection.py (res33_manga_selection.json)
T34         sensitivity of T32 and T33 to threshold and scatter             diag_selection_sensitivity.py (res34_selection_sensitivity.json)
T35         the Milky Way with T18's contracted formulas                     diag_milkyway_contraction.py
T36         contraction growing with stellar mass (all-data look, no protocol) harness_T36.py
T37         selection-aware expected halo in the formula, held out (MOGF0IES) harness_T37.py, chunk_split.py (res37_*)
T38         frozen T37 against Mandelbaum, the KiDS split, the Milky Way     check38_frozen_T37.py (res38_frozen_T37_checks.json)
T39         T37 with CLASH out of every training set (MOGF0IESY)              harness_T39.py, chunk_run_loco.py (res39_*)
T40         blind test on 216 HeCS-omnibus clusters (Sohn et al. 2020)       blind_hecs_T40.py (res40_hecs_blind.json, res40_posthoc_notes.txt; data in ext_hecs/)
T41         halo mass function sensitivity (Tinker, Despali, Watson, Bocquet)  harness_T41.py, describe24.py, check41.py (res41_42_summary.json, log41_*)
T42         concentration relation sensitivity (DM14, Diemer-Joyce, Ludlow)    describe42.py (res41_42_summary.json, log42_*)
T43         switch on for halo-selected samples (MOGF0IESG), held out          harness_T43.py, chunk_split.py (res43_*, log43_describe.txt, res43_posthoc_hecs.txt)
(look)      global concentration offset, all data, no protocol (a look before T44) explore_cscan.py (logx_cscan_*)
T44         all-halo concentrations for the KiDS stacks (MOGF0IESR), held out  harness_T44.py, chunk_split.py (res44_*, log44_describe.txt, res44_kids_split.txt)
(look)      residuals by stack, dwarf and cluster; threshold from Mandelbaum; mixture form (all data) diag_lens45.py, diag_dwarfs45.py, diag_clash45.py, explore45.py (res45_explore.json)
T45         threshold M_t set on SDSS red centrals (MOGF0IESM), held out    harness_T45.py, chunk_split.py (res45_*, parts45_*, log45_*, res45_secondary.json, check45_*.py)

RUNNING
python3 harness_Txx.py <primary|strict> <out.json> <config@0.1>      (all 60 held-out partitions)
python3 chunk_run.py <harness module> <mask> <config@0.1> <parts.jsonl> <budget s>; then
python3 chunk_run.py merge <harness module> <mask> <parts.jsonl> <res.json>   (resumable, used on the Mac)
python3 analyze_gen.py <mask> <baseline config> <result pattern with %s> ...   (scores, pass marks, bootstraps)
python3 chunk_split.py <module> <mask> <cfg> <parts.jsonl> <budget> <k> <n>   (one run split over n processes; used to
        redo T22 MOGHf on the stricter mask after the first cloud run was stopped at 14 of 15 parts)
python3 describe24.py <harness module> <model>                        (all-data refit, primary mask)

RESULTS
res<round>_<config>_<mask>.json: summaries and all 60 records. part*.jsonl: chunked parts. log*: run logs.
Machines: cloud (T17 to T22, T24 strict parts), Gary's Mac 4-core workspace (T21 MOGF and MOGF85, T24 primary and
first strict parts, T25 to T29), Gary's server (T22 backup, slower; shared; ext/server_results_T22/ holds its two finished
runs, which match the cloud's to six decimals in every score). Same code and data files everywhere.

DIAGNOSTICS (all data, not held out)
diag_*.py: residual patterns, contraction variants, DC14, Einasto, tides (Pace et al. 2022 pericenters in ext/),
Behroozi 2013, hot gas, halo-mass scatter in the lensing stacks (diag_scatter.py), IMF what-ifs (diag_imf.py,
diag_imf_etg.py). kids_split.py: the unfitted KiDS early/late lensing check. diag_milkyway.py: the Milky Way's
rotation curve predicted by the final formula from McMillan (2017) baryons, against Eilers et al. (2019) in
ext/milkyway/eilers2019_table1.txt.

EXTERNAL DATA (ext/)
pace2022_table3.csv: Pace, Erkal and Li 2022, Table 3 (pericenters). atlas3d/Cappellari2013a_Table1_16ETG_rows.txt:
sigma_e of the 16 early types from ATLAS3D XV Table 1 (page-text SHA-256 recorded in the file). lvdb_new/: the current
Local Volume Database dwarf_all.csv from GitHub (apace7/local_volume_database), read 7 Nov 2026 01:21 PST.
server_kit/: the scripts that ran the jobs on Gary's server.
