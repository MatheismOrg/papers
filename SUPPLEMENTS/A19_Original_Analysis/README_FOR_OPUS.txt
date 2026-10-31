ASTRA DM AND A42 RECOVERY HANDOFF FOR OPUS
30 October 2026

THE REQUESTED MATERIAL WAS FOUND

This ZIP supplies the dm_ directories and all A42 passes 10-18 requested in
the quoted SA72 T357 message. It also retains the other A42 passes through22
and sibling DM branches so imports, comparison inputs and later corrections
are available. The scientific files are copied without editing their bytes.

START HERE: TOPIC MAP

Feedback profiles / requested A19 Section3.4:
  DM/A43_Dark_Matter_Fits/dm_feedback_profiles/
  README.txt, PROTOCOL.txt, run.py, audit.py, input.csv, galaxy_features.csv,
  results/summary.json, complete per-fold results and frozen predictions.
  Keep sibling dm_cosmological_halo/: its NFW_Planck_size_cv_*.csv files
  provide the original fold assignments. All siblings are included here.

Early-type galaxies / requested A19 Section5:
  DM/A43_Dark_Matter_Fits/dm_broad_groups/
  validate_external.py, data/ETG_Lelli2017.mrt, data/Rotmod_ETG.zip,
  results/etg_external_pairs.csv and external_validation_results.json.
  This archived branch describes16galaxies with32inner/outer measurements.
  The same folder includes the group and KiDS/GAMA lensing checks.

Clusters / requested A19 Section6:
  DM/A43_Dark_Matter_Fits/dm_fit/data/xcop/ (cached FITS inputs)
  DM/A43_Dark_Matter_Fits/dm_fit/test_clusters.py
  DM/A43_Dark_Matter_Fits/dm_fixed_radius/
  DM/A43_Dark_Matter_Fits/dm_broad_lensing/ (CLASH, MaNGA tables and preparation)
  DM/A43_Dark_Matter_Fits/dm_cluster_repair/ (gas inputs, model, diagnostics)
  DM/A43_Dark_Matter_Fits/dm_best_audit/ (independent cluster checks)
  dm_gas_bridge/ and the broader fit branches are also retained.
  See the individual provenance files for the target estimators and assumptions.

A42 requested passes10-18 / Section8.6 verification:
  A42/Research/Pass10/ through A42/Research/Pass18/
  A42_PASS_10_TO_18_INDEX.json lists every Python script and available guide.
  Internal paths inside each pass are unchanged; a few have an extra original
  package-root folder. Earlier/later passes are supplied for dependencies and
  interpretation. There is no need to unpack a long chain of prior-pass ZIPs.

REPRODUCE FROM THE DM/A43_Dark_Matter_Fits DIRECTORY

  OPENBLAS_NUM_THREADS=1 python3 dm_feedback_profiles/run.py
  OPENBLAS_NUM_THREADS=1 python3 dm_feedback_profiles/audit.py
  python3 dm_broad_groups/validate_external.py

For cluster runs and A42 runs use the corresponding README; do not assume a
single shared environment for every experiment. The DM requirements.txt and
individual environment records are retained. Full fitting can take time and
some scripts overwrite adjacent result files: use a working copy for reruns.

VERSION AND CLAIM BOUNDARIES

The dm_ files were recovered from the saved A43_Code_and_Data.zip (version9).
They remain labelled A43 because that is their source-project name. The topic
mapping above follows Opus's descriptions; the named SA72_T357 A19/A42 revised
PDFs were not supplied here, so their exact section text was not compared.
Do not rename source data to imply an exact match to an unseen revised table.

In particular the older DM branches contain a44-dwarf benchmark and older
coefficients. Opus's quoted rebuild uses42dwarfs. Both are different snapshots;
this archive does NOT replace, validate or reverse that newer rebuild. Preserve
its selection and corrected coefficients when checking the newly recovered
branches. The reported40% improvement and other claims in the quoted message
were not re-evaluated in this packaging task.

The A42 manuscript here is the saved53-page Pass18 version. Passes19-22 remain
supplements. These older manuscripts and cover letters are reference records,
not replacements for Opus's latest A19/A42 paper or book edits.

The current DM archive contains every DM member in the latest recovered main
portfolio. The only differing shared members are its paper PDF, paper Markdown
and current manifest; the later CLAIM_REVISION.txt describes an editorial
revision. No additional missing DM script was found in that portfolio copy.
Thus the cited folders survived, even if absent from the copy Opus had.

Two earlier small DM archives are retained under DM/Historical_Archives/.
They are provenance, not preferred current results. Other portfolio books and
unrelated academic papers are excluded.

VERIFICATION AND LIMITS

Every copied file is checked against the source bytes. The final ZIP is CRC
tested and PACKAGE_SHA256SUMS.txt covers its payload. SOURCE_MANIFEST.json
records each original archive and member, and the source archive hashes.
Old internal manifests still describe their original layouts; use the new
manifest/checksums for this handoff's renamed top-level folders.

This is file recovery and packaging, not a new numerical replication. No
downloaded scientific script was executed for the task. Cached cluster and
ETG numerical inputs are present, so Perplexity's downloader was not needed
for these recovered calculations. Original research article PDFs, every
possible cluster catalog, and full-resolution field dumps are not promised;
the source-package READMEs specify their actual retained inputs and limits.
