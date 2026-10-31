# A19: original analysis code, inputs and saved results

The analysis branches (DM/A43_Dark_Matter_Fits/dm_*) that produced the numbers marked "Rerun" or "Saved" in the number map of A19_Reproduce/README.md: the curved relation, the feedback profiles, the population scatter, the kernel integrations, the cross-population fits, the cluster repair, the gas bridge and the cluster audits. Version 9 of the original code archive, copied without editing its bytes.

Each branch keeps its code, protocol, inputs and saved results. Large per-fold prediction tables and figures are left out to keep the repository small (LEFT_OUT.txt lists them, with sizes); the scripts regenerate them. dm_broad_groups/data/kids/ holds the KiDS-1000 lensing data of Brouwer et al. (2021). DM_DIRECTORY_INDEX.txt lists the branches, and README_FOR_OPUS.txt maps paper sections to branches and gives the run commands. Some scripts overwrite their saved results, so rerun them in a copy.

The independent rebuild, from public data with separate code, is A19_Reproduce.
