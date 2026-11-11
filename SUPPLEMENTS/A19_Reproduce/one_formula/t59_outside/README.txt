Checks of the T57 term beyond the benchmark (session SA73, turns T59 and T60, 10 November 2026 PST). All coefficients
frozen: T48 A, B, C = 11.7481, 0.8178, -0.7119; T57 adds t = 0.361, s = 0.146 (harness_T57d.py, one folder up). The
formula of the paper stays T48. Run from this folder with PYTHONPATH set to the folder above.

Sources, downloaded 10 November 2026 (lt59_inputs.tar.gz, sha256 ee13ce4a...; SOURCES.txt inside it):
  Iorio et al. (2017, MNRAS 466, 4159): finalrot.zip, 28,572 bytes, sha256 d4c456c6..., from
    https://www.dropbox.com/s/t4j8dacmnwgj0yb/finalrot.zip?dl=1 (linked on https://www.filippofraternali.com/downloads)
  Dale et al. (2009, ApJ 703, 517): Tables 1 and 2 from the arXiv source, https://arxiv.org/e-print/0907.4722
    (2,979,216 bytes, sha256 f0b7e629...)
  Read et al. (2017, MNRAS 467, 2019): Table 2 from the arXiv source, https://arxiv.org/e-print/1607.03127
    (2,042,016 bytes, sha256 fc65999f...)
Nothing in the downloads is executed; parse59.py reads them with python3 -I.

  parse59.py <bundle dir> lt59_inputs.json   the three sources into one file
  model59.py        gas (razor-thin ring sum from the HI surface density, times 1.33), stars (exponential disk), L[3.6];
                    run alone it checks the ring sum against Freeman's disk
  validate59.py     the inputs on the seven galaxies shared with SPARC, no formula -> res59_validate.json (medians: LVL
                    luminosity 1.172 of SPARC's, 1.678 R_* 1.424 of SPARC's R_eff, HI mass 0.97)
  test59.py         FORECAST_T59 (hash 54e8b138..., written 15:07 PST before any outside galaxy was scored) -> res59_test.json
                    As worded (benchmark cuts): 3 of 8 outside galaxies keep points; T48 0.1897, T57 0.1966 (+3.6 percent;
                    paired range -0.020 to +0.033); MOND 0.169, standard halos 0.188. FORECAST FAILED as worded.
                    Clean subsample: one galaxy (DDO 52), lower. Union with the 17 of T58: -2.0 percent (-0.016 to +0.003).
  explore59_cuts.py after the fact, not a test: the same with looser cuts on the velocity error -> res59_explore_cuts.json
                    No cut: 8 galaxies, 86 points, T48 0.4108, T57 0.3571 (-13 percent; -0.085 to -0.021; 7 of 8 lower),
                    MOND 0.414, standard halos 0.459.
  dropped60.py      FORECAST_T60 Part A (hash 30057ffe..., written before it was run): the 582 SPARC points, in 139
                    galaxies, that the benchmark's 10 percent cut removes -> res60_dropped.json
                    T48 0.4156, T57 0.3922 (-5.6 percent; -0.030 to -0.018). MET. Not an outside test: the galaxies are
                    SPARC's. MOND 0.381, standard halos 0.4257: MOND's relation is best on these points.
  describe60.py     description, no forecast: the dropped points by error, radius and speed -> res60_describe.json
                    (slower than 30 km/s: 0.565, 0.471, 0.578 for the formula, MOND, standard halos; faster: 0.323, 0.326,
                    0.326; all 3,389 points, no cut: 0.260, 0.244, 0.274)
  noise60_dropped.py 2000   truth = the formula; inputs and velocities redrawn, the cut applied to the redrawn velocities
                    -> res60_noise_dropped.json (dropped points: 0.251, 0.215 to 0.314, mean -0.01; the formula would beat
                    MOND's relation by 0.030, 0.014 to 0.045)
FORECAST_T60 Part B is the seventh forecast of Section 5.7 (data not in hand).
  drift61.py        no formula: on the seven galaxies shared with SPARC, SPARC's speeds at the points the benchmark drops,
                    brought to Iorio et al.'s inclination, over their drift-corrected circular speeds: median 1.00 (59
                    points), 1.05 below 30 km/s (30 points) -> res61_drift.json. Pressure support is not the cause.
The core that followed (PROTOCOL_T61) is in ../t61_core/.
