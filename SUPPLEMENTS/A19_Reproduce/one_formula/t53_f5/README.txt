Forecast F5 test (session SA73, turn T53, 10 November 2026 PST). No protocol for the test itself; F5 was hashed at T52
(FORECASTS_T52.txt, b674324f...), F6 at T53 (FORECASTS_T53.txt, 96c0e2e4...) before any test of it.
Input: nsa_manga_bcg.csv (sha256 5d481c17...), the NASA-Sloan Atlas columns (nsa_z, nsa_elpetro_th50_r, nsa_sersic_th50,
nsa_sersic_n, nsa_elpetro_mass, nsa_sersic_mass, nsa_elpetro_ba) for the 50 plate-IFUs, extracted from the MaNGA DR17
catalog drpall-v3_1_1.fits (https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/drpall-v3_1_1.fits,
75,360,960 bytes, sha256 9d34371c157bf7d1b98eb389afd8553243af1f37f2efd2527d02fd8698d862b4). NSA masses are h^-2; divided by 0.7^2.
Run from this folder (scripts exec jeans53.py from the working directory):
  python3 -I jeans53.py <package>/data/Clusters/manga_bcg_table1_extracted.txt   # rebuild check at the tabulated radius
  python3 -I jeans53b.py <same table>                                             # radius r* from g_bar; g_obs within 0.02 dex
  python3 -I f5_full.py <dmfit dir> nsa_manga_bcg.csv                             # F5 as worded, full rebuild, robustness, MOND/standard
Exploratory, superseded by f5_full.py: f5_test.py and f5_more.py (3D Prugniel-Simien ratios at the tabulated radius),
f5_split.py (stars/halo split), f5_contract.py (Cautun contraction with and without the IMF factor).
Results: res53_f5_full.json. F5 as worded: correlation 0.509 (fails, > 0.45). Full rebuild: 0.0872 dex, correlation 0.284.
