"""SA73 T54 exploration (no protocol; a look at forecast F2 with the data in hand): the dwarfs' residual trend with
half-light radius at fixed stellar mass, split by pericenter (Pace, Erkal and Li 2022, with the LMC) and by host.
Usage: python3 dwarfs54_peri.py <dmfit dir>"""
import sys, csv, json; sys.path.insert(0, sys.argv[1])
import numpy as np
import harness_T48 as M8
M = M8.P; H = M.H; B = M.B; B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
AB = {'Ant': 'antlia', 'Aqu': 'aquarius', 'Boo': 'bootes', 'CVn': 'canes_venatici', 'Car': 'carina', 'Col': 'columba', 'CB': 'coma_berenices',
      'Cra': 'crater', 'Dra': 'draco', 'Eri': 'eridanus', 'For': 'fornax', 'Gru': 'grus', 'Her': 'hercules', 'Hor': 'horologium', 'Hyd': 'hydra',
      'Hyi': 'hydrus', 'Leo': 'leo', 'Peg': 'pegasus', 'Phx': 'phoenix', 'Pis': 'pisces', 'Ret': 'reticulum', 'Sgr': 'sagittarius',
      'Scl': 'sculptor', 'Seg': 'segue', 'Sxt': 'sextans', 'Tri': 'triangulum', 'Tuc': 'tucana', 'UMa': 'ursa_major', 'UMi': 'ursa_minor', 'Wil': 'willman'}
ROM = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}
def key(n):
    p = n.split(); return f"{AB[p[0]]}_{1 if len(p) == 1 else (ROM[p[1]] if p[1] in ROM else int(p[1]))}"
PERI = {key(r['name']): float(r['rperi_kpc']) for r in csv.DictReader(open(H.S / 'ext/pace2022_table3.csv', encoding='utf-8-sig'))}
rows = H.dwarf_rows(); a = H.arr(rows)
rr, _ = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(H.PKG / 'blind_test')); R = {x['key']: x for x in rr}
def f(x, k, d=np.nan):
    try: return float(x[k])
    except Exception: return d
es = np.array([.5 * (f(R[n], 'vlos_sigma_em', np.nan) + f(R[n], 'vlos_sigma_ep', np.nan)) / f(R[n], 'vlos_sigma') for n in a['name']])
host = np.array([R[n].get('host', '') for n in a['name']]); peri = np.array([PERI.get(n, np.nan) for n in a['name']])
e = np.log10(a['go'] / H.gpred(M.MODEL, th, a)); lm = np.log10(a['ms']); lr = np.log10(a['re'])
print('hosts:', {h: int((host == h).sum()) for h in np.unique(host)}, '| with pericenter', int(np.isfinite(peri).sum()), '| sigma err < 10%:', int((es < .1).sum()), '< 20%:', int((es < .2).sum()))
rng = np.random.default_rng(54); out = {}
def slope(k, lab):
    k = np.flatnonzero(k)
    if len(k) < 6: print(f'{lab:46s} n {len(k)} (too few)'); return
    X = np.c_[np.ones(len(k)), lm[k], lr[k]]; s = np.linalg.lstsq(X, e[k], rcond=None)[0][2]
    bs = []
    for _ in range(4000):
        i = rng.choice(k, len(k)); Xi = np.c_[np.ones(len(i)), lm[i], lr[i]]
        try: bs.append(np.linalg.lstsq(Xi, e[i], rcond=None)[0][2])
        except Exception: pass
    lo, hi = np.percentile(bs, [2.5, 97.5]); out[lab] = dict(n=len(k), slope=round(float(s), 3), ci95=[round(float(lo), 3), round(float(hi), 3)], rms=round(float(np.sqrt(np.mean(e[k] ** 2))), 3), mean=round(float(np.mean(e[k])), 3))
    print(f'{lab:46s} n {len(k):3d}  slope {s:+.2f}  [{lo:+.2f}, {hi:+.2f}]  rms {np.sqrt(np.mean(e[k]**2)):.3f}  mean {np.mean(e[k]):+.3f}')
slope(np.ones(len(e), bool), 'all 75')
slope(np.isfinite(peri), 'Milky Way satellites with orbits')
slope(np.isfinite(peri) & (peri < 50), '  pericenter < 50 kpc')
slope(np.isfinite(peri) & (peri >= 50), '  pericenter >= 50 kpc')
slope(~np.isfinite(peri), 'no orbit (M31 satellites, field, others)')
slope(es < .2, 'dispersion error < 20%')
slope((es < .2) & ~(np.isfinite(peri) & (peri < 50)), 'dispersion error < 20%, not pericenter < 50')
slope((es < .1), 'dispersion error < 10%')
slope((es < .1) & np.isfinite(peri) & (peri >= 50), 'F2 as worded (error < 10%, pericenter >= 50)')
k = np.isfinite(peri); print('corr(residual, log pericenter) among MW satellites:', round(float(np.corrcoef(np.log10(peri[k]), e[k])[0, 1]), 3), ' partial given logM*, logRe:',
      round(float(np.corrcoef(np.log10(peri[k]) - np.c_[np.ones(k.sum()), lm[k], lr[k]] @ np.linalg.lstsq(np.c_[np.ones(k.sum()), lm[k], lr[k]], np.log10(peri[k]), rcond=None)[0],
                              e[k] - np.c_[np.ones(k.sum()), lm[k], lr[k]] @ np.linalg.lstsq(np.c_[np.ones(k.sum()), lm[k], lr[k]], e[k], rcond=None)[0])[0, 1]), 3))
json.dump(out, open('res54_dwarfs_peri.json', 'w'), indent=1)
