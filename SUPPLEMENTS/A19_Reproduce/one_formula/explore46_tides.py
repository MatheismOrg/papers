"""SA73 T46 exploration (all data, no protocol): tidal truncation of Milky Way satellite halos at their pericenters
(Pace, Erkal and Li 2022, Table 3, with the LMC) in the McMillan (2017) Milky Way (NFW rho0 = 0.00854 Msun/pc^3,
r_s = 19.6 kpc, plus 6.6e10 Msun of baryons as a point mass). Jacobi radius r_t = r_p [m(<r_t) / ((2 - dlnM/dlnr) M_MW(<r_p))]^(1/3)
(Tormen 1998 form); the subhalo's mass stops growing beyond r_t. Formula T45 otherwise. Prints the dwarf score."""
import csv, json, numpy as np
import harness_T45 as M
H = M.H; S = M.S; K = H.K; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.339])
AB = {'Ant': 'antlia', 'Aqu': 'aquarius', 'Boo': 'bootes', 'CVn': 'canes_venatici', 'Car': 'carina', 'Col': 'columba', 'CB': 'coma_berenices',
      'Cra': 'crater', 'Dra': 'draco', 'Eri': 'eridanus', 'For': 'fornax', 'Gru': 'grus', 'Her': 'hercules', 'Hor': 'horologium', 'Hyd': 'hydra',
      'Hyi': 'hydrus', 'Leo': 'leo', 'Peg': 'pegasus', 'Phx': 'phoenix', 'Pis': 'pisces', 'Ret': 'reticulum', 'Sgr': 'sagittarius',
      'Scl': 'sculptor', 'Seg': 'segue', 'Sxt': 'sextans', 'Tri': 'triangulum', 'Tuc': 'tucana', 'UMa': 'ursa_major', 'UMi': 'ursa_minor', 'Wil': 'willman'}
ROM = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}
def key(n):
    p = n.split(); base = AB[p[0]]; num = 1 if len(p) == 1 else (ROM[p[1]] if p[1] in ROM else int(p[1])); return f'{base}_{num}'
PERI = {key(r['name']): float(r['rperi_kpc']) for r in csv.DictReader(open(H.S / 'ext/pace2022_table3.csv', encoding='utf-8-sig'))}
RHO0, RS, MB = 0.00854e9, 19.6, 6.6e10
def m_mw(r): x = r / RS; return 4 * np.pi * RHO0 * RS ** 3 * (np.log1p(x) - x / (1 + x)) + MB
def slope(r): return (np.log(m_mw(r * 1.01)) - np.log(m_mw(r / 1.01))) / (2 * np.log(1.01))
rows = H.dwarf_rows(); a = H.arr(rows)
m = M.P.G.mass(th, a['ms'], a['mp'], a['re'], a['z'], S.thresholds(a)); c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
rv = (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
def menc(i, r): return m[i] * H.ff(r * c[i] / rv[i]) / H.ff(c[i])
rt = np.full(len(m), np.inf); names = a['name']
for i, n in enumerate(names):
    if n not in PERI: continue
    rp = PERI[n]; Mhost = m_mw(rp) * (2 - slope(rp)); r = rv[i]
    for _ in range(200):
        r_new = rp * (menc(i, r) / Mhost) ** (1 / 3)
        if abs(r_new - r) < 1e-6 * r: break
        r = .5 * (r + r_new)
    rt[i] = r
reff = np.minimum(a['r'], rt)
mt = m * H.ff(reff * c / rv) / H.ff(c)
g_t = a['gb'] + H.G * (1 - H.FB) * mt / a['r'] ** 2 * H.ACC
g_0 = H.gpred(M.MODEL, th, a)
e0, e1 = np.log10(a['go'] / g_0), np.log10(a['go'] / g_t)
def rms(e): return float(np.sqrt(np.mean([np.mean(e[names == n] ** 2) for n in np.unique(names)])))
print('matched', sum(n in PERI for n in np.unique(names)), 'of', len(np.unique(names)), 'unmatched Pace:', sorted(set(PERI) - set(names))[:30])
for i in np.argsort(rt / a['r']):
    if rt[i] < 3 * a['r'][i]: print(f'{names[i]:18s} peri {PERI[names[i]]:6.1f}  r_t {rt[i]:6.2f} kpc  r {a["r"][i]:5.2f}  res {e0[i]:+.2f} -> {e1[i]:+.2f}')
print(json.dumps(dict(dwarfs_T45=round(rms(e0), 4), dwarfs_tides=round(rms(e1), 4), mean0=round(float(np.mean(e0)), 4), mean1=round(float(np.mean(e1)), 4))))
