"""SA73 T62 (compare46_boot.py with the final formula): paired host bootstrap of the held-out errors of the final formula
with its core (PROTOCOL_T61's run, res61_MOGF0IESBC_primary.json) against MOND's relation and standard Moster halos
(compare46.py definitions, nothing fitted, each given its better case: the spheroid IMF trend for early types and MaNGA,
and no subhalo factor for the standard halos' dwarfs). 4000 resamples of hosts. -> res62_boot.json"""
import json, numpy as np
H = __import__('harness_T15_hashed')
gs = H.disks()
def pops(etg_rows, cluster_rows):
    cl = cluster_rows()
    return {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(etg_rows()),
            'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']),
            'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
RAW = pops(H.etg_rows, H.cluster_rows)
import compare46 as C                                   # reruns compare46 (prints its table) and patches the rows
IMF = C.IMF
best = {'MOND': {'disk': RAW, 'dwarfs': RAW, 'ETG': IMF, 'MaNGA': IMF, 'CLASH_centre': RAW, 'CLASH_outer': RAW},
        'LCDM': {'disk': RAW, 'dwarfs': RAW, 'ETG': IMF, 'MaNGA': IMF, 'CLASH_centre': IMF, 'CLASH_outer': IMF}}
gf = {'MOND': C.g_mond, 'LCDM': C.g_lcdm}
recs = json.load(open('res61_MOGF0IESBC_primary.json'))['records']
rng = np.random.default_rng(20261107); out = {}
for p in ['disk', 'dwarfs', 'ETG', 'MaNGA', 'CLASH_centre', 'CLASH_outer']:
    pooled = {}
    for r in recs:
        for h, v in r['losses'][p].items(): pooled.setdefault(h, []).append(v)
    hosts = sorted(pooled); f = np.array([np.mean(pooled[h]) for h in hosts])
    row = {'formula': round(float(np.sqrt(f.mean())), 4)}
    for m in ('MOND', 'LCDM'):
        a = best[m][p][p]; er = np.log10(a['go'] / gf[m](a))
        o = np.array([np.mean(er[a['name'] == h] ** 2) for h in hosts])
        idx = rng.integers(0, len(hosts), (4000, len(hosts)))
        d = np.sqrt(f[idx].mean(1)) - np.sqrt(o[idx].mean(1))
        row[m] = round(float(np.sqrt(o.mean())), 4); row[f'formula_minus_{m}'] = round(float(np.sqrt(f.mean()) - np.sqrt(o.mean())), 4)
        row[f'ci95_{m}'] = [round(float(x), 4) for x in np.percentile(d, [2.5, 97.5])]
    out[p] = row; print(p, json.dumps(row), flush=True)
json.dump(out, open('res62_boot.json', 'w'), indent=1)
