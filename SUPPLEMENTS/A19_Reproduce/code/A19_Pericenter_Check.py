"""A19 descriptive check for Section 3.3 (exploratory, not a protocol): do the stellar rule's dwarf residuals track
orbits (pericenters, Pace, Erkal and Li 2022), present Galactocentric distance (LVDB distance_gc) or luminosity? Rank (Spearman) correlations, and rank partial
correlations from residualized ranks. Residuals are the frozen stellar rule's, Local Volume Database values: the
original 44 (results/A19_Dwarfs_LVDB.csv) and the 25 out of the original sample (blind_test/A19_OutOfSample_v2.csv).
Usage: python3 A19_Pericenter_Check.py . out.json"""
import sys, os, csv, json, numpy as np
from scipy.stats import spearmanr, pearsonr, rankdata

def partial(a, b, c):
    ra, rb, rc = rankdata(a), rankdata(b), rankdata(c)
    ea = ra - np.polyval(np.polyfit(rc, ra, 1), rc); eb = rb - np.polyval(np.polyfit(rc, rb, 1), rc)
    r, p = pearsonr(ea, eb); return round(float(r), 3), round(float(p), 3)

def main(root, out):
    peri = {r['lvdb_key']: float(r['rperi_kpc']) for r in csv.DictReader(open(os.path.join(root, 'data/Orbits/pace2022_table3_pericenters.csv')))}
    res = {}
    for x in csv.DictReader(open(os.path.join(root, 'results/A19_Dwarfs_LVDB.csv'))): res[x['lvdb_key']] = ('original_44', float(x['res_stellar_rule']), float(x['L_V_Lsun']))
    for x in csv.DictReader(open(os.path.join(root, 'blind_test/A19_OutOfSample_v2.csv'))):
        if x['lvdb_key'] not in res: res[x['lvdb_key']] = ('out_of_sample_25', float(x['res_stellar_rule']), float(x['L_V']))
    db = {r['key']: r for r in csv.DictReader(open(os.path.join(root, 'data/LVDB/dwarf_all.csv')))}
    keys = sorted(k for k in res if k in peri)
    ld = np.log10([float(db[k]['distance_gc']) for k in keys])   # present Galactocentric distance, kpc
    y = np.array([res[k][1] for k in keys]); lL = np.log10([res[k][2] for k in keys]); lp = np.log10([peri[k] for k in keys])
    sp = lambda a, b: [round(float(v), 3) for v in spearmanr(a, b)]
    o = dict(n=len(keys), by_set={s: sum(1 for k in keys if res[k][0] == s) for s in ('original_44', 'out_of_sample_25')},
             residual_vs_logL=sp(lL, y), residual_vs_log_rperi=sp(lp, y), log_rperi_vs_logL=sp(lp, lL),
             partial_residual_rperi_given_L=partial(lp, y, lL), partial_residual_L_given_rperi=partial(lL, y, lp),
             residual_vs_log_dgc=sp(ld, y), log_dgc_vs_logL=sp(ld, lL), partial_residual_dgc_given_L=partial(ld, y, lL), dwarfs=keys)
    json.dump(o, open(out, 'w'), indent=1); print(json.dumps({k: v for k, v in o.items() if k != 'dwarfs'}, indent=1))

if __name__ == '__main__': main(*sys.argv[1:3])
