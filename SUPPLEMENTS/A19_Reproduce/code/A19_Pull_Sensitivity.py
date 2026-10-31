"""A19 descriptive check for Section 5.3 (not a protocol test): how much the halo's pull at each measured
cluster-scale radius changes per dex of halo mass, d log10 g_halo / d log10 M_200, evaluated at protocol 5's
halo mass for each system (and, for the CLASH points, at the lensing NFW mass of protocol 6, check D).
Usage: python3 A19_Pull_Sensitivity.py data/Clusters blind_test out.json"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_Clusters_v5 as V5

def slope(r, M, z, d=0.01):
    lo = V5.g_nfw_z(r, M * 10 ** -d, z)[0]; hi = V5.g_nfw_z(r, M * 10 ** d, z)[0]
    return float((np.log10(hi) - np.log10(lo)) / (2 * d))

def main(dd, bt, out):
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    pts, sha = V5.load(dd); R5, O5, halo5, capped5 = V5.predict(pts, rules)
    lens = {k: 10 ** v['log_M200_lensing'] for k, v in json.load(open(os.path.join(bt, 'A19_Cluster_Audit_v6_Results.json')))['D_independent_masses']['per_cluster'].items()}
    res = {}
    for g in ('MaNGA_BCG', 'CLASH_centre', 'CLASH_outer'):
        s = [slope(p['r'], halo5[p['host']], p['z']) for p in pts if p['group'] == g]
        res[g] = dict(n=len(s), min=round(min(s), 3), median=round(float(np.median(s)), 3), max=round(max(s), 3),
                      r_kpc=[round(min(p['r'] for p in pts if p['group'] == g), 1), round(max(p['r'] for p in pts if p['group'] == g), 1)])
    s = [slope(p['r'], lens[p['host']], p['z']) for p in pts if p['group'] in ('CLASH_centre', 'CLASH_outer')]
    res['CLASH_at_lensing_mass'] = dict(n=len(s), min=round(min(s), 3), median=round(float(np.median(s)), 3), max=round(max(s), 3))
    allv = [v for k, v in res.items() if k != 'CLASH_at_lensing_mass']
    res['all_points'] = dict(min=min(v['min'] for v in allv), max=max(v['max'] for v in allv))
    json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
