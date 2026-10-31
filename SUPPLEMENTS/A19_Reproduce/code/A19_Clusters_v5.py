"""A19 protocol version 5, implemented exactly as written in PROTOCOL_v5.txt.
The published Moster relation, unchanged and unfitted, carried to brightest cluster galaxies and clusters,
against MOND's radial acceleration relation and the disk-fitted stellar rule.
Usage: python3 A19_Clusters_v5.py data/Clusters blind_test_dir out_prefix
The halo mass is the value in the protocol's search range [10^9, 10^16.5] suns that brings the Moster
relation closest to the BCG's stellar mass (a bounded search, as A19_Dwarfs_LVDB.moster_mref does); hosts
whose stellar mass lies beyond the relation's reach at 10^16.5 are counted and reported."""
import sys, os, csv, json, re, hashlib, numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import G, ACC, FB, h, g_nfw, rar

OM, OL = 0.3175, 0.6825; SEED = 20261031; NBOOT = 5000; LO, HI = 9.0, 16.5
rms_h = lambda res, hosts: float(np.sqrt(np.mean([np.mean(np.square(res[hosts == u])) for u in np.unique(hosts)])))

def moster_mstar(Mh, z, evolve=True):
    a = z / (1 + z) if evolve else 0.0
    M1 = 10 ** (11.590 + 1.195 * a); N = 0.0351 - 0.0247 * a; be = 1.376 - 0.826 * a; ga = 0.608 + 0.329 * a
    x = Mh / M1
    return 2 * N * Mh / (x ** (-be) + x ** ga)

def moster_halo(Ms, z, evolve=True):
    f = lambda lm: (np.log10(moster_mstar(10 ** lm, z, evolve)) - np.log10(Ms)) ** 2
    o = minimize_scalar(f, bounds=(LO, HI), method='bounded', options={'xatol': 1e-8})
    return 10 ** o.x, bool(o.x > HI - 1e-3 and f(o.x) > 1e-6)

def g_nfw_z(r, M, z, evolve=True):
    if not evolve: return g_nfw(r, M)
    Hz = (100 * h * 1e-3) * np.sqrt(OM * (1 + z) ** 3 + OL)            # km/s/kpc
    rho = 3 * Hz ** 2 / (8 * np.pi * G); r200 = (3 * M / (4 * np.pi * 200 * rho)) ** (1 / 3)
    a = 0.520 + 0.385 * np.exp(-0.617 * z ** 1.21); b = -0.101 + 0.026 * z
    c = 10 ** (a + b * np.log10(M * h / 1e12)); x = np.minimum(r / (r200 / c), c)
    m = lambda y: np.log(1 + y) - y / (1 + y)
    return G * (1 - FB) * M * m(x) / m(c) / r ** 2 * ACC, c

def load(dd):
    sha = {f: hashlib.sha256(open(os.path.join(dd, f), 'rb').read()).hexdigest() for f in sorted(os.listdir(dd))}
    pts = []
    tab = {}
    for l in open(os.path.join(dd, 'manga_bcg_table1_extracted.txt')):
        if not l.strip(): continue
        ident, rest = l.strip().split(maxsplit=1)
        v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]; assert len(v) == 14
        tab['MaNGA_' + ident.replace('∗', '').replace('*', '')] = dict(z=v[0], Re=v[2], logMs=v[6])
    for r in csv.DictReader(open(os.path.join(dd, 'manga_bcg_standardized.csv'))):
        t = tab[r['system']]
        pts.append(dict(group='MaNGA_BCG', host=r['system'], r=float(r['radius_kpc']), gbar=float(r['gbar_m_s2']), gobs=float(r['gobs_m_s2']),
                        z=t['z'], Ms=10 ** t['logMs'], Re=t['Re']))
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(dd, 'clash_table1_idia.tsv')) if re.match(r'^\s*\d+\t', l)]
    ct = {x[2].strip(): dict(z=float(x[3]), Re=float(x[8]), Ms=float(x[11]) * 1e11) for x in rows}
    assert len(ct) == 20
    for r in csv.DictReader(open(os.path.join(dd, 'clash_standardized.csv'))):
        t = ct[r['canonical_cluster_name']]
        pts.append(dict(group='CLASH_centre' if r['population'] == 'CLASH_BCG' else 'CLASH_outer', host=r['canonical_cluster_name'],
                        r=float(r['radius_kpc']), gbar=float(r['gbar_m_s2']), gobs=float(r['gobs_m_s2']), z=t['z'], Ms=t['Ms'], Re=t['Re']))
    return pts, sha

def predict(pts, rules, clash_ms=1.0, manga_dex=0.0, evolve=True):
    S = rules['stellar_ABC']; A0 = rules['rar_a0']
    out = {k: np.zeros(len(pts)) for k in ('H', 'R1', 'R2', 'S')}; halo = {}; capped = set()
    for i, p in enumerate(pts):
        Ms = p['Ms'] * (clash_ms if p['group'] != 'MaNGA_BCG' else 10 ** manga_dex)
        key = p['host']
        if key not in halo:
            Mh, cap = moster_halo(Ms, p['z'], evolve); Mh = max(Mh, Ms / FB)
            if cap: capped.add(key)
            halo[key] = Mh
        Mh = halo[key]
        if evolve: gh, c = g_nfw_z(p['r'], Mh, p['z'])
        else: gh = g_nfw(p['r'], Mh)
        out['H'][i] = p['gbar'] + gh
        out['R1'][i] = rar(p['gbar'], A0); out['R2'][i] = rar(p['gbar'], 1.20e-10)
        out['S'][i] = p['gbar'] + g_nfw(p['r'], Ms / FB + 10 ** S[0] * (Ms / 1e10) ** S[1] * (p['Re'] / 3) ** S[2])
    lo = np.log10([p['gobs'] for p in pts])
    return {k: lo - np.log10(v) for k, v in out.items()}, out, halo, capped

def summarize(R, out, pts, hosts, groups):
    res = {}
    for k, v in R.items():
        d = {'all_70_hosts': round(rms_h(v, hosts), 4)}
        for gname in ('MaNGA_BCG', 'CLASH_centre', 'CLASH_outer'):
            m = groups == gname; hv = hosts[m]; vv = v[m]
            frac = np.abs(out[k][m] / np.array([p['gobs'] for p, mm in zip(pts, m) if mm]) - 1)
            d[gname] = dict(score=round(rms_h(vv, hv), 4), mean_residual=round(float(np.mean([np.mean(vv[hv == u]) for u in np.unique(hv)])), 4),
                            median_abs_frac_error_pct=round(float(np.median(frac)) * 100, 1), n_points=int(m.sum()), n_hosts=int(len(np.unique(hv))))
        res[k] = d
    return res

def main(dd, bt, out):
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    pts, sha = load(dd); hosts = np.array([p['host'] for p in pts]); groups = np.array([p['group'] for p in pts])
    assert len(np.unique(hosts)) == 70 and len(pts) == 134
    R, O, halo, capped = predict(pts, rules)
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v5.sha256')).read().split()[0], 'input_sha256': sha,
           'n_points': len(pts), 'n_hosts': 70, 'hosts_at_search_bound': sorted(capped), 'scores': summarize(R, O, pts, hosts, groups)}
    uh = np.unique(hosts); rng = np.random.default_rng(SEED); diffs = {k: [] for k in ('R1', 'R2', 'S')}
    idx = {u: np.where(hosts == u)[0] for u in uh}
    for _ in range(NBOOT):
        pick = rng.choice(uh, len(uh), replace=True)
        ms = {k: [np.mean(np.square(R[k][idx[u]])) for u in pick] for k in R}
        sH = np.sqrt(np.mean(ms['H']))
        for k in diffs: diffs[k].append(sH - np.sqrt(np.mean(ms[k])))
    res['bootstrap_95_H_minus'] = {k: [round(res['scores']['H']['all_70_hosts'] - res['scores'][k]['all_70_hosts'], 4)] + np.percentile(v, [2.5, 97.5]).round(4).tolist() for k, v in diffs.items()}
    b = res['bootstrap_95_H_minus']['R1']
    res['primary_success'] = bool(b[0] < 0 and b[2] < 0)
    res['halo_masses_log10'] = {u: round(float(np.log10(halo[u])), 3) for u in uh}
    sens = {}
    for lab, kw in (('clash_stellar_mass_doubled', dict(clash_ms=2.0)), ('manga_stellar_mass_minus_0.3dex', dict(manga_dex=-0.3)), ('no_redshift_evolution', dict(evolve=False))):
        Rs, Os, hs, cs = predict(pts, rules, **kw)
        sens[lab] = dict(scores={k: v['all_70_hosts'] for k, v in summarize(Rs, Os, pts, hosts, groups).items()},
                         groups_H={g: summarize(Rs, Os, pts, hosts, groups)['H'][g]['score'] for g in ('MaNGA_BCG', 'CLASH_centre', 'CLASH_outer')},
                         hosts_at_search_bound=len(cs))
    res['sensitivity'] = sens
    res['clash_residuals'] = {u: {k: [round(float(x), 3) for x in R[k][idx[u]]] for k in ('H', 'R1')} for u in uh if not u.startswith('MaNGA')}
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '_Points.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['group', 'host', 'redshift', 'radius_kpc', 'stellar_mass_Msun', 'log10_halo_mass_H', 'log_gobs', 'res_H', 'res_R1', 'res_R2', 'res_S'])
        for i, p in enumerate(pts):
            w.writerow([p['group'], p['host'], p['z'], p['r'], f"{p['Ms']:.4g}", f"{np.log10(halo[p['host']]):.3f}", f"{np.log10(p['gobs']):.4f}"] + [f'{R[k][i]:.4f}' for k in ('H', 'R1', 'R2', 'S')])
    print(json.dumps({k: v for k, v in res.items() if k not in ('clash_residuals', 'halo_masses_log10', 'input_sha256')}, indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
