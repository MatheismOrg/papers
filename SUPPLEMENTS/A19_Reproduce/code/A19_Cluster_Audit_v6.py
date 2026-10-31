"""A19 protocol version 6: audit of the cluster-scale test of protocol 5, implemented exactly as written in
PROTOCOL_v6.txt. Checks: (A) root status of each inversion, (B) the stellar-mass range the Moster relation
was fitted on, (C) a probabilistic inversion with the published 0.15 dex scatter and the Tinker et al. (2008)
mass function, (D) halo masses against NFW masses fitted to the CLASH accelerations, (E) h = 0.704.
Usage: python3 A19_Cluster_Audit_v6.py data/Clusters blind_test_dir out_prefix   (needs colossus)"""
import sys, os, re, csv, json, numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_Clusters_v5 as V5
from A19_Dwarfs_LVDB import G, ACC, FB, rar
from colossus.cosmology import cosmology
from colossus.lss import mass_function

H = 0.671; OM, OL = 0.3175, 0.6825; SEED = 20261031; NBOOT = 5000
GRID = np.arange(9.0, 16.5 + 1e-9, 0.005)
cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624})
cosmology.setCosmology('A19planck')
m_ = lambda y: np.log(1 + y) - y / (1 + y)

def g_nfw_h(r, M, z, h=H, dark=True):
    Hz = (100 * h * 1e-3) * np.sqrt(OM * (1 + z) ** 3 + OL); rho = 3 * Hz ** 2 / (8 * np.pi * G)
    r200 = (3 * M / (4 * np.pi * 200 * rho)) ** (1 / 3)
    c = 10 ** (0.520 + 0.385 * np.exp(-0.617 * z ** 1.21) + (-0.101 + 0.026 * z) * np.log10(M * h / 1e12))
    x = np.minimum(r / (r200 / c), c)
    return G * ((1 - FB) if dark else 1.0) * M * m_(x) / m_(c) / r ** 2 * ACC

def score(res, hosts, keep=None):
    if keep is not None: res, hosts = res[keep], hosts[keep]
    return V5.rms_h(res, hosts)

def boot(a, b, hosts):
    uh = np.unique(hosts); idx = {u: np.where(hosts == u)[0] for u in uh}; rng = np.random.default_rng(SEED); d = []
    for _ in range(NBOOT):
        pick = rng.choice(uh, len(uh), replace=True)
        d.append(np.sqrt(np.mean([np.mean(a[idx[u]] ** 2) for u in pick])) - np.sqrt(np.mean([np.mean(b[idx[u]] ** 2) for u in pick])))
    return [round(V5.rms_h(a, hosts) - V5.rms_h(b, hosts), 4)] + np.percentile(d, [2.5, 97.5]).round(4).tolist()

def main(dd, bt, out):
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    pts, sha = V5.load(dd); hosts = np.array([p['host'] for p in pts]); groups = np.array([p['group'] for p in pts])
    lo = np.log10([p['gobs'] for p in pts])
    R5, O5, halo5, capped5 = V5.predict(pts, rules)
    elm = {}
    for l in open(os.path.join(dd, 'manga_bcg_table1_extracted.txt')):
        if l.strip():
            ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
            elm['MaNGA_' + ident.replace('∗', '').replace('*', '')] = v[7]
    info = {}
    for p in pts:
        if p['host'] not in info: info[p['host']] = dict(Ms=p['Ms'], z=p['z'], group='MaNGA' if p['group'] == 'MaNGA_BCG' else 'CLASH')
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v6.sha256')).read().split()[0], 'input_sha256': sha}
    # (A) root status
    A = {}
    for u, d in info.items():
        Mh, cap = V5.moster_halo(d['Ms'], d['z'])
        gap = abs(np.log10(V5.moster_mstar(Mh, d['z'])) - np.log10(d['Ms']))
        A[u] = dict(log_Ms=round(float(np.log10(d['Ms'])), 3), log_Mh_deterministic=round(float(np.log10(Mh)), 3), root=bool(gap < 1e-3), gap_dex=round(float(gap), 4))
    res['A_root_status'] = dict(n_clipped=sum(1 for v in A.values() if not v['root']), clipped=[u for u, v in A.items() if not v['root']])
    # (B) supported stellar-mass range
    keep = np.array([info[h]['Ms'] <= 10 ** 11.7 for h in hosts])
    res['B_supported_range'] = dict(n_hosts=int(len(np.unique(hosts[keep]))), n_points=int(keep.sum()),
                                    n_hosts_MaNGA=int(len({h for h in hosts[keep] if h.startswith('MaNGA')})), n_hosts_CLASH=int(len({h for h in hosts[keep] if not h.startswith('MaNGA')})),
                                    H=round(score(R5['H'], hosts, keep), 4), R1=round(score(R5['R1'], hosts, keep), 4),
                                    bootstrap_H_minus_R1=boot(R5['H'][keep], R5['R1'][keep], hosts[keep]))
    # (C) probabilistic inversion
    pred = np.zeros(len(pts)); post_med = {}
    for u, d in info.items():
        s = np.sqrt(0.15 ** 2 + (elm[u] if d['group'] == 'MaNGA' else 0.10) ** 2)
        prior = mass_function.massFunction(10 ** GRID * H, d['z'], mdef='200c', model='tinker08', q_out='dndlnM')
        like = np.exp(-0.5 * ((np.log10(d['Ms']) - np.log10(V5.moster_mstar(10 ** GRID, d['z']))) / s) ** 2)
        w = prior * like; w = w / w.sum()
        cdf = np.cumsum(w); post_med[u] = float(GRID[np.searchsorted(cdf, 0.5)])
        Mg = np.maximum(10 ** GRID, d['Ms'] / FB)
        for i in np.where(hosts == u)[0]:
            p = pts[i]; pred[i] = np.sum(w * np.log10(p['gbar'] + g_nfw_h(p['r'], Mg, p['z'])))
    RC = lo - pred
    res['C_probabilistic'] = dict(score=round(score(RC, hosts), 4), R1=round(score(R5['R1'], hosts), 4), deterministic=round(score(R5['H'], hosts), 4),
                                  groups={g: round(V5.rms_h(RC[groups == g], hosts[groups == g]), 4) for g in ('MaNGA_BCG', 'CLASH_centre', 'CLASH_outer')},
                                  mean_residual={g: round(float(np.mean([np.mean(RC[(groups == g) & (hosts == u)]) for u in np.unique(hosts[groups == g])])), 4) for g in ('MaNGA_BCG', 'CLASH_centre', 'CLASH_outer')},
                                  bootstrap_C_minus_R1=boot(RC, R5['R1'], hosts),
                                  posterior_median_log_Mh={'CLASH_range': [round(min(v for k, v in post_med.items() if not k.startswith('MaNGA')), 2), round(max(v for k, v in post_med.items() if not k.startswith('MaNGA')), 2)],
                                                           'MaNGA_range': [round(min(v for k, v in post_med.items() if k.startswith('MaNGA')), 2), round(max(v for k, v in post_med.items() if k.startswith('MaNGA')), 2)],
                                                           'n_above_16': int(sum(1 for v in post_med.values() if v > 16))})
    b = res['C_probabilistic']['bootstrap_C_minus_R1']
    res['C_probabilistic']['robust_to_scatter'] = bool(b[0] < 0 and b[2] < 0)
    # (D) independent CLASH masses from the published accelerations at 100 to 600 kpc
    D = {}
    for u, d in info.items():
        if d['group'] != 'CLASH': continue
        m = (hosts == u) & (groups == 'CLASH_outer'); r = np.array([pts[i]['r'] for i in np.where(m)[0]]); y = lo[m]
        f = lambda q: np.log10(g_nfw_tot(r, 10 ** q[0], q[1], d['z'])) - y
        fit = least_squares(f, [15.0, 4.0], bounds=([12.0, 1.0], [17.0, 20.0]))
        D[u] = dict(log_M200_lensing=round(float(fit.x[0]), 3), c_lensing=round(float(fit.x[1]), 2), fit_rms_dex=round(float(np.sqrt(np.mean(fit.fun ** 2))), 4),
                    log_ratio_deterministic=round(float(np.log10(halo5[u]) - fit.x[0]), 3), log_ratio_probabilistic=round(float(post_med[u] - fit.x[0]), 3))
    rd = np.array([v['log_ratio_deterministic'] for v in D.values()]); rp = np.array([v['log_ratio_probabilistic'] for v in D.values()])
    res['D_independent_masses'] = dict(n=len(D), median_log_ratio_deterministic=round(float(np.median(rd)), 3), scatter_deterministic=round(float(np.std(rd)), 3),
                                       median_log_ratio_probabilistic=round(float(np.median(rp)), 3), scatter_probabilistic=round(float(np.std(rp)), 3),
                                       lensing_log_M200_range=[round(min(v['log_M200_lensing'] for v in D.values()), 2), round(max(v['log_M200_lensing'] for v in D.values()), 2)],
                                       consistent=bool(abs(np.median(rd)) <= 0.2 or abs(np.median(rp)) <= 0.2), per_cluster=D)
    # (E) h = 0.704
    RE = np.zeros(len(pts))
    for i, p in enumerate(pts):
        RE[i] = lo[i] - np.log10(p['gbar'] + g_nfw_h(p['r'], halo5[p['host']], p['z'], h=0.704))
    res['E_h_0704'] = dict(score=round(score(RE, hosts), 4), R1=round(score(R5['R1'], hosts), 4))
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '_Hosts.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['host', 'group', 'redshift', 'log_Ms', 'log_Mh_deterministic', 'root', 'log_Mh_posterior_median', 'log_M200_lensing'])
        for u, d in info.items():
            w.writerow([u, d['group'], d['z'], A[u]['log_Ms'], A[u]['log_Mh_deterministic'], int(A[u]['root']), round(post_med[u], 3), D[u]['log_M200_lensing'] if u in D else ''])
    print(json.dumps({k: v for k, v in res.items() if k not in ('input_sha256',)}, indent=1)[:6000])

def g_nfw_tot(r, M, c, z, h=H):
    Hz = (100 * h * 1e-3) * np.sqrt(OM * (1 + z) ** 3 + OL); rho = 3 * Hz ** 2 / (8 * np.pi * G)
    r200 = (3 * M / (4 * np.pi * 200 * rho)) ** (1 / 3); x = np.minimum(r / (r200 / c), c)
    return G * M * m_(x) / m_(c) / r ** 2 * ACC

if __name__ == '__main__': main(*sys.argv[1:4])
