"""A19 protocol version 7: the cluster halo masses of protocols 5 and 6 against the published lensing masses of
Umetsu et al. (2016, Table 2), implemented exactly as written in PROTOCOL_v7.txt.
Usage: python3 A19_Published_Masses_v7.py data/Clusters blind_test out_prefix"""
import sys, os, csv, json, numpy as np

H_PUB, H_USE = 0.7, 0.671

def stats(r, s):
    r, s = np.asarray(r), np.asarray(s)
    return dict(n=int(len(r)), median=round(float(np.median(r)), 3), sd=round(float(np.std(r)), 3),
                rms_expected=round(float(np.sqrt(np.mean(s ** 2))), 3), n_within_2s=int(np.sum(np.abs(r) <= 2 * s)))

def main(dd, bt, out):
    pub = {r['host']: r for r in csv.DictReader(open(os.path.join(dd, 'umetsu2016_table2.csv')))}
    hosts = {r['host']: r for r in csv.DictReader(open(os.path.join(bt, 'A19_Cluster_Audit_v6_Hosts.csv'))) if r['group'] == 'CLASH'}
    assert set(pub) == set(hosts) and len(pub) == 20
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v7.sha256')).read().split()[0]}
    rows = []
    for conv, scale in (('as_published', 1.0), ('h_0671', H_PUB / H_USE)):
        R = {'mean': [], 'probabilistic': []}; S = []; route = []; D = []
        for u in sorted(pub):
            p = pub[u]; M = float(p['M200c_1e14_Msun_h70inv']) * 1e14 * scale; e = float(p['e_M200c']) * 1e14 * scale
            lm = np.log10(M); s = e / (M * np.log(10)); h = hosts[u]
            R['mean'].append(float(h['log_Mh_deterministic']) - lm); R['probabilistic'].append(float(h['log_Mh_posterior_median']) - lm)
            S.append(s); route.append(p['route']); D.append(float(h['log_M200_lensing']) - lm)
            if conv == 'as_published':
                rows.append([u, p['route'], round(lm, 3), round(s, 3), h['log_Mh_deterministic'], h['log_Mh_posterior_median'], h['log_M200_lensing'],
                             round(R['mean'][-1], 3), round(R['probabilistic'][-1], 3), round(D[-1], 3)])
        S = np.array(S); route = np.array(route); block = {}
        for k in ('mean', 'probabilistic'):
            r = np.array(R[k]); block[k] = {'all': stats(r, S)}
            for g in ('X-ray', 'high-magnification'): block[k][g] = stats(r[route == g], S[route == g])
            block[k]['consistent'] = bool(abs(block[k]['all']['median']) <= 0.2 and block[k]['all']['sd'] <= 2 * block[k]['all']['rms_expected'])
        D = np.array(D)
        block['check_D_vs_published'] = dict(median=round(float(np.median(D)), 3), max_abs=round(float(np.max(np.abs(D))), 3), faithful=bool(np.max(np.abs(D)) < 0.1))
        res[conv] = block
    res['consistent_with_published'] = bool(res['as_published']['mean']['consistent'] or res['as_published']['probabilistic']['consistent'])
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '_Clusters.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['host', 'route', 'log_M200c_published', 'sigma_log_published', 'log_Mh_mean_inversion', 'log_Mh_posterior_median', 'log_M200_checkD',
                    'log_ratio_mean', 'log_ratio_probabilistic', 'log_ratio_checkD'])
        w.writerows(rows)
    print(json.dumps(res, indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
