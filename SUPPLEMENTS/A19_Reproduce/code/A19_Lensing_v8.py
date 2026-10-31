"""A19 protocol version 8: the halo rules against the galaxy-galaxy lensing halo masses of SDSS central galaxies
(Mandelbaum et al. 2016), implemented exactly as written in PROTOCOL_v8.txt.
Usage: python3 A19_Lensing_v8.py data blind_test out_prefix   (needs colossus)"""
import sys, os, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_Disk_Fit_SPARC as DF, A19_Clusters_v5 as V5
from colossus.cosmology import cosmology
from colossus.halo import mass_defs

H, HL = 0.671, 0.673
cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624})
cosmology.setCosmology('A19planck')
LG = np.arange(5.0, 17.5, 0.0005)

def to_m200m(m200c, c, z):
    """M_200c in suns with concentration c -> log10 M_200m in h^-1 suns for h = 0.673."""
    m, r, cn = mass_defs.changeMassDefinition(m200c * H, c, z, '200c', '200m', profile='nfw')
    return float(np.log10(m / H * HL))

def c_rule(m): return 10 ** 0.905 * (m * H / 1e12) ** -0.101
def c_dm14(m, z): return 10 ** (0.520 + 0.385 * np.exp(-0.617 * z ** 1.21) + (-0.101 + 0.026 * z) * np.log10(m * H / 1e12))

def bend(ms, z):
    a = z / (1 + z); M1 = 10 ** (11.590 + 1.195 * a); N = 0.0351 - 0.0247 * a; be = 1.376 - 0.826 * a; ga = 0.608 + 0.329 * a
    M = 10 ** LG; lms = np.log10(2 * N * M / ((M / M1) ** -be + (M / M1) ** ga)); k = np.argmax(lms)
    lm_mo = np.interp(np.log10(ms), lms[:k + 1], LG[:k + 1]); lm_pl = np.log10(M1) + (np.log10(ms) - np.log10(2 * N * M1)) / (1 + be)
    return float(10 ** (lm_mo - lm_pl))

def main(dd, bt, out):
    A, B, C = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']['stellar_ABC']
    gal = DF.load(DF.table1(os.path.join(dd, 'SPARC', 'SPARC_Lelli2016c.mrt')), os.path.join(dd, 'SPARC', 'Rotmod_LTG'))
    lms = np.log10([g['Ms'] for g in gal]); lre = np.log10([g['Re'] for g in gal]); mhi = np.array([(g['Mp'] - g['Ms']) / 1.33 for g in gal])
    p1, p0 = np.polyfit(lms, lre, 1); k = mhi > 0; q1, q0 = np.polyfit(lms[k], np.log10(mhi[k]), 1)
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v8.sha256')).read().split()[0],
           'disk_scalings': dict(n_disks=len(gal), logRe=[round(float(p0), 4), round(float(p1), 4)], n_HI=int(k.sum()), logMHI=[round(float(q0), 4), round(float(q1), 4)])}
    rows = list(csv.DictReader(open(os.path.join(dd, 'Lensing', 'mandelbaum2016_central_halos.csv'))))
    def rule(lm, z, gas=True, size_shift=0.0, size_factor=1.0):
        ms = 10 ** lm; re = 10 ** (p0 + p1 * lm + size_shift) * size_factor; mp = ms + (1.33 * 10 ** (q0 + q1 * lm) if gas else 0.0)
        m = mp / DF.FB + 10 ** A * (ms / 1e10) ** B * (re / 3.0) ** C
        return to_m200m(m, c_rule(m), z)
    def moster(lm, z):
        m, cap = V5.moster_halo(10 ** lm, z); return to_m200m(m, c_dm14(m, z), z)
    def run(mcol):
        T = []
        for r in rows:
            lm, z, obs = float(r[mcol]), float(r['z_eff']), float(r['logM200m_hinv'])
            ep, em = float(r['err_plus']), float(r['err_minus'])
            d = dict(color=r['color'], logMs=lm, z=z, obs=obs, err_plus=ep, err_minus=em, moster=round(moster(lm, z), 3))
            if r['color'] == 'blue':
                d.update(rule=round(rule(lm, z), 3), rule_nogas=round(rule(lm, z, gas=False), 3),
                         rule_size_plus=round(rule(lm, z, size_shift=0.15), 3), rule_size_minus=round(rule(lm, z, size_shift=-0.15), 3))
            else:
                d.update(rule_disk_size=round(rule(lm, z, gas=False), 3), rule_half_size=round(rule(lm, z, gas=False, size_factor=0.5), 3))
            T.append(d)
        return T
    out_v = run('logMs_eff_VAGC'); out_m = run('logMs_eff_MPA')
    def crit_A(T, key='rule'):
        chk = []
        for d in T:
            if d['color'] != 'blue' or d['logMs'] >= 11.4: continue
            dev = d[key] - d['obs']; tol = max(0.2, 2 * (d['err_plus'] if dev > 0 else d['err_minus']))
            chk.append(dict(logMs=d['logMs'], dev=round(dev, 3), tol=round(tol, 3), ok=bool(abs(dev) <= tol)))
        return dict(bins=chk, consistent=bool(all(c['ok'] for c in chk)), median_dev=round(float(np.median([c['dev'] for c in chk])), 3))
    res['A_blue_VAGC'] = crit_A(out_v); res['A_blue_MPA'] = crit_A(out_m)
    res['A_sensitivity_VAGC'] = {k: crit_A(out_v, k) for k in ('rule_nogas', 'rule_size_plus', 'rule_size_minus')}
    # (B) the turn 6 type-bend prediction against the measured red-to-blue ratio
    red = {r['bin']: r for r in rows if r['color'] == 'red'}; blue = {r['bin']: r for r in rows if r['color'] == 'blue'}
    Bt = []
    for b_ in red:
        R, Bl = red[b_], blue[b_]; lm, z = float(R['logMs_eff_VAGC']), float(R['z_eff'])
        ratio = float(R['logM200m_hinv']) - float(Bl['logM200m_hinv'])
        sig = float(np.hypot((float(R['err_plus']) + float(R['err_minus'])) / 2, (float(Bl['err_plus']) + float(Bl['err_minus'])) / 2))
        lb = float(np.log10(bend(10 ** lm, z)))
        Bt.append(dict(bin=b_, logMs=lm, measured_log_ratio=round(ratio, 3), sigma=round(sig, 3), log_b=round(lb, 3), excess_sigma=round((lb - ratio) / sig, 2)))
    test = [t for t in Bt if 10.7 <= t['logMs'] <= 11.2]
    res['B_type_bend'] = dict(bins=Bt, refuted=bool(any(t['excess_sigma'] > 3 for t in test)), tested_bins=[t['bin'] for t in test])
    res['C_red_diagnostic'] = [dict(logMs=d['logMs'], obs=d['obs'], rule_disk_size=d['rule_disk_size'], rule_half_size=d['rule_half_size']) for d in out_v if d['color'] == 'red']
    res['D_moster'] = [dict(color=d['color'], logMs=d['logMs'], obs=d['obs'], moster=d['moster']) for d in out_v]
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '_Bins.csv', 'w', newline='') as fh:
        keys = ['color', 'logMs', 'z', 'obs', 'err_plus', 'err_minus', 'moster', 'rule', 'rule_nogas', 'rule_size_plus', 'rule_size_minus', 'rule_disk_size', 'rule_half_size']
        w = csv.writer(fh); w.writerow(['stellar_masses'] + keys)
        for lab, T in (('VAGC', out_v), ('MPA', out_m)):
            for d in T: w.writerow([lab] + [d.get(k, '') for k in keys])
    print(json.dumps(res, indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
