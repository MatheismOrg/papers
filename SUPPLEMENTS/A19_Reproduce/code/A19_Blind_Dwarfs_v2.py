"""A19 blind dwarf test, implementing PROTOCOL.txt exactly.
Usage: python3 A19_Blind_Dwarfs.py dwarf_all.csv frozen_disk_rules.json out_prefix [--select-only]"""
import sys, os, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import g_nfw, rar, moster_mref, FB, G, ACC
KEYS = set(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'original44_keys.txt')).read().split())
V1_BLIND = {'aquarius_2','bootes_2','carina_2','carina_3','centaurus_1','columba_1','eridanus_4','grus_1','leo_6','pegasus_3','pegasus_4','pictor_2','pisces_2','tucana_2','tucana_4','tucana_5','andromeda_17','andromeda_19','andromeda_20','andromeda_26','andromeda_27','andromeda_28','andromeda_29','cassiopeia_2','cassiopeia_3','lacerta_1','perseus_1'}
TIDAL = {'sagittarius_1': 'Sagittarius I', 'tucana_3': 'Tucana III', 'bootes_3': 'Bootes III'}
HOSTS = {'mw', 'lmc', 'smc', 'm_031'}
def select(path):
    keep, excl = [], []
    for x in csv.DictReader(open(path)):
        k = x['key']
        def out(why): excl.append((k, x['name'], why))
        if k in KEYS: continue                              # step 1 (not listed: original sample)
        if x['confirmed_real'] != '1' or x['confirmed_galaxy'] != '1': continue   # step 2 (not galaxies)
        if x['host'] not in HOSTS: continue                 # step 3 (outside the Local Group satellite systems)
        if not x['M_V']: out('no M_V'); continue
        if float(x['M_V']) < -14: out('M_V brighter than -14'); continue
        if not x['vlos_sigma'] or x['vlos_sigma_ul']: out('no measured dispersion (absent or upper limit)'); continue
        if not x['rhalf_physical'] or not x['rhalf_sph_physical']: out('no half-light radius'); continue
        LV = 10 ** (-0.4 * (float(x['M_V']) - 4.83))
        if x['mass_HI'] and not x['mass_HI_ul'] and 10 ** float(x['mass_HI']) > 0.1 * 2 * LV: out('HI above a tenth of the stellar mass'); continue
        if k in TIDAL: out('named tidally disrupting system'); continue
        keep.append(x)
    return keep, excl
def main(path, rules, out, select_only=False):
    keep, excl = select(path)
    print('selected', len(keep)); print(', '.join(x['name'] for x in keep))
    for e in excl: print('excluded', e)
    if select_only: return
    D = json.load(open(rules))['frozen_on_all_disks']; S, B, A0 = D['stellar_ABC'], D['baryonic_ABC'], D['rar_a0']
    LV = 10 ** (-0.4 * (np.array([float(x['M_V']) for x in keep]) - 4.83)); s = np.array([float(x['vlos_sigma']) for x in keep])
    Rs = np.array([float(x['rhalf_sph_physical']) for x in keep]) / 1000; Rm = np.array([float(x['rhalf_physical']) for x in keep]) / 1000
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    def run(Re, ml=2.0):
        Ms = ml * LV; r = 4 / 3 * Re; gb = G * (Ms / 2) / r ** 2 * ACC; lo = np.log10(3 * s ** 2 / r * ACC)
        H = lambda p: gb + g_nfw(r, Ms / FB + 10 ** p[0] * (Ms / 1e10) ** p[1] * (Re / 3) ** p[2])
        mo = np.array([g + g_nfw(rr, max(moster_mref(m), m / FB)) for g, rr, m in zip(gb, r, Ms)])
        return lo, {'stellar_rule': lo - np.log10(H(S)), 'baryonic_rule': lo - np.log10(H(B)), 'rar_disk_a0': lo - np.log10(rar(gb, A0)),
                    'rar_published_a0': lo - np.log10(rar(gb, 1.2e-10)), 'moster': lo - np.log10(mo)}
    lo, R = run(Rs); n = len(keep); st = R['stellar_rule']
    def boot(a, b):
        rng = np.random.default_rng(20261030); d = []
        for _ in range(20000):
            i = rng.integers(0, n, n); d.append(rms(a[i]) - rms(b[i]))
        return [round(rms(a) - rms(b), 4)] + np.percentile(d, [2.5, 97.5]).round(4).tolist()
    res = {'protocol_sha256': open(os.path.join(os.path.dirname(os.path.abspath(rules)), 'PROTOCOL_v2.sha256')).read().split()[0],
           'n': n, 'selected': [x['name'] for x in keep], 'excluded': excl,
           'rmse_dex': {k: round(rms(v), 4) for k, v in R.items()},
           'bootstrap_rmse_difference': {'stellar_minus_rar_disk': boot(st, R['rar_disk_a0']), 'stellar_minus_baryonic': boot(st, R['baryonic_rule']),
                                          'stellar_minus_moster': boot(st, R['moster']), 'stellar_minus_rar_published': boot(st, R['rar_published_a0'])},
           'closer_under_stellar': {k: int(np.sum(np.abs(st) < np.abs(v))) for k, v in R.items() if k != 'stellar_rule'}}
    p = res['bootstrap_rmse_difference']
    res['primary_success'] = bool(res['rmse_dex']['stellar_rule'] < res['rmse_dex']['rar_disk_a0'] and p['stellar_minus_rar_disk'][2] < 0)
    res['secondary_baryonic_success'] = bool(p['stellar_minus_baryonic'][2] < 0)
    res['sensitivity'] = {lab: {k: round(rms(v), 4) for k, v in run(*a)[1].items()} for lab, a in
                          [('major_axis_radii', (Rm,)), ('ML_1', (Rs, 1.0)), ('ML_3', (Rs, 3.0))]}
    v1 = np.array([x['key'] in V1_BLIND for x in keep]); res['blind_v1_subset'] = dict(n=int(v1.sum()), keys=[x['key'] for x in keep if x['key'] in V1_BLIND], **{k: round(rms(v[v1]), 4) for k, v in R.items()})
    br = LV >= 1e5; res['subsets'] = {lab: dict(n=int(m.sum()), **{k: round(rms(v[m]), 4) for k, v in R.items()}) for lab, m in [('LV_ge_1e5', br), ('LV_lt_1e5', ~br)] if m.sum()}
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['lvdb_key', 'name', 'host', 'L_V', 'R_e_circ_kpc', 'sigma_kms', 'log_g_obs'] + ['res_' + k for k in R] + ['ref_structure', 'ref_vlos', 'ref_m_v'])
        for i, x in enumerate(keep): w.writerow([x['key'], x['name'], x['host'], f'{LV[i]:.4g}', f'{Rs[i]:.4f}', x['vlos_sigma'], f'{lo[i]:.4f}'] + [f'{R[k][i]:.4f}' for k in R] + [x['ref_structure'], x['ref_vlos'], x['ref_m_v']])
    print(json.dumps({k: res[k] for k in ['n', 'rmse_dex', 'bootstrap_rmse_difference', 'closer_under_stellar', 'primary_success', 'secondary_baryonic_success', 'sensitivity', 'subsets', 'blind_v1_subset']}, indent=1))
if __name__ == '__main__': main(sys.argv[1], sys.argv[2], sys.argv[3], '--select-only' in sys.argv)
