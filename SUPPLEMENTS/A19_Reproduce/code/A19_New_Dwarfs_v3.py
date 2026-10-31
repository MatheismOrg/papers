"""A19 protocol version 3, implemented exactly as written in PROTOCOL_v3.txt.
Scores the frozen rules on dwarfs no rule has scored (measured dispersions outside the earlier sets, and
dispersion upper limits) and writes forecasts for dwarfs whose dispersion has not been measured.
Usage: python3 A19_New_Dwarfs_v3.py dwarf_all.csv blind_test_dir out_prefix
blind_test_dir holds frozen_disk_rules.json, original44_keys.txt, A19_OutOfSample_v2.csv and
PROTOCOL_v3.sha256. Nothing is fitted; every constant is frozen or published."""
import sys, os, csv, json, hashlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import g_nfw, rar, moster_mref, FB, G, ACC, KEYS as LVDB_CHECK_KEYS
from A19_Environment_Tests import c_dm14, g_nfw_c, r200_of

DATA_SHA = 'c5c029dce9970e7b70fc329cc46dbf2ed02b308901579a3f2aa756887d542b26'
V1_27 = {'aquarius_2', 'bootes_2', 'carina_2', 'carina_3', 'centaurus_1', 'columba_1', 'eridanus_4', 'grus_1', 'leo_6',
         'pegasus_3', 'pegasus_4', 'pictor_2', 'pisces_2', 'tucana_2', 'tucana_4', 'tucana_5', 'andromeda_17', 'andromeda_19',
         'andromeda_20', 'andromeda_26', 'andromeda_27', 'andromeda_28', 'andromeda_29', 'cassiopeia_2', 'cassiopeia_3',
         'lacerta_1', 'perseus_1'}
NAMED = {'sagittarius_1': 'Sagittarius I (tidally disrupting)', 'tucana_3': 'Tucana III (tidally disrupting)',
         'bootes_3': 'Bootes III (tidally disrupting)', 'ngc_4594-sucd1': 'SUCD1 (ultracompact dwarf)'}
MW_HOSTS = {'mw', 'lmc', 'smc'}
SEED = 20261031; A0_MM13 = 1.21e-10; V_M31 = 230.0; TOL = 10 ** (0.33 / 2)
rms = lambda v: float(np.sqrt(np.mean(np.square(v))))

def num(x, c):
    v = x.get(c, '')
    if v in ('', None) or v.strip().lower() == 'nan': return None
    return float(v)

def scored_before(bt):
    k44 = set(open(os.path.join(bt, 'original44_keys.txt')).read().split())
    k25 = {r['lvdb_key'] for r in csv.DictReader(open(os.path.join(bt, 'A19_OutOfSample_v2.csv')))}
    k42 = set(LVDB_CHECK_KEYS) - {'segue_2', 'andromeda_12'}
    return k44 | k25 | k42 | V1_27, dict(original_44=len(k44), protocol_2_25=len(k25), lvdb_check_42=len(k42), protocol_1_27=len(V1_27))

def select(path, before):
    sets = {'M': [], 'U': [], 'F': []}; excl = []; n_before = 0; n_notgal = 0
    for x in csv.DictReader(open(path)):
        k = x['key']
        if k in before: n_before += 1; continue                                    # step 1
        if x['confirmed_real'] != '1' or x['confirmed_galaxy'] != '1': n_notgal += 1; continue   # step 2
        def out(why): excl.append([k, x['name'], x['host'], why])
        MV = num(x, 'M_V')
        if MV is None: out('no M_V'); continue                                     # step 3
        if MV < -14: out('M_V brighter than -14'); continue
        if num(x, 'rhalf_physical') is None or num(x, 'rhalf_sph_physical') is None: out('no half-light radius'); continue   # step 4
        LV = 10 ** (-0.4 * (MV - 4.83)); hi = num(x, 'mass_HI')
        if hi is not None and not x['mass_HI_ul'] and 10 ** hi > 0.1 * 2 * LV: out('HI above a tenth of the stellar mass'); continue   # step 5
        if k in NAMED: out('named: ' + NAMED[k]); continue                         # step 6
        if num(x, 'vlos_sigma_ul') is not None: sets['U'].append(x)
        elif num(x, 'vlos_sigma') is not None: sets['M'].append(x)
        else: sets['F'].append(x)
    return sets, excl, n_before, n_notgal

def predict(rows, rules, ml=2.0, major=False, Vmw=200.0):
    S, B, A0 = rules['stellar_ABC'], rules['baryonic_ABC'], rules['rar_a0']
    LV = np.array([10 ** (-0.4 * (num(x, 'M_V') - 4.83)) for x in rows]); Ms = ml * LV
    Re = np.array([num(x, 'rhalf_physical' if major else 'rhalf_sph_physical') for x in rows]) / 1000
    r = 4 / 3 * Re; gb = G * (Ms / 2) / r ** 2 * ACC
    Mst = Ms / FB + 10 ** S[0] * (Ms / 1e10) ** S[1] * (Re / 3) ** S[2]
    Mba = Ms / FB + 10 ** B[0] * (Ms / 1e10) ** B[1] * (Re / 3) ** B[2]
    Mmo = np.array([max(moster_mref(m), m / FB) for m in Ms])
    fac = np.ones(len(rows)); gext = np.full(len(rows), np.nan); kind = []
    for i, x in enumerate(rows):
        h, dgc, dm31 = x['host'], num(x, 'distance_gc'), num(x, 'distance_m31')
        if (h in MW_HOSTS or h == '') and dgc is None or (h == 'm_031' or h == '') and dm31 is None:
            kind.append('no distance'); continue
        if h in MW_HOSTS:
            kind.append('satellite'); xs = dgc / r200_of(1.0e12); gext[i] = Vmw ** 2 / dgc
        elif h == 'm_031':
            kind.append('satellite'); xs = dm31 / r200_of(1.5e12); gext[i] = V_M31 ** 2 / dm31
        elif h == '':
            kind.append('no host'); xs = None; gext[i] = max(Vmw ** 2 / dgc, V_M31 ** 2 / dm31)
        else:
            kind.append('other group'); xs = None
        if xs is not None and xs < 1: fac[i] = 1 - 0.54 * np.log10(xs)
    gext_ms2 = gext * ACC
    s_iso = (4 / 81 * G * Ms * A0_MM13 / ACC) ** 0.25
    with np.errstate(invalid='ignore', divide='ignore'):
        s_efe = np.sqrt(G * Ms * (A0_MM13 / gext_ms2) / (3 * r))
        s_e2 = np.where(np.isnan(gext), np.nan, np.minimum(s_iso, s_efe))
    g = {'stellar_rule': gb + g_nfw(r, Mst), 'baryonic_rule': gb + g_nfw(r, Mba), 'moster': gb + g_nfw(r, Mmo),
         'rar_disk_a0': rar(gb, A0), 'rar_published_a0': rar(gb, 1.20e-10),
         'stellar_subhalo_E1': gb + g_nfw_c(r, Mst, c_dm14(Mst) * fac), 'mond_efe_E2': 3 * s_e2 ** 2 / r * ACC}
    sig = {k: np.sqrt(v / ACC * r / 3) for k, v in g.items()}
    return dict(LV=LV, Ms=Ms, Re=Re, r=r, gb=gb, g=g, sig=sig, fac=fac, kind=kind, gext_over_a0=gext_ms2 / 1.2e-10)

def boot(a, b, n):
    rng = np.random.default_rng(SEED); d = []
    for _ in range(20000):
        i = rng.integers(0, n, n); d.append(rms(a[i]) - rms(b[i]))
    return [round(rms(a) - rms(b), 4)] + np.percentile(d, [2.5, 97.5]).round(4).tolist()

def score_M(rows, rules, **kw):
    P = predict(rows, rules, **kw); s = np.array([num(x, 'vlos_sigma') for x in rows])
    lo = np.log10(3 * s ** 2 / P['r'] * ACC)
    R = {k: lo - np.log10(v) for k, v in P['g'].items()}
    out = {}
    for k, v in R.items():
        m = ~np.isnan(v); out[k] = dict(n=int(m.sum()), rmse_dex=round(rms(v[m]), 4) if m.any() else None)
    return P, lo, R, out

def score_U(rows, rules, **kw):
    P = predict(rows, rules, **kw); ul = np.array([num(x, 'vlos_sigma_ul') for x in rows]); out = {}
    for k, v in P['sig'].items():
        m = ~np.isnan(v); q = np.log10(v[m] / ul[m])
        out[k] = dict(n=int(m.sum()), exceed=int(np.sum(v[m] > ul[m])), exceed_by_factor_1_46=int(np.sum(v[m] > TOL * ul[m])),
                      median_log10_pred_over_ul=round(float(np.median(q)), 3))
    return P, ul, out

def main(lvdb, bt, out):
    sha = hashlib.sha256(open(lvdb, 'rb').read()).hexdigest(); assert sha == DATA_SHA, sha
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    before, counts = scored_before(bt); sets, excl, n_before, n_notgal = select(lvdb, before)
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v3.sha256')).read().split()[0], 'data_sha256': sha,
           'scored_before': counts, 'rows_scored_before_skipped': n_before, 'rows_not_confirmed_galaxies': n_notgal,
           'n': {k: len(v) for k, v in sets.items()}, 'excluded': excl}
    # Set M
    M = sets['M']; P, lo, R, outM = score_M(M, rules); n = len(M); st = R['stellar_rule']
    res['set_M'] = {'dwarfs': [[x['key'], x['name'], x['host'] or 'none', x['ref_vlos']] for x in M], 'rmse': outM,
                    'bootstrap_rmse_difference': {f'stellar_minus_{k}': boot(st[~np.isnan(v)], v[~np.isnan(v)], int((~np.isnan(v)).sum()))
                                                  for k, v in R.items() if k != 'stellar_rule'},
                    'closer_under_stellar': {k: int(np.sum(np.abs(st[~np.isnan(v)]) < np.abs(v[~np.isnan(v)]))) for k, v in R.items() if k != 'stellar_rule'},
                    'residuals_dex': {x['name']: {k: (None if np.isnan(R[k][i]) else round(float(R[k][i]), 3)) for k in R} for i, x in enumerate(M)},
                    'sigma_obs_kms': {x['name']: num(x, 'vlos_sigma') for x in M},
                    'sigma_pred_kms': {x['name']: {k: (None if np.isnan(P['sig'][k][i]) else round(float(P['sig'][k][i]), 2)) for k in P['sig']} for i, x in enumerate(M)}}
    b = res['set_M']['bootstrap_rmse_difference']['stellar_minus_rar_disk_a0']
    res['set_M']['primary_success'] = bool(outM['stellar_rule']['rmse_dex'] < outM['rar_disk_a0']['rmse_dex'] and b[2] < 0)
    res['set_M']['sensitivity'] = {lab: {k: v['rmse_dex'] for k, v in score_M(M, rules, **kw)[3].items()} for lab, kw in
                                   [('major_axis_radii', dict(major=True)), ('ML_1', dict(ml=1.0)), ('ML_3', dict(ml=3.0)),
                                    ('V_MW_175', dict(Vmw=175.0)), ('V_MW_220', dict(Vmw=220.0))]}
    # Set U
    U = sets['U']; PU, ul, outU = score_U(U, rules)
    res['set_U'] = {'dwarfs': [[x['key'], x['name'], x['host'], num(x, 'vlos_sigma_ul'), x['ref_vlos']] for x in U], 'tolerance_factor': round(TOL, 3),
                    'by_model': outU,
                    'sigma_pred_over_ul': {x['name']: {k: (None if np.isnan(PU['sig'][k][i]) else round(float(PU['sig'][k][i] / ul[i]), 2)) for k in PU['sig']} for i, x in enumerate(U)},
                    'sigma_pred_kms': {x['name']: {k: (None if np.isnan(PU['sig'][k][i]) else round(float(PU['sig'][k][i]), 2)) for k in PU['sig']} for i, x in enumerate(U)},
                    'subhalo_factor': {x['name']: round(float(PU['fac'][i]), 3) for i, x in enumerate(U)}}
    res['set_U']['primary_success'] = bool(outU['stellar_rule']['exceed'] <= outU['rar_disk_a0']['exceed'])
    res['set_U']['E1_lead_passes'] = bool(outU['stellar_subhalo_E1']['exceed'] <= outU['stellar_rule']['exceed'])
    res['set_U']['sensitivity'] = {lab: {k: [v['exceed'], v['exceed_by_factor_1_46']] for k, v in score_U(U, rules, **kw)[2].items()} for lab, kw in
                                   [('major_axis_radii', dict(major=True)), ('ML_1', dict(ml=1.0)), ('ML_3', dict(ml=3.0)),
                                    ('V_MW_175', dict(Vmw=175.0)), ('V_MW_220', dict(Vmw=220.0))]}
    # Set F
    F = sets['F']; PF = predict(F, rules); sg = PF['sig']
    ref = np.where(np.isnan(sg['mond_efe_E2']), sg['rar_published_a0'], sg['mond_efe_E2'])
    ratio = sg['stellar_rule'] / ref; decisive = (ratio > 2) | (ratio < 0.5)
    in_table = np.array([(x['host'] in MW_HOSTS or x['host'] == 'm_031') or (x['host'] == '' and (num(x, 'distance_gc') or 1e9) < 1500) for x in F])
    fpath = os.path.join(os.path.dirname(os.path.abspath(out)), 'A19_Forecasts_v3.csv')
    with open(fpath, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['lvdb_key', 'name', 'host', 'kind', 'distance_gc_kpc', 'M_V', 'L_V_Lsun', 'R_e_circ_pc', 'sigma_stellar_rule', 'sigma_stellar_subhalo_E1',
                    'sigma_moster', 'sigma_rar_disk_a0', 'sigma_rar_published_a0', 'sigma_mond_efe_E2', 'mond_reference', 'ratio_stellar_over_mond', 'decisive', 'in_paper_table'])
        for i, x in enumerate(F):
            f2 = lambda v: '' if np.isnan(v) else f'{v:.2f}'
            w.writerow([x['key'], x['name'], x['host'] or 'none', PF['kind'][i], ('' if num(x, 'distance_gc') is None else f"{num(x, 'distance_gc'):.1f}"), x['M_V'], f'{PF["LV"][i]:.4g}', f'{PF["Re"][i] * 1000:.1f}',
                        f2(sg['stellar_rule'][i]), f2(sg['stellar_subhalo_E1'][i]), f2(sg['moster'][i]), f2(sg['rar_disk_a0'][i]), f2(sg['rar_published_a0'][i]),
                        f2(sg['mond_efe_E2'][i]), 'E2' if not np.isnan(sg['mond_efe_E2'][i]) else 'rar_published_a0', f'{ratio[i]:.2f}', int(decisive[i]), int(in_table[i])])
    res['set_F'] = {'n': len(F), 'n_in_paper_table': int(in_table.sum()), 'n_decisive': int(decisive.sum()), 'n_decisive_in_paper_table': int((decisive & in_table).sum()),
                    'by_kind': {k: int(sum(1 for q in PF['kind'] if q == k)) for k in ('satellite', 'no host', 'other group')},
                    'median_ratio_stellar_over_mond': round(float(np.median(ratio)), 2),
                    'median_ratio_in_paper_table': round(float(np.median(ratio[in_table])), 2) if in_table.any() else None,
                    'forecast_file_sha256': hashlib.sha256(open(fpath, 'rb').read()).hexdigest()}
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    with open(out + '_Measured_and_Limits.csv', 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['set', 'lvdb_key', 'name', 'host', 'L_V', 'R_e_circ_pc', 'sigma_or_upper_limit_kms'] + ['sigma_pred_' + k for k in P['sig']] + ['ref_structure', 'ref_vlos', 'ref_m_v'])
        for lab, rows, PP, val in (('M', M, P, [num(x, 'vlos_sigma') for x in M]), ('U', U, PU, list(ul))):
            for i, x in enumerate(rows):
                w.writerow([lab, x['key'], x['name'], x['host'] or 'none', f'{PP["LV"][i]:.4g}', f'{PP["Re"][i] * 1000:.1f}', val[i]] +
                           [('' if np.isnan(PP['sig'][k][i]) else f'{PP["sig"][k][i]:.2f}') for k in PP['sig']] + [x['ref_structure'], x['ref_vlos'], x['ref_m_v']])
    print(json.dumps({k: res[k] for k in ('n', 'scored_before')}, indent=1))
    print(json.dumps({k: v for k, v in res['set_M'].items() if k not in ('residuals_dex',)}, indent=1))
    print(json.dumps(res['set_U'], indent=1)); print(json.dumps(res['set_F'], indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
