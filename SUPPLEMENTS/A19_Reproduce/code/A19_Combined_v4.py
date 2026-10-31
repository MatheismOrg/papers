"""A19 protocol version 4, implemented exactly as written in PROTOCOL_v4.txt.
One likelihood over every dwarf scored so far: a measured dispersion counts as a value with its error, an
upper limit as a ceiling. Also writes 68 percent ranges for protocol 3's forecasts.
Usage: python3 A19_Combined_v4.py dwarf_all.csv blind_test_dir out_prefix
blind_test_dir holds frozen_disk_rules.json, original44_keys.txt, A19_OutOfSample_v2.csv,
A19_New_Dwarfs_v3_Results.json, A19_Forecasts_v3.csv and PROTOCOL_v4.sha256. Nothing is fitted except each
model's scatter."""
import sys, os, csv, json, hashlib, numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import G, ACC
import A19_New_Dwarfs_v3 as V3

DATA_SHA = 'c5c029dce9970e7b70fc329cc46dbf2ed02b308901579a3f2aa756887d542b26'
SEED = 20261031; NBOOT = 5000; A0_MM13 = 1.21e-10
MODELS = ['stellar_rule', 'stellar_subhalo_E1', 'baryonic_rule', 'moster', 'rar_disk_a0', 'rar_published_a0', 'mond_efe_E2']
num = V3.num

def load_rows(lvdb, bt):
    db = {r['key']: r for r in csv.DictReader(open(lvdb))}
    k44 = open(os.path.join(bt, 'original44_keys.txt')).read().split()
    k25 = [r['lvdb_key'] for r in csv.DictReader(open(os.path.join(bt, 'A19_OutOfSample_v2.csv')))]
    v3 = json.load(open(os.path.join(bt, 'A19_New_Dwarfs_v3_Results.json')))
    keys = k44 + k25 + [d[0] for d in v3['set_M']['dwarfs']] + ['andromeda_18'] + [d[0] for d in v3['set_U']['dwarfs']]
    assert len(keys) == len(set(keys)) == 85
    src = {k: 'original_44' for k in k44}; src.update({k: 'protocol_2' for k in k25})
    src.update({d[0]: 'protocol_3' for d in v3['set_M']['dwarfs'] + v3['set_U']['dwarfs']}); src['andromeda_18'] = 'lvdb_check'
    return [db[k] for k in keys], src

def arrays(rows, rules, ml=2.0, major=False):
    P = V3.predict(rows, rules, ml=ml, major=major)
    r = P['r']; Ms = P['Ms']
    g = dict(P['g'])
    iso = np.isnan(g['mond_efe_E2'])                      # dwarfs in other groups: isolated MOND estimate
    s_iso = (4 / 81 * G * Ms * A0_MM13 / ACC) ** 0.25
    g['mond_efe_E2'] = np.where(iso, 3 * s_iso ** 2 / r * ACC, g['mond_efe_E2'])
    is_ul = np.array([num(x, 'vlos_sigma_ul') is not None for x in rows])
    sig = np.array([num(x, 'vlos_sigma_ul') if u else num(x, 'vlos_sigma') for x, u in zip(rows, is_ul)])
    lo = np.log10(3 * sig ** 2 / r * ACC)                 # log g_obs, or log g_ul for limits
    e = np.zeros(len(rows))
    for i, x in enumerate(rows):
        if is_ul[i]: continue
        ds = np.nanmean([num(x, 'vlos_sigma_em') or np.nan, num(x, 'vlos_sigma_ep') or np.nan])
        cr = 'rhalf_physical' if major else 'rhalf_sph_physical'
        dR = np.nanmean([num(x, cr + '_em') or np.nan, num(x, cr + '_ep') or np.nan])
        R = num(x, cr)
        e[i] = np.sqrt((2 * ds / sig[i]) ** 2 + (dR / R) ** 2) / np.log(10)
    bright = P['LV'] >= 1e5
    return dict(lo=lo, e=e, ul=is_ul, bright=bright, sig=sig, logg={k: np.log10(g[k]) for k in MODELS}, n_iso=int(iso.sum()))

def lnl_terms(res, e, ul, s, b=0.0):
    rr = res - b; v = s * s + e * e
    meas = -0.5 * np.log(2 * np.pi * v) - 0.5 * rr * rr / v
    lim = norm.logcdf(rr / s)                              # rr = log g_ul - log g_m for limits
    return np.where(ul, lim, meas)

def fit_s(res, e, ul, idx=None):
    if idx is not None: res, e, ul = res[idx], e[idx], ul[idx]
    f = lambda s: -np.sum(lnl_terms(res, e, ul, s))
    o = minimize_scalar(f, bounds=(0.01, 3.0), method='bounded', options={'xatol': 1e-7})
    return float(o.x), float(-o.fun)

def fit_bs(res, e, ul):
    f = lambda p: -np.sum(lnl_terms(res, e, ul, min(max(p[1], 0.01), 3.0), p[0]))
    o = minimize(f, [0.0, 0.4], method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 4000})
    return dict(offset_dex=round(float(o.x[0]), 3), s_dex=round(float(min(max(o.x[1], 0.01), 3.0)), 3), lnL=round(float(-o.fun), 2))

def score(A, keep=None):
    out = {}
    for m in MODELS:
        res = A['lo'] - A['logg'][m]; e, ul = A['e'], A['ul']
        if keep is not None: res, e, ul = res[keep], e[keep], ul[keep]
        s, L = fit_s(res, e, ul); t = lnl_terms(res, e, ul, s)
        br = A['bright'] if keep is None else A['bright'][keep]
        out[m] = dict(s_dex=round(s, 4), lnL=round(L, 2), n=int(len(res)),
                      split=dict(measured_bright=round(float(t[~ul & br].sum()), 2), measured_faint=round(float(t[~ul & ~br].sum()), 2),
                                 limits=round(float(t[ul].sum()), 2)))
    return out

def main(lvdb, bt, out):
    sha = hashlib.sha256(open(lvdb, 'rb').read()).hexdigest(); assert sha == DATA_SHA, sha
    rules = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']
    rows, src = load_rows(lvdb, bt); A = arrays(rows, rules); n = len(rows)
    res = {'protocol_sha256': open(os.path.join(bt, 'PROTOCOL_v4.sha256')).read().split()[0], 'data_sha256': sha,
           'n': n, 'n_measured': int((~A['ul']).sum()), 'n_limits': int(A['ul'].sum()), 'n_measured_bright': int((~A['ul'] & A['bright']).sum()),
           'n_E2_isolated_estimate': A['n_iso'], 'models': score(A)}
    L = res['models']
    res['delta_lnL_stellar_minus'] = {m: round(L['stellar_rule']['lnL'] - L[m]['lnL'], 2) for m in MODELS if m != 'stellar_rule'}
    # bootstrap over whole dwarfs, scatter refitted for both models each time
    rng = np.random.default_rng(SEED); R = {m: A['lo'] - A['logg'][m] for m in MODELS}
    boots = {m: [] for m in MODELS if m != 'stellar_rule'}
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        Ls = fit_s(R['stellar_rule'], A['e'], A['ul'], i)[1]
        for m in boots: boots[m].append(Ls - fit_s(R[m], A['e'], A['ul'], i)[1])
    res['bootstrap_95'] = {m: np.percentile(v, [2.5, 97.5]).round(2).tolist() for m, v in boots.items()}
    bp = res['bootstrap_95']['rar_disk_a0']
    res['primary'] = dict(delta_lnL=res['delta_lnL_stellar_minus']['rar_disk_a0'], interval_95=bp,
                          success=bool(res['delta_lnL_stellar_minus']['rar_disk_a0'] > 0 and bp[0] > 0))
    # sensitivities
    sens = {}
    sens['fitted_offset'] = {m: fit_bs(R[m], A['e'], A['ul']) for m in MODELS}
    for lab, kw in (('ML_1', dict(ml=1.0)), ('ML_3', dict(ml=3.0)), ('major_axis_radii', dict(major=True))):
        sc = score(arrays(rows, rules, **kw)); sens[lab] = {m: dict(lnL=sc[m]['lnL'], s_dex=sc[m]['s_dex']) for m in MODELS}
    keep = A['ul'] | (A['sig'] >= 3.0)
    sc = score(A, keep); sens['measured_under_3kms_dropped'] = dict(n=int(keep.sum()), **{m: dict(lnL=sc[m]['lnL'], s_dex=sc[m]['s_dex']) for m in MODELS})
    A0 = dict(A); A0['e'] = np.zeros(n); sc = score(A0); sens['no_measurement_errors'] = {m: dict(lnL=sc[m]['lnL'], s_dex=sc[m]['s_dex']) for m in MODELS}
    for lab, v in sens.items():
        if lab != 'fitted_offset':
            v['delta_lnL_stellar_minus'] = {m: round(v['stellar_rule']['lnL'] - v[m]['lnL'], 2) for m in MODELS if m != 'stellar_rule'}
    res['sensitivity'] = sens
    # forecast ranges for protocol 3's forecasts
    s_st = L['stellar_rule']['s_dex']; s_e2 = L['mond_efe_E2']['s_dex']; s_pub = L['rar_published_a0']['s_dex']
    fc = list(csv.DictReader(open(os.path.join(bt, 'A19_Forecasts_v3.csv'))))
    fpath = os.path.join(os.path.dirname(os.path.abspath(out)), 'A19_Forecast_Ranges_v4.csv')
    zs = []
    with open(fpath, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['lvdb_key', 'name', 'host', 'kind', 'sigma_stellar_rule', 'stellar_lo68', 'stellar_hi68', 'mond_reference', 'sigma_mond', 'mond_lo68', 'mond_hi68',
                    'separation_z', 'in_paper_table'])
        for x in fc:
            ss = float(x['sigma_stellar_rule'])
            if x['sigma_mond_efe_E2']: sm, ref, sM = float(x['sigma_mond_efe_E2']), 'E2', s_e2
            else: sm, ref, sM = float(x['sigma_rar_published_a0']), 'rar_published_a0', s_pub
            z = abs(np.log10(ss) - np.log10(sm)) / np.sqrt((s_st / 2) ** 2 + (sM / 2) ** 2); zs.append((x['in_paper_table'] == '1', z))
            f = lambda v, s, k: f'{v * 10 ** (k * s / 2):.2f}'
            w.writerow([x['lvdb_key'], x['name'], x['host'], x['kind'], f'{ss:.2f}', f(ss, s_st, -1), f(ss, s_st, 1), ref, f'{sm:.2f}', f(sm, sM, -1), f(sm, sM, 1),
                        f'{z:.2f}', x['in_paper_table']])
    zt = [z for t, z in zs if t]
    res['forecast_ranges'] = dict(stellar_factor_68=round(10 ** (s_st / 2), 3), mond_E2_factor_68=round(10 ** (s_e2 / 2), 3),
                                  median_z_paper_table=round(float(np.median(zt)), 2), combined_z_paper_table=round(float(np.sqrt(np.sum(np.square(zt)))), 2),
                                  file_sha256=hashlib.sha256(open(fpath, 'rb').read()).hexdigest())
    res['per_dwarf'] = {x['name']: dict(source=src[x['key']], limit=bool(A['ul'][i]), log_g=round(float(A['lo'][i]), 3), err_dex=round(float(A['e'][i]), 3),
                                        **{m: round(float(A['lo'][i] - A['logg'][m][i]), 3) for m in MODELS}) for i, x in enumerate(rows)}
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('per_dwarf', 'sensitivity')}, indent=1))
    print(json.dumps(res['sensitivity'], indent=1))

if __name__ == '__main__': main(*sys.argv[1:4])
