"""A19 supplement: environment and error tests on the dwarf spheroidals, all with published formulas and
no fitting to dwarfs.
(1) MOND with the external field effect (Chae et al. 2020, Eqs. 5-6; host field g_ext = V_host^2 / D,
    V_M31 = 230 km/s as in McGaugh and Milgrom 2013; V_MW = 175, 200, 220 km/s).
(2) Subhalo concentration near a host (Moliné et al. 2017, Eq. 6): the stellar rule and the Moster relation
    with the radial factor [1 + b log10(R_sub/R_200,host)], b = -0.54, and with the full Moliné c(m, x).
(3) Measurement-error budget: how much of each model's scatter the quoted errors of g_obs explain.
Dwarf sets: the original 44 (McConnachie January 2021, built by the original pipeline) and the 25
out-of-sample dwarfs of protocol version 2 (Local Volume Database)."""
import sys, os, csv, json, numpy as np
from scipy.optimize import minimize_scalar, minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Dwarfs_LVDB import moster_mref, rar, FB, G, ACC
h = 0.671; H0 = 100 * h
RHOC = 3 * (H0 * 1e-3) ** 2 / (8 * np.pi * G)                  # Msun/kpc^3 (G in kpc (km/s)^2/Msun)
def r200_of(M): return (3 * M / (4 * np.pi * 200 * RHOC)) ** (1 / 3)
HOST = {'mw': dict(M200=1.0e12, V=200.0), 'm31': dict(M200=1.5e12, V=230.0)}
def c_dm14(M): return 10 ** 0.905 * (M / (1e12 / h)) ** -0.101
def c_moline(M, x):
    lm = np.log10(M * h / 1e8); a = (-0.195, 0.089, 0.089)
    return 19.9 * (1 + sum((a[i] * lm) ** (i + 1) for i in range(3))) * (1 + -0.54 * np.log10(x))
def g_nfw_c(r, Mref, c):
    r200 = r200_of(Mref); x = np.minimum(r / (r200 / c), c); m = lambda y: np.log(1 + y) - y / (1 + y)
    return G * (1 - FB) * Mref * m(x) / m(c) / r ** 2 * ACC
def nu_efe(z, e):
    A = e * (1 + e / 2) / (1 + e); B = 1 + e
    return 0.5 - A / z + np.sqrt((0.5 - A / z) ** 2 + B / z)
rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
LVDB_PATH = None
def lvdb_index():
    return {x['key']: x for x in csv.DictReader(open(LVDB_PATH))}
def load_sets(bench, oos, keys44):
    lv = lvdb_index(); sets = {}
    rows = list(csv.DictReader(open(bench))); f = lambda k: np.array([float(x[k]) for x in rows])
    ks = open(keys44).read().split()  # same order as the mapping? map by order of systems below
    import re
    rom = {'i':1,'ii':2,'iii':3,'iv':4,'v':5,'vi':6,'vii':7,'viii':8,'ix':9,'x':10,'xi':11,'xii':12,'xiii':13,'xiv':14,'xv':15,'xvi':16,'xvii':17,'xviii':18,'xix':19,'xx':20,'xxi':21,'xxii':22,'xxiii':23,'xxiv':24,'xxv':25,'xxvi':26,'xxvii':27,'xxviii':28,'xxix':29}
    def key_of(s):
        s = s.replace('(1)', '1'); m = re.match(r'^(Andromeda)([IVX]+)$', s)
        if m: return 'andromeda_%02d' % rom[m.group(2).lower()]
        m = re.match(r'^([A-Za-z]+?)(\d+)$', s); base, num = (m.group(1), int(m.group(2))) if m else (s, 1)
        return re.sub(r'(?<!^)([A-Z])', r'_\1', base).lower() + f'_{num}'
    K = [key_of(x['system']) for x in rows]
    def env(k):
        x = lv[k]; hst = x['host']
        if hst == 'm_031': return 'm31', float(x['distance_m31'])
        return 'mw', float(x['distance_gc'])
    E = [env(k) for k in K]
    sets['original_44'] = dict(names=[x['system'] for x in rows], R=f('radius_kpc'), gb=f('gbar_m_s2'), lo=np.log10(f('gobs_m_s2')),
                              Re=f('projected_Re_circularized_kpc'), Ms=f('stellar_mass_Msun'), host=[e[0] for e in E], D=np.array([e[1] for e in E]),
                              err=np.array([float(x['err_log10_gobs_dex_stat_approx']) if x['err_log10_gobs_dex_stat_approx'] else
                                  2 * 0.5 * (float(x['sigma_err_plus_km_s']) + float(x['sigma_err_minus_km_s'])) / float(x['sigma_los_km_s']) / np.log(10) for x in rows]))
    rows = list(csv.DictReader(open(oos))); g = lambda k: np.array([float(x[k]) for x in rows])
    Re = g('R_e_circ_kpc'); Ms = 2 * g('L_V'); R = 4 / 3 * Re; E = [env(x['lvdb_key']) for x in rows]
    err = []
    for x in rows:
        y = lv[x['lvdb_key']]; s = float(y['vlos_sigma']); ds = 0.5 * (float(y['vlos_sigma_em'] or 0) + float(y['vlos_sigma_ep'] or 0))
        rr = float(y['rhalf_sph_physical']); dr = 0.5 * (float(y['rhalf_sph_physical_em'] or 0) + float(y['rhalf_sph_physical_ep'] or 0))
        err.append(np.sqrt((2 * ds / s) ** 2 + (dr / rr) ** 2) / np.log(10))
    sets['out_of_sample_25'] = dict(names=[x['name'] for x in rows], R=R, gb=G * (Ms / 2) / R ** 2 * ACC, lo=g('log_g_obs'), Re=Re, Ms=Ms,
                                    host=[e[0] for e in E], D=np.array([e[1] for e in E]), err=np.array(err))
    return sets
def models(d, rules, Vmw=200.0):
    S, B, A0 = rules['stellar_ABC'], rules['baryonic_ABC'], rules['rar_a0']
    R, gb, Re, Ms, D = d['R'], d['gb'], d['Re'], d['Ms'], d['D']
    Mst = Ms / FB + 10 ** S[0] * (Ms / 1e10) ** S[1] * (Re / 3) ** S[2]
    Mmo = np.array([max(moster_mref(m), m / FB) for m in Ms])
    xsub = np.array([D[i] / r200_of(HOST[d['host'][i]]['M200']) for i in range(len(D))])
    fac = 1 + -0.54 * np.log10(xsub)
    Vh = np.array([Vmw if hh == 'mw' else HOST['m31']['V'] for hh in d['host']])
    gext = Vh ** 2 / D * ACC
    out = {'stellar_rule': gb + g_nfw_c(R, Mst, c_dm14(Mst)),
           'stellar_rule_subhalo_factor': gb + g_nfw_c(R, Mst, c_dm14(Mst) * fac),
           'stellar_rule_moline_c': gb + g_nfw_c(R, Mst, c_moline(Mst, xsub)),
           'moster': gb + g_nfw_c(R, Mmo, c_dm14(Mmo)),
           'moster_subhalo_factor': gb + g_nfw_c(R, Mmo, c_dm14(Mmo) * fac),
           'mond_isolated_a0_1.2': rar(gb, 1.2e-10), 'mond_isolated_disk_a0': rar(gb, A0),
           'mond_efe_a0_1.2': gb * nu_efe(gb / 1.2e-10, gext / 1.2e-10), 'mond_efe_disk_a0': gb * nu_efe(gb / A0, gext / A0)}
    return {k: d['lo'] - np.log10(v) for k, v in out.items()}, xsub, gext
def intrinsic(res, err):
    f = lambda p: np.sum((res - p[0]) ** 2 / (err ** 2 + p[1] ** 2) + np.log(err ** 2 + p[1] ** 2))
    o = minimize(f, [0.0, 0.3], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-9})
    return {'bias_dex': round(float(o.x[0]), 3), 'intrinsic_scatter_dex': round(abs(float(o.x[1])), 3)}
def main(bench, oos, keys44, rulesjson, out):
    rules = json.load(open(rulesjson))['frozen_on_all_disks']; sets = load_sets(bench, oos, keys44); res = {}
    for name, d in sets.items():
        R, xsub, gext = models(d, rules)
        block = {'n': len(d['lo']), 'rmse_dex': {k: round(rms(v), 4) for k, v in R.items()},
                 'median_x_sub': round(float(np.median(xsub)), 3), 'median_gext_over_a0': round(float(np.median(gext / 1.2e-10)), 3),
                 'measurement_error_rms_dex': round(rms(d['err']), 3)}
        block['error_budget'] = {k: intrinsic(R[k], d['err']) for k in ('stellar_rule', 'stellar_rule_subhalo_factor', 'moster', 'mond_efe_a0_1.2')}
        block['efe_V_MW_sensitivity'] = {f'V_MW_{v:.0f}': round(rms(models(d, rules, v)[0]['mond_efe_a0_1.2']), 4) for v in (175.0, 220.0)}
        block['residual_vs_log_host_distance_corr'] = {k: round(float(np.corrcoef(np.log10(d['D']), R[k])[0, 1]), 3) for k in ('stellar_rule', 'moster', 'mond_efe_a0_1.2')}
        res[name] = block
    json.dump(res, open(out + '_Results.json', 'w'), indent=1); print(json.dumps(res, indent=1))
def mm13(sets, rules, out):
    # MOND with the external field as McGaugh and Milgrom (2013): sigma = min(sigma_iso, sigma_efe), g = 3 sigma^2 / r_1/2
    a0 = 1.21e-10; res = {}
    for name, d in sets.items():
        R, _, _ = models(d, rules); M = d['Ms']; L = M / 2; br = L >= 1e5
        for Vmw in (175.0, 200.0, 220.0):
            Vh = np.array([Vmw if h == 'mw' else 230.0 for h in d['host']]); gex = Vh ** 2 / d['D'] * ACC
            s_iso = (4 / 81 * G * M * a0 / ACC) ** 0.25; s_efe = np.sqrt(G * M * (a0 / gex) / (3 * d['R']))
            r = d['lo'] - np.log10(3 * np.minimum(s_iso, s_efe) ** 2 / d['R'] * ACC)
            res[f'{name}_V_MW_{Vmw:.0f}'] = {'rmse_dex': round(rms(r), 4), 'mean_residual_dex': round(float(np.mean(r)), 3),
                                             'bright_rmse': round(rms(r[br]), 3), 'faint_rmse': round(rms(r[~br]), 3), 'n_bright': int(br.sum())}
        res[f'{name}_stellar_rule_bright_faint'] = [round(rms(R['stellar_rule'][br]), 3), round(rms(R['stellar_rule'][~br]), 3)]
    json.dump(res, open(out + '_MM13_Results.json', 'w'), indent=1); print(json.dumps(res, indent=1))
if __name__ == '__main__':
    # Usage: python3 A19_Environment_Tests.py dwarfs_benchmark.csv A19_OutOfSample_v2.csv original44_keys.txt frozen_disk_rules.json dwarf_all.csv out_prefix
    LVDB_PATH = sys.argv[5]
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[6])
    mm13(load_sets(sys.argv[1], sys.argv[2], sys.argv[3]), json.load(open(sys.argv[4]))['frozen_on_all_disks'], sys.argv[6])
