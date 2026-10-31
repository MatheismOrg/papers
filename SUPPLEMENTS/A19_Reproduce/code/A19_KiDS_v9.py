"""A19 protocol version 9: the frozen stellar rule forward-modelled into the KiDS-1000 lensing radial acceleration
relation of isolated galaxies (Brouwer et al. 2021), implemented as written in PROTOCOL_v9.txt and PROTOCOL_v9_NOTE.txt.
Usage:
  python3 A19_KiDS_v9.py predict data blind_test out_prefix     (model only; no lensing values are read)
  python3 A19_KiDS_v9.py compare data blind_test out_prefix     (reads the KiDS profiles and covariances)
Needs numpy and scipy; colossus only for the optional NFW cross-check (check_nfw)."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_Disk_Fit_SPARC as DF

GMSUN = 1.32712440018e20          # m^3 s^-2
MPC = 3.0856775814913673e22       # m
HR, H0, OM = 0.671, 67.1, 0.3175  # the rule's cosmology
FB = DF.FB
Z0 = 0.25
GEDGES = np.logspace(-15, np.log10(5e-12), 16)          # B21: 15 logarithmic g_bar bins
R_INNER = 0.3                                            # Mpc (0.3 h70^-1 Mpc, h = 0.7)
VDW_LM = np.array([9.25, 9.75, 10.25, 10.75])
VDW_D = np.array([-0.22, -0.33, -0.28, -0.16])           # van der Wel et al. 2014, Table 2, z = 0.25: ETG - LTG medians
VDW_ETG = (np.array([9.25, 9.75, 10.25, 10.75, 11.25]), np.array([0.27, 0.28, 0.38, 0.67, 0.76]))
VDW_LTG = (np.array([9.25, 9.75, 10.25, 10.75]), np.array([0.49, 0.61, 0.66, 0.83]))
FIG9_LIMITS = [8.5, 10.3, 10.6, 10.8, 11.0]
FIG9_MGAL = [10.14, 10.57, 10.78, 10.96]                  # B21 Section 5.5, log <M_gal>
GL_X, GL_W = np.polynomial.legendre.leggauss(32)
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(9)

def rho_crit(z):
    """Critical density in suns per Mpc^3 for H0 = 67.1, Omega_m = 0.3175, flat."""
    h0 = H0 * 1e3 / MPC
    return 3 * h0 ** 2 * (OM * (1 + z) ** 3 + 1 - OM) * MPC ** 3 / (8 * np.pi * GMSUN)

def f_cold(lm): return 10 ** (-0.69 * lm + 6.63)
def m_gal(lm): return 10 ** lm * (1 + f_cold(lm))

def g_wb(x):
    """Wright and Brainerd (2000) Delta Sigma / (r_s rho_s) for an untruncated NFW."""
    x = np.asarray(x, float); out = np.empty_like(x)
    lo, hi, eq = x < 1 - 1e-6, x > 1 + 1e-6, np.abs(x - 1) <= 1e-6
    a = x[lo]; t = np.arctanh(np.sqrt((1 - a) / (1 + a)))
    out[lo] = 8 * t / (a ** 2 * np.sqrt(1 - a ** 2)) + 4 / a ** 2 * np.log(a / 2) - 2 / (a ** 2 - 1) + 4 * t / ((a ** 2 - 1) * np.sqrt(1 - a ** 2))
    b = x[hi]; t = np.arctan(np.sqrt((b - 1) / (1 + b)))
    out[hi] = 8 * t / (b ** 2 * np.sqrt(b ** 2 - 1)) + 4 / b ** 2 * np.log(b / 2) - 2 / (b ** 2 - 1) + 4 * t / (b ** 2 - 1) ** 1.5
    out[eq] = 10 / 3 + 4 * np.log(0.5)
    return out

def halo_params(mref, z, conc):
    r200 = (3 * mref / (4 * np.pi * 200 * rho_crit(z))) ** (1 / 3)
    c = conc(mref, z); rs = r200 / c
    dc = 200 / 3 * c ** 3 / (np.log(1 + c) - c / (1 + c))
    return r200, c, rs, dc * rho_crit(z)

def esd_nfw(R, mref, z, conc, dark=True, truncated=False):
    """Delta Sigma in suns per Mpc^2 at projected radius R (Mpc)."""
    r200, c, rs, rhos = halo_params(mref, z, conc)
    amp = (1 - FB) if dark else 1.0
    if not truncated:
        return amp * rs * rhos * g_wb(R / rs)
    # sharply truncated at r200: numerical projection
    def sigma(Rp):
        Rp = np.atleast_1d(Rp); out = np.zeros_like(Rp)
        for i, r in enumerate(Rp):
            if r >= r200: continue
            L = np.sqrt(r200 ** 2 - r ** 2); u = 0.5 * L * (GL_X + 1); rr = np.sqrt(r ** 2 + u ** 2) / rs
            out[i] = 2 * np.sum(0.5 * L * GL_W * rhos / (rr * (1 + rr) ** 2))
        return out
    R = np.atleast_1d(R); res = np.zeros_like(R)
    for i, r in enumerate(R):
        # mean inside R: (2/R^2) int_0^R Sigma R' dR', integrated in ln R' from 1e-5 Mpc
        lo = np.log(1e-5); hi = np.log(r); t = 0.5 * (hi - lo) * (GL_X + 1) + lo; rp = np.exp(t)
        inner = np.sum(0.5 * (hi - lo) * GL_W * sigma(rp) * rp ** 2)
        mean_in = 2 * inner / r ** 2 + 1e-10 ** 2 * 0
        # enclosed core below 1e-5 Mpc is negligible for these halos
        res[i] = mean_in - sigma(np.array([r]))[0]
    return amp * res

def annulus_mean(func, rin, rout):
    lo, hi = np.log(rin), np.log(rout); t = 0.5 * (hi - lo) * (GL_X + 1) + lo; r = np.exp(t)
    return np.sum(0.5 * (hi - lo) * GL_W * func(r) * r ** 2) / ((rout ** 2 - rin ** 2) / 2)

def c_rule(m, z): return 10 ** 0.905 * (m * HR / 1e12) ** -0.101
def c_dm14(m, z): return 10 ** (0.520 + 0.385 * np.exp(-0.617 * z ** 1.21) + (-0.101 + 0.026 * z) * np.log10(m * HR / 1e12))

def interp_hold(x, xs, ys): return float(np.interp(x, xs, ys))   # np.interp holds the end values

class Rule:
    def __init__(self, dd, bt):
        self.A, self.B, self.C = json.load(open(os.path.join(bt, 'frozen_disk_rules.json')))['frozen_on_all_disks']['stellar_ABC']
        gal = DF.load(DF.table1(os.path.join(dd, 'SPARC', 'SPARC_Lelli2016c.mrt')), os.path.join(dd, 'SPARC', 'Rotmod_LTG'))
        lms = np.log10([g['Ms'] for g in gal]); lre = np.log10([g['Re'] for g in gal]); mhi = np.array([(g['Mp'] - g['Ms']) / 1.33 for g in gal])
        self.p1, self.p0 = np.polyfit(lms, lre, 1); k = mhi > 0; self.q1, self.q0 = np.polyfit(lms[k], np.log10(mhi[k]), 1)
        assert abs(self.p0 + 1.6721) < 1e-4 and abs(self.p1 - 0.2245) < 1e-4 and abs(self.q0 - 4.1156) < 1e-4 and abs(self.q1 - 0.5396) < 1e-4
    def mref(self, lm, early, case):
        ms = 10 ** lm
        re_disk = 10 ** (self.p0 + self.p1 * lm); mhi = 10 ** (self.q0 + self.q1 * lm)
        if case == 'vdw_sizes':
            re = 10 ** interp_hold(lm, *VDW_ETG) if early else 10 ** interp_hold(lm, *VDW_LTG)
        elif early:
            if case == 'early_half_disk': re = 0.5 * re_disk
            elif case == 'early_disk_size': re = re_disk
            else: re = re_disk * 10 ** interp_hold(lm, VDW_LM, VDW_D)
        else:
            re = re_disk
        gas = (not early) or case == 'early_gas'
        mp = ms + (1.33 * mhi if gas else 0.0)
        return mp / FB + 10 ** self.A * (ms / 1e10) ** self.B * (re / 3.0) ** self.C, mp

def lens_model(rule, lm_pub, early, case):
    """Annulus-averaged model ESD (suns per pc^2) of one lens in each of the 15 g_bar bins."""
    shift = {'ms_plus01': 0.1, 'ms_minus01': -0.1}.get(case, 0.0)
    z = {'z015': 0.15, 'z035': 0.35}.get(case, Z0)
    conc = c_dm14 if case == 'c_dm14' else c_rule
    dark = case != 'nfw_full_mass'; trunc = case == 'truncated'
    mref, mp = rule.mref(lm_pub + shift, early, case)
    mg = m_gal(lm_pub)
    r = np.sqrt(GMSUN * mg / GEDGES) / MPC          # radius at each edge, Mpc (decreasing with g)
    out = np.zeros(15)
    if case == 'scatter015':
        masses = [(mref * 10 ** (0.15 * x), w / np.sqrt(2 * np.pi)) for x, w in zip(GH_X, GH_W)]
    else:
        masses = [(mref, 1.0)]
    for k in range(15):
        rout, rin = r[k], r[k + 1]
        val = 0.0
        for m, w in masses:
            val += w * annulus_mean(lambda R: esd_nfw(R, m, z, conc, dark, trunc), rin, rout)
        val += mp / np.pi * np.log(rout / rin) * 2 / (rout ** 2 - rin ** 2)
        out[k] = val / 1e12
    return out, mref, mg

def populations(fig7, case):
    shift = {'dist_plus01': 0.1, 'dist_minus01': -0.1}.get(case, 0.0)
    bins = fig7['bins']
    def rows(kind, lo, hi):
        R = []
        for b in bins:
            c = 0.5 * (b['lo'] + b['hi'])
            if not (lo <= c < hi): continue
            if kind == 'matched' and b['matched'] > 0:
                R.append((c, b['matched']))
            elif kind in ('red', 'blue') and b[kind] > 0:
                R.append((c, b[kind]))
        return [(c + shift, n, c) for c, n in R]
    P = {'fig8_late': [('late', rows('matched', 0, 99))], 'fig8_early': [('early', rows('matched', 0, 99))]}
    for i in range(4):
        lo, hi = FIG9_LIMITS[i], FIG9_LIMITS[i + 1]
        P['fig9_bin%d' % (i + 1)] = [('late', rows('blue', lo, hi)), ('early', rows('red', lo, hi))]
    P['full_isolated'] = [('late', rows('blue', 8.5, 11.0)), ('early', rows('red', 8.5, 11.0))]
    return P

def stack(rule, pop, case, cache):
    num = np.zeros(15); den = 0.0; mg_sum = 0.0; n_sum = 0.0
    for kind, rows in pop:
        for lm, n, lm0 in rows:
            key = (round(lm, 4), kind, case)
            if key not in cache: cache[key] = lens_model(rule, lm, kind == 'early', case)
            esd, mref, mg = cache[key]
            w = n * mg; num += w * esd; den += w; mg_sum += n * mg; n_sum += n
    return num / den, float(np.log10(mg_sum / n_sum))

CASES = ['primary', 'nfw_full_mass', 'truncated', 'early_half_disk', 'early_disk_size', 'vdw_sizes', 'early_gas', 'z015', 'z035',
         'ms_plus01', 'ms_minus01', 'c_dm14', 'scatter015', 'dist_plus01', 'dist_minus01']

def inner_mask(log_mgal):
    return GEDGES[:-1] >= GMSUN * 10 ** log_mgal / (R_INNER * MPC) ** 2

def predict(dd, bt, out):
    rule = Rule(dd, bt)
    fig7 = json.load(open(os.path.join(dd, 'Lensing', 'KiDS', 'brouwer2021_fig7_digitized.json')))
    res = dict(protocol_sha256=open(os.path.join(bt, 'PROTOCOL_v9.sha256')).read().split()[0],
               note_sha256=open(os.path.join(bt, 'PROTOCOL_v9_NOTE.sha256')).read().split()[0],
               gbar_edges=[float(g) for g in GEDGES], cases={})
    cache = {}
    for case in CASES:
        P = populations(fig7, case); C = {}
        for name, pop in P.items():
            esd, lmg = stack(rule, pop, case, cache)
            C[name] = dict(esd=[float('%.6g' % v) for v in esd], log_mean_mgal_digitized=round(lmg, 4))
        res['cases'][case] = C
    prim = res['cases']['primary']
    inner = {}
    for name in prim:
        lmg = FIG9_MGAL[int(name[-1]) - 1] if name.startswith('fig9') else prim[name]['log_mean_mgal_digitized']
        inner[name] = dict(log_mean_mgal_used=lmg, inner=[bool(b) for b in inner_mask(lmg)])
    res['inner_bins'] = inner
    res['fig9_mgal_check'] = {('fig9_bin%d' % (i + 1)): dict(B21=FIG9_MGAL[i], digitized=prim['fig9_bin%d' % (i + 1)]['log_mean_mgal_digitized']) for i in range(4)}
    # predicted split, primary model, inner and all bins (same model for both splits)
    E, L = np.array(prim['fig8_early']['esd']), np.array(prim['fig8_late']['esd']); k = np.array(inner['fig8_late']['inner'])
    res['predicted_split'] = dict(inner=round(float(np.log10(np.mean(E[k] / L[k]))), 4), all=round(float(np.log10(np.mean(E / L))), 4))
    res['predicted_split_by_case'] = {cs: dict(inner=round(float(np.log10(np.mean(np.array(res['cases'][cs]['fig8_early']['esd'])[k] / np.array(res['cases'][cs]['fig8_late']['esd'])[k]))), 4))
                                      for cs in CASES}
    # typical halo masses, primary
    res['typical_halos'] = {lab: dict(logMs=lm, log_Mref_late=round(float(np.log10(rule.mref(lm, False, 'primary')[0])), 3),
                                      log_Mref_early=round(float(np.log10(rule.mref(lm, True, 'primary')[0])), 3)) for lab, lm in (('10.3', 10.3), ('10.6', 10.6), ('10.9', 10.9))}
    json.dump(res, open(out + '_Predictions.json', 'w'), indent=1)
    print(json.dumps(dict(predicted_split=res['predicted_split'], by_case=res['predicted_split_by_case'], fig9_mgal_check=res['fig9_mgal_check'],
                          inner_counts={n: sum(v['inner']) for n, v in inner.items()}, typical=res['typical_halos']), indent=1))

def check_nfw():
    from colossus.cosmology import cosmology
    from colossus.halo import profile_nfw
    cosmology.addCosmology('A19planck', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624})
    cosmology.setCosmology('A19planck')
    m = 1.5e12; z = 0.25
    p = profile_nfw.NFWProfile(M=m * HR, c=c_rule(m, z), z=z, mdef='200c')
    R = np.array([0.03, 0.1, 0.3, 1.0, 3.0])
    mine = esd_nfw(R, m, z, c_rule, dark=False) / 1e12            # suns / pc^2
    col = p.deltaSigma(R * 1e3 * HR) * HR / 1e6                   # colossus: M_sun h / kpc^2 with R in kpc/h
    trunc = esd_nfw(R, m, z, c_rule, dark=False, truncated=True) / 1e12
    print('R', R.tolist()); print('mine', mine.round(3).tolist()); print('colossus', np.round(col, 3).tolist()); print('truncated', trunc.round(3).tolist())

def load_profile(path):
    d = np.loadtxt(path, comments='#')
    return dict(gbar=d[:, 0], esd=d[:, 1] / d[:, 4], err=d[:, 3] / d[:, 4], bias=d[:, 4])

def load_cov(path_json):
    """Compact covariance moved from the browser: each element already divided by its own bias column there,
    6 significant digits; blocks (a <= b) in order, upper triangle inside diagonal blocks; 'blocks_only' keeps
    only the diagonal blocks (the cross blocks are then left as zero and never used)."""
    J = json.load(open(path_json)); nb, nobs = J['nb'], len(J['obs_min']); n = nb * nobs
    C = np.zeros((n, n)); v = iter(J['v'])
    for a in range(nobs):
        for b in (range(a, a + 1) if J['blocks_only'] else range(a, nobs)):
            for i in range(nb):
                for j in (range(i, nb) if a == b else range(nb)):
                    x = next(v); C[a * nb + i, b * nb + j] = x; C[b * nb + j, a * nb + i] = x
    assert next(v, None) is None
    return C, J

def compare(dd, bt, out):
    pred = json.load(open(out + '_Predictions.json'))
    kd = os.path.join(dd, 'Lensing', 'KiDS')
    prof = lambda f: load_profile(os.path.join(kd, f))
    P = {'blue': prof('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'), 'red': prof('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'),
         'disc': prof('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'), 'bulge': prof('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt'),
         'full_isolated': prof('Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt')}
    for i in range(1, 5): P['fig9_bin%d' % i] = prof('Fig-9_RAR-KiDS-isolated_Massbin-%d.txt' % i)
    Cc, _ = load_cov(os.path.join(kd, 'Fig-8_Colorbins_covmatrix_compact.json'))
    Cs, _ = load_cov(os.path.join(kd, 'Fig-8_Sersicbins_covmatrix_compact.json'))
    C9, _ = load_cov(os.path.join(kd, 'Fig-9_Massbins_covmatrix_compact.json'))
    Cf, _ = load_cov(os.path.join(kd, 'Fig-4-5-C1_covmatrix_compact.json'))
    cov = {'blue': Cc[:15, :15], 'red': Cc[15:, 15:], 'disc': Cs[:15, :15], 'bulge': Cs[15:, 15:], 'full_isolated': Cf}
    for i in range(4): cov['fig9_bin%d' % (i + 1)] = C9[15 * i:15 * i + 15, 15 * i:15 * i + 15]
    model_of = {'blue': 'fig8_late', 'red': 'fig8_early', 'disc': 'fig8_late', 'bulge': 'fig8_early', 'full_isolated': 'full_isolated'}
    for i in range(1, 5): model_of['fig9_bin%d' % i] = 'fig9_bin%d' % i
    inner_of = lambda s: np.array(pred['inner_bins'][model_of[s]]['inner'])
    res = dict(protocol_sha256=pred['protocol_sha256'], predictions_sha256=open(out + '_Predictions.sha256').read().split()[0])
    # bins agree with the protocol's edges
    res['gbar_in_edges'] = {s: bool(np.all((P[s]['gbar'] > GEDGES[:-1]) & (P[s]['gbar'] < GEDGES[1:]))) for s in P}
    # (A) the split
    rng = np.random.default_rng(9)
    def split(dE, dL, Cfull, k):
        """Cfull is ordered [late, early] (B21's bin 1, then bin 2), so the draws' mean vector is [late, early] too.
        Post-run correction (turn 13): the code first run passed [early, late], which paired each profile with the
        other's covariance block (found by an independent review). The central value S_data and all four verdicts are
        unchanged; sigma_S for colour moves from 0.040 to 0.048 (inner bins). See blind_test/PROTOCOL_v9_CORRECTION.txt.
        Also returned: the delta-method sigma and the number of draws whose mean ratio is not positive (left out)."""
        S = lambda e, l: np.log10(np.mean(e[k] / l[k]))
        draws = rng.multivariate_normal(np.concatenate([dL, dE]), Cfull, size=20000)
        r = np.array([np.mean(d[15:][k] / d[:15][k]) for d in draws])
        ok = r > 0; Sd = np.log10(r[ok])
        idx = np.where(k)[0]; n = len(idx); rbar = np.mean(dE[k] / dL[k]); g = np.zeros(30)
        g[15 + idx] = 1 / (np.log(10) * rbar * n * dL[idx])                  # dS / d(early_j)
        g[idx] = -dE[idx] / (np.log(10) * rbar * n * dL[idx] ** 2)            # dS / d(late_j)
        return (float(S(dE, dL)), float(np.std(Sd)), float(0.5 * (np.percentile(Sd, 84) - np.percentile(Sd, 16))),
                float(np.sqrt(g @ Cfull @ g)), int((~ok).sum()))
    A = {}
    for lab, lo, hi, Cfull in (('colour', 'blue', 'red', Cc), ('sersic', 'disc', 'bulge', Cs)):
        k_in = inner_of(lo); A[lab] = {}
        for rng_lab, k in (('inner', k_in), ('all', np.ones(15, bool))):
            sd, sig, sig68, sig_dm, n_bad = split(P[hi]['esd'], P[lo]['esd'], Cfull, k)
            E, L = np.array(pred['cases']['primary']['fig8_early']['esd']), np.array(pred['cases']['primary']['fig8_late']['esd'])
            sm = float(np.log10(np.mean(E[k] / L[k])))
            T = max(2 * sig, 0.25 * sd)
            verdict = ('ACCOUNTS FOR' if abs(sm - sd) <= T else 'OVERPREDICTS' if sm > sd + T else 'ACCOUNTS FOR PART' if sm >= 0.25 * sd else 'DOES NOT')
            A[lab][rng_lab] = dict(S_data=round(sd, 4), sigma_S=round(sig, 4), sigma_S_68=round(sig68, 4), sigma_S_delta_method=round(sig_dm, 4),
                                   draws_left_out=n_bad, S_model=round(sm, 4),
                                   fraction=round(sm / sd, 3), T=round(T, 4), verdict=verdict, n_bins=int(k.sum()))
            A[lab][rng_lab]['S_model_by_case'] = {cs: round(float(np.log10(np.mean(np.array(pred['cases'][cs]['fig8_early']['esd'])[k] / np.array(pred['cases'][cs]['fig8_late']['esd'])[k]))), 4) for cs in pred['cases']}
        # chi2 of the measured difference, null and against the model difference (all bins and inner)
        dE, dL = P[hi]['esd'], P[lo]['esd']; Cd = Cfull[15:, 15:] + Cfull[:15, :15] - Cfull[15:, :15] - Cfull[:15, 15:]
        E, L = np.array(pred['cases']['primary']['fig8_early']['esd']), np.array(pred['cases']['primary']['fig8_late']['esd'])
        for rng_lab, k in (('inner', k_in), ('all', np.ones(15, bool))):
            Ci = np.linalg.inv(Cd[np.ix_(k, k)]); d = (dE - dL)[k]; m = (E - L)[k]
            A[lab][rng_lab]['chi2_no_difference'] = round(float(d @ Ci @ d), 2)
            A[lab][rng_lab]['chi2_vs_model_difference'] = round(float((d - m) @ Ci @ (d - m)), 2)
    res['A_split'] = A
    # (B) amplitudes
    Bres = {}
    for s in P:
        k = inner_of(s); d = P[s]['esd']; Cs_ = cov[s]
        Bres[s] = {}
        for cs in pred['cases']:
            m = np.array(pred['cases'][cs][model_of[s]]['esd'])
            Ci = np.linalg.inv(Cs_[np.ix_(k, k)]); mk, dk = m[k], d[k]
            a = float(mk @ Ci @ dk / (mk @ Ci @ mk)); sa = float((mk @ Ci @ mk) ** -0.5)
            la, sla = np.log10(a), sa / (a * np.log(10))
            ok = abs(la) <= 0.1 or (sla > 0.05 and abs(la) <= 2 * sla)
            ent = dict(log10_a=round(float(la), 4), err=round(float(sla), 4), consistent=bool(ok),
                       # added after the run (turn 13; not a protocol verdict): the radial-profile test, chi2 of the inner bins
                       # after the fitted amplitude a, with n_inner - 1 degrees of freedom
                       chi2_inner_after_amplitude=round(float((dk - a * mk) @ Ci @ (dk - a * mk)), 2))
            if cs != 'primary': ent['chi2_inner'] = round(float((dk - mk) @ Ci @ (dk - mk)), 2)
            if cs == 'primary':
                Cia = np.linalg.inv(Cs_)
                ent.update(chi2_inner=round(float((dk - mk) @ Ci @ (dk - mk)), 2), n_inner=int(k.sum()),
                           chi2_all=round(float((d - m) @ Cia @ (d - m)), 2), n_all=15,
                           log_ratio_by_bin=[round(float(np.log10(d[j] / m[j])), 3) if d[j] > 0 else None for j in range(15)],
                           data_esd=[float('%.5g' % v) for v in d], data_err=[float('%.4g' % v) for v in P[s]['err']], model_esd=[float('%.5g' % v) for v in m],
                           gbar=[float('%.4g' % v) for v in P[s]['gbar']], inner=[bool(v) for v in k])
            Bres[s][cs] = ent
    res['B_amplitudes'] = Bres
    res['B_late_types_consistent'] = dict(blue=Bres['blue']['primary']['consistent'], disc=Bres['disc']['primary']['consistent'])
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
    summ = dict(A={l: {r: {k: v for k, v in A[l][r].items() if k != 'S_model_by_case'} for r in A[l]} for l in A},
                B={s: {k: v for k, v in Bres[s]['primary'].items() if k in ('log10_a', 'err', 'consistent', 'chi2_inner', 'chi2_inner_after_amplitude', 'n_inner', 'chi2_all')} for s in Bres},
                B_truncated={s: {k: v for k, v in Bres[s]['truncated'].items() if k in ('log10_a', 'err', 'consistent', 'chi2_inner', 'chi2_inner_after_amplitude')} for s in Bres},
                gbar_in_edges=res['gbar_in_edges'])
    print(json.dumps(summ, indent=1))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'predict': predict(*sys.argv[2:5])
    elif mode == 'compare': compare(*sys.argv[2:5])
    elif mode == 'check_nfw': check_nfw()
