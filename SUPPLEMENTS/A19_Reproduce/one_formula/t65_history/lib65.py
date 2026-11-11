"""SA73 T65 library: the formula's halos with four readings of the core, every coefficient frozen at its all-data value.
  N  no core (T48).
  G  the final formula's core, gated by cold gas and diffuseness (T61: t = 0.1501, kappa = 0.6229): none without gas.
  H1 the published history recipe, coreNFW of Read, Agertz and Collins (2016, MNRAS 459, 2573) exactly:
     M(r) = M_NFW(r) tanh(r / r_c)^n, r_c = 1.75 R_1/2, n = tanh(0.04 t_SF / t_dyn),
     t_dyn = 2 pi sqrt(r_s^3 / (G M_NFW(r_s))), in the formula's halo, no mass term. t_SF = 13.8 Gyr minus the lookback
     time at which star formation ended (13.8 for a galaxy still forming stars).
  H2 the formula's own core with the history gate: w = s L, s = t_SF / 13.8, with t and kappa from all 149 benchmark disks
     under the gate of harness_T65.py (res65_alldata_fit.json).
  P  self-interaction of 1 cm^2/g at low speed (T64's choice), the Jeans model of sidm46.py, baryons as in t64/sidm64.py.
Dark-matter mass and density of the cored profile follow from the mass formula by differentiation:
  rho(r) = f^n rho_NFW + n f^(n-1) (1 - f^2) M_NFW / (4 pi r^2 r_c), f = tanh(r / r_c)."""
import json, os, numpy as np
import harness_T65 as M5
import harness_T64 as M
import sidm46 as SD
H = M.H; B = M.B; M61 = M.M61; S = M.S; K = M.K; G43 = M.G43
ABC = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339]); T_G, K_G = 0.1501, 0.6229; ETA = 1.75; KR = 0.04; TU = 13.8
GYR = 0.977792            # 1 kpc / (km/s) in Gyr
def halo(a, mult=1.):
    B.FROZEN['abc'] = ABC; thr = S.thresholds(a); m = G43.mass(TH, a['ms'], a['mp'], a['re'], a['z'], thr) * mult; c = K.c_dm14(m, a['z']) * a.get('fac', 1.)
    return m, c, (H.G * m / (100 * H.hz(a['z']) ** 2)) ** (1 / 3)
def nfw_m(h, r): m, c, rv = h; return (1 - H.FB) * m * H.ff(r * c / rv) / H.ff(c)
def nfw_rho(h, r): m, c, rv = h; rs = rv / c; x = r / rs; return (1 - H.FB) * m / (4 * np.pi * rs ** 3 * H.ff(c)) / (x * (1 + x) ** 2)
def core_m(h, r, n, rc): return nfw_m(h, r) * np.tanh(r / rc) ** n
def core_rho(h, r, n, rc):
    f = np.tanh(r / rc); return f ** n * nfw_rho(h, r) + n * f ** (n - 1) * (1 - f ** 2) * nfw_m(h, r) / (4 * np.pi * r ** 2 * rc)
def tdyn(h): m, c, rv = h; rs = rv / c; return 2 * np.pi * np.sqrt(rs ** 3 / (H.G * m * H.ff(1.) / H.ff(c))) * GYR
def n_read(h, tsf): return np.tanh(KR * tsf / tdyn(h))
def w_gas(a): return M61.gl(a['ms'], a['mp'], a['re'], S.thresholds(a))
def w_hist(a, tsf): return np.where(S.thresholds(a) > 0, 0., tsf / TU * M5.diffuse(a['ms'], a['re']))
def reading(a, r, kind, tsf=None, tk=None):
    """dark mass inside r and dark density at r (Msun, Msun/kpc^3) for a one-galaxy array a (first row's inputs)"""
    a1 = {k: np.asarray(v)[:1] for k, v in a.items()}
    if kind == 'N': h = halo(a1); return nfw_m(h, r), nfw_rho(h, r)
    if kind == 'G': w = w_gas(a1); t, kap = T_G, K_G
    elif kind == 'H2': w = w_hist(a1, tsf); t, kap = tk
    elif kind == 'H1': h = halo(a1); n = n_read(h, tsf); rc = ETA * a1['re']; return core_m(h, r, n, rc), core_rho(h, r, n, rc)
    else: raise ValueError(kind)
    h = halo(a1, 10 ** (t * w)); n = np.tanh(kap * w); rc = ETA * a1['re']
    return (core_m(h, r, n, rc), core_rho(h, r, n, rc)) if float(w[0]) > 0 else (nfw_m(h, r), nfw_rho(h, r))
def mb_from_points(r, mb):
    o = np.argsort(r); r, mb = r[o], np.maximum.accumulate(mb[o]); lr, lm = np.log(r), np.log(mb)
    def f(x):
        x = np.asarray(x, float); out = np.exp(np.interp(np.log(np.maximum(x, 1e-12)), lr, lm)); return np.where(x < r[0], mb[0] * (x / r[0]) ** 2, out)
    return f
def baryons(a, dsph):
    r = a['r']; mb = a['gb'] / H.ACC * r ** 2 / H.G
    if dsph:
        ap = a['re'][0]; norm = mb[0] / (r[0] ** 3 / (r[0] ** 2 + ap ** 2) ** 1.5); return lambda x: norm * np.asarray(x, float) ** 3 / (np.asarray(x, float) ** 2 + ap ** 2) ** 1.5
    return mb_from_points(r, mb)
def sidm(a, r, dsph, s0=1.0):
    a1 = {k: np.asarray(v)[:1] for k, v in a.items()}; m, c, rv = (float(x[0]) for x in halo(a1))
    sol = SD.solve((1 - H.FB) * m, rv, c, baryons(a, dsph), lambda v: SD.xs_vdep(v, s0), 10.); r = np.atleast_1d(np.asarray(r, float)); e = 0.02
    mm = sol['mdm'](r); rho = (sol['mdm'](r * (1 + e)) - sol['mdm'](r * (1 - e))) / (4 * np.pi * r ** 2 * 2 * e * r)
    return mm, rho, sol['flag'], sol['r1']
# ---------------------------------------------------------------- name maps and inputs
DSPH = {'UMi': 'ursa_minor_1', 'Draco': 'draco_1', 'Sculptor': 'sculptor_1', 'Sextans': 'sextans_1', 'Leo I': 'leo_1', 'Leo II': 'leo_2', 'Carina': 'carina_1', 'Fornax': 'fornax_1'}
DIRR = {'WLM': ('sparc', 'UGCA444'), 'DDO 52': ('lt', 'DDO 52'), 'DDO 87': ('sparc', 'UGC05918'), 'DDO 154': ('sparc', 'DDO154'), 'Aquarius': ('lt', 'DDO 210'), 'NGC 2366': ('sparc', 'NGC2366'), 'CVnIdwA': ('lt', 'CVnIdwA'), 'DDO 168': ('sparc', 'DDO168')}
WEISZ = {'Andromeda I': 'andromeda_01', 'Andromeda II': 'andromeda_02', 'Andromeda III': 'andromeda_03', 'Andromeda V': 'andromeda_05', 'Andromeda VI': 'andromeda_06', 'Andromeda XI': 'andromeda_11',
         'Andromeda XIII': 'andromeda_13', 'Carina': 'carina_1', 'Canes Venatici I': 'canes_venatici_1', 'Canes Venatici II': 'canes_venatici_2', 'Draco': 'draco_1', 'Fornax': 'fornax_1', 'Hercules': 'hercules_1',
         'Leo I': 'leo_1', 'Leo II': 'leo_2', 'Leo IV': 'leo_4', 'Phoenix': 'phoenix_1', 'Sculptor': 'sculptor_1', 'Tucana': 'tucana_1', 'Ursa Minor': 'ursa_minor_1'}
def lt_rows(t59, name):
    import model59 as L
    g = json.load(open(os.path.join(t59, 'lt59_inputs.json')))[name]; KL, KRR = 1.172, 1.424
    R = np.array(g['R']); V = np.array(g['Vc']); Ms = 0.5 * L.l36(g['f36'][0], g['D']) / KL; Re = 1.678 * g['read']['Rstar'] * g['D'] / g['read']['D'] / KRR
    vg2, mhi, mcut = L.gas(R, g['Sd'], 'exp'); vb2 = vg2 + L.stars_v2(R, Ms, Re / 1.678); k = (R > 0) & (V > 0) & (vb2 > 0); Mp = Ms + 1.33 * mhi
    return H.arr([dict(name=name, kind='disk', r=x, gb=y, go=w, ms=Ms, mp=Mp, re=Re, z=0.) for x, y, w in zip(R[k], vb2[k] / R[k] * H.ACC, V[k] ** 2 / R[k] * H.ACC)])
def galaxy(t59, label):
    """A19's own inputs for one of Read et al.'s sixteen dwarfs"""
    if label in DSPH:
        d = H.arr(H.dwarf_rows()); k = d['name'] == DSPH[label]; assert k.sum() == 1; return {key: np.asarray(v)[k] for key, v in d.items()}, True
    src, nm = DIRR[label]
    if src == 'sparc': k = M.ALL['name'] == nm; assert k.sum() > 3; return {key: np.asarray(v)[k] for key, v in M.ALL.items()}, False
    return lt_rows(t59, nm), False
