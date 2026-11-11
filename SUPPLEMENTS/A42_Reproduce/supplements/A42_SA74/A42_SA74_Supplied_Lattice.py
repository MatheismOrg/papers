#!/usr/bin/env python3
"""A42 SA74 T5: the known laws on a supplied four-dimensional lattice.
1. Bits: the A38 symbolic statement of the Standard Model with general relativity, plus a lattice form and two completions
   (declared data, or the SMASH mechanisms), compressed exactly as A38 does (raw LZMA2, preset 9 extreme); plus the
   fine-tuning a physical lattice needs for Lambda and mu^2, estimated from the statement's own inputs.
2. Chirality on a lattice: naive, Wilson and overlap Dirac operators on a 2D U(1) lattice with topological charge Q:
   doublers, the Ginsparg-Wilson relation, and the index of the overlap operator against Q.
3. Anomalies of the statement's fermions, per family: every gauge anomaly, Witten's SU(2) anomaly, and the mod-16 anomaly
   of X = 5(B-L) - 4Y, for 16 and for 15 Weyl fermions.
4. Lorentz symmetry: the photon's dispersion on a 4D hypercubic lattice, and the largest spacing LHAASO allows.
Python 3 with NumPy. Usage: python3 A42_SA74_Supplied_Lattice.py SPEC.txt OUT.json"""
import sys, json, lzma, math, time
from fractions import Fraction as F
import numpy as np

SPEC = sys.argv[1] if len(sys.argv) > 1 else 'A38_Symbolic_Benchmark_Spec.txt'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'A42_SA74_Supplied_Lattice_Results.json'
res = {'started': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}

# ---------------------------------------------------------------- 1. bits
def lz(b): return 8 * len(lzma.compress(b, format=lzma.FORMAT_RAW, filters=[{'id': lzma.FILTER_LZMA2, 'preset': 9 | lzma.PRESET_EXTREME}]))
spec = open(SPEC, 'rb').read()
LATTICE = b"""
Lattice form (the 4-dimensional lattice supplied as data):
Lattice: periodic hypercubic in 4 dimensions, spacing a = sqrt(GN) (the Planck length), physical rather than a regulator.
Gauge fields: link variables U_m(x) = exp(i a k X_m(x)) for each gauge field X with coupling k, Wilson plaquette action.
Metric: the lengths of the edges of the Freudenthal triangulation of each cell; the term (R - 2 Lambda)/(16 pi GN) by the Regge action.
Fermions: overlap operator D = (1 + gamma5 sign(gamma5 (a D_W - 1)))/a with D_W the Wilson operator;
left-handed fields f = (1 - gamma5 (1 - a D)) f/2 and fbar = fbar (1 + gamma5)/2 (Luscher).
Amplitudes: the Euclidean lattice path integral; one time step is the transfer matrix along a lattice direction.
Bare Lambda and mu^2 at the lattice scale are tuned so that the measured values above result.
"""
DECLARED = b"""
Completion by declared data:
Dark matter: a real scalar S (1, 1, 0) with action - (1/2) (d_m S)(d^m S) - (1/2) mS^2 S^2, coupled only through the metric; mS = 100 (declared).
Initial state: a spatially flat expanding universe at temperature 1e10 in thermal equilibrium, with adiabatic Gaussian perturbations
of amplitude A_s = 2.100e-9 and tilt n_s = 0.9649 at k = 0.05/Mpc, baryon-to-photon ratio 6.12e-10, dark-to-baryon mass ratio 5.36.
"""
SMASH = b"""
Completion by mechanisms (SMASH):
Add a complex scalar s (1, 1, 0) and Weyl fermions q (3, 1, -1/3), qc (3bar, 1, 1/3).
Global U(1) charges: s 1, q -1/2, qc -1/2, L_n 1/2, ec_n -1/2, nc_n -1/2, all others 0.
Add - ( (1/2) yN_jk s nc_j nc_k + yq s q qc + h.c. ), replace the Higgs potential by
lam (H^dagger H - v^2/2)^2 + ls (s* s - vs^2/2)^2 + 2 lhs (H^dagger H - v^2/2)(s* s - vs^2/2),
and the gravity term by (1/(16 pi GN) + xs s* s) R - Lambda/(8 pi GN).
Inputs: vs = 1.0e11, ls = 1.0e-10, lhs = 1.0e-6, xs = 0.10, yq = 0.010, yN = diag(1.00e-4, 1.04e-4, 2.00e-4);
yn = (sqrt(2)/v) U diag(sqrt(m_1 M_1), sqrt(m_2 M_2), sqrt(m_3 M_3)) with M = yN vs/sqrt(2); th dropped (the phase of s absorbs it).
Initial state: inflation driven by |s| in a generic patch.
"""
bits = {'A38 statement': lz(spec), 'plus lattice form': lz(spec + LATTICE),
        'plus declared completion': lz(spec + DECLARED), 'plus SMASH completion': lz(spec + SMASH),
        'plus lattice and declared': lz(spec + LATTICE + DECLARED), 'plus lattice and SMASH': lz(spec + LATTICE + SMASH)}
# fine-tuning on a physical lattice of spacing a = sqrt(GN), from the statement's own inputs
GN, Lam, mH = 6.7088e-39, 4.25e-84, 125.20
MP = GN ** -0.5
rho_lat = Lam / (8 * math.pi * MP ** 2)          # vacuum energy density in lattice units a^-4 = MP^4
higgs = mH ** 2 * 16 * math.pi ** 2 / MP ** 2      # m_H^2 against MP^2/(16 pi^2)
tune = {'Planck mass GeV': MP, 'vacuum energy in lattice units': rho_lat, 'Lambda extra bits': math.log2(1 / rho_lat),
        'Higgs ratio': higgs, 'Higgs extra bits': math.log2(1 / higgs)}
tune['total extra bits'] = tune['Lambda extra bits'] + tune['Higgs extra bits']
res['bits'] = bits; res['tuning'] = tune
res['bits_cap'] = {'lattice + declared + tuning': bits['plus lattice and declared'] + tune['total extra bits'],
                   'lattice + SMASH + tuning': bits['plus lattice and SMASH'] + tune['total extra bits']}
res['texts'] = {'lattice': LATTICE.decode(), 'declared': DECLARED.decode(), 'smash': SMASH.decode()}

# ---------------------------------------------------------------- 2. chirality on a 2D U(1) lattice
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1., -1.]).astype(complex)
GAM = [s1, s2]; G5 = s3

def links(L, Q, eps=0.0, seed=0):
    x0 = np.arange(L)[:, None] * np.ones((1, L))
    th1 = 2 * np.pi * Q * x0 / L ** 2
    th0 = np.zeros((L, L)); th0[L - 1, :] = -2 * np.pi * Q * np.arange(L) / L
    if eps:
        rng = np.random.default_rng(seed)
        th0 = th0 + eps * rng.standard_normal((L, L)); th1 = th1 + eps * rng.standard_normal((L, L))
    return [np.exp(1j * th0), np.exp(1j * th1)]

def charge(U):
    P = U[0] * np.roll(U[1], -1, axis=0) * np.conj(np.roll(U[0], -1, axis=1)) * np.conj(U[1])
    return float(np.angle(P).sum() / (2 * np.pi)), float(np.max(np.abs(np.angle(P))))

def dirac(U, r):
    L = U[0].shape[0]; n = L * L; D = np.zeros((2 * n, 2 * n), complex); I2 = np.eye(2)
    for a in range(L):
        for b in range(L):
            i = a * L + b
            D[2*i:2*i+2, 2*i:2*i+2] += 2 * r * I2
            for mu, (da, db) in enumerate([(1, 0), (0, 1)]):
                j = ((a + da) % L) * L + (b + db) % L; k = ((a - da) % L) * L + (b - db) % L
                D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (r * I2 - GAM[mu]) * U[mu][a, b]
                D[2*i:2*i+2, 2*k:2*k+2] += -0.5 * (r * I2 + GAM[mu]) * np.conj(U[mu][(a - da) % L, (b - db) % L])
    return D

def g5(n2): return np.kron(np.eye(n2 // 2), G5)

def overlap(DW):
    n2 = DW.shape[0]; g = g5(n2)
    H = g @ (DW - np.eye(n2)); H = (H + H.conj().T) / 2
    w, V = np.linalg.eigh(H)
    return np.eye(n2) + g @ ((V * np.sign(w)) @ V.conj().T), float(np.min(np.abs(w)))

def zero_modes(D, tol=1e-8):
    g = g5(D.shape[0]); u, s, vh = np.linalg.svd(D)
    V0 = vh[s < tol].conj().T
    if V0.shape[1] == 0: return 0, 0, float(np.min(s))
    c = np.linalg.eigvalsh(V0.conj().T @ g @ V0)
    return int(np.sum(c > 0.5)), int(np.sum(c < -0.5)), float(np.min(s))

rows = []
for L, Q, eps, seed in [(16, 0, 0, 0), (16, 1, 0, 0), (16, -1, 0, 0), (16, 2, 0, 0), (16, -2, 0, 0), (16, 3, 0, 0),
                        (16, 1, 0.3, 1), (16, 2, 0.3, 2), (16, -1, 0.5, 3), (24, 1, 0, 0), (24, 3, 0, 0), (24, -2, 0.3, 4)]:
    U = links(L, Q, eps, seed); q, pmax = charge(U); n2 = 2 * L * L; g = g5(n2)
    DW = dirac(U, 1.0); Dov, hgap = overlap(DW)
    gw = float(np.linalg.norm(g @ Dov + Dov @ g - Dov @ g @ Dov) / np.linalg.norm(Dov))
    nplus, nminus, smin = zero_modes(Dov)
    ev = np.linalg.eigvals(Dov); n_two = int(np.sum(np.abs(ev - 2) < 1e-8))
    evW = np.linalg.eigvals(DW); wmin = float(np.min(np.abs(evW)))
    rows.append({'L': L, 'Q_set': Q, 'noise': eps, 'Q_measured': round(q, 6), 'max_plaquette_angle': round(pmax, 3),
                 'overlap_n_plus': nplus, 'overlap_n_minus': nminus, 'index_n_minus_minus_n_plus': nminus - nplus,
                 'overlap_modes_at_2': n_two, 'GW_violation_rel': gw, 'H_W_gap': round(hgap, 4),
                 'wilson_smallest_abs_eig': round(wmin, 4)})
    print(rows[-1], flush=True)
# doublers in the free field
U0 = links(16, 0)
naive = dirac(U0, 0.0); wil = dirac(U0, 1.0); ov, _ = overlap(wil)
cnt = lambda D: int(np.sum(np.abs(np.linalg.eigvals(D)) < 1e-8))
res['chirality'] = {'rows': rows, 'free_zero_modes': {'naive': cnt(naive), 'wilson': cnt(wil), 'overlap': cnt(ov), 'continuum': 2},
                    'wilson_GW_violation_rel_free': float(np.linalg.norm(g5(512) @ wil + wil @ g5(512) - wil @ g5(512) @ wil) / np.linalg.norm(wil))}
print(res['chirality']['free_zero_modes'], flush=True)

# ---------------------------------------------------------------- 3. anomalies (left-handed Weyl fields; (SU3 dim, SU2 dim, Y, B-L))
fam = {'Q': (3, 2, F(1, 6), F(1, 3)), 'uc': (-3, 1, F(-2, 3), F(-1, 3)), 'dc': (-3, 1, F(1, 3), F(-1, 3)),
       'L': (1, 2, F(-1, 2), F(-1)), 'ec': (1, 1, F(1), F(1)), 'nc': (1, 1, F(0), F(1))}
def anomalies(fields):
    A = {'grav^2 Y': F(0), 'Y^3': F(0), 'SU2^2 Y': F(0), 'SU3^2 Y': F(0), 'SU3^3': 0, 'SU2 doublets': 0, 'mod16': 0}
    for c, w, Y, BL in fields.values():
        nc = abs(c); mult = nc * w
        A['grav^2 Y'] += mult * Y; A['Y^3'] += mult * Y ** 3
        if w == 2: A['SU2^2 Y'] += nc * Y; A['SU2 doublets'] += nc
        if nc == 3: A['SU3^2 Y'] += w * Y; A['SU3^3'] += w * (1 if c > 0 else -1)
        X = 5 * BL - 4 * Y; assert X.denominator == 1 and X.numerator % 2 == 1
        A['mod16'] += mult * (1 if X.numerator % 4 == 1 else -1)
    return {k: str(v) for k, v in A.items()}
f15 = {k: v for k, v in fam.items() if k != 'nc'}
res['anomalies'] = {'16 per family, one family': anomalies(fam), '15 per family, one family': anomalies(f15),
                    'mod16 for 3 families': {'16 per family': (3 * int(anomalies(fam)['mod16'])) % 16, '15 per family': (3 * int(anomalies(f15)['mod16'])) % 16}}
smash_extra = {'q': (3, 1, F(-1, 3), F(1, 3)), 'qc': (-3, 1, F(1, 3), F(-1, 3))}
res['anomalies']['SMASH extra quarks'] = anomalies(smash_extra)
print(res['anomalies'], flush=True)

# ---------------------------------------------------------------- 4. Lorentz on the 4D hypercubic lattice (a = 1)
def E_of_p(p, n):
    n = np.asarray(n, float) / np.linalg.norm(n); s = np.sum(np.sin(p * n / 2) ** 2)
    return 2 * np.arcsinh(np.sqrt(s))
lor = {}
for name, n in [('axis', (1, 0, 0)), ('face diagonal', (1, 1, 0)), ('body diagonal', (1, 1, 1))]:
    p = 1e-2; c = (p - E_of_p(p, n)) / p ** 3       # E = p - c a^2 p^3
    EQG2 = 1 / math.sqrt(2 * c)                     # E^2 = p^2 (1 - (p/EQG2)^2) with EQG2 in units of 1/a
    lor[name] = {'c': c, 'EQG2_times_a': EQG2}
LHAASO = 6.9e11   # GeV, subluminal quadratic limit, LHAASO 2024
hbarc, lP = 1.97327e-16, 1.616255e-35  # GeV m, m
for k in lor:
    amax = lor[k]['EQG2_times_a'] / LHAASO
    lor[k]['a_max_m'] = amax * hbarc; lor[k]['a_max_planck'] = amax * hbarc / lP
res['lorentz'] = lor
print(json.dumps({'bits': bits, 'tuning': tune, 'cap': res['bits_cap'], 'lorentz': lor}, indent=1), flush=True)
res['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
json.dump(res, open(OUT, 'w'), indent=1, default=str)
