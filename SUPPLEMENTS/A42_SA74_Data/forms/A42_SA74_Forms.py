#!/usr/bin/env python3
"""A42 SA74 T7: fermions as differential forms (Kahler-Dirac fields) on causal triangulations and on the supplied lattice.
A Kahler-Dirac field is a sum of 0- to d-forms; on a simplicial complex its operator is K = d - d^dagger on cochains, so it
needs no local frames and lives on any triangulation. K^2 = -(Hodge Laplacian), so its spectrum is the union of the Hodge
Laplacians' spectra on k-forms.
1. Two-dimensional causal triangulations of a torus: slice lengths sampled by Metropolis with the strip counts
   C(l + l' - 1, l - 1) of 2D causal triangulations (Ambjorn and Loll, 1998), strips drawn uniformly given the lengths,
   equilateral triangles of side 1 with circumcentric (DEC) Hodge stars. Reports exact zero modes against the Betti
   numbers (1, 2, 1) and the low spectrum against Weyl's law, N(lam) = 4 A lam / (4 pi): doubled modes would push the
   ratio far above 1. A regular triangular torus is the control.
2. The supplied four-dimensional lattice: the Freudenthal triangulation of a periodic 3^4 hypercubic lattice, with random
   positive Hodge weights standing in for a bent geometry. Reports the zero modes of each Hodge Laplacian against the
   Betti numbers of the 4-torus (1, 4, 6, 4, 1).
Python 3 with NumPy and SciPy. Usage: python3 A42_SA74_Forms.py OUT.json"""
import sys, json, math, time, itertools
import numpy as np, scipy.sparse as sp
from scipy.special import gammaln
OUT = sys.argv[1] if len(sys.argv) > 1 else 'A42_SA74_Forms_Results.json'
rng = np.random.default_rng(7)
res = {'started': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}

def logC(n, k): return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

def complex_from_top(tops, top_dim):
    S = [dict() for _ in range(top_dim + 1)]
    for s in tops:
        s = tuple(sorted(s))
        for k in range(1, len(s) + 1):
            for f in itertools.combinations(s, k):
                if f not in S[k - 1]: S[k - 1][f] = len(S[k - 1])
    D = []
    for k in range(top_dim):
        r, c, v = [], [], []
        for s, j in S[k + 1].items():
            for i in range(len(s)):
                r.append(j); c.append(S[k][s[:i] + s[i + 1:]]); v.append((-1) ** i)
        D.append(sp.csr_matrix((v, (r, c)), shape=(len(S[k + 1]), len(S[k]))).toarray())
    return S, D

def hodge_spectra(D, W):
    # symmetric form of the Hodge Laplacian on k-forms with diagonal inner products W[k]
    A = [np.sqrt(W[k + 1])[:, None] * D[k] / np.sqrt(W[k])[None, :] for k in range(len(D))]
    spec = []
    for k in range(len(W)):
        M = np.zeros((len(W[k]), len(W[k])))
        if k > 0: M += A[k - 1] @ A[k - 1].T
        if k < len(D): M += A[k].T @ A[k]
        spec.append(np.linalg.eigvalsh(M))
    return spec

def zeros(spec):
    top = max(float(s.max()) for s in spec)
    return [int(np.sum(s < 1e-9 * top)) for s in spec]

# ---------------------------------------------------------------- 1. 2D causal triangulations
def sample_slices(T, N, sweeps, lmin=3):
    L = np.full(T, N // (2 * T))
    lw = lambda a, b: logC(a + b - 1, a - 1)
    for _ in range(sweeps * T):
        s, t = rng.integers(T, size=2)
        if s == t or L[s] - 1 < lmin: continue
        pairs = {(s - 1) % T, s, (t - 1) % T, t}
        old = sum(lw(L[i], L[(i + 1) % T]) for i in pairs)
        L[s] -= 1; L[t] += 1
        if math.log(rng.random()) > sum(lw(L[i], L[(i + 1) % T]) for i in pairs) - old: L[s] += 1; L[t] -= 1
    return L

def causal_torus(L, regular=False):
    T = len(L); vid = {}; n = 0
    for t in range(T):
        for i in range(L[t]): vid[(t, i)] = n; n += 1
    tris = []
    for t in range(T):
        l, lp, t1 = L[t], L[(t + 1) % T], (t + 1) % T
        if regular: word = [0, 1] * l; off = (t * 0)
        else:
            word = [0] * l + [1] * lp; rng.shuffle(word); off = int(rng.integers(lp))
        a, b = 0, off
        for w in word:
            if w == 0: tris.append((vid[(t, a % l)], vid[(t, (a + 1) % l)], vid[(t1, b % lp)])); a += 1
            else: tris.append((vid[(t1, b % lp)], vid[(t1, (b + 1) % lp)], vid[(t, a % l)])); b += 1
    return n, tris

def run2d(L, regular=False):
    n, tris = causal_torus(L, regular)
    S, D = complex_from_top(tris, 2)
    V, E, F = (len(S[k]) for k in range(3))
    assert np.allclose(D[1] @ D[0], 0) and np.all(np.abs(D[1]).sum(axis=0) == 2), 'not a closed surface'
    deg = np.zeros(V)
    for tr in S[2]:
        for v in tr: deg[S[0][(v,)]] += 1
    area = math.sqrt(3) / 4
    W = [deg * area / 3, np.full(E, 1 / math.sqrt(3)), np.full(F, 1 / area)]
    spec = hodge_spectra(D, W); A = F * area
    allev = np.concatenate(spec)
    weyl = {}
    for lam in (0.1, 0.25, 0.5):
        nz = int(np.sum((allev > 1e-9 * allev.max()) & (allev < lam)))
        weyl[str(lam)] = round(nz / (4 * A * lam / (4 * math.pi)), 3)
    return {'triangles': F, 'slices': len(L), 'slice_min': int(min(L)), 'slice_max': int(max(L)),
            'euler': V - E + F, 'zero_modes_by_degree': zeros(spec), 'zero_modes_total': sum(zeros(spec)),
            'low_modes_over_weyl': weyl}

rows2d = []
for T, N in [(16, 512), (24, 1024), (32, 2048)]:
    for rep in range(2):
        L = sample_slices(T, N, sweeps=400)
        r = run2d(L); r['kind'] = 'causal, sampled'; rows2d.append(r); print(r, flush=True)
    r = run2d(np.full(T, N // (2 * T)), regular=True); r['kind'] = 'regular control'; rows2d.append(r); print(r, flush=True)
res['two_d'] = rows2d
json.dump(res, open(OUT, 'w'), indent=1)

# ---------------------------------------------------------------- 2. the supplied 4D lattice, randomly bent
def freudenthal(Lsz, d=4):
    idx = lambda x: int(np.ravel_multi_index(tuple(np.mod(x, Lsz)), (Lsz,) * d))
    tops = []
    for x in itertools.product(range(Lsz), repeat=d):
        for perm in itertools.permutations(range(d)):
            y = np.array(x); vs = [idx(y)]
            for ax in perm: y = y.copy(); y[ax] += 1; vs.append(idx(y))
            tops.append(vs)
    return tops

import scipy.sparse.linalg as sla
def hodge_low(Ds, W, k_eig=24):
    A = [sp.diags(np.sqrt(W[k + 1])) @ Ds[k] @ sp.diags(1 / np.sqrt(W[k])) for k in range(len(Ds))]
    out = []
    for k in range(len(W)):
        M = sp.csr_matrix((len(W[k]), len(W[k])))
        if k > 0: M = M + A[k - 1] @ A[k - 1].T
        if k < len(Ds): M = M + A[k].T @ A[k]
        ev = sla.eigsh(M.tocsc(), k=min(k_eig, M.shape[0] - 2), sigma=-1e-3, which='LM', return_eigenvectors=False)
        out.append(np.sort(ev))
    return out
S4, D4 = complex_from_top(freudenthal(3), 4)
Ds = [sp.csr_matrix(d) for d in D4]; del D4
counts = [len(S4[k]) for k in range(5)]
assert all(abs(Ds[k + 1] @ Ds[k]).sum() == 0 for k in range(3))
rows4d = []
for sigma in (0.0, 1.0):
    W = [np.exp(sigma * rng.standard_normal(c)) for c in counts]
    low = hodge_low(Ds, W)
    z = [int(np.sum(ev < 1e-8)) for ev in low]
    gap = [round(float(ev[ev >= 1e-8].min()), 4) for ev in low]
    rows4d.append({'lattice': '3^4 Freudenthal', 'cells': counts, 'euler': counts[0] - counts[1] + counts[2] - counts[3] + counts[4],
                   'weight_spread': sigma, 'zero_modes_by_degree': z, 'zero_modes_total': sum(z), 'smallest_nonzero_by_degree': gap})
    print(rows4d[-1], flush=True)
res['four_d'] = rows4d
res['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
json.dump(res, open(OUT, 'w'), indent=1)
