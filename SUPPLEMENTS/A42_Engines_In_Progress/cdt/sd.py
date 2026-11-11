"""SA74: spectral dimension of a cdt4 configuration by diffusion on the dual graph (4-simplices as nodes, 5 neighbours each).
P(sigma) = probability that a walk started at simplex i is back at i after sigma steps, averaged over starts in the slices that hold the
blob; D_s(sigma) = -2 d ln P / d ln sigma (even sigma). Published for (k0, Delta) = (2.2, 0.6), N4 = 181,000 (Ambjorn, Jurkiewicz, Loll 2005):
D_s = 4.02 - 119/(54 + sigma) on sigma = 40 to 400. Usage: sd.py config [nstart] [sigma_max] [seed]"""
import sys, numpy as np, scipy.sparse as sp
def load(fn):
    f = open(fn); h = f.readline().split(); T, n0, n4 = int(h[1]), int(h[2]), int(h[3])
    vt = np.array([int(f.readline()) for _ in range(n0)]); S = np.loadtxt(f, dtype=np.int64); return T, vt, S[:, :5], S[:, 5]
def dual(V):
    n = len(V); Vs = np.sort(V, 1); faces = np.concatenate([np.delete(Vs, i, 1) for i in range(5)]); owner = np.tile(np.arange(n), 5)
    o = np.lexsort(faces.T[::-1]); f = faces[o]; ow = owner[o]; same = (f[1:] == f[:-1]).all(1); a = ow[:-1][same]; b = ow[1:][same]
    assert len(a) * 2 == 5 * n, 'every tetrahedron must be shared by exactly two simplices'
    return sp.csr_matrix((np.full(2 * len(a), 0.2), (np.concatenate([a, b]), np.concatenate([b, a]))), shape=(n, n))
def spectral(fn, nstart=200, smax=400, seed=1):
    T, vt, V, tt = load(fn); M = dual(V); n = len(V); rng = np.random.default_rng(seed)
    prof = np.bincount(tt, minlength=T); thr = 0.5 * prof.max(); cand = np.where(prof[tt] >= thr)[0]      # starts in slabs with at least half the peak volume
    starts = rng.choice(cand, min(nstart, len(cand)), replace=False); K = np.zeros((n, len(starts))); K[starts, np.arange(len(starts))] = 1.0
    P = np.zeros(smax + 1); P[0] = 1
    for s in range(1, smax + 1):
        K = M @ K; P[s] = K[starts, np.arange(len(starts))].mean()
    return n, P
if __name__ == '__main__':
    fn = sys.argv[1]; ns = int(sys.argv[2]) if len(sys.argv) > 2 else 200; smax = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    n, P = spectral(fn, ns, smax, int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    s = np.arange(2, smax + 1, 2); lp = np.log(P[s]); ls = np.log(s); Ds = -2 * np.gradient(lp, ls)
    print(f'{fn}: N4 = {n}; sigma, D_s, published fit 4.02 - 119/(54 + sigma)')
    for x in (20, 40, 60, 100, 150, 200, 300, 400):
        if x <= smax: i = list(s).index(x); print(f'  {x:4d}  {Ds[i]:.2f}  {4.02 - 119 / (54 + x):.2f}')
    np.save(fn + '.P.npy', P)
