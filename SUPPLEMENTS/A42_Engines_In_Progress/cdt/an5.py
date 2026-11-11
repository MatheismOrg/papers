"""SA74 plan B analysis of a cdt4 run: the measurement rows split into NB blocks (default 8); in each block the volume profile N3(t),
centered row by row on its center of volume (periodic), is averaged and fitted to c0 + A cos^3((t - t0)/B). Reports B per block, its mean
and standard error over blocks, s0 = B / N4^(1/4) with N4 the total number of four-simplices (the published convention, arXiv:0807.4481
eq. 12: s0 = 0.59 at (kappa0, Delta) = (2.2, 0.6), verified there for N4 = 45,500 to 362,000; with N41 in place of N4 the published
value is 0.59 (N4/N41)^(1/4) = 0.72), N4/N41, N0/N41, k4, and the fit of the whole-run average. PHASE=0 FROM=sweep analyzes the
thermalization rows (tuner on) from that sweep. Usage: an5.py run.dat [NB] [skip_rows]"""
import os
import sys, numpy as np
from scipy.optimize import curve_fit
f = sys.argv[1]; NB = int(sys.argv[2]) if len(sys.argv) > 2 else 8; skip = int(sys.argv[3]) if len(sys.argv) > 3 else 0
L = [l.split('|') for l in open(f) if not l.startswith('#') and '|' in l]
T = len(L[0][1].split()); L = [x for x in L if len(x[1].split()) == T and len(x[0].split()) == 6]
h = np.array([x[0].split() for x in L], float); pr = np.array([x[1].split() for x in L], float)
m = (h[:, 1] == int(os.environ.get('PHASE', 1))) & (h[:, 0] >= int(os.environ.get('FROM', 0))); h = h[m][skip:]; pr = pr[m][skip:]
if len(h) < 40: print(f, 'measurement rows', len(h), '(too few)'); sys.exit()
def center(P):
    ang = np.angle((P * np.exp(2j * np.pi * np.arange(T) / T)).sum(1)); sh = np.rint(ang / (2 * np.pi) * T).astype(int) % T
    return np.array([np.roll(p, T // 2 - s) for p, s in zip(P, sh)])
def f3(t, c0, A, B, t0):
    x = (t - t0) / B; return c0 + A * np.where(abs(x) < np.pi / 2, np.cos(x) ** 3, 0)
C = center(pr); t = np.arange(T) - T // 2; N41 = h[:, 4].mean(); N4 = h[:, 3].mean()
def fit(a):
    w = (a * (t - (a * t).sum() / a.sum()) ** 2).sum() / a.sum()
    p, _ = curve_fit(f3, t, a, p0=[a.min(), a.max(), max(2.0, np.sqrt(w)), 0], maxfev=20000); r = a - f3(t, *p)
    return p, np.sqrt((r ** 2).mean()) / p[1]
print(f'{f}: T {T}, measurement rows {len(h)} (sweeps {int(h[0, 0])} to {int(h[-1, 0])}), N41 {N41:.0f}, N4/N41 {(h[:, 3] / h[:, 4]).mean():.3f}, N0/N41 {(h[:, 2] / h[:, 4]).mean():.4f}, k4 {h[:, 5].mean():.4f}')
Bs = []
for i, idx in enumerate(np.array_split(np.arange(len(C)), NB)):
    p, rr = fit(C[idx].mean(0)); Bs.append(p[2]); print(f'  block {i + 1}: stalk {p[0]:.1f} height {p[1]:.0f} B {p[2]:.3f} residual {100 * rr:.1f}% N41 {h[idx, 4].mean():.0f}')
Bs = np.array(Bs); p, rr = fit(C.mean(0)); q = N4 ** 0.25
print(f'  B over blocks {Bs.mean():.3f}({Bs.std(ddof=1) / np.sqrt(NB):.3f}); whole run: stalk {p[0]:.1f} height {p[1]:.0f} B {p[2]:.3f} residual {100 * rr:.1f}%; N4 {N4:.0f}, N4^(1/4) {q:.3f}; s0 = B/N4^(1/4) {Bs.mean() / q:.3f}({Bs.std(ddof=1) / np.sqrt(NB) / q:.3f}) (published 0.59); with N41: {Bs.mean() / N41 ** 0.25:.3f} (published 0.72)')
a = C.mean(0); print('  profile', ' '.join('%d' % x for x in a))
