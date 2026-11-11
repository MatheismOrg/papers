"""SA74 plan C analysis, as fixed in jobs4/srv/C0_build.job (Nov 10 13:42 PST, before the runs) and its addendum more_runs.txt:
the first fifth of each run's measured sweeps is dropped; the rest is cut into blocks of 100,000 sweeps; in each block the profile,
centered row by row, is fitted to c0 + A cos^3((t - t0)/B). A run's B is the mean over its blocks, its error the largest batch-means
error over batches of 1, 2, 4, ... blocks (at least 4 batches). Runs at one size are combined by inverse variance, the error not below
the spread between runs. Reported: B; s0 = B / N4^(1/4) with the run's mean N4 (published 0.59); N4/N41 (published 2.26 to 2.275); the
first against the second half of the kept blocks; the power of B in N4 over the sizes (four-dimensional de Sitter space gives 1/4).
A run still going is analyzed on the sweeps it has and marked partial. Usage: planC.py OUT.json run.dat:planned_measured_sweeps ..."""
import sys, json, numpy as np
from scipy.optimize import curve_fit
def load(f):
    rows = [l.split('|') for l in open(f) if not l.startswith('#') and '|' in l]
    T = len(rows[0][1].split()); rows = [x for x in rows if len(x[1].split()) == T and len(x[0].split()) == 6]
    return np.array([x[0].split() for x in rows], float), np.array([x[1].split() for x in rows], float), T
def center(P, T):
    ang = np.angle((P * np.exp(2j * np.pi * np.arange(T) / T)).sum(1)); sh = np.rint(ang / (2 * np.pi) * T).astype(int) % T
    return np.array([np.roll(p, T // 2 - s) for p, s in zip(P, sh)])
def f3(t, c0, A, B, t0):
    x = (t - t0) / B; return c0 + A * np.where(abs(x) < np.pi / 2, np.cos(x) ** 3, 0)
def fit(a, t):
    w = (a * (t - (a * t).sum() / a.sum()) ** 2).sum() / a.sum()
    p, _ = curve_fit(f3, t, a, p0=[a.min(), a.max(), max(2.0, np.sqrt(w)), 0], maxfev=20000); r = a - f3(t, *p)
    return p, float(np.sqrt((r ** 2).mean()) / p[1])
def batch_se(v):
    v = np.asarray(v, float); best = 0.0; m = 1; used = None
    while len(v) // m >= 4:
        nb = len(v) // m; bm = v[:nb * m].reshape(nb, m).mean(1); se = bm.std(ddof=1) / np.sqrt(nb)
        if se > best: best, used = se, m
        m *= 2
    return float(best), used
out = {}
for arg in sys.argv[2:]:
    f, planned = arg.split(':'); planned = int(planned); name = f.split('/')[-1][:-4]
    h, pr, T = load(f); m = h[:, 1] == 1; h, pr = h[m], pr[m]
    if len(h) < 10: out[name] = dict(note='no measured sweeps yet'); print(name, 'no measured sweeps yet'); continue
    step = int(round(h[1, 0] - h[0, 0])); done = len(h) * step; s0 = h[0, 0]
    keep = h[:, 0] >= s0 + done / 5; h, pr = h[keep], pr[keep]
    per = 100000 // step; nblk = len(h) // per
    if nblk < 4: out[name] = dict(note='fewer than four blocks of 100,000 sweeps after the first fifth', measured=done, planned=planned); print(name, 'too short:', done, 'measured sweeps'); continue
    t = np.arange(T) - T // 2; Bs, res = [], []
    for b in range(nblk):
        C = center(pr[b * per:(b + 1) * per], T); p, rr = fit(C.mean(0), t); Bs.append(abs(p[2])); res.append(rr)
    Bs = np.array(Bs); hk = h[:nblk * per]; N4 = hk[:, 3].mean(); N41 = hk[:, 4].mean(); se, bsz = batch_se(Bs)
    h1, h2 = Bs[:nblk // 2], Bs[nblk // 2:]
    r = dict(run=name, T=T, measured_sweeps=done, planned=planned, partial=bool(done < planned), blocks=int(nblk), B=float(Bs.mean()), B_se=se, batch_blocks=bsz,
             N4=float(N4), N41=float(N41), N4_over_N41=float((hk[:, 3] / hk[:, 4]).mean()), N0_over_N41=float((hk[:, 2] / hk[:, 4]).mean()), k4=float(hk[:, 5].mean()),
             s0=float(Bs.mean() / N4 ** 0.25), s0_se=float(se / N4 ** 0.25), B_first_half=float(h1.mean()), B_second_half=float(h2.mean()),
             halves_z=float((h1.mean() - h2.mean()) / np.hypot(h1.std(ddof=1) / np.sqrt(len(h1)), h2.std(ddof=1) / np.sqrt(len(h2)))) if len(h1) > 1 and len(h2) > 1 else None,
             B_min=float(Bs.min()), B_max=float(Bs.max()), residual_percent=float(100 * np.mean(res)), B_blocks=[round(float(x), 3) for x in Bs])
    out[name] = r
    print('%-10s %s measured %9d of %9d, blocks %2d, N41 %6.0f, N4 %6.0f, N4/N41 %.3f, k4 %.5f, B %.3f(%.3f) [batches of %s], range %.2f to %.2f, halves %.3f | %.3f, s0 %.3f(%.3f), fit residual %.1f%%' % (
        name, 'PARTIAL' if r['partial'] else 'complete', done, planned, nblk, N41, N4, r['N4_over_N41'], r['k4'], r['B'], se, bsz, Bs.min(), Bs.max(), h1.mean(), h2.mean(), r['s0'], r['s0_se'], r['residual_percent']))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
