"""SA74, exploratory (Nov 10, 23:55 PST): the causal order of a causal triangulation, read with A42's three dimension tests.
Elements: the vertices. Direct causes of a vertex at slice t+1: the vertices at slice t joined to it by a timelike edge (an edge of a
4-simplex between the two slices). The order is the transitive closure. Reported: vertices per slice; links per vertex; how fast one vertex's
future cone fills the slices above it; and, for intervals of growing height in the fat part of the universe, their size, the share of each
slice they hold, and the three readings (ordering fraction, midpoint, longest chain) against the bands of sprinkled flat intervals of the
nearest size (A42 Table 2). Usage: cdt_gate.py config [seed]"""
import sys, math, json, numpy as np
sys.path.insert(0, '/home/claude/sa74/lorentz'); sys.path.insert(0, '/home/claude/sa74/cdt')
import A42_Rule_Search as A
from lat2 import band
from sd import load
fn = sys.argv[1]; rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
T, vt, V, tt = load(fn); n0 = len(vt)
cnt = np.bincount(vt, minlength=T); t0 = int(np.argmin(cnt))            # unroll time from the thinnest slice
tau = (vt - t0) % T
order = np.argsort(tau, kind='stable'); rank = np.empty(n0, int); rank[order] = np.arange(n0)   # ids grow with time
par = [set() for _ in range(n0)]
for s, t in zip(V, tt):
    lo = [v for v in s if vt[v] == t]; hi = [v for v in s if vt[v] == (t + 1) % T]
    if (t - t0) % T == T - 1: continue                                    # the slab that closes the periodic time is cut
    for h in hi:
        for l in lo: par[rank[h]].add(int(rank[l]))
parents = [tuple(sorted(p)) for p in par]; anc, desc = A.from_parents(parents)
tr = tau[order]; slice_ids = [np.where(tr == k)[0] for k in range(T)]; size = np.array([len(x) for x in slice_ids])
nch = np.zeros(n0, int)
for e in range(n0):
    for p in parents[e]: nch[p] += 1
fat = [k for k in range(T) if size[k] >= 0.5 * size.max()]
print(f'{fn}: T {T}, vertices {n0}, 4-simplices {len(V)}; vertices per slice (unrolled): ' + ' '.join(map(str, size)))
print(f'fat slices {fat[0]} to {fat[-1]} ({len(fat)} slices, {size[fat].min()} to {size[fat].max()} vertices each); direct causes per vertex there: mean {np.mean([len(parents[e]) for k in fat for e in slice_ids[k]]):.1f}, '
      f'max {max(len(parents[e]) for k in fat for e in slice_ids[k])}; direct effects: mean {np.mean([nch[e] for k in fat[:-1] for e in slice_ids[k]]):.1f}')
# cone growth: share of slice a + k inside the future of a vertex at slice a
a = fat[0]; cov = []
for p in rng.choice(slice_ids[a], min(40, size[a]), replace=False):
    d = desc[int(p)]; cov.append([sum((d >> int(e)) & 1 for e in slice_ids[a + k]) / size[a + k] for k in range(1, min(9, T - a))])
cov = np.array(cov); print('share of the slice k steps later inside one vertex\'s future, k = 1..: ' + ' '.join(f'{x:.2f}' for x in cov.mean(0)))
# intervals of height H between a vertex at slice a and one at slice a + H, the largest of 12 tries per bottom vertex
res = {}
for H in range(2, 11):
    out = []
    for a in fat:
        if a + H > fat[-1] + 2 or a + H >= T: continue
        for p in rng.choice(slice_ids[a], min(6, size[a]), replace=False):
            p = int(p); tops = [int(e) for e in slice_ids[a + H] if (desc[p] >> int(e)) & 1]
            if not tops: continue
            best = max((int(q) for q in rng.choice(tops, min(12, len(tops)), replace=False)), key=lambda q: (desc[p] & anc[q]).bit_count())
            inner = desc[p] & anc[best]; members = [p] + [e for k in range(a + 1, a + H) for e in map(int, slice_ids[k]) if (inner >> e) & 1] + [best]
            fill = max(sum((inner >> int(e)) & 1 for e in slice_ids[k]) / size[k] for k in range(a + 1, a + H))
            if len(members) < 30: out.append((len(members), fill, None)); continue
            n, f, d_mm, d_mid, L, d_ch = A.diamond_dims(members, anc, desc, parents); out.append((n, fill, (d_mm, d_mid, d_ch, L)))
    if not out: continue
    ns = np.array([o[0] for o in out]); fl = np.array([o[1] for o in out]); rd = np.array([o[2][:3] for o in out if o[2]]); s = f'height {H:2d}: {len(out):3d} intervals, size median {np.median(ns):6.0f} (max {ns.max()}), largest share of a slice inside: median {np.median(fl):.2f}'
    if len(rd):
        m = np.median(rd, 0); b, k = band(float(np.median([o[0] for o in out if o[2]]))); z = [(m[i] - b[i][0]) / b[i][1] for i in range(3)]
        s += f'; readings {m[0]:.2f} / {m[1]:.2f} / {m[2]:.2f} (band at {k}: {b[0][0]:.2f}, {b[1][0]:.2f}, {b[2][0]:.2f}; z {z[0]:+.1f} {z[1]:+.1f} {z[2]:+.1f})'
        res[H] = dict(n=float(np.median(ns)), fill=float(np.median(fl)), mm=float(m[0]), mid=float(m[1]), chain=float(m[2]), z=z, band=k, intervals=len(out))
    print(s)
json.dump(res, open(fn + '.gate.json', 'w'), indent=1)
