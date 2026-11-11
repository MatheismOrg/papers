#!/usr/bin/env python3
"""A42 search: the shortest local rewriting rules, the causal order of their updates, and its dimension.

Rules rewrite binary relations (directed edges between nodes), as in Wolfram's models. A rule replaces a connected pattern of
relations (its left side) by new relations (its right side), which may use new nodes. Updating is generational: in each
generation the alive relations are scanned oldest first, each is matched with the oldest compatible partner, non-overlapping
matches are collected, and all of them are applied at once. Every application is an event (an update). An event's direct
causes are the events that created the relations it consumed; a root event creates the initial condition, two self-loops on
one node. The direct causes, closed under transitivity, form a causal order, the structure a causal set is made of.

Dimension is read as causal-set theory reads it, on causal intervals (all events between two related events x < y):
  * Myrheim-Meyer: the fraction f of related pairs in an interval of d-dimensional flat spacetime is
    f(d) = Gamma(d+1) Gamma(d/2) / (2 Gamma(3d/2)); invert the measured f.
  * Midpoint scaling: the largest interval [x,z] or [z,y] with z between x and y holds about N / 2^d elements.
  * Longest chain: the longest chain grows as m_d N^(1/d); m_d is calibrated on sprinkled flat intervals (validate).
  * Coarse-grained chain (a check): keep each element with probability p and fit how the longest chain shrinks.
A manifold-like order gives the same d from all of them; orders that are not manifold-like generally do not.
validate() also counts links (relations with nothing between) in sprinkled four-dimensional intervals and reads two grids
as controls: the light-cone lattice of 1+1 dimensions and a three-dimensional cellular automaton.

Python 3.10 or later (int.bit_count), standard library only. Results are deterministic, apart from the time limits.

usage:
  python3 -I a42_search.py validate [OUT.json]   # sprinkled flat diamonds, d = 2 to 6, link counts, two grid controls
  python3 -I a42_search.py enumerate OUT.json    # the 2,935 canonical rules rewriting two relations into three
  python3 -I a42_search.py sample OUT.json [count]  # a random sample of rules rewriting two relations into four
  python3 -I a42_search.py run RULES.json OUT.jsonl [--events N] [--start i --stop j] [--sec S] [--series]
  python3 -I a42_search.py one "RULE" [events]   # run one rule, e.g. "{{x,y},{x,z}}->{{x,y},{x,w},{z,x}}"
  python3 -I a42_search.py sample2 OUT.json n_left n_right arity count seed   # second pass: any signature
  python3 -I a42_search.py run2 RULES.json OUT.jsonl [--events N] [--start i --stop j] [--sec S]
  python3 -I a42_search.py long RULES.json OUT.jsonl --events 120000 --sec 600   # read inside a top update's past
  python3 -I a42_search.py desitter [OUT.json]    # sprinkled four-dimensional de Sitter intervals, 0 to 4 e-folds
The paper's runs: enumerate, then run with --events 3000; sample, then run with --events 3000; the 27 rules run long with
--events 20000 --sec 60 --series. Second pass: sample2 3 4 2 1000 42 and sample2 2 3 3 1000 42, each run with run2
--events 3000; slow rules that stop below 300 updates rerun with --sec 30; ten rules run with long --events 120000
--sec 600. run_general reproduces run exactly for two-relation left sides.
"""
import json, math, random, sys, time, itertools
from collections import defaultdict

# ---------------------------------------------------------------- dimension estimators
def f_of_d(d):
    return math.exp(math.lgamma(d + 1) + math.lgamma(d / 2) - math.log(2) - math.lgamma(3 * d / 2))

def mm_dimension(f):
    """Invert f(d) on [1, 20] by bisection; f(1) = 1, f decreases with d."""
    if f >= 1.0: return 1.0
    if f <= f_of_d(20.0): return 20.0
    lo, hi = 1.0, 20.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if f_of_d(mid) > f: lo = mid
        else: hi = mid
    return (lo + hi) / 2

# Longest chain ~ m_d N^(1/d) in a sprinkled diamond of N elements. These constants are the calibration that validate()
# measures (three seeds at N = 400 and 800 for each d) and are used as they stand; validate() recomputes and replaces them.
CHAIN_M = {2: 1.9544, 3: 2.0567, 4: 2.1443, 5: 2.2492, 6: 2.3}

def chain_dimension(n, length):
    """Solve length = m(d) n^(1/d) for d, with m(d) interpolated from the calibration table."""
    ds = sorted(CHAIN_M)
    def m(d):
        if d <= ds[0]: return CHAIN_M[ds[0]]
        if d >= ds[-1]: return CHAIN_M[ds[-1]]
        lo = max(x for x in ds if x <= d); hi = min(x for x in ds if x >= d)
        return CHAIN_M[lo] if lo == hi else CHAIN_M[lo] + (CHAIN_M[hi] - CHAIN_M[lo]) * (d - lo) / (hi - lo)
    g = lambda d: m(d) * n ** (1 / d) - length          # decreasing in d
    lo, hi = 1.0, 20.0
    if g(lo) <= 0: return 1.0
    if g(hi) >= 0: return 20.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if g(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def longest_chain(sub, anc):
    """Longest chain (number of elements) in the order induced on the ids in sub, which must be listed so that every
    element comes after its ancestors (ids here grow with causality). Levels are antichains, and if an element's past
    meets level k it meets every level below k, so each element finds its level by bisection."""
    levels = []                                   # levels[k] = bitset of elements whose longest chain ending there is k + 1
    for e in sub:
        a = anc[e]; lo, hi = 0, len(levels)       # find the number of levels a meets
        while lo < hi:
            mid = (lo + hi) // 2
            if a & levels[mid]: lo = mid + 1
            else: hi = mid
        if lo == len(levels): levels.append(0)
        levels[lo] |= 1 << e
    return len(levels)

def subsample_slope(members, anc, ps=(1.0, 0.5, 0.25, 0.125), reps=4, seed=7):
    """Coarse-graining by random selection: keep each element with probability p (always keeping the two tips) and
    measure the longest chain. In a d-dimensional manifold-like order the chain grows as (number kept)^(1/d), so the
    slope of log chain against log size gives d with no fitted constant. Returns (d, [(p, mean size, mean chain)])."""
    rnd = random.Random(seed)
    mask = 0
    for e in members: mask |= 1 << e
    lo = min(members, key=lambda e: (anc[e] & mask).bit_count())
    hi = max(members, key=lambda e: (anc[e] & mask).bit_count())
    order = sorted(members); pts = []
    for p in ps:
        sizes = []; chains = []
        for r in range(1 if p == 1.0 else reps):
            sub = [e for e in order if e in (lo, hi) or p == 1.0 or rnd.random() < p]
            sizes.append(len(sub)); chains.append(longest_chain(sub, anc))
        pts.append((p, sum(sizes) / len(sizes), sum(chains) / len(chains)))
    xs = [math.log(n) for _, n, _ in pts]; ys = [math.log(c) for _, _, c in pts]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys); vx = sum((u - mx) ** 2 for u in xs)
    slope = sum((u - mx) * (v - my) for u, v in zip(xs, ys)) / vx
    return (1 / slope if slope > 0 else float('inf')), [(p, round(n, 1), round(c, 2)) for p, n, c in pts]

def diamond_dims(members, anc, desc, parents=None):
    """members: element ids in the diamond (both tips included), in an order compatible with causality;
    anc/desc: strict ancestor/descendant bitsets; parents: direct causes (optional, for a fast longest chain).
    Returns (N, ordering fraction, Myrheim-Meyer d, midpoint d, longest chain, chain d)."""
    mask = 0
    for e in members: mask |= 1 << e
    n = len(members)
    rel = sum((anc[e] & mask).bit_count() for e in members)
    f = rel / (n * (n - 1) / 2)
    d_mm = mm_dimension(f)
    lo = min(members, key=lambda e: (anc[e] & mask).bit_count())   # bottom tip
    hi = max(members, key=lambda e: (anc[e] & mask).bit_count())   # top tip
    best = 0
    for z in members:
        if z in (lo, hi): continue
        a = (anc[z] & mask).bit_count() + 1
        b = (desc[z] & mask).bit_count() + 1
        best = max(best, min(a, b))
    d_mid = math.log2(n / best) if best > 0 else float('inf')
    length = longest_chain(sorted(members), anc)   # an interval is convex, so its longest chain runs along direct causes
    return n, f, d_mm, d_mid, length, chain_dimension(n, length)

# ---------------------------------------------------------------- validation on sprinkled flat diamonds
def sprinkle(d, n, seed=1):
    rnd = random.Random(seed); pts = []
    while len(pts) < n:
        t = rnd.uniform(-0.5, 0.5); x = [rnd.uniform(-0.5, 0.5) for _ in range(d - 1)]
        if math.sqrt(sum(v * v for v in x)) <= 0.5 - abs(t): pts.append((t, x))
    pts.sort(key=lambda p: p[0])
    pts = [(-0.5, [0.0] * (d - 1))] + pts + [(0.5, [0.0] * (d - 1))]      # the two tips
    m = len(pts); anc = [0] * m; desc = [0] * m
    for j in range(m):
        tj, xj = pts[j]
        for i in range(j):
            ti, xi = pts[i]
            if tj - ti >= math.sqrt(sum((a - b) ** 2 for a, b in zip(xi, xj))):
                anc[j] |= 1 << i; desc[i] |= 1 << j
    return list(range(m)), anc, desc

def from_parents(parents):
    """Ancestor and descendant bitsets of a graph given by direct causes; ids must grow with causality."""
    E = len(parents); anc = [0] * E; children = [[] for _ in range(E)]
    for e in range(E):
        a = 0
        for p in parents[e]: a |= anc[p] | (1 << p); children[p].append(e)
        anc[e] = a
    desc = [0] * E
    for e in range(E - 1, -1, -1):
        d = 0
        for c in children[e]: d |= desc[c] | (1 << c)
        desc[e] = d
    return anc, desc

def lattice2(n):
    """Control: the light-cone lattice of 1+1 dimensions, the causal order of a one-dimensional cellular automaton.
    Element (i, j) precedes (i', j') when i <= i' and j <= j' (i and j count steps along the two light rays). The
    diamond from (0, 0) to (n, n) holds (n + 1)^2 elements, and every element has at most two direct causes."""
    pts = sorted(((i, j) for i in range(n + 1) for j in range(n + 1)), key=lambda p: (p[0] + p[1], p[0]))
    idx = {p: k for k, p in enumerate(pts)}
    parents = [tuple(sorted(idx[q] for q in ((i - 1, j), (i, j - 1)) if q in idx)) for (i, j) in pts]
    anc, desc = from_parents(parents)
    return list(range(len(pts))), anc, desc, parents

def ca3(T):
    """Control: the causal order of a three-dimensional cellular automaton with the 27-cell neighbourhood, a 3+1 grid.
    Cell x at step t has as direct causes the cells within one step in every axis at step t - 1, so signals spread
    in a square-sided cone. The diamond from (0, origin) to (T, origin) is the set with max|x_i| <= min(t, T - t)."""
    pts = []
    for t in range(T + 1):
        r = min(t, T - t)
        for x in range(-r, r + 1):
            for y in range(-r, r + 1):
                for z in range(-r, r + 1): pts.append((t, x, y, z))
    idx = {p: k for k, p in enumerate(pts)}
    parents = []
    for (t, x, y, z) in pts:
        ps = [idx[q] for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
              for q in [(t - 1, x + dx, y + dy, z + dz)] if q in idx]
        parents.append(tuple(sorted(ps)))
    anc, desc = from_parents(parents)
    return list(range(len(pts))), anc, desc, parents

def link_counts(members, anc, desc):
    """Links are the relations with nothing between them. Returns (mean links per element, largest number)."""
    mask = 0
    for e in members: mask |= 1 << e
    deg = {e: 0 for e in members}
    for f in members:
        x = anc[f] & mask
        while x:
            lb = x & -x; e = lb.bit_length() - 1; x ^= lb
            if not (desc[e] & anc[f] & mask): deg[e] += 1; deg[f] += 1
    inner = [deg[e] for e in members[1:-1]]          # the two tips sit on the boundary and are left out
    return sum(inner) / len(inner), max(inner)

def validate():
    print('Sprinkled causal diamonds in d-dimensional Minkowski spacetime:')
    out = []; consts = {}
    for d in (2, 3, 4, 5, 6):
        for n in (400, 800):
            vals = []; slopes = []
            for seed in (1, 2, 3):
                members, anc, desc = sprinkle(d, n, seed)
                vals.append(diamond_dims(members, anc, desc))
                slopes.append(subsample_slope(members, anc, seed=seed)[0])
            mm = sum(v[2] for v in vals) / 3; mid = sum(v[3] for v in vals) / 3; L = sum(v[4] for v in vals) / 3
            N = sum(v[0] for v in vals) / 3
            consts.setdefault(d, []).append(L / N ** (1 / d))
            out.append({'d': d, 'n': n, 'mm': round(mm, 3), 'mid': round(mid, 3), 'chain': round(L, 2),
                        'coarse_d': round(sum(slopes) / 3, 3)})
    for d in consts: CHAIN_M[d] = sum(consts[d]) / len(consts[d])
    for o in out:
        o['chain_d'] = round(chain_dimension(o['n'] + 2, o['chain']), 3)
        print(f"  d={o['d']} N={o['n']}: Myrheim-Meyer {o['mm']:.2f}, midpoint {o['mid']:.2f}, longest chain {o['chain']:.1f} -> {o['chain_d']:.2f}, coarse-grained chain {o['coarse_d']:.2f}")
    print('  chain constants m_d:', {d: round(v, 3) for d, v in CHAIN_M.items()})
    links = []
    for n in (200, 400, 800, 1600):
        members, anc, desc = sprinkle(4, n, 1)
        mean, mx = link_counts(members, anc, desc)
        links.append({'d': 4, 'n': n, 'mean_links': round(mean, 2), 'max_links': mx})
        print(f'  links per element, d=4 N={n}: mean {mean:.1f}, largest {mx}')
    controls = []
    for name, (members, anc, desc, parents) in (('lattice 1+1, n=29', lattice2(29)), ('cellular automaton 3+1, T=10', ca3(10))):
        n, f, d_mm, d_mid, L, d_ch = diamond_dims(members, anc, desc, parents)
        d_co = subsample_slope(members, anc)[0]
        mean, mx = link_counts(members, anc, desc)
        controls.append({'order': name, 'n': n, 'mm': round(d_mm, 3), 'mid': round(d_mid, 3), 'chain': L, 'chain_d': round(d_ch, 3),
                         'coarse_d': round(d_co, 3), 'mean_links': round(mean, 2), 'max_links': mx})
        print(f'  control {name}: N={n}, Myrheim-Meyer {d_mm:.2f}, midpoint {d_mid:.2f}, chain {L} -> {d_ch:.2f}, coarse-grained {d_co:.2f}, links mean {mean:.1f} max {mx}')
    return {'diamonds': out, 'chain_m': {str(d): round(v, 4) for d, v in CHAIN_M.items()}, 'links': links, 'controls': controls}

# ---------------------------------------------------------------- rules
def parse_rule(s):
    s = s.replace(' ', '')
    l, r = s.split('->')
    def rels(x):
        x = x.strip('{}')
        return [tuple(p.split(',')) for p in x.split('},{')]
    L, R = rels(l), rels(r)
    names = []
    for rel in L + R:
        for v in rel:
            if v not in names: names.append(v)
    idx = {v: i for i, v in enumerate(names)}
    k = len({v for rel in L for v in rel})
    return tuple(tuple(idx[v] for v in rel) for rel in L), tuple(tuple(idx[v] for v in rel) for rel in R), k

def rule_str(L, R):
    letters = 'xyzwuvst'
    f = lambda rels: '{' + ','.join('{' + ','.join(letters[v] for v in rel) + '}' for rel in rels) + '}'
    return f(L) + '->' + f(R)

def connected(rels):
    verts = {v for rel in rels for v in rel}
    if not verts: return False
    adj = defaultdict(set)
    for a, b in rels: adj[a].add(b); adj[b].add(a)
    start = next(iter(verts)); seen = {start}; stack = [start]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen: seen.add(w); stack.append(w)
    return seen == verts

def canonical(L, R, k, m):
    best = None
    for perm in itertools.permutations(range(k)):
        for nperm in itertools.permutations(range(k, k + m)):
            mp = {i: perm[i] for i in range(k)}; mp.update({k + j: nperm[j] for j in range(m)})
            key = (tuple(sorted(tuple(mp[v] for v in rel) for rel in L)), tuple(sorted(tuple(mp[v] for v in rel) for rel in R)))
            if best is None or key < best: best = key
    return best

def first_appearance_ok(seq, start, count):
    """New symbols start..start+count-1 must first appear in increasing order."""
    nxt = start
    for v in seq:
        if v >= start:
            if v > nxt: return False
            if v == nxt: nxt += 1
    return nxt == start + count

def enumerate_rules(n_left=2, n_right=3, max_new=2):
    """Connected left sides of n_left binary relations; right sides of n_right binary relations over the left-side
    variables and up to max_new new nodes, connected and using every new node. Returns canonical rules, shortest first."""
    lefts = set()
    for vals in itertools.product(range(2 * n_left), repeat=2 * n_left):
        k = len(set(vals))
        if not first_appearance_ok(vals, 0, k): continue
        L = tuple(zip(vals[0::2], vals[1::2]))
        if connected(L): lefts.add((L, k))
    rules = {}
    for L, k in sorted(lefts):
        for m in range(0, max_new + 1):
            alpha = k + m
            for vals in itertools.product(range(alpha), repeat=2 * n_right):
                if m and not first_appearance_ok(vals, k, m): continue
                R = tuple(zip(vals[0::2], vals[1::2]))
                if not connected(R): continue
                key = canonical(L, R, k, m)
                if key not in rules: rules[key] = (k, m)
    out = []
    for (L, R), (k, m) in rules.items():
        bits = (2 * (len(L) + len(R))) * math.log2(k + m)       # symbols times bits per symbol
        out.append({'rule': rule_str(L, R), 'L': L, 'R': R, 'k': k, 'm': m, 'bits': round(bits, 2)})
    out.sort(key=lambda r: (r['bits'], r['rule']))
    return out


def sample_rules(n_left=2, n_right=4, max_new=2, count=1500, seed=42):
    """Random canonical rules from a larger class, for a first look beyond the shortest one."""
    rnd = random.Random(seed)
    lefts = []
    for vals in itertools.product(range(2 * n_left), repeat=2 * n_left):
        k = len(set(vals))
        if not first_appearance_ok(vals, 0, k): continue
        L = tuple(zip(vals[0::2], vals[1::2]))
        if connected(L): lefts.append((L, k))
    seen = {}
    tries = 0
    while len(seen) < count and tries < count * 200:
        tries += 1
        L, k = rnd.choice(lefts); m = rnd.randint(0, max_new)
        vals = [rnd.randrange(k + m) for _ in range(2 * n_right)]
        if m and not first_appearance_ok(vals, k, m): continue
        R = tuple(zip(vals[0::2], vals[1::2]))
        if not connected(R): continue
        key = canonical(L, R, k, m)
        if key not in seen: seen[key] = (k, m)
    out = []
    for (L, R), (k, m) in seen.items():
        bits = (2 * (len(L) + len(R))) * math.log2(k + m)
        out.append({'rule': rule_str(L, R), 'L': L, 'R': R, 'k': k, 'm': m, 'bits': round(bits, 2)})
    out.sort(key=lambda r: (r['bits'], r['rule']))
    return out

# ---------------------------------------------------------------- the rewriting engine
def run(L, R, k, max_events=2000, max_gen=400, max_rel=60000, init=((0, 0), (0, 0)), max_sec=8.0):
    t_start = time.time()
    rel = {}; creator = {}
    out_idx = defaultdict(set); in_idx = defaultdict(set)
    nid = 0; node = 1
    def add(a, b, ev):
        nonlocal nid
        rel[nid] = (a, b); creator[nid] = ev; out_idx[a].add(nid); in_idx[b].add(nid); nid += 1
    def remove(i):
        a, b = rel.pop(i); out_idx[a].discard(i); in_idx[b].discard(i); creator.pop(i, None)
    for a, b in init: add(a, b, 0)
    parents = [()]; gens = [0]
    P = L
    def extend(binding, pat, a, b):
        p, q = pat
        bd = dict(binding)
        for var, val in ((p, a), (q, b)):
            if var in bd:
                if bd[var] != val: return None
            else: bd[var] = val
        return bd
    def candidates(binding, pat):
        p, q = pat
        if p in binding: return out_idx[binding[p]]
        if q in binding: return in_idx[binding[q]]
        return None
    for g in range(1, max_gen + 1):
        used = set(); events = []
        for r in sorted(rel):
            if r in used: continue
            best = None
            for first in range(len(P)):
                bd = extend({}, P[first], *rel[r])
                if bd is None: continue
                if len(P) == 1:
                    best = ((r,), bd); break
                other = 1 - first
                cand = candidates(bd, P[other])
                for c in sorted(cand):
                    if c == r or c in used: continue
                    bd2 = extend(bd, P[other], *rel[c])
                    if bd2 is None: continue
                    ids = (r, c) if first == 0 else (c, r)
                    if best is None or max(ids) < max(best[0]): best = (ids, bd2)
                    break
            if best is not None:
                used.update(best[0]); events.append(best)
        if not events: break
        for ids, bd in events:
            ev = len(parents)
            parents.append(tuple(sorted({creator[i] for i in ids}))); gens.append(g)
            for i in ids: remove(i)
            fresh = {}
            def node_of(v):
                nonlocal node
                if v < k: return bd[v]
                if v not in fresh: fresh[v] = node; node += 1
                return fresh[v]
            for a, b in R: add(node_of(a), node_of(b), ev)
        if len(parents) - 1 >= max_events or len(rel) > max_rel or time.time() - t_start > max_sec: break
    return parents, gens, len(rel), node

def pick_x(y, frac, anc, desc, depth, E):
    """The event at depth round(frac * depth[y]) in y's past whose interval with y is largest (0 = the root event)."""
    if frac == 0.0: return 0
    target = max(1, round(depth[y] * frac))
    cands = [x for x in range(E) if (anc[y] >> x) & 1 and depth[x] == target]
    return max(cands, key=lambda x: (desc[x] & anc[y]).bit_count()) if cands else None

def measure(parents, n_diamonds=6, min_size=120, series=False):
    E = len(parents)
    anc, desc = from_parents(parents)
    nchild = [0] * E
    for e in range(E):
        for p in parents[e]: nchild[p] += 1
    depth = [0] * E
    for e in range(1, E): depth[e] = 1 + max(depth[p] for p in parents[e])
    maxdepth = max(depth)
    tree = all(len(p) <= 1 for p in parents)
    tops = sorted(range(1, E), key=lambda e: (anc[e].bit_count(), depth[e]), reverse=True)[:n_diamonds * 4]
    out = {'events': E - 1, 'depth': maxdepth, 'tree': tree, 'past_max': anc[tops[0]].bit_count() if tops else 0,
           'max_causes': max(len(p) for p in parents[1:]), 'max_effects': max(nchild[1:])}
    # cone dimension: the fullest past cone at each depth level grows as depth^D in a D-dimensional spacetime
    best = defaultdict(int)
    for e in range(1, E): best[depth[e]] = max(best[depth[e]], anc[e].bit_count() + 1)
    lv = [(L, best[L]) for L in range(max(2, maxdepth // 3), maxdepth + 1) if best[L] > 1]
    if len(lv) >= 3:
        xs = [math.log(L) for L, _ in lv]; ys = [math.log(v) for _, v in lv]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys); vx = sum((u - mx) ** 2 for u in xs)
        out['cone_dim'] = round(sum((u - mx) * (v - my) for u, v in zip(xs, ys)) / vx, 3) if vx else None
    med = lambda v: sorted(v)[len(v) // 2]
    for label, frac in (('root', 0.0), ('mid', 0.3)):
        res = []; seen = set()
        for y in tops:
            x = pick_x(y, frac, anc, desc, depth, E)
            if x is None: continue
            inner = desc[x] & anc[y]
            members = [x] + [e for e in range(E) if (inner >> e) & 1] + [y]
            if len(members) < min_size or (x, y) in seen: continue
            seen.add((x, y))
            res.append(diamond_dims(members, anc, desc, parents) + (subsample_slope(members, anc)[0],))
            if len(res) >= n_diamonds: break
        if res:
            out[label] = {'diamonds': len(res), 'n': round(sum(r[0] for r in res) / len(res)),
                          'f': round(sum(r[1] for r in res) / len(res), 4), 'mm': round(med([r[2] for r in res]), 3),
                          'mid': round(med([r[3] for r in res]), 3), 'chain': round(sum(r[4] for r in res) / len(res), 1),
                          'chain_d': round(med([r[5] for r in res]), 3), 'coarse_d': round(med([r[6] for r in res]), 3)}
    if series and tops:
        # how the readings change with the size of the interval: one top event, bottom tips at growing depth
        y = tops[0]; ser = []
        for frac in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
            x = pick_x(y, frac, anc, desc, depth, E)
            if x is None: continue
            inner = desc[x] & anc[y]
            members = [x] + [e for e in range(E) if (inner >> e) & 1] + [y]
            if len(members) < 60: break
            n, f, d_mm, d_mid, L, d_ch = diamond_dims(members, anc, desc, parents)
            d_co, pts = subsample_slope(members, anc)
            ser.append({'frac': frac, 'n': n, 'f': round(f, 4), 'mm': round(d_mm, 3), 'mid': round(d_mid, 3), 'chain': L,
                        'chain_d': round(d_ch, 3), 'coarse_d': round(d_co, 3), 'height': depth[y] - depth[x]})
        out['series'] = ser
    return out

def growth(gens):
    cnt = defaultdict(int)
    for g in gens[1:]: cnt[g] += 1
    G = max(cnt) if cnt else 0
    c = [cnt[g] for g in range(1, G + 1)]
    fit = {}
    if G >= 6:
        half = [(g, c[g - 1]) for g in range(max(2, G // 2), G) if c[g - 1] > 0]   # last generation may be cut off
        if len(half) >= 3:
            xs = [g for g, _ in half]; ls = [math.log(v) for _, v in half]; lg = [math.log(g) for g in xs]
            def slope(a, b):
                ma, mb = sum(a) / len(a), sum(b) / len(b)
                va = sum((u - ma) ** 2 for u in a)
                return sum((u - ma) * (v - mb) for u, v in zip(a, b)) / va if va else 0.0
            fit = {'exp_rate': round(slope(xs, ls), 3), 'power': round(slope(lg, ls), 3)}
    return c, fit

# ---------------------------------------------------------------- second pass (SO70 T214): any signature, long runs, de Sitter

def connected_general(rels):
    verts = {v for rel in rels for v in rel}
    if not verts: return False
    parent = {v: v for v in verts}
    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for rel in rels:
        for a, b in zip(rel, rel[1:]): parent[find(a)] = find(b)
    return len({find(v) for v in verts}) == 1

def run_general(L, R, k, max_events=2000, max_gen=400, max_rel=60000, init=None, max_sec=8.0):
    """Generational updating as in run(): each generation scans the relations present at its start, oldest first, and
    matches each free one, as every pattern it fits, with the oldest compatible free partners (the remaining patterns are
    filled in pattern order, each from the oldest candidate that fits, backtracking when a later pattern finds none);
    of the ways a relation can enter a match, the one whose newest relation is oldest wins. Matches in a generation share
    no relation and are applied together. For two-relation left sides this is exactly run()."""
    t_start = time.time()
    P = [tuple(p) for p in L]; n = len(P); arity = len(P[0])
    if init is None: init = [tuple([0] * arity)] * n
    rel = {}; creator = {}; pos_idx = defaultdict(set)
    nid = 0; node = 1
    def add(t, ev):
        nonlocal nid
        rel[nid] = t; creator[nid] = ev
        for i, v in enumerate(t): pos_idx[(i, v)].add(nid)
        nid += 1
    def remove(r):
        t = rel.pop(r); creator.pop(r, None)
        for i, v in enumerate(t): pos_idx[(i, v)].discard(r)
    for t in init: add(tuple(t), 0)
    parents = [()]; gens = [0]
    def extend(bd, pat, t):
        if len(pat) != len(t): return None
        b = dict(bd)
        for var, val in zip(pat, t):
            if var in b:
                if b[var] != val: return None
            else: b[var] = val
        return b
    def cands(bd, pat):
        for i, var in enumerate(pat):                      # first bound variable in pattern order, as run() does
            if var in bd: return pos_idx[(i, bd[var])]
        return None
    def fill(bd, todo, chosen, used):
        """Assign relations to the patterns in todo (a list of pattern indices), oldest first, with backtracking."""
        if not todo: return chosen, bd
        # next pattern: the first in order that shares a bound variable (left sides are connected)
        j = next((j for j in todo if any(v in bd for v in P[j])), todo[0])
        rest = [x for x in todo if x != j]
        c = cands(bd, P[j])
        pool = sorted(c) if c is not None else sorted(rel)
        taken = set(chosen.values())
        for r in pool:
            if r in used or r in taken: continue
            b2 = extend(bd, P[j], rel[r])
            if b2 is None: continue
            ch = dict(chosen); ch[j] = r
            res = fill(b2, rest, ch, used)
            if res is not None: return res
        return None
    for g in range(1, max_gen + 1):
        used = set(); events = []
        for r in sorted(rel):
            if r in used: continue
            best = None
            for first in range(n):
                bd = extend({}, P[first], rel[r])
                if bd is None: continue
                res = fill(bd, [j for j in range(n) if j != first], {first: r}, used)
                if res is None: continue
                chosen, b2 = res
                ids = tuple(chosen[j] for j in range(n))
                if best is None or max(ids) < max(best[0]): best = (ids, b2)
            if best is not None:
                used.update(best[0]); events.append(best)
        if not events: break
        for ids, bd in events:
            ev = len(parents)
            parents.append(tuple(sorted({creator[i] for i in ids}))); gens.append(g)
            for i in ids: remove(i)
            fresh = {}
            def node_of(v):
                nonlocal node
                if v < k: return bd[v]
                if v not in fresh: fresh[v] = node; node += 1
                return fresh[v]
            for t in R: add(tuple(node_of(v) for v in t), ev)
        if len(parents) - 1 >= max_events or len(rel) > max_rel or time.time() - t_start > max_sec: break
    return parents, gens, len(rel), node

def left_sides(n_left, arity):
    out = []
    size = n_left * arity
    for vals in itertools.product(range(size), repeat=size):
        k = len(set(vals))
        if not first_appearance_ok(vals, 0, k): continue
        L = tuple(tuple(vals[i * arity:(i + 1) * arity]) for i in range(n_left))
        if connected_general(L): out.append((L, k))
    return out

def sample_general(n_left, n_right, arity, max_new=2, count=1000, seed=42, lefts=None):
    """Random canonical rules: a connected left side of n_left relations of the given arity, a connected right side of
    n_right relations over the left side's variables and up to max_new new nodes that it uses. Uniform over left sides,
    then over the number of new nodes, then over right sides, before renamings are merged."""
    rnd = random.Random(seed)
    lefts = lefts or left_sides(n_left, arity)
    seen = {}; tries = 0
    while len(seen) < count and tries < count * 500:
        tries += 1
        L, k = rnd.choice(lefts); m = rnd.randint(0, max_new)
        vals = [rnd.randrange(k + m) for _ in range(n_right * arity)]
        if m and not first_appearance_ok(vals, k, m): continue
        R = tuple(tuple(vals[i * arity:(i + 1) * arity]) for i in range(n_right))
        if not connected_general(R): continue
        key = canonical(L, R, k, m)
        if key not in seen: seen[key] = (k, m)
    out = []
    for (L, R), (k, m) in seen.items():
        bits = (sum(len(r) for r in L) + sum(len(r) for r in R)) * math.log2(k + m)
        out.append({'rule': rule_str(L, R), 'L': L, 'R': R, 'k': k, 'm': m, 'bits': round(bits, 2)})
    out.sort(key=lambda r: (r['bits'], r['rule']))
    return out

def restrict_to_past(parents, y):
    """The sub-graph of y's past and y, relabeled in id order (which keeps causes first)."""
    keep = set(); stack = [y]
    while stack:
        e = stack.pop()
        if e in keep: continue
        keep.add(e); stack.extend(parents[e])
    order = sorted(keep); idx = {e: i for i, e in enumerate(order)}
    return [tuple(idx[p] for p in parents[e]) for e in order], order

def past_size(parents, y):
    seen = set(); stack = [y]
    while stack:
        e = stack.pop()
        if e in seen: continue
        seen.add(e); stack.extend(parents[e])
    return len(seen)

def sprinkle_ds(s, n, seed=1, d=4):
    """Conformal time drawn from its exact density (cross-section volume times eta^-4) by inverting a fine cumulative
    table, then a uniform point in that slice's ball."""
    import bisect
    rnd = random.Random(seed); e1, e2 = -s, -1.0; pts = []
    K = 20000; grid = [e1 + (e2 - e1) * (i + 0.5) / K for i in range(K)]
    w = [min(g - e1, e2 - g) ** (d - 1) / g ** 4 for g in grid]
    cum = []; acc = 0.0
    for v in w: acc += v; cum.append(acc)
    while len(pts) < n:
        i = bisect.bisect_left(cum, rnd.random() * acc)
        eta = grid[i] + (rnd.random() - 0.5) * (e2 - e1) / K
        r = min(eta - e1, e2 - eta)
        if r <= 0: continue
        while True:
            x = [rnd.uniform(-r, r) for _ in range(d - 1)]
            if sum(v * v for v in x) <= r * r: break
        pts.append((eta, x))
    pts.sort(key=lambda p: p[0])
    pts = [(e1, [0.0] * (d - 1))] + pts + [(e2, [0.0] * (d - 1))]
    m = len(pts); anc = [0] * m; desc = [0] * m
    for j in range(m):
        tj, xj = pts[j]
        for i in range(j):
            ti, xi = pts[i]
            if tj - ti >= math.sqrt(sum((a - b) ** 2 for a, b in zip(xi, xj))):
                anc[j] |= 1 << i; desc[i] |= 1 << j
    return list(range(m)), anc, desc


def desitter(out_path=None):
    """What the equations give for an expanding universe: sprinkled intervals of four-dimensional de Sitter space spanning
    0 to 4 e-folds of expansion, read by the three tests (N = 800, three sprinklings each)."""
    out = []
    for s in (1.01, math.e, math.e ** 2, math.e ** 3, math.e ** 4):
        vals = []
        for seed in (1, 2, 3):
            members, anc, desc = sprinkle_ds(s, 800, seed)
            n, f, dmm, dmid, L, dch = diamond_dims(members, anc, desc)
            vals.append((dmm, dmid, L, dch))
        mean = lambda i: sum(v[i] for v in vals) / 3
        out.append({'expansion': round(s, 3), 'efolds': round(math.log(s), 2), 'n': 800, 'mm': round(mean(0), 3),
                    'mid': round(mean(1), 3), 'chain': round(mean(2), 2), 'chain_d': round(mean(3), 3)})
        print(out[-1], flush=True)
    if out_path: json.dump(out, open(out_path, 'w'), indent=1)
    return out

def run_measure_general(rule, events=3000, sec=8.0, long=False):
    """Run a rule of any signature with run_general (start: one self-loop per left-side relation) and read it. A long run
    is read inside the past of a top update, the deepest one whose past fits in memory (70,000 updates at most)."""
    L = tuple(map(tuple, rule['L'])); R = tuple(map(tuple, rule['R'])); k = rule['k']
    t0 = time.time()
    parents, gens, nrel, nnodes = run_general(L, R, k, max_events=events, max_sec=sec,
                                              max_rel=10 ** 7 if long else 60000, max_gen=4000 if long else 400)
    c, fit = growth(gens)
    out = {'rule': rule['rule'], 'bits': rule['bits'], 'gens': max(gens), 'relations': nrel, 'nodes': nnodes,
           'per_gen': c[:60], 'fit': fit, 'run_sec': round(time.time() - t0, 2)}
    if len(parents) - 1 < 300:
        out.update({'events': len(parents) - 1, 'status': 'stops or stays small'})
    elif not long:
        out.update(measure(parents)); out['status'] = 'measured'
    else:
        E = len(parents); depth = [0] * E
        for e in range(1, E): depth[e] = 1 + max(depth[p] for p in parents[e])
        md = max(depth); cands = [e for e in range(E) if depth[e] >= md - 2][-40:]
        sizes = {e: past_size(parents, e) for e in cands}; CAP = 70000
        ok = [e for e in cands if sizes[e] <= CAP]
        if ok: y = max(ok, key=lambda e: sizes[e])
        else:
            y = None
            for d0 in range(md - 3, 0, -1):
                fit_ = [(past_size(parents, e), e) for e in [e for e in range(E) if depth[e] == d0][-20:]]
                fit_ = [x for x in fit_ if x[0] <= CAP]
                if fit_: y = max(fit_)[1]; break
        sub, order = restrict_to_past(parents, y)
        out.update({'events': E - 1, 'top': y, 'past_of_top': len(sub)})
        global subsample_slope
        full = subsample_slope
        subsample_slope = lambda members, anc, **kw: full(members, anc, ps=(1.0, 0.5, 0.25), reps=1)   # lighter check
        try:
            m = measure(sub, n_diamonds=1, series=True)
        finally:
            subsample_slope = full
        out.update({k2: v for k2, v in m.items() if k2 != 'events'}); out['status'] = 'measured'
    out['sec'] = round(time.time() - t0, 2)
    return out

def run_and_measure(rule, events=2000, series=False, max_sec=8.0):
    L, R, k = tuple(map(tuple, rule['L'])), tuple(map(tuple, rule['R'])), rule['k']
    t0 = time.time()
    parents, gens, nrel, nnodes = run(L, R, k, max_events=events, max_sec=max_sec)
    c, fit = growth(gens)
    out = {'rule': rule['rule'], 'bits': rule['bits'], 'gens': max(gens), 'relations': nrel, 'nodes': nnodes,
           'per_gen': c[:60], 'fit': fit}
    if len(parents) - 1 < 300:
        out.update({'events': len(parents) - 1, 'status': 'stops or stays small'})
    else:
        out.update(measure(parents, series=series)); out['status'] = 'measured'
    out['sec'] = round(time.time() - t0, 2)
    return out

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'validate':
        res = validate()
        if len(sys.argv) > 2: json.dump(res, open(sys.argv[2], 'w'), indent=1)
    elif cmd == 'enumerate':
        rules = enumerate_rules()
        json.dump(rules, open(sys.argv[2], 'w'))
        print(len(rules), 'canonical rules; bits from', rules[0]['bits'], 'to', rules[-1]['bits'])
    elif cmd == 'sample':
        rules = sample_rules(count=int(sys.argv[3]) if len(sys.argv) > 3 else 1500)
        json.dump(rules, open(sys.argv[2], 'w')); print(len(rules), 'sampled rules')
    elif cmd == 'run':
        rules = json.load(open(sys.argv[2])); outp = sys.argv[3]
        ev = int(sys.argv[sys.argv.index('--events') + 1]) if '--events' in sys.argv else 2000
        a = int(sys.argv[sys.argv.index('--start') + 1]) if '--start' in sys.argv else 0
        b = int(sys.argv[sys.argv.index('--stop') + 1]) if '--stop' in sys.argv else len(rules)
        sec = float(sys.argv[sys.argv.index('--sec') + 1]) if '--sec' in sys.argv else 8.0
        with open(outp, 'a') as fh:
            for r in rules[a:b]:
                fh.write(json.dumps(run_and_measure(r, ev, series='--series' in sys.argv, max_sec=sec)) + '\n'); fh.flush()
    elif cmd == 'sample2':          # sample2 OUT.json n_left n_right arity count seed
        nl, nr, ar, cnt, seed = map(int, sys.argv[3:8])
        rules = sample_general(nl, nr, ar, count=cnt, seed=seed)
        json.dump(rules, open(sys.argv[2], 'w')); print(len(rules), 'sampled rules')
    elif cmd in ('run2', 'long'):   # run2|long RULES.json OUT.jsonl [--events N] [--start i --stop j] [--sec S]
        rules = json.load(open(sys.argv[2])); outp = sys.argv[3]
        opt = lambda name, d, f=int: f(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else d
        ev = opt('--events', 3000); a = opt('--start', 0); b = opt('--stop', len(rules)); sec = opt('--sec', 8.0, float)
        with open(outp, 'a') as fh:
            for r in rules[a:b]:
                fh.write(json.dumps(run_measure_general(r, ev, sec, long=(cmd == 'long'))) + '\n'); fh.flush()
    elif cmd == 'desitter':
        desitter(sys.argv[2] if len(sys.argv) > 2 else None)
    elif cmd == 'one':
        L, R, k = parse_rule(sys.argv[2])
        m = len({v for rel in R for v in rel} - set(range(k)))
        rule = {'rule': rule_str(L, R), 'L': L, 'R': R, 'k': k, 'bits': round(2 * (len(L) + len(R)) * math.log2(k + m), 2)}
        ev = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
        print(json.dumps(run_and_measure(rule, ev, series=True, max_sec=120.0)))
