#!/usr/bin/env python3
"""A42 SA74 extensions: noise bands from matched-size sprinklings, other update orders, transitive percolation,
links per element, and a critical-pair check of causal invariance. Builds on a42_rule_search.py (the A42 supplement,
unchanged) and uses its interval choice and its three tests, so every reading here is made the same way as in the
first passes. Requires numpy.

usage:
  python3 a42x.py bands OUT.json                      # sprinkled flat intervals, d = 2, 3, 4, the protocol's sizes
  python3 a42x.py orders RULES.json OUT.jsonl ORDER [--start i --stop j]   # ORDER: oldest | newest | random1 | random2
  python3 a42x.py tp OUT.jsonl N p seed               # transitive percolation, one causal set, read as a rule's graph
  python3 a42x.py ci RULES.json OUT.jsonl [--steps 2] # critical-pair check of causal invariance
  python3 a42x.py links2d OUT.json                    # links per element in sprinkled two-dimensional intervals
"""
import sys, os, json, math, random, time, itertools
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a42_rule_search as A

# ---------------------------------------------------------------- sprinklings with numpy
def sprinkle_np(d, n, seed=1, return_rel=False):
    """n points uniform in the causal interval |x| <= 1/2 - |t| of d-dimensional Minkowski space, plus the two tips.
    Returns members, anc, desc (Python int bitsets, as a42_rule_search uses), and optionally the relation matrix."""
    rng = np.random.default_rng(seed)
    t_list = []; x_list = []; have = 0
    while have < n:
        m = max(1000, 3 * (n - have) * (2 ** (d - 1)))
        t = rng.uniform(-0.5, 0.5, m); x = rng.uniform(-0.5, 0.5, (m, d - 1))
        ok = np.sqrt((x * x).sum(axis=1)) <= 0.5 - np.abs(t)
        t_list.append(t[ok]); x_list.append(x[ok]); have += int(ok.sum())
    t = np.concatenate(t_list)[:n]; x = np.concatenate(x_list)[:n]
    o = np.argsort(t, kind='stable'); t = t[o]; x = x[o]
    t = np.concatenate([[-0.5], t, [0.5]]); x = np.vstack([np.zeros((1, d - 1)), x, np.zeros((1, d - 1))])
    M = n + 2
    rel = np.zeros((M, M), dtype=bool)
    for a in range(0, M, 512):
        b = min(M, a + 512)
        dt = t[None, :] - t[a:b, None]
        dx2 = ((x[None, :, :] - x[a:b, None, :]) ** 2).sum(axis=2)
        rel[a:b] = (dt > 0) & (dt * dt >= dx2)
    anc = bits_from_rows(rel.T); desc = bits_from_rows(rel)
    if return_rel: return list(range(M)), anc, desc, rel
    return list(range(M)), anc, desc

def bits_from_rows(mat):
    mat = np.ascontiguousarray(mat)
    packed = np.packbits(mat, axis=1, bitorder='little')
    return [int.from_bytes(packed[i].tobytes(), 'little') for i in range(mat.shape[0])]

def read_interval(members, anc, desc, coarse=True):
    n, f, d_mm, d_mid, L, d_ch = A.diamond_dims(members, anc, desc)
    out = {'n': n, 'f': f, 'mm': d_mm, 'mid': d_mid, 'chain': L, 'chain_d': d_ch}
    if coarse: out['coarse_d'] = A.subsample_slope(members, anc)[0]
    return out


def sprinkle_ds_np(efolds, n, seed=1, d=4):
    """Sprinkled interval of d-dimensional de Sitter space across which the universe grows by `efolds` e-folds:
    conformal time from -s to -1 with s = e^efolds, density proportional to the cross-section volume times eta^-d,
    uniform in each slice's ball; conformal flatness keeps the light cones of flat space. Numpy version of
    a42_rule_search.sprinkle_ds."""
    rng = np.random.default_rng(seed)
    s = math.exp(efolds) if efolds > 0 else 1.0001
    e1, e2 = -s, -1.0
    K = 200000
    grid = e1 + (e2 - e1) * (np.arange(K) + 0.5) / K
    w = np.minimum(grid - e1, e2 - grid) ** (d - 1) / np.abs(grid) ** d
    cum = np.cumsum(w); cum /= cum[-1]
    idx = np.searchsorted(cum, rng.random(n))
    eta = grid[idx] + (rng.random(n) - 0.5) * (e2 - e1) / K
    r = np.minimum(eta - e1, e2 - eta)
    v = rng.normal(size=(n, d - 1)); v /= np.linalg.norm(v, axis=1)[:, None]
    rad = r * rng.random(n) ** (1.0 / (d - 1))
    x = v * rad[:, None]
    o = np.argsort(eta, kind='stable'); eta = eta[o]; x = x[o]
    t = np.concatenate([[e1], eta, [e2]]); x = np.vstack([np.zeros((1, d - 1)), x, np.zeros((1, d - 1))])
    M = n + 2
    rel = np.zeros((M, M), dtype=bool)
    for a in range(0, M, 512):
        b = min(M, a + 512)
        dt = t[None, :] - t[a:b, None]
        dx2 = ((x[None, :, :] - x[a:b, None, :]) ** 2).sum(axis=2)
        rel[a:b] = (dt > 0) & (dt * dt >= dx2)
    return list(range(M)), bits_from_rows(rel.T), bits_from_rows(rel)

def ds_table(out_path):
    res = []
    for n in (1000, 3000, 10000):
        for k in (0, 1, 2, 3, 4, 6, 8):
            vals = []
            for seed in (1, 2):
                m, a, de = sprinkle_ds_np(k, n, seed)
                vals.append(read_interval(m, a, de, coarse=False))
            row = {'n': n, 'efolds': k}
            for key in ('mm', 'mid', 'chain_d'): row[key] = round(sum(v[key] for v in vals) / len(vals), 3)
            row['chain'] = round(sum(v['chain'] for v in vals) / len(vals), 1)
            res.append(row); print(json.dumps(row), flush=True)
            json.dump(res, open(out_path, 'w'), indent=1)
    return res

SIZES = (120, 250, 500, 1000, 2000, 4000)
REPS = {120: 40, 250: 40, 500: 40, 1000: 20, 2000: 12, 4000: 6}

def bands(out_path):
    res = {}
    for d in (4, 2, 3):
        for n in SIZES:
            vals = []
            t0 = time.time()
            for seed in range(1, REPS[n] + 1):
                members, anc, desc = sprinkle_np(d, n, seed)
                vals.append(read_interval(members, anc, desc, coarse=(n <= 2000)))
            row = {'d': d, 'n': n, 'reps': len(vals)}
            for key in ('mm', 'mid', 'chain_d', 'coarse_d'):
                v = [x[key] for x in vals if key in x and math.isfinite(x[key])]
                if v:
                    row[key + '_mean'] = round(float(np.mean(v)), 4)
                    row[key + '_sd'] = round(float(np.std(v, ddof=1)), 4)
            row['chain_mean_len'] = round(float(np.mean([x['chain'] for x in vals])), 2)
            row['sec'] = round(time.time() - t0, 1)
            res[f'{d}_{n}'] = row
            print(json.dumps(row), flush=True)
            json.dump(res, open(out_path, 'w'), indent=1)
    return res

def links2d(out_path):
    """Links per element (mean over interior elements) in sprinkled intervals, two and four dimensions."""
    res = []
    for d in (2, 4):
        for n in (500, 1000, 2000, 4000):
            reps = 3 if n <= 2000 else 1
            ms = []; mx = []
            for seed in range(1, reps + 1):
                members, anc, desc, rel = sprinkle_np(d, n, seed, return_rel=True)
                R = rel.astype(np.float32)
                two = (R @ R) > 0.5
                link = rel & ~two
                deg = link.sum(axis=0) + link.sum(axis=1)
                inner = deg[1:-1]
                ms.append(float(inner.mean())); mx.append(int(inner.max()))
            row = {'d': d, 'n': n, 'reps': reps, 'mean_links': round(float(np.mean(ms)), 2), 'max_links': max(mx)}
            res.append(row); print(json.dumps(row), flush=True)
            json.dump(res, open(out_path, 'w'), indent=1)
    return res

# ---------------------------------------------------------------- other update orders
def run_order(L, R, k, order='oldest', max_events=20000, max_gen=4000, max_rel=10 ** 7, init=((0, 0), (0, 0)), max_sec=900.0):
    """As a42_rule_search.run, with the scan order and the choice of partner set by `order`:
    oldest: the first passes (oldest relation first, oldest compatible partner, the match whose newest relation is oldest);
    newest: newest relation first, newest compatible partner, the match whose oldest relation is newest;
    random1, random2: relations scanned in a random order, partners tried in a random order (seeded)."""
    t_start = time.time()
    rnd = random.Random({'random1': 1, 'random2': 2}.get(order, 0))
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
    def ordered(ids):
        if order == 'oldest': return sorted(ids)
        if order == 'newest': return sorted(ids, reverse=True)
        lst = list(ids); rnd.shuffle(lst); return lst
    for g in range(1, max_gen + 1):
        used = set(); events = []
        for r in ordered(rel):
            if r in used: continue
            best = None
            firsts = list(range(len(P)))
            if order.startswith('random'): rnd.shuffle(firsts)
            for first in firsts:
                bd = extend({}, P[first], *rel[r])
                if bd is None: continue
                other = 1 - first
                cand = candidates(bd, P[other])
                for c in ordered(cand):
                    if c == r or c in used: continue
                    bd2 = extend(bd, P[other], *rel[c])
                    if bd2 is None: continue
                    ids = (r, c) if first == 0 else (c, r)
                    if best is None: best = (ids, bd2)
                    elif order == 'oldest' and max(ids) < max(best[0]): best = (ids, bd2)
                    elif order == 'newest' and min(ids) > min(best[0]): best = (ids, bd2)
                    break
                if best is not None and order.startswith('random'): break
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

def link_stats(parents):
    """Links per update in a rule's causal graph. A direct cause p of e is a link unless p precedes another direct
    cause q of e. Returns mean and largest number of links per update (root excluded), and the in/out bounds."""
    E = len(parents)
    anc, _ = A.from_parents(parents)
    deg = [0] * E
    for e in range(1, E):
        ps = parents[e]
        for p in ps:
            if any(q != p and (anc[q] >> p) & 1 for q in ps): continue
            deg[e] += 1; deg[p] += 1
    inner = deg[1:]
    return {'mean_links': round(sum(inner) / len(inner), 3), 'max_links': max(inner),
            'max_causes': max(len(p) for p in parents[1:])}

def orders_run(rules_path, out_path, order, a=0, b=None, events=20000):
    rules = json.load(open(rules_path))
    b = len(rules) if b is None else b
    with open(out_path, 'a') as fh:
        for rule in rules[a:b]:
            L, Rr, k = A.parse_rule(rule['rule'])
            t0 = time.time()
            parents, gens, nrel, nnodes = run_order(L, Rr, k, order=order, max_events=events)
            c, fit = A.growth(gens)
            out = {'rule': rule['rule'], 'group': rule.get('group'), 'order': order, 'gens': max(gens),
                   'relations': nrel, 'per_gen': c[:60], 'fit': fit, 'run_sec': round(time.time() - t0, 2)}
            if len(parents) - 1 < 300:
                out.update({'events': len(parents) - 1, 'status': 'stops or stays small'})
            else:
                out.update(A.measure(parents, series=True)); out['status'] = 'measured'
                out.update(link_stats(parents))
                E = len(parents); depth = [0] * E
                for e in range(1, E): depth[e] = 1 + max(depth[p] for p in parents[e])
                hist = defaultdict(int)
                for e in range(1, E): hist[depth[e]] += 1
                out['depth_hist'] = [hist[i] for i in range(1, min(max(depth), 400) + 1)]
            out['sec'] = round(time.time() - t0, 2)
            fh.write(json.dumps(out) + '\n'); fh.flush()
            print(order, rule['rule'], out.get('events'), out['sec'], flush=True)

# ---------------------------------------------------------------- transitive percolation
def transitive_percolation(N, p, seed=1):
    """Transitive percolation (Rideout and Sorkin's simplest classical sequential growth): elements 1..N are born in
    order, and each earlier element becomes a direct cause of the newborn with probability p; the order is the
    transitive closure. Element 0 is a root added below all of them, as the bottom of the largest intervals."""
    rng = np.random.default_rng(seed)
    parents = [()]
    for j in range(N):
        if j:
            hits = np.nonzero(rng.random(j) < p)[0] + 1
            parents.append(tuple([0] + hits.tolist()))
        else:
            parents.append((0,))
    return parents

def originary_percolation(N, p, seed=1):
    """Originary percolation (Ahmed and Rideout's model): element 0 is the origin; each newborn element picks each
    earlier element as a direct cause with probability p and tries again if it picked none, so every element lies
    in the origin's future."""
    rng = np.random.default_rng(seed)
    parents = [()]
    for j in range(1, N + 1):
        while True:
            hits = np.nonzero(rng.random(j) < p)[0]
            if len(hits): break
        parents.append(tuple(hits.tolist()))
    return parents

def tp_run(out_path, N, p, seed, model='tp'):
    t0 = time.time()
    if model == 'op':
        parents = originary_percolation(N, p, seed)
        out = {'model': 'originary percolation', 'N': N, 'p': p, 'seed': seed, 'pN': round(p * N, 3)}
    else:
        parents = transitive_percolation(N, p, seed)
        out = {'model': 'transitive percolation', 'N': N, 'p': p, 'seed': seed, 'pN': round(p * N, 3)}
    m = A.measure(parents, series=True)
    out.update({k: v for k, v in m.items()})
    E = len(parents)
    anc, desc = A.from_parents(parents)
    tot = sum(a.bit_count() for a in anc[1:]) - (E - 1)        # related pairs among the N elements (root excluded)
    out['whole_f'] = round(tot / (N * (N - 1) / 2), 5)
    out['whole_mm'] = round(A.mm_dimension(out['whole_f']), 3)
    # posts: elements related to every other element
    posts = [e for e in range(1, E) if (anc[e].bit_count() - 1) + desc[e].bit_count() == N - 1]
    out['posts'] = len(posts)
    out['sec'] = round(time.time() - t0, 2)
    with open(out_path, 'a') as fh: fh.write(json.dumps(out) + '\n')
    return out


def root_mm(parents):
    """Myrheim-Meyer reading of the largest interval: from the root to the element with the largest past."""
    anc, desc = A.from_parents(parents)
    E = len(parents)
    y = max(range(1, E), key=lambda e: anc[e].bit_count())
    members = [e for e in range(E) if (anc[y] >> e) & 1] + [y]
    mask = 0
    for e in members: mask |= 1 << e
    n = len(members)
    rel = sum((anc[e] & mask).bit_count() for e in members)
    return A.mm_dimension(rel / (n * (n - 1) / 2)), n

def tune(out_path, N, target=4.0, model='tp', seed=101):
    """Exploratory, outside the protocol: find p so that the largest interval reads `target` on the ordering
    fraction, as Rideout and Sorkin's trajectories hold the ordering fraction fixed, by bisection in log p."""
    gen = originary_percolation if model == 'op' else transitive_percolation
    lo, hi = math.log(0.5 / N), math.log(200.0 / N)        # sparse (high reading) to dense (low reading)
    hist = []
    for _ in range(12):
        mid = (lo + hi) / 2; p = math.exp(mid)
        d, n = root_mm(gen(N, p, seed))
        hist.append((round(p, 7), round(d, 3), n))
        if d > target: lo = mid
        else: hi = mid
    p = math.exp((lo + hi) / 2)
    res = {'model': model, 'N': N, 'target': target, 'p': p, 'pN': p * N, 'steps': hist}
    with open(out_path, 'a') as fh: fh.write(json.dumps(res) + '\n')
    return res

# ---------------------------------------------------------------- causal invariance by critical pairs
def matches(state, L):
    """All matches of a two-relation left side in a state (a tuple of relations, indexed by position).
    Returns a list of (ids, binding); ids are distinct positions; variables may share a node."""
    out = []
    n = len(state)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            bd = {}
            ok = True
            for pat, rel in ((L[0], state[i]), (L[1], state[j])):
                for var, val in zip(pat, rel):
                    if var in bd and bd[var] != val: ok = False; break
                    bd[var] = val
                if not ok: break
            if ok: out.append(((i, j), bd))
    return out

def apply(state, m, R, k, next_node):
    ids, bd = m
    rest = [r for x, r in enumerate(state) if x not in ids]
    fresh = {}
    def node_of(v):
        nonlocal next_node
        if v < k: return bd[v]
        if v not in fresh: fresh[v] = next_node; next_node += 1
        return fresh[v]
    new = [(node_of(a), node_of(b)) for a, b in R]
    return tuple(rest + new), next_node

def canon(state, fixed):
    """Canonical form of a state up to renaming the nodes not in `fixed` (the nodes of the overlap state)."""
    nodes = sorted({v for r in state for v in r})
    free = [v for v in nodes if v not in fixed]
    best = None
    if len(free) > 7:                                   # too many to permute; fall back to a coarse invariant
        return ('coarse', tuple(sorted(state)))
    for perm in itertools.permutations(range(len(free))):
        mp = {v: 10 ** 6 + perm[i] for i, v in enumerate(free)}
        key = tuple(sorted((mp.get(a, a), mp.get(b, b)) for a, b in state))
        if best is None or key < best: best = key
    return best

def overlap_states(L):
    """Minimal overlap states: two matches of the left side sharing at least one relation. Built by taking the
    left side as a pattern twice, identifying one or both relations of the second copy with relations of the first,
    and then every way of identifying the remaining variables (set partitions), which covers non-injective matches."""
    k = len({v for r in L for v in r})
    states = set()
    for a_ids in itertools.product(range(2), repeat=2):
        pass
    # variables of copy 1: 0..k-1; of copy 2: k..2k-1
    for share in (((0, 0),), ((0, 1),), ((1, 0),), ((1, 1),), ((0, 0), (1, 1)), ((0, 1), (1, 0))):
        # share: pairs (relation index in copy 2, relation index in copy 1) that are the same relation
        eq = []
        for r2, r1 in share:
            for u, v in zip(L[r2], L[r1]): eq.append((k + u, v))
        nvars = 2 * k
        for part in set_partitions(nvars):
            cls = {}
            for ci, block in enumerate(part):
                for v in block: cls[v] = ci
            if any(cls[u] != cls[v] for u, v in eq): continue
            rels1 = [(cls[a], cls[b]) for a, b in L]
            rels2 = [(cls[k + a], cls[k + b]) for a, b in L]
            # the state: copy 1's relations plus copy 2's relations that are not identified with copy 1's
            shared2 = {r2 for r2, _ in share}
            state = tuple(rels1 + [rels2[i] for i in range(len(L)) if i not in shared2])
            m1 = ((0, 1), {v: cls[v] for v in range(k)})
            pos2 = []
            extra = len(L)
            for i in range(len(L)):
                if i in shared2: pos2.append(dict(share)[i])
                else: pos2.append(extra); extra += 1
            m2 = (tuple(pos2), {v: cls[k + v] for v in range(k)})
            if m1[0] == m2[0] and m1[1] == m2[1]: continue
            states.add((state, m1[0], tuple(sorted(m1[1].items())), m2[0], tuple(sorted(m2[1].items()))))
    return states

def set_partitions(n):
    """All set partitions of range(n) (restricted growth strings)."""
    def rec(i, a, mx):
        if i == n:
            blocks = defaultdict(list)
            for v, c in enumerate(a): blocks[c].append(v)
            yield list(blocks.values()); return
        for c in range(mx + 2):
            a.append(c)
            yield from rec(i + 1, a, max(mx, c))
            a.pop()
    yield from rec(0, [], -1)

def reach(state, L, R, k, next_node, steps, fixed, limit=4000):
    """Canonical forms of all states reachable in 0..steps single updates."""
    seen = {canon(state, fixed)}
    frontier = [(state, next_node)]
    for _ in range(steps):
        nxt = []
        for s, nn in frontier:
            for m in matches(s, L):
                s2, nn2 = apply(s, m, R, k, nn)
                c = canon(s2, fixed)
                if c not in seen:
                    seen.add(c); nxt.append((s2, nn2))
                if len(seen) > limit: return seen, True
        frontier = nxt
    return seen, False

def critical_pairs(rule, steps=2):
    L, R, k = A.parse_rule(rule)
    res = {'rule': rule, 'pairs': 0, 'resolved': 0, 'unresolved_example': None, 'capped': 0}
    for state, ids1, b1, ids2, b2 in overlap_states(L):
        fixed = {v for r in state for v in r}
        nn = max(fixed) + 1
        s1, nn1 = apply(state, (ids1, dict(b1)), R, k, nn)
        s2, nn2 = apply(state, (ids2, dict(b2)), R, k, nn1)
        res['pairs'] += 1
        if canon(s1, fixed) == canon(s2, fixed): res['resolved'] += 1; continue
        A1, c1 = reach(s1, L, R, k, nn2, steps, fixed)
        A2, c2 = reach(s2, L, R, k, nn2 + 50, steps, fixed)
        if c1 or c2: res['capped'] += 1
        if A1 & A2: res['resolved'] += 1
        elif res['unresolved_example'] is None:
            res['unresolved_example'] = {'state': state, 'm1': ids1, 'm2': ids2}
    res['causally_invariant_within_steps'] = (res['resolved'] == res['pairs'])
    res['steps'] = steps
    return res

def ci_run(rules_path, out_path, steps=2, a=0, b=None):
    rules = json.load(open(rules_path))
    b = len(rules) if b is None else b
    with open(out_path, 'a') as fh:
        for rule in rules[a:b]:
            t0 = time.time()
            r = critical_pairs(rule['rule'], steps)
            r['sec'] = round(time.time() - t0, 2)
            fh.write(json.dumps(r) + '\n'); fh.flush()


# ---------------------------------------------------------------- strict critical pairs: relation instances tracked
def matches_ids(state, L):
    out = []
    n = len(state)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            bd = {}; ok = True
            for pat, rel in ((L[0], state[i]), (L[1], state[j])):
                for var, val in zip(pat, rel[:2]):
                    if var in bd and bd[var] != val: ok = False; break
                    bd[var] = val
                if not ok: break
            if ok: out.append(((i, j), bd))
    return out

def apply_ids(state, m, R, k, next_node):
    ids, bd = m
    rest = [r for x, r in enumerate(state) if x not in ids]
    fresh = {}
    def node_of(v):
        nonlocal next_node
        if v < k: return bd[v]
        if v not in fresh: fresh[v] = next_node; next_node += 1
        return fresh[v]
    new = [(node_of(a), node_of(b), -1) for a, b in R]
    return tuple(rest + new), next_node

def canon_ids(state, fixed):
    nodes = sorted({v for r in state for v in r[:2]})
    free = [v for v in nodes if v not in fixed]
    if len(free) > 7: return ('coarse', tuple(sorted(state)))
    best = None
    for perm in itertools.permutations(range(len(free))):
        mp = {v: 10 ** 6 + perm[i] for i, v in enumerate(free)}
        key = tuple(sorted((mp.get(a, a), mp.get(b, b), rid) for a, b, rid in state))
        if best is None or key < best: best = key
    return best

def reach_ids(state, L, R, k, next_node, steps, fixed, limit=6000):
    seen = {canon_ids(state, fixed)}; frontier = [(state, next_node)]
    for _ in range(steps):
        nxt = []
        for s, nn in frontier:
            for m in matches_ids(s, L):
                s2, nn2 = apply_ids(s, m, R, k, nn)
                c = canon_ids(s2, fixed)
                if c not in seen:
                    seen.add(c); nxt.append((s2, nn2))
                if len(seen) > limit: return seen, True
        frontier = nxt
    return seen, False

def critical_pairs_strict(rule, steps=2):
    """As critical_pairs, but each relation of the overlap state keeps its identity (the event that made it): two
    branches resolve only if they reach the same state with the same surviving original relations, so consuming a
    different one of two parallel relations counts as a different history, as it gives a different causal graph."""
    L, R, k = A.parse_rule(rule)
    res = {'rule': rule, 'pairs': 0, 'resolved': 0, 'capped': 0, 'unresolved_example': None}
    for state, ids1, b1, ids2, b2 in overlap_states(L):
        st = tuple((a, b, i) for i, (a, b) in enumerate(state))
        fixed = {v for r in state for v in r}
        nn = max(fixed) + 1
        s1, nn1 = apply_ids(st, (ids1, dict(b1)), R, k, nn)
        s2, nn2 = apply_ids(st, (ids2, dict(b2)), R, k, nn1)
        res['pairs'] += 1
        if canon_ids(s1, fixed) == canon_ids(s2, fixed): res['resolved'] += 1; continue
        A1, c1 = reach_ids(s1, L, R, k, nn2, steps, fixed)
        A2, c2 = reach_ids(s2, L, R, k, nn2 + 50, steps, fixed)
        if c1 or c2: res['capped'] += 1
        if A1 & A2: res['resolved'] += 1
        elif res['unresolved_example'] is None: res['unresolved_example'] = {'state': state, 'm1': ids1, 'm2': ids2}
    res['causally_invariant_within_steps'] = res['resolved'] == res['pairs']
    res['steps'] = steps
    return res

def ci_strict_run(rules_path, out_path, steps=2, a=0, b=None):
    rules = json.load(open(rules_path)); b = len(rules) if b is None else b
    with open(out_path, 'a') as fh:
        for rule in rules[a:b]:
            t0 = time.time(); r = critical_pairs_strict(rule['rule'], steps); r['sec'] = round(time.time() - t0, 2)
            fh.write(json.dumps(r) + '\n'); fh.flush()

if __name__ == '__main__':
    cmd = sys.argv[1]
    opt = lambda name, d, f=int: f(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else d
    if cmd == 'bands': bands(sys.argv[2])
    elif cmd == 'links2d': links2d(sys.argv[2])
    elif cmd == 'orders': orders_run(sys.argv[2], sys.argv[3], sys.argv[4], opt('--start', 0), opt('--stop', None),
                                     opt('--events', 20000))
    elif cmd == 'tp':
        N, p, seed = int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
        print(json.dumps({k: v for k, v in tp_run(sys.argv[2], N, p, seed).items() if k != 'series'}))
    elif cmd == 'op':
        N, p, seed = int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
        print(json.dumps({k: v for k, v in tp_run(sys.argv[2], N, p, seed, model='op').items() if k != 'series'}))
    elif cmd == 'tune':
        print(json.dumps(tune(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5])))
    elif cmd == 'ds': ds_table(sys.argv[2])
    elif cmd == 'cistrict': ci_strict_run(sys.argv[2], sys.argv[3], opt('--steps', 2), opt('--start', 0), opt('--stop', None))
    elif cmd == 'ci': ci_run(sys.argv[2], sys.argv[3], opt('--steps', 2), opt('--start', 0), opt('--stop', None))
