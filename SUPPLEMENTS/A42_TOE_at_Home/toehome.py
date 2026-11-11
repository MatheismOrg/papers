#!/usr/bin/env python3
"""TOE@home: a volunteer search for a short rule that grows four-dimensional spacetime.

What it does. A rule rewrites a few relations between points into a few new ones, again and again. Every rewrite is an
update, and the order of cause and effect among the updates is a candidate spacetime. This program takes a block of rules,
runs each rule for up to 3,000 updates and reads the dimension of what grew with the three tests of causal-set theory
(ordering fraction, midpoint, longest chain). Our spacetime reads four on all three. No rule tested so far does.

How it shares the work. A coordinator hands out blocks. Every block goes to two volunteers, and a block counts only when
both results agree. The results are the same on every machine: there is no time limit inside a run, only a fixed budget
of steps. The coordinator sends nothing but the label of a block (for example b2_4:3:1:400); it cannot send code, and this
program runs no code but its own two files.

Requirements: Python 3.10 or later. Standard library only; nothing to install.

  python3 toehome.py                 run until stopped (Ctrl-C); progress is kept and resumed
  python3 toehome.py --cores 4       use four processes
  python3 toehome.py --nick Ada      a label for the list of contributors (optional; no other personal data is sent)
  python3 toehome.py --server URL    use this coordinator (or put {"server": "URL"} in toehome.json beside this file)
  python3 toehome.py --offline       no network: take blocks at random and keep the result files in toehome_data/done
  python3 toehome.py --send          upload result files kept from offline runs or failed uploads
  python3 toehome.py --selftest      check that this machine reproduces a known block (about a minute)
  python3 toehome.py --unit b2_3:5:1:0   run one block and print its result

The engine is A42_Rule_Search.py, the search program published with the paper, unchanged; this file checks its SHA-256
before use. The only difference from the paper's runs: a rule that needs more than MAX_OPS matching steps is cut off at
that count, where the paper cut off at 8 seconds of computer time.
"""
import sys, os, json, math, time, hashlib, itertools, argparse, gzip, uuid, platform, random, importlib.util
import urllib.request, urllib.error
from collections import defaultdict

VERSION = 1                      # fixes the computation below; results of different versions are not compared
ENGINE = 'A42_Rule_Search.py'
ENGINE_SHA256 = '4739063aceab87cecb1054c12c603c2ea9d433a2773962bfcbc7ab05a4a2418e'
EVENTS, MAX_GEN, MAX_REL, MIN_EVENTS = 3000, 400, 60000, 300     # as in the paper's scan
MAX_OPS = 6_000_000              # matching steps per rule (replaces the paper's 8-second limit)
MAX_NEW = 2                      # new nodes a right side may use
# class: (relations on the left, relations on the right, nodes per relation)
CLASSES = {'b2_3': (2, 3, 2),    # the paper's complete class, 2,935 rules: kept as a check of every volunteer machine
           'b2_4': (2, 4, 2),    # two relations into four (the paper sampled 1,500 of these)
           'b2_5': (2, 5, 2),    # two relations into five
           'b3_4': (3, 4, 2),    # three relations into four (the paper sampled 1,000)
           'b3_5': (3, 5, 2),    # three relations into five
           't2_3': (2, 3, 3),    # two three-node relations into three (the paper sampled 1,000)
           'b4_5': (4, 5, 2),    # four relations into five
           't2_4': (2, 4, 3)}    # two three-node relations into four
ORDER = ('b2_3', 'b2_4', 'b2_5', 'b3_4', 'b3_5', 't2_3', 'b4_5', 't2_4')          # shortest rules first
SPAN = {'b2_3': 3000, 'b2_4': 3000, 'b2_5': 6000, 'b3_4': 9000, 'b3_5': 1500, 't2_3': 4500, 'b4_5': 2400, 't2_4': 1500}   # right sides per block (about five minutes)
LONG_EVENTS, LONG_GEN, LONG_REL, LONG_OPS = 20000, 4000, 1000000, 60_000_000   # second stage: one rule, run long
# four-dimensional reference: mean and standard deviation of each test on random point sets in flat spacetime, by size
BANDS4 = {120: (3.6608, 0.1554, 3.9093, 0.238, 3.9762, 0.414), 250: (3.8432, 0.1383, 4.0887, 0.2063, 3.9941, 0.226),
          500: (3.916, 0.0754, 4.1272, 0.1599, 4.0161, 0.2249), 1000: (3.9599, 0.0642, 4.12, 0.1305, 4.0214, 0.1242),
          2000: (3.9911, 0.064, 4.1358, 0.1208, 4.0812, 0.1385), 4000: (3.9949, 0.0356, 4.0745, 0.0537, 4.0323, 0.0957)}
SELFTEST = ('b2_3:5:1:0', '64d8cbdf534b7c990e7cdd7a408584a457aa0283b8dade948bc274694d289fd1')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'toehome_data')

# ------------------------------------------------------------------ the engine and the classes
_A = None
def engine():
    """The paper's search program, loaded from the file beside this one after its hash is checked."""
    global _A
    if _A is None:
        path = os.path.join(HERE, ENGINE)
        if not os.path.exists(path): sys.exit('missing %s: keep it in the same folder as toehome.py' % ENGINE)
        if hashlib.sha256(open(path, 'rb').read()).hexdigest() != ENGINE_SHA256:
            sys.exit('%s differs from the published file (SHA-256 mismatch). Download the package again.' % ENGINE)
        spec = importlib.util.spec_from_file_location('A42_Rule_Search', path)
        _A = importlib.util.module_from_spec(spec); spec.loader.exec_module(_A)
    return _A

def left_sides(nl, ar):
    """The same list as left_sides() of the paper's program (connected left sides whose node labels first appear in
    order), built directly instead of by filtering every labeling, which is slow for four relations."""
    A = engine(); size = nl * ar; out = []
    def grow(vals, k):
        if len(vals) == size:
            L = tuple(tuple(vals[i * ar:(i + 1) * ar]) for i in range(nl))
            if A.connected_general(L): out.append((L, k))
            return
        for v in range(k + 1): grow(vals + [v], max(k, v + 1))
    grow([], 0)
    return out

def unrank(n, r, index):
    """The index-th non-decreasing sequence of r numbers below n, in dictionary order (the order in which
    itertools.combinations_with_replacement lists them)."""
    seq = []; lo = 0
    for j in range(r):
        v = lo
        while True:
            cnt = math.comb(n - v + r - j - 2, r - j - 1)
            if index < cnt: break
            index -= cnt; v += 1
        seq.append(v); lo = v
    return seq

def right_sides(rels, nr, a, b):
    """(index, right side) for the indices a to b - 1 of the list of right sides."""
    n = len(rels); seq = unrank(n, nr, a)
    for i in range(a, b):
        yield i, tuple(rels[c] for c in seq)
        j = nr - 1
        while j >= 0 and seq[j] == n - 1: j -= 1
        if j < 0: return
        seq[j:] = [seq[j] + 1] * (nr - j)

_LEFTS = {}
def lefts(cls):
    """The left sides of a class, one per renaming class, each with the renamings that leave it unchanged."""
    if cls not in _LEFTS:
        nl, nr, ar = CLASSES[cls]; out = []
        for L, k in sorted(left_sides(nl, ar)):
            auts = []; ok = True
            for perm in itertools.permutations(range(k)):
                L2 = tuple(sorted(tuple(perm[v] for v in rel) for rel in L))
                if L2 < L: ok = False; break
                if L2 == L: auts.append(perm)
            if ok: out.append((L, k, tuple(auts)))
        _LEFTS[cls] = out
    return _LEFTS[cls]

def domain(cls, k, m):
    """Number of right sides over k left-side nodes and m new ones, each listed once (its relations in sorted order)."""
    nl, nr, ar = CLASSES[cls]
    return math.comb((k + m) ** ar + nr - 1, nr)

def unit_list(cls):
    out = []
    for li, (L, k, auts) in enumerate(lefts(cls)):
        for m in range(MAX_NEW + 1):
            size = domain(cls, k, m)
            out.extend((cls, li, m, a, min(size, a + SPAN[cls])) for a in range(0, size, SPAN[cls]))
    return out

def uid(u): return 'long:%s:%d:%d:%d' % tuple(u[1:5]) if u[0] == 'long' else '%s:%d:%d:%d' % tuple(u[:4])

def parse_uid(s):
    """A block label, checked. class:left:new_nodes:start is a block of rules run to 3,000 updates;
    long:class:left:new_nodes:index is one rule run to 20,000. Anything else is refused."""
    try:
        parts = str(s).split(':'); long = parts[0] == 'long'
        if long: parts = parts[1:]
        cls, li, m, a = parts; li, m, a = int(li), int(m), int(a)
    except Exception:
        raise ValueError('bad block label')
    if cls not in CLASSES or not 0 <= m <= MAX_NEW or li < 0 or a < 0: raise ValueError('bad block label')
    ls = lefts(cls)
    if li >= len(ls): raise ValueError('bad block label')
    size = domain(cls, ls[li][1], m)
    if a >= size: raise ValueError('bad block label')
    if long:
        L, k, auts = ls[li]
        i, R = next(right_sides(list(itertools.product(range(k + m), repeat=CLASSES[cls][2])), CLASSES[cls][1], a, a + 1))
        if not accepted(R, k, m, auts): raise ValueError('bad block label')
        return ('long', cls, li, m, a)
    if a % SPAN[cls]: raise ValueError('bad block label')
    return (cls, li, m, a, min(size, a + SPAN[cls]))

def accepted(R, k, m, auts):
    """True if this right side, with the block's left side, is the one listed form of its rule: it uses every new node,
    is connected, and no renaming that keeps the left side gives a smaller listing."""
    if len({v for rel in R for v in rel if v >= k}) != m: return False
    if not engine().connected_general(R): return False
    for perm in auts:
        for nperm in itertools.permutations(range(k, k + m)):
            mp = perm + nperm
            if tuple(sorted(tuple(mp[v] for v in rel) for rel in R)) < R: return False
    return True

def run_ops(L, R, k, max_events=EVENTS, max_gen=MAX_GEN, max_rel=MAX_REL, max_ops=MAX_OPS):
    """run_general() of the paper's program, with its clock replaced by a count of matching steps. A generation that
    passes the budget is dropped whole, so the result does not depend on the machine."""
    P = [tuple(p) for p in L]; n = len(P); arity = len(P[0])
    init = [tuple([0] * arity)] * n
    rel = {}; creator = {}; pos_idx = defaultdict(set)
    nid = 0; node = 1; ops = 0
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
        nonlocal ops
        ops += 1
        if len(pat) != len(t): return None
        b = dict(bd)
        for var, val in zip(pat, t):
            if var in b:
                if b[var] != val: return None
            else: b[var] = val
        return b
    def cands(bd, pat):
        for i, var in enumerate(pat):
            if var in bd: return pos_idx[(i, bd[var])]
        return None
    def fill(bd, todo, chosen, used):
        nonlocal ops
        if not todo: return chosen, bd
        j = next((j for j in todo if any(v in bd for v in P[j])), todo[0])
        rest = [x for x in todo if x != j]
        c = cands(bd, P[j])
        pool = sorted(c) if c is not None else sorted(rel)
        ops += len(pool) >> 4                      # sorting the candidates is work too
        taken = set(chosen.values())
        for r in pool:
            if r in used or r in taken: continue
            b2 = extend(bd, P[j], rel[r])
            if b2 is None: continue
            ch = dict(chosen); ch[j] = r
            res = fill(b2, rest, ch, used)
            if res is not None: return res
        return None
    cut = False
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
            if ops > max_ops: cut = True; break
        if cut or not events: break
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
        if len(parents) - 1 >= max_events or len(rel) > max_rel: break
    return parents, gens, len(rel), node, ops, cut

def clean(x):
    """JSON-safe: a reading that is not a finite number becomes null."""
    if isinstance(x, float): return x if math.isfinite(x) else None
    if isinstance(x, dict): return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [clean(v) for v in x]
    return x

def run_rule(i, L, R, k, long=False):
    A = engine()
    if long: parents, gens, nrel, nnodes, ops, cut = run_ops(L, R, k, LONG_EVENTS, LONG_GEN, LONG_REL, LONG_OPS)
    else: parents, gens, nrel, nnodes, ops, cut = run_ops(L, R, k)
    c, fit = A.growth(gens)
    rec = {'i': i, 'rule': A.rule_str(L, R), 'ev': len(parents) - 1, 'g': max(gens), 'rel': nrel, 'nodes': nnodes,
           'ops': ops, 'cut': int(cut), 'fit': fit}
    if rec['ev'] >= MIN_EVENTS:
        m = A.measure(parents, series=long)
        for key in ('depth', 'tree', 'past_max', 'max_causes', 'max_effects', 'cone_dim', 'root', 'mid', 'series'):
            if key in m: rec[key] = m[key]
    return clean(rec)

def rules_of(u):
    """(index, right side) of every rule in a block, without running anything."""
    if u[0] == 'long': u = (u[1], u[2], u[3], u[4], u[4] + 1)
    cls, li, m, a, b = u
    nl, nr, ar = CLASSES[cls]; L, k, auts = lefts(cls)[li]
    rels = list(itertools.product(range(k + m), repeat=ar))
    return [(i, R) for i, R in right_sides(rels, nr, a, b) if accepted(R, k, m, auts)]

def run_unit(u, resume=None, say=None):
    """Run every rule of a block. With a resume file, records are appended as they finish and a restart continues."""
    long = u[0] == 'long'
    L, k, auts = lefts(u[1] if long else u[0])[u[2] if long else u[1]]; recs = []
    if resume and os.path.exists(resume):
        for line in open(resume):
            try: recs.append(json.loads(line))
            except ValueError: break                     # a line cut off by a crash: redo from there
    done = {r['i'] for r in recs}
    fh = open(resume, 'a') if resume else None
    todo = [(i, R) for i, R in rules_of(u) if i not in done]
    for n, (i, R) in enumerate(todo):
        rec = run_rule(i, L, R, k, long); recs.append(rec)
        if fh: fh.write(json.dumps(rec) + '\n'); fh.flush()
        if say: say(rec, n + 1, len(todo))
    if fh: fh.close()
    recs.sort(key=lambda r: r['i'])
    return {'v': VERSION, 'uid': uid(u), 'rules': recs}

# ------------------------------------------------------------------ comparing results
def int_view(rec):
    """The whole-number part of a record: equal on every machine, compared exactly."""
    iv = [rec['i'], rec['rule'], rec['ev'], rec['g'], rec['rel'], rec['nodes'], rec['ops'], rec['cut'],
          rec.get('depth'), rec.get('tree'), rec.get('past_max'), rec.get('max_causes'), rec.get('max_effects')]
    for lab in ('root', 'mid'):
        r = rec.get(lab)
        iv.append([r['diamonds'], r['n'], int(round(r['chain'] * 10))] if r else None)
    iv.append([[s['n'], s['chain'], s['height']] for s in rec.get('series', [])])
    return iv

def float_view(rec):
    """The readings with decimals: compared to within 0.002, since logarithms may differ in the last digit by machine."""
    fv = [rec.get('fit', {}).get('exp_rate'), rec.get('fit', {}).get('power'), rec.get('cone_dim')]
    for lab in ('root', 'mid'):
        r = rec.get(lab) or {}
        fv.extend(r.get(key) for key in ('f', 'mm', 'mid', 'chain_d', 'coarse_d'))
    for s in rec.get('series', []): fv.extend(s.get(key) for key in ('frac', 'f', 'mm', 'mid', 'chain_d', 'coarse_d'))
    return fv

def digest(result):
    body = json.dumps([result['v'], result['uid'], [int_view(r) for r in result['rules']]], separators=(',', ':'))
    return hashlib.sha256(body.encode()).hexdigest()

def floats_agree(r1, r2, tol=0.002):
    for a, b in zip(r1['rules'], r2['rules']):
        for x, y in zip(float_view(a), float_view(b)):
            if (x is None) != (y is None): return False
            if x is not None and abs(x - y) > tol: return False
    return True

def band4(n):
    sizes = sorted(BANDS4)
    if n <= sizes[0]: return BANDS4[sizes[0]]
    if n >= sizes[-1]: return BANDS4[sizes[-1]]
    for a, b in zip(sizes, sizes[1:]):
        if a <= n <= b:
            w = (math.log(n) - math.log(a)) / (math.log(b) - math.log(a))
            return tuple(x * (1 - w) + y * w for x, y in zip(BANDS4[a], BANDS4[b]))

def readings(rec):
    """Per interval: (label, size, three readings, largest distance from four dimensions in standard deviations)."""
    out = []
    for lab in ('root', 'mid'):
        r = rec.get(lab)
        if not r or None in (r.get('mm'), r.get('mid'), r.get('chain_d')): continue
        bd = band4(r['n']); v = (r['mm'], r['mid'], r['chain_d'])
        z = max(abs(v[j] - bd[2 * j]) / bd[2 * j + 1] for j in range(3))
        out.append((lab, r['n'], v, z))
    return out

def near_four(rec):
    """The paper's window: all three tests between 3.5 and 4.5 on the same interval."""
    return [x for x in readings(rec) if all(3.5 <= t <= 4.5 for t in x[2])]

def reads_four(iv):
    """Every test within three standard deviations of the four-dimensional reference at this size."""
    if None in (iv.get('mm'), iv.get('mid'), iv.get('chain_d')): return False
    bd = band4(iv['n'])
    return all(abs(iv[key] - bd[2 * j]) <= 3 * bd[2 * j + 1] for j, key in enumerate(('mm', 'mid', 'chain_d')))

def candidate(rec):
    """The paper's criterion fixed in advance (protocol item 5), on a long run: the largest interval of at least 1,000
    elements reads four, and so do the two largest intervals of the size series with at least 1,000 elements."""
    ivs = [rec[lab] for lab in ('root', 'mid') if rec.get(lab)]
    ser = [s for s in rec.get('series', []) if s['n'] >= 1000]
    big = [iv for iv in ivs + ser if iv['n'] >= 1000]
    if not big: return False, 'no interval of 1,000 elements'
    if not reads_four(max(big, key=lambda iv: iv['n'])): return False, 'largest interval fails'
    ser = sorted(ser, key=lambda s: s['n'], reverse=True)[:2]
    if len(ser) < 2: return False, 'fewer than two series intervals of 1,000'
    if not all(reads_four(s) for s in ser): return False, 'series drifts'
    return True, 'candidate'

# ------------------------------------------------------------------ the volunteer side
def settings():
    path = os.path.join(HERE, 'toehome.json')
    try: return json.load(open(path)) if os.path.exists(path) else {}
    except ValueError: sys.exit('toehome.json is not valid JSON')

def client_id():
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'client.json')
    if os.path.exists(path):
        try: return json.load(open(path))
        except ValueError: pass
    me = {'id': uuid.uuid4().hex}                     # a random label; it says nothing about the machine or its owner
    json.dump(me, open(path, 'w'))
    return me

def call(server, path, obj, timeout=90):
    if not server.startswith(('https://', 'http://')): raise ValueError('the coordinator address must start with https:// or http://')
    body = gzip.compress(json.dumps(obj, separators=(',', ':')).encode())
    req = urllib.request.Request(server.rstrip('/') + path, data=body, method='POST',
                                 headers={'Content-Type': 'application/json', 'Content-Encoding': 'gzip',
                                          'User-Agent': 'toehome/%d' % VERSION})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        raw = r.read(1 << 20)
    return json.loads(raw.decode())

def low_priority():
    try:
        if hasattr(os, 'nice'): os.nice(19)
        elif os.name == 'nt':
            import ctypes
            ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)   # idle
    except Exception:
        pass

def describe(rec):
    rd = readings(rec)
    if not rd: return 'stops or stays small' if rec['ev'] < MIN_EVENTS else 'no interval large enough to read'
    lab, n, v, z = min(rd, key=lambda x: x[3])
    return 'reads %.2f / %.2f / %.2f on %d elements' % (v[0], v[1], v[2], n)

def keep(result, extra):
    """Store a finished block (and what is needed to upload it) in toehome_data/done."""
    d = os.path.join(DATA, 'done'); os.makedirs(d, exist_ok=True)
    path = os.path.join(d, result['uid'].replace(':', '_') + '.json.gz')
    with gzip.open(path, 'wt') as fh: json.dump(dict(extra, result=result), fh, separators=(',', ':'))
    return path

def upload(server, path, me, nick):
    pack = json.load(gzip.open(path, 'rt'))
    msg = {'v': VERSION, 'client': me['id'], 'nick': nick, 'uid': pack['result']['uid'], 'token': pack.get('token'),
           'result': pack['result'], 'sec': pack.get('sec'), 'py': platform.python_version()}
    ans = call(server, '/result', msg)
    if ans.get('status') in ('validated', 'waiting', 'conflict', 'duplicate', 'refused'):
        os.replace(path, path.replace('.json.gz', '.sent.json.gz'))
    return ans

def pending():
    d = os.path.join(DATA, 'done')
    return sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith('.json.gz') and '.sent.' not in f) if os.path.isdir(d) else []

def worker(w, server, nick, offline, classes, limit):
    low_priority()
    me = client_id(); rnd = random.Random(); done = 0
    wdir = os.path.join(DATA, 'work'); os.makedirs(wdir, exist_ok=True)
    lease_path = os.path.join(wdir, 'w%d.unit.json' % w); resume = os.path.join(wdir, 'w%d.jsonl' % w)
    tag = '[%d] ' % w
    try:
        while limit is None or done < limit:
            lease = None
            if os.path.exists(lease_path):                                   # a block left unfinished by the last run
                try: lease = json.load(open(lease_path))
                except ValueError: lease = None
            if lease is None:
                if os.path.exists(resume): os.remove(resume)
                if offline:
                    cls = rnd.choice(classes); lease = {'uid': uid(rnd.choice(unit_list(cls)))}
                else:
                    try: ans = call(server, '/work', {'v': VERSION, 'client': me['id'], 'nick': nick})
                    except (urllib.error.URLError, OSError, ValueError) as e:
                        print(tag + 'coordinator not reachable (%s); trying again in a minute' % e, flush=True); time.sleep(60); continue
                    if ans.get('done'): print(tag + 'every block is assigned or finished. Thank you.', flush=True); return
                    if 'uid' not in ans:
                        print(tag + str(ans.get('msg', 'no work right now'))[:200], flush=True); time.sleep(min(600, max(10, int(ans.get('wait', 60))))); continue
                    lease = {'uid': ans['uid'], 'token': ans.get('token')}
                json.dump(lease, open(lease_path, 'w'))
            try: u = parse_uid(lease['uid'])
            except ValueError:
                print(tag + 'the coordinator sent a block label this version does not know; skipped', flush=True)
                os.remove(lease_path); time.sleep(30); continue
            t0 = time.time(); print(tag + 'block %s' % lease['uid'], flush=True)
            def say(rec, n, total):
                hit = near_four(rec)
                if candidate(rec)[0]: print(tag + '*** CANDIDATE: %s passes the four-dimension criterion on this run. It will be rerun and checked before anything is claimed.' % rec['rule'], flush=True)
                elif hit: print(tag + '*** all three tests between 3.5 and 4.5: %s (%s). Reported for a long run.' % (rec['rule'], describe(rec)), flush=True)
                elif n % 10 == 0 or n == total: print(tag + '  %d of %d rules; last: %s' % (n, total, describe(rec)), flush=True)
            result = run_unit(u, resume, say)
            path = keep(result, {'token': lease.get('token'), 'sec': round(time.time() - t0, 1), 'client': me['id']})
            os.remove(lease_path)
            if os.path.exists(resume): os.remove(resume)
            done += 1
            print(tag + 'block %s done: %d rules in %.0f s' % (lease['uid'], len(result['rules']), time.time() - t0), flush=True)
            if not offline:
                for p in pending():
                    try:
                        ans = upload(server, p, me, nick)
                        print(tag + 'uploaded %s: %s' % (os.path.basename(p), ans.get('status') or ans.get('msg')), flush=True)
                    except (urllib.error.URLError, OSError, ValueError) as e:
                        print(tag + 'upload failed (%s); the file is kept and sent later' % e, flush=True); break
    except KeyboardInterrupt:
        print(tag + 'stopped; the unfinished block continues on the next start', flush=True)

def selftest(quiet=False):
    label, want = SELFTEST
    t0 = time.time(); got = digest(run_unit(parse_uid(label)))
    ok = got == want
    if not quiet or not ok:
        print('self-test %s in %.0f s (block %s)' % ('passed' if ok else 'FAILED', time.time() - t0, label))
        if not ok: print('  expected %s\n  got      %s\nThis machine would return results that do not match. Please report it; nothing was sent.' % (want, got))
    return ok

def main():
    ap = argparse.ArgumentParser(description='TOE@home volunteer client')
    ap.add_argument('--server'); ap.add_argument('--nick'); ap.add_argument('--cores', type=int, default=1)
    ap.add_argument('--offline', action='store_true'); ap.add_argument('--send', action='store_true')
    ap.add_argument('--selftest', action='store_true'); ap.add_argument('--unit')
    ap.add_argument('--classes', default='b2_4'); ap.add_argument('--blocks', type=int, help='stop after this many blocks per process')
    a = ap.parse_args()
    if sys.version_info < (3, 10): sys.exit('Python 3.10 or later is needed (this is %s)' % platform.python_version())
    engine()
    if a.unit:
        res = run_unit(parse_uid(a.unit), say=lambda rec, n, total: print('  %s  %s' % (rec['rule'], describe(rec)), flush=True))
        print(json.dumps({'uid': res['uid'], 'rules': len(res['rules']), 'digest': digest(res)})); return
    if a.selftest: sys.exit(0 if selftest() else 1)
    cfg = settings(); server = a.server or cfg.get('server') or ''
    nick = ''.join(ch for ch in (a.nick or cfg.get('nick') or '') if ch.isalnum() or ch in ' _-')[:24]
    me = client_id()
    if me.get('selftest') != VERSION:                       # once per machine and version
        if not selftest(): sys.exit(1)
        me['selftest'] = VERSION; json.dump(me, open(os.path.join(DATA, 'client.json'), 'w'))
    if a.send:
        if not server: sys.exit('no coordinator address: use --server URL or toehome.json')
        for p in pending(): print(os.path.basename(p), upload(server, p, me, nick))
        return
    offline = a.offline or not server
    if offline and not a.offline: print('No coordinator address is set, so this runs offline: blocks are taken at random and results are kept in toehome_data/done.')
    classes = [c for c in a.classes.split(',') if c in CLASSES] or ['b2_4']
    print('TOE@home version %d, %d process(es), %s. Stop with Ctrl-C; progress is kept.' % (VERSION, a.cores, 'offline' if offline else server), flush=True)
    if a.cores <= 1:
        worker(0, server, nick, offline, classes, a.blocks)
    else:
        import multiprocessing as mp
        ps = [mp.Process(target=worker, args=(w, server, nick, offline, classes, a.blocks)) for w in range(a.cores)]
        for p in ps: p.start()
        try:
            for p in ps: p.join()
        except KeyboardInterrupt:
            for p in ps: p.join()

if __name__ == '__main__':
    main()
