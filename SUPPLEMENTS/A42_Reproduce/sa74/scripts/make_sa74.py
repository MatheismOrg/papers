#!/usr/bin/env python3
"""Apply the SA74 changes to A42 (source: the SA72 T372 Markdown). Every number comes from sa74_numbers.json and the
run files. Edits are exact replacements with a count check, so nothing else in the paper moves."""
import json, sys, re, math
sys.path.insert(0, '/home/claude/sa74/scripts')
import a42_classify as C
BASE = '/home/claude/sa74/src/A42_SA72_T372_package.md'   # the SA72 T372 source, unchanged
OUT = '/home/claude/sa74/out/A42_Shortest_Rules_Search.md'
R = '/home/claude/sa74/runs/'
N = json.load(open(R + 'sa74_numbers.json'))
md = open(BASE, encoding='utf-8').read()
CHANGES = []

def rep(old, new, count=1, label=''):
    global md
    n = md.count(old)
    assert n == count, (label, n, old[:80])
    md = md.replace(old, new)
    CHANGES.append((label, old, new))

def arrow(rule):
    l, r = rule.split('->'); return f'{l} → {r}'

def f2(x): return f'{x:.2f}'

# ------------------------------------------------------------------ numbers used in the text
ci = N['ci']; ci_total = ci['pass2']; pct = round(100 * ci_total / 2935)
beh = dict(ci['pass2_behavior'])
n_re3 = N['ci_steps3_done']; n_p3 = ci['pass3_extra']
tuned = sorted(((int(k), v) for k, v in N['tp_tuned'].items()), key=lambda x: x[0])
mids = [s['mid'] for _, t in tuned for s in t['seeds']]; chains = [s['chain_d'] for _, t in tuned for s in t['seeds']]
zch = [s['z']['chain_d'] for _, t in tuned for s in t['seeds']]
transfer = json.load(open('/home/claude/sa74/transfer/transfer_expanding_results.json'))
tr_last = transfer['runs'][-1]; tr_ratios = transfer['error_ratios']
dsb = N['ds_big']['8.0']; op32 = N['op']['32000_0.001']['root']; op32b = N['op']['32000_0.0005']['root']
big = {}
for l in open('/home/claude/sa74/big/trajectory2.jsonl'):
    r = json.loads(l); big.setdefault(r['N'], []).append(r)
bigm = {Nn: {k: sum(r[k] for r in v) / len(v) for k in ('p', 'n', 'mm', 'mid', 'chain_d')} for Nn, v in big.items()}
spr = [json.loads(l) for l in open('/home/claude/sa74/big/sprinkle4.jsonl') if l.startswith('{')]
spr = [x for x in spr if x['N'] >= 8000]
sr = {k: (min(x[k] for x in spr), max(x[k] for x in spr)) for k in ('mm', 'mid', 'chain_d')}

# ------------------------------------------------------------------ abstract and keywords
if False: rep("If spacetime is the causal order of discrete updates, a theory of everything is an update rule, and the search can run shortest first, since each extra bit halves a rule's prior weight.",
    "Spacetime is the causal order of discrete updates, so a theory of everything is an update rule, and the search can run shortest first, since each extra bit halves a rule's prior weight.", label='Abstract, first sentence')
if False: rep("A mirror pair of rules grows two-dimensional orders, on which all three tests agree up to 49,000 events.",
    "Run oldest relations first, a mirror pair of rules grows two-dimensional orders, on which all three tests agree up to 49,000 events.", label='Abstract, two dimensions')
if False: rep("A rewriting rule can match four-dimensional spacetime only at large scales. Run from the other end,",
    f"A rewriting rule can match four-dimensional spacetime only at large scales. A criterion fixed in advance, with noise bands from random point sets of matched size, puts the near misses' chains 4.7 to 21 standard deviations short of four dimensions. The update order matters: under random orders the two-dimensional pair grows orders that read near one on the first two tests, although it passes a critical-pair test of causal invariance, as {ci_total:,} of the 2,935 rules do, and no rule reads four under any order. Sequential growth, which escapes the link bound, also fails: transitive percolation held at four on the first test reads {min(mids):.1f} to {max(mids):.1f} on the second and {min(chains):.1f} to {max(chains):.1f} on the chain at up to 40,000 elements, though its chain reading climbs as it grows. Run from the other end,",
    label='Abstract, new results')
if False: rep("a first upper bound from the top for matter on a given spacetime.",
    "a first upper bound from the top for matter on a given spacetime, and run unchanged on an expanding universe it was not built for, it creates fermion pairs as the continuum does, to within 3 x 10^-7.",
    label='Abstract, transfer test')
rep("Schwinger model; semiclassical gravity\n", "Schwinger model; semiclassical gravity; causal invariance; transitive percolation; sequential growth; preregistration\n", label='Keywords')

# ------------------------------------------------------------------ introduction
rep("and supplies the code and data. The result is negative",
    "and supplies the code and data. It fixes a four-dimension criterion in advance, with noise bands from random point sets of matched size, tests every rule of the smallest class for causal invariance, reruns the long runs under four update orders, and reads sequential growth with the same tests. The result is negative",
    label='Introduction, what the paper does')
rep("The rules that come closest read four on two tests and between 1.7 and 3.0 on the third.",
    "The rules that come closest read four on two tests and between 1.7 and 3.0 on the third. The two-dimensional orders hold under the first passes' update order only, and sequential growth, the way around the obstacle below, does not read four on all three tests anywhere it was run.",
    label='Introduction, results')

# ------------------------------------------------------------------ section 2
rep("Rules that differ only by renaming variables count once, and so do rules whose relations are listed in a different order.",
    "Rules that differ only by renaming variables count once, and so do rules whose relations are listed in a different order. Under the fixed update order of Section 2.2 the listing is not neutral, since new relations are aged in the order the right side lists them: listed as {{x,y},{x,z}} → {{z,x},{x,w},{x,y}}, the two-dimensional rule of Section 4.2 grows a single chain. The census runs one listing of each rule, its canonical one.",
    label='Section 2.1, listing')
rep("all of them are applied at once. Every",
    "all of them are applied at once. Section 4.6 varies this order. Every", label='Section 2.2')

# ------------------------------------------------------------------ tables renumbered in order of appearance
rep("**Table 2.** What the rules do in 3,000 updates.", "**Table 3.** What the rules do in 3,000 updates.", label='Table 2 -> 3')
rep("Table 2 sorts the rules", "Table 3 sorts the rules", label='Table 2 -> 3, text')
rep("**Table 3.** Long runs to 20,000 updates", "**Table 4.** Long runs to 20,000 updates", label='Table 3 -> 4')
rep("Table 3 shows the two cleanest", "Table 4 shows the two cleanest", label='Table 3 -> 4, text')
rep("Table 2 counts them that way", "Table 3 counts them that way", label='Table 2 -> 3, Section 6')
rep("**Table 4.** Links per element: sprinkled four-dimensional intervals and two grids.",
    "**Table 9.** Links per element: sprinkled intervals, two grids and the two-dimensional pair.", label='Table 4 -> 9')
rep("(Table 4)", "(Table 9)", label='Table 4 -> 9, text')

# ------------------------------------------------------------------ section 3: noise bands and the fixed criterion
b = N['bands']
def brow(d, n):
    r = b[f'{d}_{n}']
    return f"| {d} | {n:,} | {r['reps']} | {r['mm_mean']:.2f} ± {r['mm_sd']:.2f} | {r['mid_mean']:.2f} ± {r['mid_sd']:.2f} | {r['chain_d_mean']:.2f} ± {r['chain_d_sd']:.2f} |"
rows = [brow(4, n) for n in (120, 250, 500, 1000, 2000, 4000)] + [brow(3, n) for n in (1000, 4000)] + [brow(2, n) for n in (1000, 4000)]
r1000 = b['4_1000']
sec3 = ("\n\nA fixed window ignores how much each test scatters at a given size, and the scatter shrinks as intervals grow. "
        "Table 2 measures it on sprinkled intervals of 120 to 4,000 points, with 40 sprinklings at the smaller sizes and 6 at the largest. "
        f"At 1,000 points a four-dimensional interval reads {r1000['mm_mean']:.2f} ± {r1000['mm_sd']:.2f} on the ordering fraction, {r1000['mid_mean']:.2f} ± {r1000['mid_sd']:.2f} on the midpoint and {r1000['chain_d_mean']:.2f} ± {r1000['chain_d_sd']:.2f} on the longest chain. "
        "From Section 4.4 on, an interval reads d-dimensional when each of its three readings lies within three standard deviations of the d-dimensional mean at its size. "
        "A rule is a four-dimensional candidate when its largest interval of at least 1,000 elements reads four-dimensional and so do the two largest intervals of its size series, so that its readings do not drift away from four as it grows. "
        "The size floor matters, since the bands are wide at 120 points.\n\n"
        "This criterion, and the runs of Sections 4.6 and 5.1, were written down and fingerprinted with SHA-256 before they ran (04eca2085b16830e…; the full value and the protocol are in the supplements), so a near miss cannot be reread after the fact. Runs added after the protocol are marked exploratory.\n\n"
        "**Table 2.** The scatter of the three tests on sprinkled flat intervals: mean ± standard deviation.\n\n"
        "| Dimension | Points | Sprinklings | Ordering fraction | Midpoint | Longest chain |\n|---:|---:|---:|---:|---:|---:|\n" + '\n'.join(rows))
rep("and as three-dimensional if all read between 2.5 and 3.5.",
    "and as three-dimensional if all read between 2.5 and 3.5." + sec3, label='Section 3, noise bands and Table 2')

# ------------------------------------------------------------------ section 4.2: the pair against the bands
p2 = N['pair2d']
rep("So the tests find a dimension when a short rule has one, and two dimensions are within reach of the shortest rules.",
    "So the tests find a dimension when a short rule has one, and two dimensions are within reach of the shortest rules.\n\n"
    f"The noise bands of Table 2 sharpen this. The pair's orders are two-dimensional and are not two-dimensional sprinklings: their chains run {round(100 * (p2['chain_longer'][0] - 1))} to {round(100 * (p2['chain_longer'][1] - 1))} percent longer than a sprinkled interval's of the same size, {abs(p2['z_chain'][1]):.1f} to {abs(p2['z_chain'][0]):.0f} standard deviations below the band, and at their largest sizes the first two readings fall to {p2['first_two_largest'][0]:.2f} to {p2['first_two_largest'][1]:.2f}. "
    "Their updates have few links (Section 5), so at the scale of single updates the order is closer to a lattice than to a sprinkling. Section 4.6 shows that the two-dimensional readings also depend on the update order.",
    label='Section 4.2, the pair against the bands')

# ------------------------------------------------------------------ section 4.4: distances
nf, ns = N['near_first'], N['near_second']
rep("would have to rise by 0.65 to 1.8 to reach the four-dimensional range.",
    "would have to rise by 0.65 to 1.8 to reach the four-dimensional range.\n\n"
    f"The noise bands of Table 2 turn the miss into a distance. On the {nf['count']} intervals, the chain readings lie {abs(nf['z_chain'][1]):.1f} to {abs(nf['z_chain'][0]):.0f} standard deviations below the four-dimensional mean at their sizes. "
    f"The runs past 100,000 updates add {ns['count']} intervals that read four on the first two tests; their chains lie {abs(ns['z_chain'][1]):.1f} to {abs(ns['z_chain'][0]):.0f} standard deviations below and run {ns['ratio'][0]:.1f} to {ns['ratio'][1]:.0f} times longer than four dimensions allow. "
    "No interval of the long runs, at any size, reads four-dimensional within the bands.",
    label='Section 4.4, distances')

# ------------------------------------------------------------------ section 4.6: orders and causal invariance
ords = N['orders']
pair = ['{{x,y},{x,z}}->{{x,y},{x,w},{z,x}}', '{{x,y},{z,y}}->{{x,y},{y,z},{w,y}}']
other2d = ['{{x,y},{x,z}}->{{y,x},{y,w},{z,x}}', '{{x,y},{z,y}}->{{y,x},{y,z},{w,x}}']
near = ['{{x,y},{y,z}}->{{x,y},{y,w},{z,y}}', '{{x,y},{y,z}}->{{x,y},{y,w},{z,y},{w,y}}', '{{x,y},{y,z}}->{{x,y},{y,x},{z,y},{w,y}}',
        '{{x,y},{x,z}}->{{y,y},{z,w},{w,y}}', '{{x,y},{z,y}}->{{x,x},{x,w},{w,z}}', '{{x,y},{y,z}}->{{y,x},{y,w},{z,x}}',
        '{{x,y},{y,z}}->{{z,x},{z,w},{w,x}}', '{{x,y},{z,y}}->{{x,z},{x,w},{y,x},{w,y}}', '{{x,y},{x,z}}->{{x,y},{z,y},{z,w},{w,z}}']
def cell(rule, o):
    c = ords[rule].get(o)
    if not c: return ''
    if c['kind'] == 'single chain': return 'single chain'
    L = c['largest']
    if not L: return c['kind']
    return f"{L['mm']:.2f}, {L['mid']:.2f}, {L['chain_d']:.2f}"
def ci_word(rule):
    s = ords[rule]['ci']; return {'2': 'passes', '3': 'passes (3)', 'no': 'fails'}[s]
t8 = ["| Rule | Critical pairs | Oldest first | Newest first | Random 1 | Random 2 |", "|---|---|---:|---:|---:|---:|"]
for rule in pair + other2d + near:
    t8.append(f"| {arrow(rule)} | {ci_word(rule)} | " + ' | '.join(cell(rule, o) for o in ('oldest', 'newest', 'random1', 'random2')) + ' |')
nchain = sum(1 for k, v in ords.items() if not (v['group'] or '').startswith('two-dimensional') and v.get('newest', {}).get('kind') == 'single chain')
mc = N['orders_maxchain']; mcz = C.zscores({'n': mc[1][2], 'mm': mc[1][3], 'mid': mc[1][4], 'chain_d': mc[0]}, 4)['chain_d']
# the pair under the orders
det = [ords[r][o] for r in pair for o in ('oldest', 'newest')]
rnd = [ords[r][o] for r in pair for o in ('random1', 'random2')]
lo_det = min(min(c['range'][k][0] for k in ('mm', 'mid', 'chain_d')) for c in det); hi_det = max(max(c['range'][k][1] for k in ('mm', 'mid', 'chain_d')) for c in det)
fr = [json.loads(l) for o in ('random1', 'random2') for l in open(R + f'orders_{o}.jsonl')]
fr = [x for x in fr if x['rule'] in pair]
pastfrac = [x['past_max'] / x['events'] for x in fr]
rnd12 = (min(min(c['range']['mm'][0], c['range']['mid'][0]) for c in rnd), max(max(c['range']['mm'][1], c['range']['mid'][1]) for c in rnd))
rndch = (min(c['range']['chain_d'][0] for c in rnd), max(c['range']['chain_d'][1] for c in rnd))
links_pair = [ords[r][o]['links_mean'] for r in pair for o in ('oldest', 'newest', 'random1', 'random2')]
sec46 = (
 "### 4.6 Other update orders, and causal invariance\n\n"
 "A rule fixes its spacetime by itself only if its causal graph does not depend on the order of its updates; such a rule is causally invariant (Wolfram, 2020; Gorard, 2020b). "
 "The local test is the critical pair. Two updates compete for a relation; apply each in the smallest state that holds both, and ask whether further updates bring the two results to the same state, up to the names of new nodes. "
 f"Every rule of the full class was tested with up to two further updates on each side. {ci_total:,} of the 2,935 pass, {pct} percent: {beh.get('polynomial', 0)} grow polynomially, {beh.get('exponential', 0)} exponentially, {beh.get('tree', 0)} as single chains or trees, and {beh.get('stops or stays small', 0)} stop. "
 f"With three further updates, {n_p3} of the first {n_re3:,} failing rules rechecked also pass, so the count is a floor. "
 "Of the 18 rules that read two, only the cleanest pair passes; 15 of the other 16 still fail with three further updates.\n\n"
 "Table 8 reruns the 27 rules of the first pass's long runs to 20,000 updates under four orders: the first passes' order, which takes the oldest relations first and pairs each with its oldest partner; its reverse, newest first; and two random orders. "
 f"No rule reads four on all three tests under any order. The highest chain reading under any order, {mc[0]:.2f} at {mc[1][2]:,} elements, belongs to {arrow(mc[1][0])} under newest first, still {abs(mcz):.1f} standard deviations below four dimensions. "
 f"Under newest first, {nchain} of the 23 candidates grow one-dimensional orders, one or two updates per generation, and the exponential near misses lose their readings near four.\n\n"
 f"The two-dimensional pair passes the critical-pair test, and its readings still change with the order. Under the two deterministic orders it reads {lo_det:.2f} to {hi_det:.2f} on all three tests. "
 f"Under the two random orders the past of its latest update holds {round(100 * min(pastfrac))} to {round(100 * max(pastfrac))} percent of all updates, and its intervals read {rnd12[0]:.2f} to {rnd12[1]:.2f} on the first two tests and {rndch[0]:.2f} to {rndch[1]:.2f} on the chain: nearly one-dimensional on two tests. "
 "Listing its right side in another order makes the oldest-first engine grow a single chain (Section 2.1); under the random orders that listing grows the same nearly one-dimensional orders.\n\n"
 "So a finite run shows the rule together with its order. Passing the critical-pair test is a local condition. When the rewriting runs without end, local resolution does not guarantee that every run grows the same causal graph, and here the runs do not: each grows only part of the graph, and different orders grow different parts. "
 "A rule's update order, and under a fixed order the listing of its right side, therefore belong in its description, and their bits count toward its length. "
 "The first passes' order costs few bits, so the pair's two-dimensional orders stand as the output of a short program, the rule run oldest first. A four-dimensional candidate must read four under every order, or carry its order in its length.\n\n"
 "**Table 8.** The rules that read two and the near misses, run to 20,000 updates under four orders: the readings of the largest interval of each size series (ordering fraction, midpoint, longest chain), or the run's behavior when it grows a single chain. The second column gives the critical-pair test, with three further updates where two did not resolve.\n\n"
 + '\n'.join(t8) + "\n\n")
rep("## 5. Why a rule can match four dimensions only at large scales", sec46 + "## 5. Why a rule can match four dimensions only at large scales", label='Section 4.6 and Table 8')

# ------------------------------------------------------------------ section 5: links of the pair; Table 9 rows
ls = {(x['d'], x['n']): x for x in N['links_sprinkled']}
rep("at the scale of single updates.\n",
    "at the scale of single updates.\n\n"
    f"The two-dimensional pair shows the cost in numbers. Its updates have {min(links_pair):.1f} to {max(links_pair):.1f} links on average under the four orders, and at most five, while sprinkled two-dimensional intervals of 1,000 to 4,000 points give their elements {ls[(2, 1000)]['mean_links']:.1f} to {ls[(2, 4000)]['mean_links']:.1f} on average and up to {max(ls[(2, 1000)]['max_links'], ls[(2, 4000)]['max_links'])} (Table 9).\n",
    label='Section 5, links of the pair')
rep("| Sprinkled, d = 4 | 1,600 | 47.6 | 171 |\n",
    "| Sprinkled, d = 4 | 1,600 | 47.6 | 171 |\n"
    f"| Sprinkled, d = 4 | 4,000 | {ls[(4, 4000)]['mean_links']:.1f} | {ls[(4, 4000)]['max_links']} |\n"
    f"| Sprinkled, d = 2 | 1,000 | {ls[(2, 1000)]['mean_links']:.1f} | {ls[(2, 1000)]['max_links']} |\n"
    f"| Sprinkled, d = 2 | 4,000 | {ls[(2, 4000)]['mean_links']:.1f} | {ls[(2, 4000)]['max_links']} |\n", label='Table 9 rows')
rep("| Any rule, two into three | any | at most 5 | at most 5 |",
    "| Any rule, two into three | any | at most 5 | at most 5 |\n"
    f"| Two-dimensional pair, four orders | 20,000 updates | {min(links_pair):.1f} to {max(links_pair):.1f} | 5 |", label='Table 9, pair row')

# ------------------------------------------------------------------ section 5.1: sequential growth
t10 = ["| Model | Elements | p | Largest interval | Ordering fraction | Midpoint | Longest chain |", "|---|---:|---:|---:|---:|---:|---:|"]
for n, t in tuned:
    t10.append(f"| Transitive, held at four | {n:,} | {t['p']:.5f} | {round(t['n']):,} | {t['mm']:.2f} | {t['mid']:.2f} | {t['chain_d']:.2f} |")
for Nn in (128000, 256000, 512000, 1024000, 2048000):
    b_ = bigm[Nn]; t10.append(f"| Transitive, held at four | {Nn:,} | {b_['p']:.6f} | {round(b_['n']):,} | {b_['mm']:.2f} | {b_['mid']:.2f} | {b_['chain_d']:.2f} |")
t10.append(f"| Sprinkled, flat, d = 4 | 8,000 to 128,000 | | 8,002 to 128,002 | {sr['mm'][0]:.2f} to {sr['mm'][1]:.2f} | {sr['mid'][0]:.2f} to {sr['mid'][1]:.2f} | {sr['chain_d'][0]:.2f} to {sr['chain_d'][1]:.2f} |")
t10.append(f"| Originary | 32,000 | 0.001 | {op32['n']:,} | {op32['mm']:.2f} | {op32['mid']:.2f} | {op32['chain_d']:.2f} |")
t10.append(f"| Originary | 32,000 | 0.0005 | {op32b['n']:,} | {op32b['mm']:.2f} | {op32b['mid']:.2f} | {op32b['chain_d']:.2f} |")
t10.append(f"| de Sitter, 8 e-folds | {dsb['n'] - 2:,} | | {dsb['n']:,} | {dsb['mm'][2]:.2f} | {dsb['mid'][2]:.2f} | {dsb['chain_d'][2]:.2f} |")
sec51 = (
 "### 5.1 Sequential growth under the same tests\n\n"
 "Sequential growth escapes the bound, so it gets the same tests. In transitive percolation, the simplest such dynamics, elements are born one at a time, and each earlier element becomes a direct cause of the newborn with probability p; the order is the transitive closure (Rideout and Sorkin, 2000). "
 "It has one constant, so it is short, and its elements can have as many links as they need. A root added below all elements serves as the bottom of the largest intervals.\n\n"
 "Over the protocol's grid, from 1,000 to 8,000 elements and p from 0.0003 to 0.3, no causal set is a candidate. Dense ones, with pN of 30 or more, are nearly one long chain: their largest intervals read between 1.0 and 2.2. "
 "Sparse ones, with pN of 3 or less, fall apart into pieces too small to test, and at pN of 4 to 8 they read 3.4 to 6.5 on the first two tests and 2.2 to 3.4 on the chain. "
 "In between, at pN of 8 to 24, the intervals that read near four on the first two tests read 2.3 to 3.1 on the chain: wide and long at once, like the rewriting near misses.\n\n"
 "Rideout and Sorkin (2001) found evidence that transitive percolation has a continuum limit when p falls as N grows so that the fraction of related pairs stays fixed, including at its four-dimensional value. "
 f"Table 10 follows that trajectory as an exploratory run outside the protocol: p is tuned so that the largest interval reads four on the ordering fraction, and the causal set grows from 4,000 to 64,000 elements. "
 f"The midpoint reads {min(mids):.2f} to {max(mids):.2f}, above four, and the chain {min(chains):.2f} to {max(chains):.2f}, below it; the chain lies {abs(max(zch)):.0f} to {abs(min(zch)):.0f} standard deviations below four dimensions. "
 "The chain reading climbs by about 0.1 each time the causal set doubles, while the midpoint drifts slowly up.\n\n"
 "A compiled engine followed the trajectory further, to two million elements, reading the same largest interval: the chain exactly, the ordering fraction from 3,000 sampled ancestor counts, and the midpoint from the 600 elements that bottom-k sketches of their pasts and futures rank highest, each then counted exactly. At 16,000 and 128,000 elements this search finds the same midpoint as checking every element. "
 f"The chain reading reaches {bigm[1024000]['chain_d']:.2f} at one million elements and {bigm[2048000]['chain_d']:.2f} at two million, through four, while the midpoint drifts on to {bigm[2048000]['mid']:.2f} (Table 10). "
 f"Sprinkled four-dimensional intervals of 8,000 to 128,000 points read {sr['mm'][0]:.2f} to {sr['mm'][1]:.2f} on the ordering fraction, {sr['mid'][0]:.2f} to {sr['mid'][1]:.2f} on the midpoint and {sr['chain_d'][0]:.2f} to {sr['chain_d'][1]:.2f} on the chain. "
 "So the trajectory reads four on the ordering fraction and on the chain at about a million elements, and misses the midpoint by more than one and a half: no element of its largest intervals has both a large past and a large future, as an element near the middle of a four-dimensional interval does. Transitive percolation is not a four-dimensional spacetime at any size reached, and the lead it offered closes on the midpoint.\n\n"
 "Ahmed and Rideout (2010) fitted the chains and largest intervals of originary percolation, the variant in which every element descends from a first one, to four-dimensional de Sitter space for p from 0.0001 to 0.005. "
 "On the three tests it comes close to de Sitter space on two and misses on the third. "
 f"At p = 0.001 and 32,000 elements its largest intervals read {op32['mm']:.2f} on the ordering fraction and {op32['mid']:.2f} on the midpoint, as a de Sitter interval of the same size does across about eight e-folds of expansion, {dsb['mm'][2]:.2f} and {dsb['mid'][2]:.2f}, "
 f"yet their longest chains run {round(op32['chain'])} steps where de Sitter's run {dsb['chain'][0]} to {dsb['chain'][1]}, a chain reading of {op32['chain_d']:.2f} against {dsb['chain_d'][0]:.2f} to {dsb['chain_d'][1]:.2f} in three sprinklings. "
 "Wherever originary percolation reads near four on the first two tests, its chain reads 2.3 to 3.1, while de Sitter space with the same ordering fraction reads above 3.7. "
 "Glaser and Surya (2013) reached the same verdict on transitive percolation by another route: it looks manifoldlike on coarse measures and fails a count of small intervals that sprinklings pass.\n\n"
 "**Table 10.** Sequential growth: the medians of the largest intervals, means of two seeds for transitive percolation held at four on the ordering fraction up to 64,000 elements and of three seeds from 128,000 on (compiled engine), one seed for originary percolation, sprinkled flat intervals, and de Sitter space for comparison (means of three sprinklings).\n\n"
 + '\n'.join(t10) + "\n\n")
rep("## 6. Limits of this pass", sec51 + "## 6. Limits of this pass", label='Section 5.1 and Table 10')

# ------------------------------------------------------------------ section 6
rep("One update order and one start were used. A rule whose causal graph does not depend on the update order is causally invariant (Wolfram, 2020; Gorard, 2020b). This pass did not test for it, and a different order can change the causal graph of a rule that lacks it.",
    "One start was used. The first passes used one update order and one listing of each rule; Section 4.6 adds three orders for the long runs and a critical-pair test, a local condition, for the full class. The census itself was not rerun under other orders.",
    label='Section 6, orders')
rep("No rule run long shows such a climb on all three tests; the chain readings of the near misses peak between 1.9 and 3.0 and fall back.",
    "No rule run long shows such a climb on all three tests; the chain readings of the near misses peak between 1.9 and 3.0 and fall back. Transitive percolation held at four on the first test climbs on the chain, through four at about a million elements, while its midpoint moves away (Section 5.1).",
    label='Section 6, sizes')

# ------------------------------------------------------------------ section 7
rep("It should run several update orders and test causal invariance.",
    "It should keep the critical-pair test as a first filter, read every candidate under several update orders, and count the order a rule needs in its length.", label='Section 7, orders')
rep("and rules of the sequential-growth kind, whose updates draw on unbounded parts of their past.",
    "and rules of the sequential-growth kind, whose updates draw on unbounded parts of their past. Transitive percolation, the lead this pass followed to two million elements, reads four on two tests and fails the midpoint. A growth law that favors manifold-like orders is the next candidate: weighting each causal set by the causal-set action of Benincasa and Dowker (2010) gives, in two dimensions, a phase whose orders look like a continuum (Surya, 2012).",
    label='Section 7, leads')
rep("It reproduces the steps behind those tests to 10^-15.",
    f"It reproduces the steps behind those tests to 10^-15. Run unchanged on a background it was not built or checked for, an expanding 1+1 universe, with only its mass angle changing from step to step, the same program creates fermion pairs as the continuum does: {tr_last['pair_probability']:.7f} at {tr_last['sites']:,} sites against {transfer['continuum_pair_probability']:.7f}, with errors falling {min(tr_ratios):.2f} to {max(tr_ratios):.2f} times per halving.",
    label='Section 7, transfer test')

rep("growing the spacetime is what the search from the bottom must add, and two-dimensional orders already appear at 20 bits.",
    "growing the spacetime is what the search from the bottom must add, and two-dimensional orders already appear at 20 bits. Supplying a four-dimensional lattice as data would cost few bits, since the number four is one short symbol, and the statement's general relativity can bend it, as numerical relativity bends the metric's values on a grid; Section 8.6 runs a homogeneous version driven by quantum matter. What a fixed lattice leaves open is the bending at the scale of its own spacing, where the infinities of quantum gravity are only cut off, the lattice picks out a frame (Section 5), and no limit independent of the lattice has been established. A rule that grows its own spacetime would settle that scale too.",
    label='Section 7, what a supplied 4D lattice leaves open (T3, corrected T4)')
# ------------------------------------------------------------------ section 8
rep("cancelling", "canceling", count=2, label='Section 8.6, spelling')
rep("### 8.2 Results\n\n| Situation | Result | Scope |", "### 8.2 Results\n\n**Table 11.** Local steps against continuum references.\n\n| Situation | Result | Scope |", label='Table 11 caption')

# ------------------------------------------------------------------ conclusion
rep("Among the shortest rewriting rules, under one update order, none grows a four-dimensional causal order.",
    "Among the shortest rewriting rules, under the update orders tried, none grows a four-dimensional causal order.", label='Conclusion, orders')
rep("Two-dimensional orders do appear, so the tests can find a dimension when a short rule has one.",
    "Two-dimensional orders do appear under the first passes' order, so the tests can find a dimension when a short rule has one; under random orders the same rules grow orders that are nearly one-dimensional on two tests, so a rule's update order belongs in its length.",
    label='Conclusion, two dimensions')
rep("Such a rule can therefore match our spacetime only at large scales, where its small-scale frame must stay hidden.",
    "Such a rule can therefore match our spacetime only at large scales, where its small-scale frame must stay hidden. A criterion fixed in advance with noise bands puts the near misses' chains 4.7 to 21 standard deviations from four dimensions. Sequential growth escapes the link bound and still misses: transitive percolation held at four on the first test, followed to two million elements, reaches four on the third while its midpoint drifts to 5.8, and originary percolation, read against de Sitter space, runs its chains too long.",
    label='Conclusion, new results')
rep("and the same steps create pairs in an expanding universe",
    "and the same steps, run as that program unchanged, create pairs in an expanding universe", label='Conclusion, transfer')


# ------------------------------------------------------------------ SA74 T5: the known laws on a supplied lattice (Gary's route)
LT = json.load(open('/home/claude/sa74/lattice/results.json'))
lb = LT['bits']; ch5 = LT['chirality']; lo = LT['lorentz']['axis']; tu = LT['tuning']
cap = lb['plus lattice and declared']
gwmax = max(r['GW_violation_rel'] for r in ch5['rows'])
gwexp = int(f'{gwmax:.0e}'.split('e')[1]); gwman = math.ceil(gwmax / 10 ** gwexp)
assert all(r['index_n_minus_minus_n_plus'] == round(r['Q_measured']) for r in ch5['rows'])
amax_e = math.floor(math.log10(lo['a_max_m'])); amax_c = lo['a_max_m'] / 10 ** amax_e
vals = dict(naive=ch5['free_zero_modes']['naive'], nrows=len(ch5['rows']),
            pmax=f"{max(r['max_plaquette_angle'] for r in ch5['rows']):.1f}",
            gw=f'${gwman} \\times 10^{{{gwexp}}}$', m15=LT['anomalies']['mod16 for 3 families']['15 per family'],
            amax_m=f'${amax_c:.1f} \\times 10^{{{amax_e}}}$', amax_t=f'{amax_c:.1f} × 10' + str(amax_e).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')),
            amax_lp=round(lo['a_max_planck'] / 1e6), dsmash=f"{lb['plus lattice and SMASH'] - cap:,}",
            dlam=round(math.log10(1 / tu['vacuum energy in lattice units'])), dhig=round(math.log10(1 / tu['Higgs ratio'])),
            blam=round(tu['Lambda extra bits']), bhig=round(tu['Higgs extra bits']), btune=round(tu['total extra bits']),
            cap=f'{cap:,}', b0=f"{lb['A38 statement']:,}", blat=f"{lb['plus lattice form'] - lb['A38 statement']:,}",
            bdec=f"{cap - lb['plus lattice form']:,}")
assert vals['amax_lp'] == 43 and cap == 16272, vals
sec87 = open('/home/claude/sa74/t5/sec87_template.md').read().format(**vals).strip()
rep("### 8.7 Scope\n", sec87 + "\n\n### 8.8 Scope\n", label='New Section 8.7: the known laws on a supplied lattice, with Table 12 (T5)')
rep("would cost few bits, since the number four is one short symbol, and",
    "would cost few bits, since the number four is one short symbol, though the lattice form of each field costs more, and",
    label='Section 7, lattice form costs more (T5)')
rep("A rule that grows its own spacetime would settle that scale too.",
    f"A rule that grows its own spacetime would settle that scale too. Section 8.7 writes this candidate out in full, with the lattice form of each field and the completion the Standard Model needs: {cap:,} bits.",
    label='Section 7, pointer to Section 8.7 (T5)')
rep("and it exposes a lattice that passes every wave test yet doubles the source of gravity.",
    f"and it exposes a lattice that passes every wave test yet doubles the source of gravity. Section 8.7 sets the known laws in full on a supplied four-dimensional lattice, completed with dark matter and an initial state, at {cap:,} bits.",
    label='Introduction, the cap from the top (T5)')
rep("while a separate 3+1 lattice model drives its expansion.",
    f"while a separate 3+1 lattice model drives its expansion. Written out in full on a supplied four-dimensional lattice, with dark matter and an initial state, the known laws compress to {cap:,} bits, enough for present data within its known tensions; they join general relativity with quantum theory if the lattice's own quantum gravity yields smooth four-dimensional spacetime, the test that remains.",
    label='Conclusion, the cap (T5)')
rep("sequential growth; preregistration", "sequential growth; preregistration; Ginsparg-Wilson relation; Regge calculus", label='Keywords (T5)')


rep("Supplying a four-dimensional lattice as data would cost few bits, since the number four is one short symbol, though the lattice form of each field costs more, and the statement's general relativity can bend it, as numerical relativity bends the metric's values on a grid; Section 8.6 runs a homogeneous version driven by quantum matter. What a fixed lattice leaves open is the bending at the scale of its own spacing, where the infinities of quantum gravity are only cut off, the lattice picks out a frame (Section 5), and no limit independent of the lattice has been established. A rule that grows its own spacetime would settle that scale too.",
    "Supplying a four-dimensional lattice as data costs few bits, and the statement's general relativity can bend it, as numerical relativity bends the metric's values on a grid; Section 8.6 runs a homogeneous version driven by quantum matter. What a fixed lattice leaves open is gravity at the scale of its own spacing, which a rule that grows its own spacetime would settle too.",
    label='Section 7, trimmed (T13)')
# ------------------------------------------------------------------ references
newrefs = [
 "Ahmed, M., & Rideout, D. (2010). Indications of de Sitter spacetime from classical sequential growth dynamics of causal sets. *Physical Review D*, 81, 083528. https://arxiv.org/abs/0909.4771",
 "Benincasa, D. M. T., & Dowker, F. (2010). The scalar curvature of a causal set. *Physical Review Letters*, 104, 181301. https://arxiv.org/abs/1001.2725",
 "Glaser, L., & Surya, S. (2013). Towards a definition of locality in a manifoldlike causal set. *Physical Review D*, 88, 124026. https://arxiv.org/abs/1309.3403",
 "Rideout, D. P., & Sorkin, R. D. (2001). Evidence for a continuum limit in causal set dynamics. *Physical Review D*, 63, 104011. https://arxiv.org/abs/gr-qc/0003117",
 "Surya, S. (2012). Evidence for the continuum in 2D causal set quantum gravity. *Classical and Quantum Gravity*, 29, 132001. https://arxiv.org/abs/1110.6244",
 "Ambjørn, J., Görlich, A., Jurkiewicz, J., & Loll, R. (2008). Planckian birth of a quantum de Sitter universe. *Physical Review Letters*, 100, 091304. https://arxiv.org/abs/0712.2485",
 "Ambjørn, J., Jurkiewicz, J., & Loll, R. (2004). Emergence of a 4D world from causal quantum gravity. *Physical Review Letters*, 93, 131301. https://arxiv.org/abs/hep-th/0404156",
 "Ballesteros, G., Redondo, J., Ringwald, A., & Tamarit, C. (2017). Unifying inflation with the axion, dark matter, baryogenesis, and the seesaw mechanism. *Physical Review Letters*, 118, 071802. https://arxiv.org/abs/1608.05414",
 "Cao, Z., et al. (LHAASO Collaboration). (2024). Stringent tests of Lorentz invariance violation from LHAASO observations of GRB 221009A. *Physical Review Letters*, 133, 071501. https://arxiv.org/abs/2402.06009",
 "Chaitin, G. J. (1974). Information-theoretic limitations of formal systems. *Journal of the ACM*, 21, 403-424.",
 "García-Etxebarria, I., & Montero, M. (2019). Dai-Freed anomalies in particle physics. *Journal of High Energy Physics*, 2019(8), 003. https://arxiv.org/abs/1808.00009",
 "Gibbons, G. W., Hawking, S. W., & Perry, M. J. (1978). Path integrals and the indefiniteness of the gravitational action. *Nuclear Physics B*, 138, 141-150.",
 "Ginsparg, P. H., & Wilson, K. G. (1982). A remnant of chiral symmetry on the lattice. *Physical Review D*, 25, 2649-2657.",
 "Golterman, M., & Shamir, Y. (2023). Propagator zeros and lattice chiral gauge theories. arXiv:2311.12790. https://arxiv.org/abs/2311.12790",
 "Golterman, M., & Shamir, Y. (2026). Symmetric mass generation and the Nielsen-Ninomiya theorem. arXiv:2603.15985. https://arxiv.org/abs/2603.15985",
 "Lüscher, M. (1998). Exact chiral symmetry on the lattice and the Ginsparg-Wilson relation. *Physics Letters B*, 428, 342-345. https://arxiv.org/abs/hep-lat/9802011",
 "Lüscher, M. (1999). Abelian chiral gauge theories on the lattice with exact gauge invariance. *Nuclear Physics B*, 549, 295-334. https://arxiv.org/abs/hep-lat/9811032",
 "Neuberger, H. (1998). Exactly massless quarks on the lattice. *Physics Letters B*, 417, 141-144. https://arxiv.org/abs/hep-lat/9707022",
 "Regge, T. (1961). General relativity without coordinates. *Il Nuovo Cimento*, 19, 558-571.",
 "Wang, J., & Wen, X.-G. (2020). A non-perturbative definition of the Standard Models. *Physical Review Research*, 2, 023356. https://arxiv.org/abs/1809.11171",
 "Wilson, K. G. (1974). Confinement of quarks. *Physical Review D*, 10, 2445-2459.",
 "Ambjørn, J., & Loll, R. (1998). Non-perturbative Lorentzian quantum gravity, causality and topology change. *Nuclear Physics B*, 536, 407-434. https://arxiv.org/abs/hep-th/9805108",
 "Becher, P., & Joos, H. (1982). The Dirac-Kähler equation and fermions on the lattice. *Zeitschrift für Physik C*, 15, 343-365.",
 "Catterall, S., Laiho, J., & Unmuth-Yockey, J. (2018). Kähler-Dirac fermions on Euclidean dynamical triangulations. *Physical Review D*, 98, 114503. https://arxiv.org/abs/1810.10626",
 "Catterall, S. (2021). Chiral lattice fermions from staggered fields. *Physical Review D*, 104, 014503. https://arxiv.org/abs/2010.02290",
 "Butt, N., Catterall, S., & Hasenfratz, A. (2025). Symmetric mass generation with four SU(2) doublet fermions. *Physical Review Letters*, 134, 031602. https://arxiv.org/abs/2409.02062",
 "Hamber, H. W. (1992). Phases of four-dimensional simplicial quantum gravity. *Physical Review D*, 45, 507-512.",
 "Hamber, H. W. (2019). Vacuum condensate picture of quantum gravity. *Symmetry*, 11(1), 87. https://arxiv.org/abs/1707.08188",
 "Hamber, H. W. (2015). Scaling exponents for lattice quantum gravity in four dimensions. *Physical Review D*, 92, 064017. https://arxiv.org/abs/1506.07795",
 "Shaposhnikov, M., & Wetterich, C. (2010). Asymptotic safety of gravity and the Higgs boson mass. *Physics Letters B*, 683, 196-200. https://arxiv.org/abs/0912.0208",
 "Hasenfratz, A. (2022). Emergent strongly coupled ultraviolet fixed point in four dimensions with eight Kähler-Dirac fermions. *Physical Review D*, 106, 014513. https://arxiv.org/abs/2204.04801",
 "Wang, J. (2022). CT or P problem and symmetric gapped fermion solution. *Physical Review D*, 106, 125007.",
 "Zeng, M., Zhu, Z., Wang, J., & You, Y.-Z. (2022). Symmetric mass generation in the 1+1 dimensional chiral fermion 3-4-5-0 model. *Physical Review Letters*, 128, 185301. https://arxiv.org/abs/2202.12355",
 "Lee, T. D., & Yang, C. N. (1956). Question of parity conservation in weak interactions. *Physical Review*, 104, 254-258.",
 "Kähler, E. (1962). Der innere Differentialkalkül. *Rendiconti di Matematica e delle sue Applicazioni*, 21, 425-523.",
 "Dai, M., Laiho, J., Schiffer, M., & Unmuth-Yockey, J. (2021). Newtonian binding from lattice quantum gravity. arXiv:2102.04492. https://arxiv.org/abs/2102.04492",
 "Green, S. R., & Wald, R. M. (2014). How well is our universe described by an FLRW model? *Classical and Quantum Gravity*, 31, 234003.",
 "Hamber, H. W., & Williams, R. M. (1995). Newtonian potential in quantum Regge gravity. *Nuclear Physics B*, 435, 361-397. https://arxiv.org/abs/hep-th/9406163",
 "Seifert, A., Lane, Z. G., Galoppo, M., Ridden-Harper, R., & Wiltshire, D. L. (2025). Supernovae evidence for foundational change to cosmological models. *Monthly Notices of the Royal Astronomical Society: Letters*, 537, L55-L60.",
 "Wiltshire, D. L. (2007). Cosmic clocks, cosmic variance and cosmic averages. *New Journal of Physics*, 9, 377.",
 "Donoghue, J. F. (1994). General relativity as an effective field theory: The leading quantum corrections. *Physical Review D*, 50, 3874-3888. https://arxiv.org/abs/gr-qc/9405057",
]
head, refs = md.split('## References\n')
items = [x for x in refs.strip().split('\n\n') if x.strip()]
import unicodedata
def key(r): return unicodedata.normalize('NFKD', r).encode('ascii', 'ignore').decode().lower()
items = sorted(items + newrefs, key=key)
md = head + '## References\n\n' + '\n\n'.join(items) + '\n'
CHANGES.append(('References', '', '\n'.join(newrefs)))

# ------------------------------------------------------------------ the abstract cut to 250 words (SA74 T2), replacing the whole paragraph
abs_new = open('/home/claude/sa74/t5/abstract_t6.txt').read().strip()
i0 = md.index('## Abstract\n\n') + len('## Abstract\n\n'); i1 = md.index('\n\nKeywords:')
CHANGES.append(('Abstract cut to 250 words or fewer', md[i0:i1], abs_new))
md = md[:i0] + abs_new + md[i1:]


# ------------------------------------------------------------------ SA74 T19: merge of the SO74 (OpenAI) audit revision, then Perplexity's framing fixes
SO = json.load(open('/home/claude/sa74/so74/pkg/A42_T19_Lean/PAPER/CHANGES.json'))
for i, c in enumerate(SO):
    new = c['new'].split('\n')[0] if i == 0 else c['new']   # title only; the build adds the byline
    rep(c['old'], new, label=f'SO74 T19 merge {i}')
rep("![Regge audit trajectories. Each curve averages 100 sweeps for display; all diagnostics use the unaveraged measurement histories. Gray shading marks equilibration sweeps.](regge_audit.png)\n\n", "",
    label='Section 8.9, figure kept in the audit archive (T19)')
rep("The original histories remain unavailable; this audit supplies its own complete histories rather than claiming to reconstruct them.",
    "The original SA74 histories are kept in its history archive, and this audit supplies its own complete histories and trajectory figure in the audit archive.",
    label='Section 8.9, histories (T19)')
rep("in 4 of four parameter groups", "in all four parameter groups", label='Conclusion, consistency (T19)')
rep("Spacetime is the causal order of discrete updates, so a theory of everything is an update rule,",
    "If spacetime is the causal order of discrete updates, a theory of everything is an update rule,", label='Abstract, premise as hypothesis (T19)')
rep("Both approaches take change to come in steps, for a reason that holds in any universe: a cause produces its effect as the next state, and a continuum has no next state, since between any two points lies a third. The reason is an a priori proof in the logician's sense. Its premise is that change is produced, each state making the next;",
    "Both approaches take change to come in steps. This paper adopts that as its working hypothesis, motivated by an argument rather than a measurement: a cause produces its effect as the next state, and a continuum has no next state, since between any two points lies a third. The argument rests on the premise that change is produced, each state making the next;",
    label='Introduction, discreteness as hypothesis (T19)')
rep("### 8.7 The known laws on a supplied lattice", "### 8.7 A proposed specification on a supplied lattice", label='Section 8.7 heading (T19)')
rep("The target is the rule both routes reach.",
    "The target is the rule both routes reach. The paper's results are the census and its failures (Sections 3 to 6), the local-step benchmarks (Sections 8.1 to 8.6) and the Regge audit (Section 8.9); Section 8.7 is a proposal, a specification and its bit count, whose dynamics are untested as a whole.",
    label='Introduction, results and proposals (T19)')
# ------------------------------------------------------------------ SA74 T25 proposals, applied at T26 (Gary 21:50: "apply if good")
# Texts come from t5/props_T25.json, built from the rerun results; figures stay in the supplement, as for Section 8.9 (T19).
P25 = {p['name']: p for p in json.load(open('/home/claude/sa74/t5/props_T25.json'))}
sec = [b for b in P25['A42 Section 8.10 (new subsection after 8.9)']['plain_new'].split('\n\n') if not b.startswith('![')]
i2 = next(i for i, b in enumerate(sec) if b.startswith('The rerun repeats the Section 8.9 protocol'))
assert sec[i2].endswith('measurement sweeps.')
sec[i2] += ' Its code, fixed plans, per-sweep histories and figures are in the supplement.'
rep('## 9. Conclusion', '\n\n'.join(sec) + '\n\n\n## 9. Conclusion', label='Section 8.10, a curvature-squared term (T25)')
for name, label in [('A42 Section 8.7 table, gravity row', 'Section 8.7 table, gravity row (T25)'),
                    ('A42 Conclusion', 'Conclusion, curvature-squared rerun (T25)'),
                    ('A42 Abstract', 'Abstract, curvature-squared rerun (T25)')]:
    rep(P25[name]['old'], P25[name]['plain_new'], label=label)
rep('and the Regge audit (Section 8.9); Section 8.7 is a proposal', 'and the Regge audit and its rerun (Sections 8.9 and 8.10); Section 8.7 is a proposal',
    label='Introduction, results include the rerun (T25)')
rep('The independent execution in Section 8.9 tests its dependence on initial conditions and sampling history.',
    'The independent execution in Section 8.9 tests its dependence on initial conditions and sampling history, and Section 8.10 repeats it with a curvature-squared term.',
    label='Section 8.7, open questions, pointer to the rerun (T25)')
refs25 = P25['A42 References (three added)']['plain_new'].split('\n\n')
head, refs = md.split('## References\n')
items = [x for x in refs.strip().split('\n\n') if x.strip()]
for r in refs25: assert r not in items
items = sorted(items + refs25, key=key)
md = head + '## References\n\n' + '\n\n'.join(items) + '\n'
CHANGES.append(('References (T25)', '', '\n'.join(refs25)))

# ------------------------------------------------------------------ SA74 T27 (Gary 01:30, Nov 7: "u can apply whatever u think is best")
# The T23 proposals (empty-set growth, the semiorder theorem, large runs), with the abstract re-trimmed to 250 words; Hamber (1993)
# added beside Hamber (1992) in Section 8.10, from the paper itself (Table 1, Section 3.4); one "names" -> "labels" style fix.
rep("A rule's causal graph can therefore match four-dimensional spacetime only at large scales. Section 5 states the bound",
    "A rule's causal graph can therefore match four-dimensional spacetime only at large scales. Large scales mean large runs. A short rule may have no shortcut to what it grows, so its dimension is read by running it, and in four dimensions an interval 200 steps tall holds about $10^8$ elements; small runs can still rule a rule out, and a theorem about a rule's limit can stand in for size (Section 5.1). Section 5 states the bound",
    label='Introduction, large runs (T27)')
rep("Transitive percolation is not a four-dimensional spacetime at any size reached, and the lead it offered closes on the midpoint.\n\n",
    "Transitive percolation is not a four-dimensional spacetime at any size reached, and the lead it offered closes on the midpoint.\n\n"
    "A theorem explains the miss: any continuum limit of classical sequential growth, the class transitive percolation belongs to, is essentially a semiorder, so none is Minkowski space in any dimension of two or more (Brightwell and Georgiou, 2010). Growth from the empty set, each new set made of earlier ones (Ackermann, 1937), reaches the two extremes. With each earlier set joining at even odds it is transitive percolation at p = 1/2, and at 2,000 sets 99.8 percent of its pairs are related: nearly one chain. In Ackermann's own order, where set n holds the sets its binary digits mark, the longest chain stays at five sets while the ordering-fraction reading climbs from 5.0 at 256 sets to 10.5 at 65,536: nearly all width.\n\n",
    label='Section 5.1, the semiorder theorem and growth from the empty set (T27)')
rep("A growth law that favors manifold-like orders is the next candidate: weighting",
    "A growth law that favors manifold-like orders is the next candidate, and it must lie outside classical sequential growth (Section 5.1): weighting",
    label='Section 7, outside classical sequential growth (T27)')
rep("bring the two results to the same state, up to the names of new nodes.", "bring the two results to the same state, up to the labels of new nodes.",
    label='Section 4.6, style: labels (T27)')
rep("Hamber (1992) found his estimated critical exponents independent of a, and that the theories with and without the term may share a phase, which would give them a common unitary continuum limit.",
    "Hamber (1992) found his estimated critical exponents independent of a, and that the theories with and without the term may share a phase, which would give them a common unitary continuum limit. His 1993 runs separated the two cases (Hamber, 1993). With the term, on lattices up to 16⁴, the curvature exponent is δ = 0.624(6) at a = 0.005 and 0.628(20) at a = 0.02, which rules out a first-order transition there; without it, δ = 0.30(3), the transition seemed close to first order, and he placed a tricritical point between a = 0 and a = 0.005. His a = 0 lattices up to 32⁴ later gave ν = 0.334(4) (Section 8.7), close to the 0.325 that δ = 0.30 implies through δ = 4ν - 1, so the transitions with and without the term differ in exponent.",
    label='Section 8.10, Hamber (1993) beside Hamber (1992) (T27)')
abs27 = open('/home/claude/sa74/t5/abstract_t27.txt').read().strip()
i0 = md.index('## Abstract\n\n') + len('## Abstract\n\n'); i1 = md.index('\n\nKeywords:')
CHANGES.append(('Abstract, large runs, re-trimmed to 250 words (T27)', md[i0:i1], abs27))
md = md[:i0] + abs27 + md[i1:]
refs27 = ["Ackermann, W. (1937). Die Widerspruchsfreiheit der allgemeinen Mengenlehre. *Mathematische Annalen*, 114, 305-315.",
          "Brightwell, G., & Georgiou, N. (2010). Continuum limits for classical sequential growth models. *Random Structures and Algorithms*, 36, 218-250. https://doi.org/10.1002/rsa.20278"]
head, refs = md.split('## References\n')
items = [x for x in refs.strip().split('\n\n') if x.strip()]
for r in refs27: assert r not in items
items = sorted(items + refs27, key=key)
md = head + '## References\n\n' + '\n\n'.join(items) + '\n'
CHANGES.append(('References (T27)', '', '\n'.join(refs27)))

# ------------------------------------------------------------------ SA74 T28 (Gary 01:30, Nov 7: "u can apply whatever u think is best"; 03:24: "u can continue past 4am without limit")
# The night's runs on larger grids (Section 8.10, two new paragraphs), the costs clause, Table 12's gravity row, the conclusion, and
# one reference. Numbers from regge/r2/night: multihist.py, reweight.py, analyze_fss.py, hysteresis.py, defects.py, analyze_au.py, AS3/AS4/AS5.
T28_HYST = open('/home/claude/sa74/t5/T28/fill_hyst.txt').read().strip()
_MH = json.load(open('/home/claude/sa74/regge/r2/night/multihist_T28.json')); _MH2 = json.load(open('/home/claude/sa74/regge/r2/night/multihist_skip50_T28.json'))
_pk = [v['peak_chi'] for v in _MH['per_L'].values()]; _cr = [v['crossing'] for v in _MH['per_L'].values()]
def _pw(p): return ('%.2f' % p['power']).replace('0.', '0.', 1) + '(%d)' % round(100 * p['se'])
_F = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results_T28.json')); _tau = {}
for _key, _rows in _F['per_chain'].items():
    _L, _k = _key.split('_')
    if 4 <= int(_L) <= 10 and 0.2425 <= float(_k) <= 0.2575:
        for _r in _rows:
            if _r['n'] >= 800: _tau.setdefault(int(_L), []).append(_r['n'] / max(_r['ess'], 1e-9))
_tm = [sorted(v)[len(v) // 2] if len(v) % 2 else 0.5 * (sorted(v)[len(v) // 2 - 1] + sorted(v)[len(v) // 2]) for v in _tau.values()]
T28_N = dict(chi_lo='%.1f' % min(_pk), chi_hi='%.1f' % max(_pk), pw=_pw(_MH['peak_power']), pw2=_pw(_MH2['peak_power']),
             cr_lo='%.4f' % min(_cr), cr_hi='%.4f' % max(_cr), tau_lo='%d' % (10 * math.floor(min(_tm) / 10)), tau_hi='%d' % (10 * math.ceil(max(_tm) / 10)))
assert _MH['L_range'] == [4, 10] and _MH2['skipfrac'] == 0.5
T28_P5 = ("Larger grids test whether the transition is critical. At a critical point the curvature fluctuation $χ = 4\\,\\mathrm{var}(Σ_h δ_h A_h)/\\langle Σ_s V_s \\rangle$ "
          "grows with the grid as $L^{2/ν - 4}$ (Hamber, 2015), as $L^{0.94}$ for the ν = 0.405 that δ = 0.62 implies, and the chains slow down as the grid grows. "
          "Flat-start chains on grids of 4⁴ to 10⁴ sites at k = 0.245 to 0.255 show neither. Combined across couplings by reweighting (Ferrenberg and Swendsen, 1989), "
          "the largest fluctuation between k = 0.2415 and 0.258 is <<chi_lo>> to <<chi_hi>> on every grid, a fitted power of <<pw>> in L, or <<pw2>> from the second halves "
          "of the chains, and autocorrelation times measured within chains near the crossing stay between <<tau_lo>> and <<tau_hi>> sweeps, with no rise from 4⁴ to 10⁴. "
          "The average curvature crosses zero between k = <<cr_lo>> and <<cr_hi>> on grids of 4⁴ to 10⁴, a little above Hamber's fitted $k_c$, and its slope, "
          "which the fluctuation sets, does not steepen with size. Hamber (1993) fitted his fluctuation at k = 0.244 on grids up to 16⁴ as $L^{0.64(11)}$ with a large "
          "correction for small grids; that fitted curve rises by only about a quarter from L = 4 to 16, so these grids cannot yet tell its growth from none. "
          "Grids below 4⁴, which Hamber (1993) left out because the term couples next-nearest neighbors, are left out here too, and so are chains on 12⁴ "
          "and 14⁴, which still drift after 2,000 sweeps.")
for _k, _v in T28_N.items(): T28_P5 = T28_P5.replace('<<' + _k + '>>', _v)
assert '<<' not in T28_P5
T28_P6 = ("The slowest motion near the crossing is local. Now and then the edges of a small cluster grow to two or three times the typical length, "
          "and while the cluster lasts, often thousands of sweeps, the curvature sits higher; so on small grids chains at one coupling can "
          "differ by several times their internal errors. Within flat-start chains the longest edge and the curvature rise and fall together, and the "
          "curvature shift per cluster falls as the volume grows, as a local defect's share of an average should. Rough starts carry many such clusters. On 8⁴ "
          "they reach the flat-start values at k = 0.250 and 0.255 within 3,000 sweeps" + T28_HYST + ". The mean edge length matches the flat starts "
          "throughout, so the gap between starts comes from these local defects, with no sign of a second bulk phase. Toward smaller a, mixing worsens: "
          "near the crossing, flat and rough starts on 4⁴ still pass the screen at a = 0.003 (split R-hat 1.04), agree but drift slowly at a = 0.002 (1.14), "
          "and disagree at a = 0.001 (1.21), as Hamber's tricritical point between a = 0 and 0.005 would lead one to expect. At a = 0.02, where Hamber (1993) "
          "found $k_c = 0.4244(36)$, the fluctuation is a tenth as large and grows by at most a fifth from 4⁴ to 6⁴, where critical growth with his δ = 0.628 "
          "would add nearly half, and at k = 0.44 the curvature is still negative on 4⁴ and 6⁴. The term buys an ensemble that mixes; on these grids its "
          "transition is a smooth crossing whose slowest motion is local, and the critical point a continuum limit needs has not appeared. Plans, fixed before "
          "their runs, per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T28.")
rep("Past $k_c$, at k = 0.30 on 3⁴, the two starts begin to differ.\n\n",
    "Past $k_c$, at k = 0.30 on 3⁴, the two starts begin to differ.\n\n" + T28_P5 + "\n\n" + T28_P6 + "\n\n",
    label='Section 8.10, larger grids and the rough starts (T28, two new paragraphs)')
rep("lattices of 3⁴ and 4⁴ sites leave finite-size scaling, the critical exponent and a Newtonian potential to larger runs;",
    "grids up to 10⁴ show no critical growth yet, which leaves the critical exponent, the continuum limit and a Newtonian potential to larger runs;",
    label='Section 8.10, costs (T28)')
rep("with a curvature-squared term both starts settle at the audit couplings (Section 8.10); larger-lattice critical behavior in cited work |",
    "with a curvature-squared term both starts settle at the audit couplings (Section 8.10); larger-lattice critical behavior in cited work; no critical growth up to 10⁴ here (Section 8.10) |",
    label='Table 12, gravity row (T28)')
rep("and the average curvature rises to zero beside the known critical coupling, where mixing slows.",
    "and the average curvature rises to zero beside the known critical coupling, where mixing slows. On grids up to 10⁴ that transition shows no critical growth, so the continuum limit the gravity sector needs remains open.",
    label='Conclusion (T28)')
refs28 = ["Ferrenberg, A. M., & Swendsen, R. H. (1989). Optimized Monte Carlo data analysis. *Physical Review Letters*, 63, 1195-1198. https://doi.org/10.1103/PhysRevLett.63.1195"]
head, refs = md.split('## References\n')
items = [x for x in refs.strip().split('\n\n') if x.strip()]
for r in refs28: assert r not in items
items = sorted(items + refs28, key=key)
md = head + '## References\n\n' + '\n\n'.join(items) + '\n'
CHANGES.append(('References (T28)', '', '\n'.join(refs28)))

# ------------------------------------------------------------------ SA74 T29 (Nov 7 morning; Gary 01:30: "u can apply whatever u think is best"; 03:24: "u can continue past 4am without limit")
# Folds in the runs that finished after T28: plan A1's 12^4 chains (Mac runner 2, fixed 03:27 before its runs), AU2 (a = 0.02 up to 8^4),
# FSS4 and the rest of FSS3. Numbers from regge/r2/night: multihist_T29.json and multihist_skip50_T29.json (A1INC=1, LMAX=12, KPAD=0.002),
# a1_results.json, au_results.json; regge/r2/fss/fss_results.json (rerun at T29).
_ND = '/home/claude/sa74/regge/r2/night/'
_M9 = json.load(open(_ND + 'multihist_T29frozen.json')); _M9b = json.load(open(_ND + 'multihist_skip50_T29frozen.json'))   # T29 inputs, rebuilt at T30 from the chains T29 used
_A1 = json.load(open(_ND + 'a1_results_T29frozen.json')); _AU = json.load(open(_ND + 'au_results_T29.json'))   # a = 0.02 as at T29 (frozen at T31)
assert _M9['L_range'] == [4, 12] and _M9b['skipfrac'] == 0.5 and sorted(int(L) for L in _M9['per_L']) == [4, 5, 6, 8, 10, 12]
_a1_12 = [r for r in _A1['chains'] if r['L'] == 12 and r['finished']]
assert len(_a1_12) >= 6 and _M9['per_L']['12']['chains'] == len(_a1_12), (len(_a1_12), _M9['per_L']['12']['chains'])
def _p9(p): return '%.2f' % p['power'] + '(%d)' % round(100 * p['se'])
def _fk(k): return [v['fixed_k'][k]['chi'] for v in _M9['per_L'].values()]
_P = _A1['fixed_k_power']
assert _P['0.2480']['L'] == [4, 5, 6, 8, 10, 12] and _P['0.2520']['L'] == [4, 5, 6, 8, 10, 12]
_F9 = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results_T29.json')); _tau9 = {}
for _key, _rows in _F9['per_chain'].items():
    _L, _k = _key.split('_')
    if 4 <= int(_L) <= 10 and 0.2425 <= float(_k) <= 0.2575:
        for _r in _rows:
            if _r['n'] >= 800: _tau9.setdefault(int(_L), []).append(_r['n'] / max(_r['ess'], 1e-9))
_tau9[12] = [r['tau'] for r in _a1_12 if 0.2425 <= r['k'] <= 0.2575]
_tm9 = [float(sorted(v)[len(v) // 2]) if len(v) % 2 else 0.5 * (sorted(v)[len(v) // 2 - 1] + sorted(v)[len(v) // 2]) for v in _tau9.values()]
_c14 = {r['k']: r['chi'] for r in _A1['chains'] if r['L'] == 14 and r['finished']}
_pk9 = [v['peak_chi'] for v in _M9['per_L'].values()]; _cr9 = [v['crossing'] for v in _M9['per_L'].values()]
T29_N = dict(f248_lo='%.1f' % min(_fk('0.2480')), f248_hi='%.1f' % max(_fk('0.2480')), p248=_p9(_P['0.2480']),
             f252_lo='%.1f' % min(_fk('0.2520')), f252_hi='%.1f' % max(_fk('0.2520')), p252=_p9(_P['0.2520']),
             chi_lo='%.1f' % min(_pk9), chi_hi='%.1f' % max(_pk9), pw=_p9(_M9['peak_power']), pw2=_p9(_M9b['peak_power']),
             tau_lo='%d' % (10 * math.floor(min(_tm9) / 10)), tau_hi='%d' % (10 * math.ceil(max(_tm9) / 10)),
             cr_lo='%.4f' % min(_cr9), cr_hi='%.4f' % max(_cr9), c14a='%.1f' % _c14[0.2475], c14b='%.1f' % _c14[0.2525])
T29_OLD5 = md[md.index("Larger grids test whether the transition is critical."):]; T29_OLD5 = T29_OLD5[:T29_OLD5.index('\n\n')]
T29_P5 = ("Larger grids test whether the transition is critical. At a critical point the curvature fluctuation $χ = 4\\,\\mathrm{var}(Σ_h δ_h A_h)/\\langle Σ_s V_s \\rangle$ "
          "grows with the grid as $L^{2/ν - 4}$ (Hamber, 2015), as $L^{0.94}$ for the ν = 0.405 that δ = 0.62 implies, and the chains slow down as the grid grows. "
          "Flat-start chains on grids of 4⁴ to 12⁴ sites at k = 0.244 to 0.255 show neither. Combined across couplings by reweighting (Ferrenberg and Swendsen, 1989), "
          "the fluctuation at k = 0.248 is <<f248_lo>> to <<f248_hi>> on every grid, a fitted power of <<p248>> in L, and at k = 0.252 it is <<f252_lo>> to <<f252_hi>>, "
          "a power of <<p252>>; its largest value near the crossing, <<chi_lo>> to <<chi_hi>>, gives a power of <<pw>>, or <<pw2>> from the second halves of the chains. "
          "Autocorrelation times measured within chains near the crossing stay between <<tau_lo>> and <<tau_hi>> sweeps, with no rise from 4⁴ to 12⁴. "
          "The average curvature crosses zero between k = <<cr_lo>> and <<cr_hi>> on grids of 4⁴ to 12⁴, a little above Hamber's fitted $k_c$, and its slope, "
          "which the fluctuation sets, does not steepen with size. Hamber (1993) fitted his fluctuation at k = 0.244 on grids up to 16⁴ as $L^{0.64(11)}$ with a large "
          "correction for small grids; that fitted curve rises by only about a sixth from L = 4 to 12 and a quarter to 16, so these grids cannot yet tell its growth from none. "
          "Grids below 4⁴, which Hamber (1993) left out because the term couples next-nearest neighbors, are left out here too. On 12⁴, where chains of 2,000 sweeps "
          "still drift, each chain runs 3,000 sweeps at its coupling before 4,000 measured ones; the two 14⁴ chains, of 1,600 measured sweeps each, are too short "
          "for the fits, and their fluctuations, <<c14a>> at k = 0.2475 and <<c14b>> at 0.2525, show no growth either.")
for _k, _v in T29_N.items(): T29_P5 = T29_P5.replace('<<' + _k + '>>', _v)
assert '<<' not in T29_P5
rep(T29_OLD5, T29_P5, label='Section 8.10, larger grids (T29: 12⁴ from plan A1; fluctuation at fixed couplings; the 14⁴ chains)')
# a = 0.02 (AU1 + AU2): chi at k = 0.46 on 4^4, 6^4, 8^4; zero crossings on 4^4 and 6^4 by interpolation, on 8^4 from R at k = 0.46 and the fluctuation slope
_au = {(r['L'], round(r['k'], 3)): r for r in _AU['chains']}
_x8 = _au[(8, 0.46)]; _cr8 = 0.46 - _x8['R'] / _x8['slope']
assert _AU['crossing']['8'] is None and abs(_x8['R']) < 2 * _x8['R_se'] and _AU['crossing']['4'] > _AU['crossing']['6'] > _cr8 > 0.4244
_g = max(_au[(L, k)]['chi'] / _au[(4, k)]['chi'] for (L, k) in _au if L > 4 and (4, k) in _au)
assert 1.0 < _g < 1.2, _g
T29_AU_OLD = ("the fluctuation is a tenth as large and grows by at most a fifth from 4⁴ to 6⁴, where critical growth with his δ = 0.628 "
              "would add nearly half, and at k = 0.44 the curvature is still negative on 4⁴ and 6⁴.")
T29_AU_NEW = ("the fluctuation is a tenth as large and grows by at most a fifth on grids from 4⁴ to 8⁴ (at k = 0.46 it is %.2f, %.2f and %.2f on 4⁴, 6⁴ and 8⁴), "
              "where critical growth with his δ = 0.628 would nearly double it, and the curvature crosses zero at k = %.3f, %.3f and about %.3f "
              "on 4⁴, 6⁴ and 8⁴, drifting down with size but still above his $k_c$."
              % (_au[(4, 0.46)]['chi'], _au[(6, 0.46)]['chi'], _au[(8, 0.46)]['chi'], _AU['crossing']['4'], _AU['crossing']['6'], _cr8))
rep(T29_AU_OLD, T29_AU_NEW, label='Section 8.10, a = 0.02 (T29: AU2, grids up to 8⁴)')
rep("The mean edge length matches the flat starts throughout, so the gap between starts comes from these local defects, with no sign of a second bulk phase.",
    "The mean edge length starts up to 0.6% above the flat starts and settles to them as the long edges shrink, so the gap between starts comes from these local defects, with no sign of a second bulk phase.",
    label='Section 8.10, rough starts: mean edge length (T29 precision fix)')
rep("Plans, fixed before their runs, per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T28.",
    "Plans, fixed before their runs, per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T29.",
    label='Section 8.10, supplement (T29)')
rep("grids up to 10⁴ show no critical growth yet, which leaves the critical exponent,",
    "grids up to 12⁴ show no critical growth yet, which leaves the critical exponent,", label='Section 8.10, costs (T29)')
rep("no critical growth up to 10⁴ here (Section 8.10) |", "no critical growth up to 12⁴ here (Section 8.10) |", label='Table 12, gravity row (T29)')
rep("On grids up to 10⁴ that transition shows no critical growth,", "On grids up to 12⁴ that transition shows no critical growth,", label='Conclusion (T29)')

# ------------------------------------------------------------------ SA74 T30 (Nov 7, Gary 10:15: "can connect to mac and server as needed? continue.")
# Plan A1 finished at 10:08: 16^4 at k = 0.248 and 0.252 (one chain each), the 12^4 flat-start replicas (k = 0.246 to 0.254) and the 8^4 chains.
# 16^4 joins the fixed-coupling comparison; "up to 12^4" becomes "up to 16^4". Numbers from regge/r2/night: multihist_T30.json and
# multihist_skip50_T30.json (A1INC=1, LMAX=12, KPAD=0.002; all finished A1 chains), a1_results_T30.json (single-chain chi errors: the larger of
# the five-block and the effective-sample-size estimates); regge/r2/fss/fss_results_T30.json.
import numpy as np
_MT = json.load(open(_ND + 'multihist_T30.json')); _MTb = json.load(open(_ND + 'multihist_skip50_T30.json')); _AT = json.load(open(_ND + 'a1_results_T30.json'))
assert _MT['L_range'] == [4, 12] and _MTb['skipfrac'] == 0.5 and sorted(int(L) for L in _MT['per_L']) == [4, 5, 6, 8, 10, 12]
_PT = _AT['fixed_k_power']
assert _PT['0.2480']['L'] == [4, 5, 6, 8, 10, 12, 16] and _PT['0.2520']['L'] == [4, 5, 6, 8, 10, 12, 16]
_c16 = {r['k']: r for r in _AT['chains'] if r['L'] == 16 and r['finished']}
assert sorted(_c16) == [0.248, 0.252]
_fk30 = lambda k: [v['chi'] for v in _AT['fixed_k'][k].values()]
_FT = json.load(open('/home/claude/sa74/regge/r2/fss/fss_results_T30.json')); _tauT = {}
for _key, _rows in _FT['per_chain'].items():
    _L, _k = _key.split('_')
    if 4 <= int(_L) <= 10 and 0.2425 <= float(_k) <= 0.2575:
        for _r in _rows:
            if _r['n'] >= 800: _tauT.setdefault(int(_L), []).append(_r['n'] / max(_r['ess'], 1e-9))
for _L in (12, 16): _tauT[_L] = [r['tau'] for r in _AT['chains'] if r['L'] == _L and r['finished'] and r['chain'].startswith('A1_') and 0.2425 <= r['k'] <= 0.2575]
_tmT = [float(np.median(v)) for v in _tauT.values()]
_cr16 = 0.248 + (0.252 - 0.248) * (0 - _c16[0.248]['R']) / (_c16[0.252]['R'] - _c16[0.248]['R'])
_crT = [v['crossing'] for v in _MT['per_L'].values()] + [_cr16]
_pkT = [v['peak_chi'] for v in _MT['per_L'].values()]
_r12 = {x['chain']: x for x in _AT['chains'] if x['L'] == 12 and x['finished']}
_d12 = [(_r12['A1_L12_k%.3f_r1' % k]['R'] - _r12['A1_L12_k%.3f_r2' % k]['R'], (_r12['A1_L12_k%.3f_r1' % k]['R'] - _r12['A1_L12_k%.3f_r2' % k]['R']) / math.hypot(_r12['A1_L12_k%.3f_r1' % k]['R_se'], _r12['A1_L12_k%.3f_r2' % k]['R_se'])) for k in (0.244, 0.246, 0.248, 0.250, 0.252, 0.254)]
assert all(d < 0 for d, z in _d12) and 0.0009 < min(-d for d, z in _d12) and max(-d for d, z in _d12) < 0.0026 and 0.9 < min(-z for d, z in _d12) and max(-z for d, z in _d12) < 2.0, _d12
T30_N = dict(f248_lo='%.1f' % min(_fk30('0.2480')), f248_hi='%.1f' % max(_fk30('0.2480')), p248=_p9(_PT['0.2480']),
             f252_lo='%.1f' % min(_fk30('0.2520')), f252_hi='%.1f' % max(_fk30('0.2520')), p252=_p9(_PT['0.2520']),
             chi_lo='%.1f' % min(_pkT), chi_hi='%.1f' % max(_pkT), pw=_p9(_MT['peak_power']), pw2=_p9(_MTb['peak_power']),
             tau_lo='%d' % (10 * math.floor(min(_tmT) / 10)), tau_hi='%d' % (10 * math.ceil(max(_tmT) / 10)),
             cr_lo='%.4f' % min(_crT), cr_hi='%.4f' % max(_crT), c14a=T29_N['c14a'], c14b=T29_N['c14b'])
T30_OLD5 = md[md.index("Larger grids test whether the transition is critical."):]; T30_OLD5 = T30_OLD5[:T30_OLD5.index('\n\n')]
T30_P5 = ("Larger grids test whether the transition is critical. At a critical point the curvature fluctuation $χ = 4\\,\\mathrm{var}(Σ_h δ_h A_h)/\\langle Σ_s V_s \\rangle$ "
          "grows with the grid as $L^{2/ν - 4}$ (Hamber, 2015), as $L^{0.94}$ for the ν = 0.405 that δ = 0.62 implies, and the chains slow down as the grid grows. "
          "Flat-start chains on grids of 4⁴ to 16⁴ sites at k = 0.244 to 0.255 show neither. Combined across couplings by reweighting (Ferrenberg and Swendsen, 1989) "
          "on grids up to 12⁴, with one chain per coupling on 16⁴, the fluctuation at k = 0.248 is <<f248_lo>> to <<f248_hi>> on every grid, a fitted power of <<p248>> in L, "
          "and at k = 0.252 it is <<f252_lo>> to <<f252_hi>>, a power of <<p252>>. Its largest value near the crossing, in a range that includes Hamber's $k_c$, is "
          "<<chi_lo>> to <<chi_hi>> on grids up to 12⁴, a power of <<pw>>, or <<pw2>> from the second halves of the chains. "
          "Autocorrelation times measured within chains near the crossing stay between <<tau_lo>> and <<tau_hi>> sweeps, with no rise from 4⁴ to 16⁴. "
          "The average curvature crosses zero between k = <<cr_lo>> and <<cr_hi>> on grids of 4⁴ to 16⁴, a little above Hamber's fitted $k_c$, and its slope, "
          "which the fluctuation sets, does not steepen with size. Hamber (1993) fitted his fluctuation at k = 0.244 on grids up to 16⁴ as $L^{0.64(11)}$ with a large "
          "correction for small grids; that fitted curve rises by only about a quarter from L = 4 to 16, so these grids cannot yet tell its growth from none. "
          "Grids below 4⁴, which Hamber (1993) left out because the term couples next-nearest neighbors, are left out here too. On 12⁴ and 16⁴, where shorter chains "
          "still drift, each chain runs 3,000 sweeps at its coupling before 3,000 or 4,000 measured ones. Even so, on 12⁴ the chains continued from earlier runs sit "
          "0.001 to 0.0025 below fresh flat starts in the average curvature at all six couplings, one to two standard errors each, a small memory of the start; "
          "their fluctuations show no such pattern. The two 14⁴ chains, of 1,600 measured sweeps each, are too short for the fits, and their fluctuations, "
          "<<c14a>> at k = 0.2475 and <<c14b>> at 0.2525, show no growth either.")
for _k, _v in T30_N.items(): T30_P5 = T30_P5.replace('<<' + _k + '>>', _v)
assert '<<' not in T30_P5
rep(T30_OLD5, T30_P5, label='Section 8.10, larger grids (T30: 16⁴ from plan A1)')
rep("per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T29.",
    "per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T30.", label='Section 8.10, supplement (T30)')
rep("grids up to 12⁴ show no critical growth yet, which leaves the critical exponent,",
    "grids up to 16⁴ show no critical growth yet, which leaves the critical exponent,", label='Section 8.10, costs (T30)')
rep("no critical growth up to 12⁴ here (Section 8.10) |", "no critical growth up to 16⁴ here (Section 8.10) |", label='Table 12, gravity row (T30)')
rep("On grids up to 12⁴ that transition shows no critical growth,", "On grids up to 16⁴ that transition shows no critical growth,", label='Conclusion (T30)')

# ------------------------------------------------------------------ SA74 T31 (Nov 9, Gary 08:57: "update and package all for a42 and code to verify claims in a zip")
# Plans A2 (16^4: k = 0.244 from the 0.248 and 0.252 configurations, second chains at 0.248 and 0.252 from the other coupling's configuration,
# 0.250; a = 0.02 on 10^4 and 12^4 at k = 0.45, 0.46), A3 (two 16^4 forks at 0.244), A4 (a = 0.02 at 0.47 on 10^4 and 12^4) and A5 (20^4 at
# 0.244, 0.248, 0.252 with the multi-core engine) finished on the Mac on Nov 7 by 19:33. Numbers from regge/r2/night: multihist_T30.json and
# multihist_skip50_T30.json (4^4 to 12^4, unchanged since T30), a1_results_T31.json (16^4 pooled per coupling by lineage, 20^4 single chains;
# fixed_k_power_lowbias; start_check), au_results_T31.json; regge/r2/fss/fss_results_T30.json (autocorrelation on 4^4 to 10^4, as at T30).
_A31 = json.load(open(_ND + 'a1_results_T31.json')); _U31 = json.load(open(_ND + 'au_results_T31.json'))
_P31 = _A31['fixed_k_power']; _P31b = _A31['fixed_k_power_lowbias']
for _k in ('0.2440', '0.2480', '0.2520'):
    assert _P31[_k]['L'] == [4, 5, 6, 8, 10, 12, 16, 20], (_k, _P31[_k]['L'])
    assert _A31['fixed_k'][_k]['16']['lineages'] == 2 and _A31['fixed_k'][_k]['20']['lineages'] == 1
_c31 = [r for r in _A31['chains'] if r['L'] >= 16 and r['finished']]
assert len(_c31) == 12 and all(r['chi'] is not None for r in _c31)
_fk31 = lambda k: [v['chi'] for v in _A31['fixed_k'][k].values()]
# autocorrelation near the crossing: medians per grid (4^4 to 12^4 as at T30; 16^4 and 20^4 from all their finished chains)
_tau31 = {L: v for L, v in _tauT.items() if L <= 12}
for _L in (16, 20): _tau31[_L] = [r['tau'] for r in _c31 if r['L'] == _L and 0.2425 <= r['k'] <= 0.2575]
_tm31 = [float(np.median(v)) for v in _tau31.values()]
# zero crossing of the average curvature on 16^4 and 20^4 from the chain averages per coupling (20^4: a short extrapolation above 0.252)
def _cross31(L):
    by = {}
    for r in _c31:
        if r['L'] == L: by.setdefault(r['k'], []).append(r['R'])
    p = sorted((k, float(np.mean(v))) for k, v in by.items())
    for (k1, r1), (k2, r2) in zip(p, p[1:]):
        if r1 < 0 <= r2: return k1 + (k2 - k1) * (-r1) / (r2 - r1)
    (k1, r1), (k2, r2) = p[-2], p[-1]; c = k1 + (k2 - k1) * (-r1) / (r2 - r1)
    assert r2 < 0 and c - k2 < 0.0005, (L, p, c); return c
_cr31 = [v['crossing'] for v in _MT['per_L'].values()] + [_cross31(16), _cross31(20)]
assert '%.4f' % min(_cr31) == T30_N['cr_lo'] and '%.4f' % max(_cr31) == T30_N['cr_hi']   # the 16^4 and 20^4 crossings fall inside the earlier range
# start check on 16^4: chains from another coupling's final configuration against flat starts
_s31 = _A31['start_check']; _z31 = [_s31['A1_L16_k0.248 - A2_L16_k0.248_r2']['z'], _s31['A1_L16_k0.252 - A2_L16_k0.252_r2']['z']]
assert max(abs(z) for z in _z31) < 1.2, _z31
# Hamber's (1993) fitted curve at k = 0.244: its rise from L = 4 to 20 as a power, and the gap to the fitted power here
_h31 = lambda L: (0.16 + 1.11 / L) * L ** 0.64
_hp31 = math.log(_h31(20) / _h31(4)) / math.log(5)
_g31a = (_hp31 - _P31['0.2440']['power']) / _P31['0.2440']['se']; _g31b = (_hp31 - _P31b['0.2440']['power']) / _P31b['0.2440']['se']
assert 0.35 < _h31(20) / _h31(4) - 1 < 0.42 and 2.0 < _g31b < _g31a < 3.5, (_h31(20) / _h31(4), _g31a, _g31b)
T31_N = dict(crit='%.1f' % 5 ** 0.94,
             f244_lo='%.1f' % min(_fk31('0.2440')), f244_hi='%.1f' % max(_fk31('0.2440')), p244=_p9(_P31['0.2440']),
             f248_lo='%.1f' % min(_fk31('0.2480')), f248_hi='%.1f' % max(_fk31('0.2480')), p248=_p9(_P31['0.2480']),
             f252_lo='%.1f' % min(_fk31('0.2520')), f252_hi='%.1f' % max(_fk31('0.2520')), p252=_p9(_P31['0.2520']),
             chi_lo=T30_N['chi_lo'], chi_hi=T30_N['chi_hi'], pw=T30_N['pw'], pw2=T30_N['pw2'],
             tau_lo='%d' % (10 * math.floor(min(_tm31) / 10)), tau_hi='%d' % (10 * math.ceil(max(_tm31) / 10)),
             cr_lo='%.4f' % min(_cr31), cr_hi='%.4f' % max(_cr31), hp='%.2f' % _hp31, g1='%.1f' % _g31a, g2='%.1f' % _g31b,
             c14a=T29_N['c14a'], c14b=T29_N['c14b'])
T31_OLD5 = md[md.index("Larger grids test whether the transition is critical."):]; T31_OLD5 = T31_OLD5[:T31_OLD5.index('\n\n')]
T31_P5 = ("Larger grids test whether the transition is critical. At a critical point the curvature fluctuation $χ = 4\\,\\mathrm{var}(Σ_h δ_h A_h)/\\langle Σ_s V_s \\rangle$ "
          "grows with the grid as $L^{2/ν - 4}$ (Hamber, 2015), as $L^{0.94}$ for the ν = 0.405 that δ = 0.62 implies, a factor of <<crit>> from 4⁴ to 20⁴, "
          "and the chains slow down as the grid grows. Chains on grids of 4⁴ to 20⁴ sites at k = 0.244 to 0.255 show neither. Combined across couplings "
          "by reweighting (Ferrenberg and Swendsen, 1989) on grids up to 12⁴, with two chains from independent starts per coupling on 16⁴ (and two forks "
          "at k = 0.244) and one chain on 20⁴, the fluctuation at k = 0.244, Hamber's $k_c$, is <<f244_lo>> to <<f244_hi>> on every grid, a fitted power "
          "of <<p244>> in L; at k = 0.248 it is <<f248_lo>> to <<f248_hi>>, a power of <<p248>>, and at k = 0.252 it is <<f252_lo>> to <<f252_hi>>, "
          "a power of <<p252>>. Its largest value near the crossing, in a range that includes Hamber's $k_c$, is <<chi_lo>> to <<chi_hi>> on grids up to 12⁴, "
          "a power of <<pw>>, or <<pw2>> from the second halves of the chains. Autocorrelation times measured within chains near the crossing stay between "
          "<<tau_lo>> and <<tau_hi>> sweeps, with no rise from 4⁴ to 20⁴. The average curvature crosses zero between k = <<cr_lo>> and <<cr_hi>> on grids "
          "of 4⁴ to 20⁴, a little above Hamber's fitted $k_c$, and its slope, which the fluctuation sets, does not steepen with size. Hamber (1993) fitted "
          "his fluctuation at k = 0.244 on grids up to 16⁴ as $L^{0.64(11)}$ with a large correction for small grids; that fitted curve rises by about two "
          "fifths from L = 4 to 20, a power of <<hp>> over that range, <<g1>> standard errors above the power found here at his $k_c$, or <<g2>> after "
          "correcting the 16⁴ and 20⁴ values for the low bias of a variance from few independent sweeps, so these grids lean against his fitted growth "
          "without ruling it out. Grids below 4⁴, which Hamber (1993) left out because the term couples next-nearest neighbors, are left out here too. "
          "On 12⁴, 16⁴ and 20⁴, where shorter chains still drift, each chain runs at least 1,500 sweeps at its coupling before 1,500 to 4,000 measured "
          "ones. Even so, on 12⁴ the chains continued from earlier runs sit 0.001 to 0.0025 below fresh flat starts in the average curvature at all six "
          "couplings, one to two standard errors each, a small memory of the start; their fluctuations show no such pattern. On 16⁴, chains started from "
          "the final configuration of the coupling above or below agree with flat starts within about one standard error. The two 14⁴ chains, of 1,600 "
          "measured sweeps each, are too short for the fits, and their fluctuations, <<c14a>> at k = 0.2475 and <<c14b>> at 0.2525, show no growth either.")
for _k, _v in T31_N.items(): T31_P5 = T31_P5.replace('<<' + _k + '>>', _v)
assert '<<' not in T31_P5
rep(T31_OLD5, T31_P5, label='Section 8.10, larger grids (T31: 16⁴ at Hamber\'s k_c and start check; 20⁴)')
# a = 0.02 (plans A2 and A4 add 10^4 and 12^4): the crossing stops moving; T29's "drifting down with size" rested on 4^4 to 8^4 and is corrected
_au31 = {(r['L'], round(r['k'], 3)): r for r in _U31['chains']}
_x8b = _au31[(8, 0.46)]; _cr8b = 0.46 - _x8b['R'] / _x8b['slope']; _cu = {int(L): v for L, v in _U31['crossing'].items()}
assert _cu[8] is None and abs(_cr8b - _cr8) < 1e-12 and _cu[4] > 0.475 and all(0.466 < _cu[L] < 0.469 for L in (6, 10, 12))
_g31 = max(_au31[(L, k)]['chi'] / _au31[(4, k)]['chi'] for (L, k) in _au31 if L > 4 and (4, k) in _au31)
_nu02 = (1 + 0.628) / 4; _tr = 3 ** (2 / _nu02 - 4)
assert 1.0 < _g31 < 1.2 and 2.5 < _tr < 3.0, (_g31, _tr)
T31_AU_NEW = ("the fluctuation is a tenth as large and grows by at most a fifth on grids from 4⁴ to 12⁴ (at k = 0.46 it is %.2f, %.2f, %.2f, %.2f and %.2f "
              "on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴), where critical growth with his δ = 0.628 would nearly triple it, and the curvature crosses zero at k = %.3f on 4⁴, "
              "%.3f on 6⁴, about %.3f on 8⁴ (from one coupling) and %.3f and %.3f on 10⁴ and 12⁴, with no drift toward his $k_c$ beyond 6⁴."
              % tuple([_au31[(L, 0.46)]['chi'] for L in (4, 6, 8, 10, 12)] + [_cu[4], _cu[6], _cr8b, _cu[10], _cu[12]]))
rep(T29_AU_NEW, T31_AU_NEW, label='Section 8.10, a = 0.02 (T31: 10⁴ and 12⁴; the crossing stops moving, correcting "drifting down with size")')
rep("per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T30.",
    "per-sweep histories and analysis code are in the supplement A42_Regge_R2_Night_T31.", label='Section 8.10, supplement (T31)')
rep("grids up to 16⁴ show no critical growth yet, which leaves the critical exponent,",
    "grids up to 20⁴ show no critical growth yet, which leaves the critical exponent,", label='Section 8.10, costs (T31)')
rep("no critical growth up to 16⁴ here (Section 8.10) |", "no critical growth up to 20⁴ here (Section 8.10) |", label='Table 12, gravity row (T31)')
rep("On grids up to 16⁴ that transition shows no critical growth,", "On grids up to 20⁴ that transition shows no critical growth,", label='Conclusion (T31)')

# ------------------------------------------------------------------ SA74 T32 (Nov 10; Gary 19:03 forwarded the portfolio session's review of the T31 package)
# Editorial, no result changed: (1) the one authorial "we" goes; (2) the abstract states the larger-grid result and stays at 250 words;
# (3) every paragraph over 250 words, and Section 8.10's 248-word one, is broken at its topic changes; in Section 8.10's larger-grids
# paragraph the reweighting clause becomes a sentence of its own so that the new paragraph opens on its numbers; (4) the autocorrelation
# sentence says "median", which is what the 120 to 240 sweeps are (single chains on 12^4 run from 76 to 354); (5) the a = 0.02 values carry
# their standard errors (au_results_T31.json; crossings: the two bracketing curvature errors carried through the interpolation).
rep("None of these, as far as we can find, reports a short rule", "None of these appears to report a short rule", label='Introduction: the authorial "we" removed (T32)')
rep("transitive percolation reaches four on the chain at a million elements while its midpoint drifts to 5.8. From the top,",
    "transitive percolation's chain reaches four at a million elements while its midpoint drifts to 5.8. From the top,", label='Abstract: two words saved (T32)')
rep("All four groups in a new 24-chain Regge audit fail the mixing screen; a curvature-squared term, which quantum matter requires, makes all four pass.",
    "All four groups in a Regge audit fail the mixing screen; a curvature-squared term, which quantum matter requires, makes all pass, without critical growth through 20⁴.",
    label='Abstract: the larger-grid result (T32)')
_abs32 = md[md.index('## Abstract'):md.index('Keywords:')].replace('## Abstract', '').strip()
assert len(_abs32.split()) == 250, len(_abs32.split())
# Section 8.10, larger grids: the paragraph on the numbers opens with them; the reweighting clause follows as its own sentence
_rw = ("by reweighting (Ferrenberg and Swendsen, 1989) on grids up to 12⁴, with two chains from independent starts per coupling on 16⁴ (and two forks at k = 0.244) and one chain on 20⁴")
rep("Combined across couplings " + _rw + ", the fluctuation at k = 0.244, Hamber's $k_c$, is ", "At k = 0.244, Hamber's $k_c$, the fluctuation is ", label='Section 8.10, larger grids: the numbers first (T32)')
_m = re.search(r"a power of -?\d\.\d\d\(\d+\)\. Its largest value near the crossing,", md); assert _m and md.count(_m.group(0)) == 1
rep(_m.group(0), _m.group(0).replace(" Its largest value near the crossing,", " These values combine the chains across couplings " + _rw + ". The fluctuation's largest value near the crossing,"),
    label='Section 8.10, larger grids: the reweighting sentence (T32)')
rep("Autocorrelation times measured within chains near the crossing stay between ", "Median autocorrelation times of chains near the crossing stay between ", label='Section 8.10: "median" (T32)')
_mt = re.search(r"Median autocorrelation times of chains near the crossing stay between (\d+) and (\d+) sweeps, with no rise from 4⁴ to 20⁴\.", md); assert _mt
rep(_mt.group(0), _mt.group(0).replace(" sweeps, with no rise", " sweeps on every grid, with no rise"), label='Section 8.10: per grid (T32)')
assert all(int(_mt.group(1)) <= m <= int(_mt.group(2)) for m in _tm31)
# a = 0.02: standard errors
_e = lambda x, s, n=2: ('%.' + str(n) + 'f(%d)') % (x, max(1, round(s * 10 ** n)))
_au32 = {(r['L'], round(r['k'], 3)): r for r in _U31['chains']}
def _cross32(L):
    q = sorted((r['k'], r['R'], r['R_se']) for (l, k), r in _au32.items() if l == L)
    (k1, r1, s1), (k2, r2, s2) = next((a, b) for a, b in zip(q, q[1:]) if a[1] < 0 <= b[1])
    return k1 + (k2 - k1) * (-r1) / (r2 - r1), abs(k2 - k1) / (r2 - r1) ** 2 * math.hypot(r2 * s1, r1 * s2)
_c32 = {L: _cross32(L) for L in (4, 6, 10, 12)}
assert all(abs(_c32[L][0] - _cu[L]) < 1e-9 for L in _c32)
_old_chi = "(at k = 0.46 it is %.2f, %.2f, %.2f, %.2f and %.2f on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴)" % tuple(_au32[(L, 0.46)]['chi'] for L in (4, 6, 8, 10, 12))
_new_chi = "(at k = 0.46 it is %s, %s, %s, %s and %s on 4⁴, 6⁴, 8⁴, 10⁴ and 12⁴)" % tuple(_e(_au32[(L, 0.46)]['chi'], _au32[(L, 0.46)]['chi_se']) for L in (4, 6, 8, 10, 12))
rep(_old_chi, _new_chi, label='Section 8.10, a = 0.02: standard errors of the fluctuation (T32)')
_old_cr = "the curvature crosses zero at k = %.3f on 4⁴, %.3f on 6⁴, about %.3f on 8⁴ (from one coupling) and %.3f and %.3f on 10⁴ and 12⁴," % (_cu[4], _cu[6], _cr8b, _cu[10], _cu[12])
_new_cr = "the curvature crosses zero at k = %s on 4⁴, %s on 6⁴, about %.3f on 8⁴ (from one coupling) and %s and %s on 10⁴ and 12⁴," % (
    _e(*_c32[4], 3), _e(*_c32[6], 3), _cr8b, _e(*_c32[10], 3), _e(*_c32[12], 3))
rep(_old_cr, _new_cr, label='Section 8.10, a = 0.02: standard errors of the crossings (T32)')
# paragraph breaks at topic changes (no word changed)
T32_BREAKS = [
    ("yet doubles the source of gravity.", "Section 8.7 records a proposed specification"),
    ("The target is the rule both routes reach.", "The paper's results are the census and its failures"),
    ("needs about $10^8$.", "Dimension can depend on scale:"),
    ("could read four at large scales.", "The test that counts is the trend with size:"),
    ("while its midpoint moves away (Section 5.1).", "A rule that stays three-dimensional at large scales,"),
    ("Section 8 runs the first tests.", "The hard parts of this route"),
    ("finite sizes it tests.", "The first cap is in hand for part of the known laws."),
    ("times per halving.", "The geometry and the gauge fields are supplied, not grown,"),
    ("What matter does at the center must come from the rule.", "Discrete quantum geometry shows what such dynamics can do."),
    ("instead of reaching a singularity.", "A bounce can survive the quantum matter it holds."),
    ("The fields are free and massless.", "The coupling to curvature squared, which the quantum fields leave free, decides the outcome."),
    ("it does not test four-dimensional SU(3) dynamics.", "The weak force acts on left-handed fields only."),
    ("that Lee and Yang (1956) proposed.", "The anomalies must cancel,"),
    ("has not been performed here.", "Joining general relativity with quantum theory"),
    ("Jurkiewicz and Loll, 2008).", "On the supplied lattice itself,"),
    ("show neither.", "At k = 0.244, Hamber's $k_c$, the fluctuation is"),
    ("from the second halves of the chains.", "Median autocorrelation times of chains near the crossing"),
    ("does not steepen with size.", "Hamber (1993) fitted his fluctuation at k = 0.244"),
    ("are left out here too.", "On 12⁴, 16⁴ and 20⁴, where shorter chains still drift,"),
    ("as a local defect's share of an average should.", "Rough starts carry many such clusters."),
    ("with no sign of a second bulk phase.", "Toward smaller a, mixing worsens:"),
    ("would lead one to expect.", "At a = 0.02, where Hamber (1993) found"),
    ("with no drift toward his $k_c$ beyond 6⁴.", "The term buys an ensemble that mixes;"),
    ("differ in exponent.", "The term also bounds the action below."),
    ("the bound stops short of $k_c$.", "The term has costs."),
    ("the opposite of an expanding universe.", "A rule that rewrites a few relations at a time bounds the links of every update,"),
    ("runs its chains too long.", "The search goes on with longer rules,"),
    ("remain separate requirements.", "The independent Regge audit fails the prespecified mixing screen"),
    ("the continuum limit the gravity sector needs remains open.", "Harder tests sharpen the search:"),
]
for _a, _b in T32_BREAKS:
    rep(_a + ' ' + _b, _a + '\n\n' + _b, label='Paragraph break (T32): ' + _b[:40])
_long32 = [len(l.split()) for l in md.split('\n') if not l.startswith(('|', '#')) and len(l.split()) > 250]
assert not _long32, _long32

open(OUT, 'w', encoding='utf-8').write(md)
json.dump([{'label': l, 'old': o, 'new': n} for l, o, n in CHANGES], open('/home/claude/sa74/out/changes_sa74.json', 'w'), indent=1)
words = sum(len(l.split()) for l in md.split('\n') if not l.startswith('#'))
print('written', OUT, 'words (non-heading lines):', words, 'changes:', len(CHANGES))
