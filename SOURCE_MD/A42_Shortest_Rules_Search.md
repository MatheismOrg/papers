# Searching the Shortest Rewriting Rules for Four-Dimensional Spacetime: First Passes and an Obstacle

## Abstract

If spacetime is the causal order of discrete updates, a theory of everything is an update rule, and the search can run shortest first, since each extra bit halves a rule's prior weight. This paper runs the first passes. It takes all 2,935 distinct rules that rewrite two connected directed relations into three, with at most two new nodes. It adds random samples of 1,500 that rewrite two into four, 1,000 that rewrite three into four and 1,000 that rewrite two three-node relations into three. The record of which update used which update's output is a rule's causal graph. Three tests from causal-set theory read its dimension: the fraction of related pairs in a causal interval, the size of half-intervals, and the longest chain. On random points in flat four-dimensional spacetime they read 3.9, 4.0 and 4.0. A mirror pair of rules grows two-dimensional orders, on which all three tests agree up to 49,000 events. No rule reads four on all three, in runs of up to 144,000 updates. The closest read four on the first two, but their chains run 1.7 to 11 times too long, the opposite of what expanding de Sitter space gives. A bound explains part of the gap: each update has at most as many direct causes as its rule consumes relations, so its links stay few, while a random four-dimensional causal set gives its elements ever more. A rewriting rule can match four-dimensional spacetime only at large scales. Code and data are supplied.

Keywords: causal sets; rewriting systems; Wolfram model; dimension estimators; Myrheim-Meyer dimension; discrete spacetime; Lorentz invariance; algorithmic probability

## 1. Introduction

Two discrete approaches to spacetime point at the same object. Causal-set theory builds spacetime from events and their order: which events can influence which, and how many events there are, give the geometry (Bombelli, Lee, Meyer and Sorkin, 1987; Surya, 2019). Rewriting systems build the world from a rule applied again and again to its own output (Wolfram, 2002, 2020). Each application of the rule is an update, and each update uses relations that earlier updates made. The record of which update used which update's output is a causal graph, and its transitive closure, the relation "can influence through a chain of updates", is a causal order. So every rewriting rule grows a candidate spacetime. The question is which rule grows ours.

Both approaches take change to come in steps, for a reason that holds in any universe: a cause produces its effect as the next state, and a continuum has no next state, since between any two points lies a third. The reason is an a priori proof in the logician's sense. Russell's continuum, in which motion is only being at different places at different times (Russell, 1903), gives the emperor's new points: positions everywhere and nothing that moves them.

The rules can be listed in order of length. Write each rule as a string of bits, and let each extra bit halve its prior weight, as algorithmic probability does (Solomonoff, 1964). Then the short rules carry most of the weight, and a search should start with them. Others have looked at the causal orders of short rules. Wolfram (2020) surveyed the behavior of many short rules, including the growth of their causal graphs. Bolognesi (2010) built causal sets from Turing machines, string rewriting and other simple models and found growth dimensions up to about 3.5. Gorard (2020a, 2020b) connected Wolfram's models to causal sets and to relativity and related their dimension estimates to causal-set ones, with examples in two dimensions. None of these, as far as we can find, reports a short rule whose causal order reads four on the tests causal-set theory uses to recognize a four-dimensional spacetime.

This paper runs those tests on the smallest class of rules with room for more than a single line of causes, all 2,935 of its rules with at most two new nodes, and on samples of three larger classes. It calibrates three tests on random scatterings of points in flat spacetime and classifies 6,435 rules. It runs the most promising to 20,000 updates and ten of them past 100,000, compares them with what the equations give for an expanding universe, and supplies the code and data. The result is negative for four dimensions and positive for two. Several rules grow two-dimensional causal orders, so the tests can find a dimension when a rule has one. No rule grows a four-dimensional order. The rules that come closest read four on two tests and between 1.7 and 3.0 on the third.

The search also turned up an obstacle that holds for every rule that rewrites a fixed number of relations at a time. An update has at most as many direct causes as its rule consumes relations, so the causal graph of such a rule gives each update a bounded number of links, the related updates with nothing between it and them. A random four-dimensional causal set gives its elements ever more links as the region grows. No graph with a bounded number of edges per element can sit on one without picking out a preferred frame (Bombelli, Henson and Sorkin, 2009). A rule's causal graph can therefore match four-dimensional spacetime only at large scales. Section 5 states the bound and what it asks of the next pass.

## 2. Rules and their causal graphs

### 2.1 Rules

A state is a collection of relations, each an ordered pair of nodes; it is a directed graph that may contain loops and repeated edges. A rule has a left side, a pattern of relations written with variables, and a right side, relations written with the left side's variables and possibly new ones. An update finds relations that match the left side, removes them, and adds the right side, creating a fresh node for each new variable. Different variables may match the same node. For example, the rule

{{x,y},{x,z}} → {{x,y},{x,w},{z,x}}

takes two relations leaving the same node x. It keeps the first, adds a relation from x to a new node w, and reverses the second.

The search starts with rules that consume two relations and make three, for two reasons. A rule that consumes one relation gives each update one direct cause, the update that made that relation, so its causal order is a tree. A rule that consumes two and makes two keeps the number of relations fixed, so from a small start its updates can only follow one another in a line. Two into three is the smallest class with room for more.

The class searched in full has left sides of two relations, connected through a shared node, and right sides of three relations, connected, built from the left side's variables and at most two new nodes. Rules that differ only by renaming variables count once, and so do rules whose relations are listed in a different order. The class then holds 2,935 distinct rules. Reversing every relation turns a rule into its mirror image, which usually grows the same kind of causal graph; the count keeps both. A second set holds 1,500 distinct rules that rewrite two relations into four, drawn at random. The draw picks right sides at random before renamings are merged, so it is not uniform over distinct rules.

A second pass samples two larger classes, 1,000 distinct rules each, drawn the same way: rules that rewrite three binary relations into four, and rules that rewrite two relations of three nodes into three. A relation of three nodes, such as {x,y,z}, links three nodes at once, and a structure built from such links is a hypergraph, the kind many of Wolfram's examples use (Wolfram, 2020). The rule stays fixed; only the hypergraph it rewrites changes. These runs start from as many self-loops as the left side has relations. They use an engine for relations of any size and left sides of any length, which reproduces the first engine's causal graphs exactly on 277 rules of the first pass.

A rule's length in bits is the number of symbols it uses, one for each node of each relation, times the bits each symbol needs, the base-two logarithm of the number of distinct variables. The full class runs from 0 to 23.2 bits, the two-into-four sample to 27.9, and the second pass's samples to 36.2 and 42.1. By this count the classes overlap: 903 rules in the sample, at 19 bits or fewer, are shorter than most rules in the full class. A strict shortest-first search would interleave the classes; this pass takes whole classes, the smallest first.

### 2.2 Updating and the causal graph

The update order is fixed. Updates happen in generations. Each generation lists the relations present at its start, oldest first, and matches each one that is still free with its oldest compatible partner. The matches found in a generation share no relation, and all of them are applied at once. Every run of the first pass starts from two self-loops on one node, a minimal start of the kind Wolfram (2020) uses. It ends with the generation that passes 3,000 updates, or earlier if no match is left, 400 generations pass, 60,000 relations exist or 8 seconds elapse. The first pass's long runs go to 20,000 updates or 60 seconds; the second pass's go to about 120,000 updates or 600 seconds, without the limits on generations and relations.

An update's direct causes are the updates that made the relations it consumed; the starting state counts as a root update. One update precedes another when a chain of direct causes runs from the first to the second. This order is a causal set in the technical sense, a partial order in which only finitely many elements lie between any two (Bombelli et al., 1987).

## 3. Reading a dimension

A causal interval is a pair of related elements together with every element between them. In flat spacetime of d dimensions, counting time, scatter points at random with uniform density, a sprinkling, and keep the order light cones give them. The intervals of a sprinkling have known statistics, and three of them give tests of dimension.

The first test is the ordering fraction (Myrheim, 1978; Meyer, 1988). In an interval of N elements, let f be the fraction of the N(N - 1)/2 pairs that are related. For a sprinkled interval in d dimensions,

$$f(d)=\frac{\Gamma(d+1)\,\Gamma(d/2)}{2\,\Gamma(3d/2)},$$

where Γ is the gamma function, the factorial extended beyond whole numbers. So f is 1/2 in two dimensions and 1/10 in four, and the measured f gives d.

The second test is midpoint scaling (Reid, 2003). Among the elements z between the two ends x and y, take the one that makes the smaller of the intervals [x, z] and [z, y] as large as possible. In d dimensions that smaller interval holds about $N/2^d$ elements.

The third test is the longest chain, the largest set of elements each preceding the next. In a sprinkled interval it grows as $m_d N^{1/d}$, where $m_d$ is a constant (Brightwell and Gregory, 1991). Three sprinklings of 400 and of 800 points for each d give $m_d$ = 1.95, 2.06, 2.14 and 2.25 for d = 2, 3, 4 and 5. A chain of length L in an interval of N elements then gives d.

A fourth reading serves as a check that needs no fitted constant: keep each element with probability 1/2, 1/4 and 1/8, four times each, keeping both ends, and see how fast the longest chain shrinks. In d dimensions it shrinks as the number kept to the power 1/d. This check is noisier, and on sprinklings of 400 and 800 points in four dimensions it read 4.45 and 4.15.

Table 1 shows the three tests on sprinklings. They agree within about 0.2 up to four dimensions and within 0.4 at five; at 800 points in four dimensions they read 3.88, 4.00 and 4.02.

**Table 1.** The three tests on sprinkled intervals (means of three sprinklings) and on two grids.

| Order | Elements | Ordering fraction | Midpoint | Longest chain |
|---|---:|---:|---:|---:|
| Sprinkled, d = 2 | 400 / 800 | 1.96 / 1.99 | 1.95 / 2.00 | 1.99 / 2.01 |
| Sprinkled, d = 3 | 400 / 800 | 2.93 / 2.98 | 3.05 / 2.98 | 3.02 / 2.98 |
| Sprinkled, d = 4 | 400 / 800 | 3.87 / 3.88 | 4.06 / 4.00 | 3.98 / 4.02 |
| Sprinkled, d = 5 | 400 / 800 | 4.66 / 4.83 | 4.92 / 5.19 | 4.86 / 5.10 |
| Light-cone lattice, 1+1 | 900 | 1.92 | 1.91 | 2.00 |
| Cellular automaton, 3+1 | 3,781 | 3.69 | 3.63 | 5.18 |

The two grids show what the tests catch. The light-cone lattice in one dimension of space and one of time, the causal order of a one-dimensional cellular automaton, reads two on all three tests. The second grid is a three-dimensional cellular automaton, a cube of cells each updated from the 27 cells around it, itself included, over ten steps. Its signals spread in square-sided pyramids, and its interval reads 3.69 and 3.63 on the first two tests and 5.18 on the third: its chains are too short for its size. Agreement of all three tests at four is therefore a demanding standard, which a grid with a preferred frame can fail.

The tests are applied to intervals of each rule's causal graph. Take the six updates with the largest pasts. The largest intervals run from the root to each of them. The inner intervals run to each of them from an earlier update at 30 percent of its depth, depth being the longest chain of causes from the root, chosen to make the interval largest. Each must hold at least 120 elements, and a rule's reading on each test is the median over its intervals. For the long runs a size series adds single intervals, of at least 60 elements, from one top update down to bottom ends at 0, 10, 20 and up to 70 percent of its depth. A rule counts as four-dimensional on an interval if all three tests read between 3.5 and 4.5 there, and as three-dimensional if all read between 2.5 and 3.5.

## 4. Results

### 4.1 A census of the shortest rules

Table 2 sorts the rules by what they do in 3,000 updates. A third of the two-into-three rules stop or stay below 300 updates: the median one makes a single update and then finds no match. More than a quarter grow a single chain, one update per generation, which is one-dimensional. Another 28 percent grow exponentially, about 10 percent or more updates in each generation than in the one before, over the second half of the run; their median rate multiplies the updates by about 1.47 per generation. These rules spread so fast that pasts stay thin: the largest past in a run holds a median 1 percent of its updates, and only 42 of 825 make an interval of 120 elements. The remaining 374 grow polynomially. If the updates per generation grow as the generation number to a power p, the total grows like a volume of about p + 1 dimensions. Of these rules, 291 have p near 0, 75 near 1, 4 near 2 and 4 near 3.

**Table 2.** What the rules do in 3,000 updates.

| Behavior | Two into three, all | Two into four, sample |
|---|---:|---:|
| Rules | 2,935 | 1,500 |
| Stop or stay small | 928 (31.6%) | 279 (18.6%) |
| Single chain or tree | 808 (27.5%) | 760 (50.7%) |
| Exponential growth | 825 (28.1%) | 357 (23.8%) |
| Polynomial growth | 374 (12.7%) | 104 (6.9%) |
| Volume growth near 1, 2, 3, 4 dimensions | 291, 75, 4, 4 | 92, 9, 3, 0 |
| Polynomial, reading two on the first two tests in the largest intervals | 18 | 5 |
| Read four on the first two tests, any interval | 2 | 0 |

The rules whose volume grows like three dimensions or more do not read three or four on the tests. Their pasts stay thin: the eight such rules in the two-into-three class, run to 20,000 updates, read between 1.0 and 2.1 on all three tests in every interval. The tests read the order of the updates, and here that order stays thin.

### 4.2 Two dimensions

Eighteen of the polynomially growing two-into-three rules read near two on the first two tests in their largest intervals, twelve of them in intervals of 500 elements or more, and five rules in the sample do. Six ran to 20,000 updates; Table 3 shows the two cleanest,

{{x,y},{x,z}} → {{x,y},{x,w},{z,x}}

and its mirror image {{x,y},{z,y}} → {{x,y},{y,z},{w,y}}, at 20 bits each. In about 250 generations they make intervals of about 1,000 to 9,400 elements, which read 1.86 to 2.09 on the first two tests and 1.86 to 1.93 on the longest chain; the check reads 1.83 to 2.04. Their updates per generation grow in proportion to the generation, as a two-dimensional volume does. These rules grow two-dimensional causal orders, flat at the sizes reached: spacetimes of one space dimension and time. Two other rules read two, with a bump to between 2.5 and 3.0 on the first two tests at middle sizes, and one rule in the sample reads between 1.8 and 2.5.

So the tests find a dimension when a short rule has one, and two dimensions are within reach of the shortest rules.

**Table 3.** Long runs to 20,000 updates: the range of each test over a rule's size series.

| Rule | Interval sizes | Ordering fraction | Midpoint | Longest chain |
|---|---:|---:|---:|---:|
| {{x,y},{x,z}} → {{x,y},{x,w},{z,x}} | 1,008 to 9,355 | 1.86 to 2.07 | 1.89 to 2.08 | 1.86 to 1.90 |
| {{x,y},{z,y}} → {{x,y},{y,z},{w,y}} | 974 to 9,220 | 1.87 to 2.07 | 1.95 to 2.09 | 1.89 to 1.93 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y}} | 78 to 5,064 | 1.79 to 4.33 | 1.89 to 5.05 | 1.43 to 2.04 |
| {{x,y},{x,z}} → {{y,y},{z,w},{w,y}} | 98 to 918 | 2.87 to 3.97 | 3.15 to 4.67 | 2.07 to 2.54 |
| {{x,y},{y,z}} → {{y,x},{y,w},{z,x}} | 79 to 964 | 2.92 to 3.85 | 3.13 to 4.15 | 2.02 to 2.50 |
| {{x,y},{x,z}} → {{x,y},{z,y},{z,w},{w,z}} | 61 to 616 | 2.80 to 4.10 | 3.12 to 4.79 | 2.28 to 2.83 |
| {{x,y},{z,y}} → {{x,z},{x,w},{y,x},{w,y}} | 86 to 401 | 3.08 to 3.64 | 3.43 to 4.14 | 2.49 to 2.85 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y},{w,y}} | 72 to 3,100 | 1.00 to 4.02 | 1.00 to 4.63 | 1.19 to 1.87 |

### 4.3 Near four

In 3,000 updates, two rules read between 3.5 and 4.5 on both of the first two tests, both in inner intervals, and none in the sample does. The long runs took 23 candidates further. They include every rule that came near four on either of the first two tests, the four exponential rules that read highest on the ordering fraction, the rules whose volume grew fastest, and the closest rules in the sample. In the size series of the long runs, 12 intervals in 7 rules read between 3.5 and 4.5 on both of the first two tests; two of the rules are a mirror pair with identical graphs. On every one of these intervals the longest chain reads between 1.7 and 2.85: the chains run 1.7 to 9.6 times longer than in a four-dimensional interval of the same size. Across all intervals of all 27 long runs, the chain test reads 2.85 at most.

The near misses come in two families. The first adds updates in proportion to the generation, so its volume grows like two dimensions. In {{x,y},{y,z}} → {{x,y},{y,w},{z,y}} and three related rules, the readings swing with size: the ordering fraction gives 1.8 at 80 elements, 4.3 at 1,300 to 3,500, and 2.7 at 5,000. Where it reads 4.3, the midpoint overshoots to between 4.5 and 5.1, and the longest chain stays between 1.4 and 2.0 throughout. A rule of the same kind in the sample, {{x,y},{y,z}} → {{x,y},{y,w},{z,y},{w,y}}, reads 4.0 and 3.85 on the first two tests at 966 elements, with a chain 9.6 times too long. The second family grows exponentially: {{x,y},{x,z}} → {{y,y},{z,w},{w,y}}, its mirror image, {{x,y},{y,z}} → {{y,x},{y,w},{z,x}}, {{x,y},{y,z}} → {{z,x},{z,w},{w,x}} and two rules in the sample. It supplies 11 of the 12 intervals. Away from the start its intervals read up to 4.1 on the ordering fraction and 4.8 on the midpoint, while the chain reads 2.0 to 2.85. The readings peak at middle sizes, and the largest intervals, which reach back to the start, when the graph was small, read lower.

These orders are wide, since most pairs of updates in an interval are unrelated, as in four dimensions, and long, since chains of causes run on for many steps, as in two or three. No sprinkled flat spacetime, of any one dimension, combines the two in these proportions, so these orders are not four-dimensional spacetimes at the sizes reached. No interval in these runs reads three on all three tests either; in the closest the chain test reads 2.49 while the midpoint reads 3.43.

### 4.4 How firm the negative result is

The pattern holds across the choices this pass made. The largest intervals, the inner intervals and the size series agree. Mirror images mostly agree: 56 of the 1,420 mirror pairs in the full class land in different classes, since the fixed update order treats the two members differently. Taking candidates from 3,000 to 20,000 updates moved no rule into the four-dimensional range. On the 12 intervals that read four on the first two tests, the check that needs no fitted constant reads 1.0 to 2.5, with one interval at 3.3, against 4.15 to 4.45 on four-dimensional sprinklings. The chain constants were fitted at 400 and 800 elements and drift slowly with size, but the chain readings on those intervals would have to rise by 0.65 to 1.8 to reach the four-dimensional range.

### 4.5 The second pass

Table 5 shows the two larger classes. They split differently from the first, with fewer single chains and more exponential growth: about a third stop, a tenth grow single chains, half grow exponentially, and the rest grow polynomially, mostly with one-dimensional volume. They end the same way. Neither class has a rule that reads four on the first two tests in any interval. Each has one polynomially growing rule that reads two, among them {{x,y,x},{z,w,x}} → {{x,y,y},{x,u,x},{w,z,x}}, whose chain test reads 1.9 as well. The exponential rules again come closest: their inner intervals reach 3.3 and 3.6 on the first two tests, with the chain at 2.4 to 2.5.

**Table 5.** The second pass's samples in 3,000 updates.

| Behavior | Three into four, sample | Two three-node into three, sample |
|---|---:|---:|
| Rules | 1,000 | 1,000 |
| Stop or stay small | 344 (34.4%) | 348 (34.8%) |
| Single chain or tree | 95 (9.5%) | 104 (10.4%) |
| Exponential growth | 496 (49.6%) | 494 (49.4%) |
| Polynomial growth | 65 (6.5%) | 54 (5.4%) |
| Volume growth near 1, 2, 3, 4 or more dimensions | 60, 3, 1, 1 | 46, 7, 1, 0 |
| Polynomial, reading two on the first two tests in the largest intervals | 1 | 1 |
| Read four on the first two tests, any interval | 0 | 0 |

The samples ran to 3,000 updates. Ten rules from the first pass then ran six times longer, to between 100,000 and 144,000 updates: the two-dimensional pair, a rule from the sample that reads two, and seven near misses (Table 6). The tests were applied inside the past of the top update, so that the stored pasts fit in memory. The two-dimensional rules stay two: intervals of 5,500 to 59,000 elements read 1.7 to 2.4 on all three tests, with no drift toward three or four. The linear near misses overshoot: at middle sizes their ordering fraction reaches 4.7 to 5.4 and their midpoint 5.0 to 5.7, while the longest chain stays at about 2.0 or below. The exponential near misses' pasts stay thin, so their intervals grow only to between 790 and 2,700 elements. They read up to 4.6 and 5.7 on the first two tests and 2.0 to 3.0 on the chain. One small interval, 144 elements in {{x,y},{z,y}} → {{x,z},{x,w},{y,x},{w,y}}, falls inside the three-dimensional window, at 3.24, 3.47 and 2.56; the same rule's larger intervals read 3.8 to 4.1, 4.4 to 4.5 and 2.7 to 3.0. In the exponential near misses all three readings rise from small intervals to middle ones and fall back at the largest. Their chain test peaks between 2.4 and 3.0, and in the linear near misses it peaks at 1.9 to 2.0. No interval in these runs reads four on all three tests.

**Table 6.** Runs past 100,000 updates: the range of each test over a rule's size series.

| Rule | Updates | Interval sizes | Ordering fraction | Midpoint | Longest chain |
|---|---:|---:|---:|---:|---:|
| {{x,y},{x,z}} → {{x,y},{x,w},{z,x}} | 120,269 | 5,483 to 44,991 | 1.71 to 1.91 | 1.71 to 1.91 | 1.84 to 1.87 |
| {{x,y},{z,y}} → {{x,y},{y,z},{w,y}} | 120,001 | 5,581 to 49,220 | 1.76 to 1.97 | 1.78 to 1.98 | 1.89 to 1.92 |
| {{x,y},{z,y}} → {{x,y},{y,w},{u,x},{u,z}} | 120,077 | 5,700 to 58,939 | 1.78 to 2.35 | 1.95 to 2.38 | 1.98 to 2.05 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y}} | 120,067 | 106 to 30,157 | 1.00 to 5.13 | 1.00 to 5.46 | 1.17 to 2.03 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y},{w,y}} | 101,763 | 97 to 14,867 | 1.00 to 4.74 | 0.98 to 5.02 | 1.15 to 1.91 |
| {{x,y},{y,z}} → {{x,y},{y,x},{z,y},{w,y}} | 117,026 | 61 to 12,376 | 1.00 to 5.36 | 0.98 to 5.73 | 1.18 to 2.04 |
| {{x,y},{x,z}} → {{y,y},{z,w},{w,y}} | 123,647 | 194 to 2,725 | 3.37 to 4.64 | 3.51 to 5.66 | 2.19 to 2.76 |
| {{x,y},{y,z}} → {{y,x},{y,w},{z,x}} | 123,261 | 134 to 1,099 | 2.82 to 3.66 | 2.74 to 3.96 | 2.02 to 2.36 |
| {{x,y},{x,z}} → {{x,y},{z,y},{z,w},{w,z}} | 144,021 | 63 to 1,102 | 2.50 to 4.44 | 2.39 to 4.88 | 2.11 to 2.85 |
| {{x,y},{z,y}} → {{x,z},{x,w},{y,x},{w,y}} | 121,899 | 67 to 786 | 2.76 to 4.07 | 2.74 to 4.49 | 2.34 to 2.96 |

The comparison so far is with flat spacetime, which is what the equations give for empty space. A rule that grows exponentially should also be compared with what they give for an expanding universe: de Sitter space, the vacuum with a positive cosmological constant. Table 7 reads sprinkled intervals of four-dimensional de Sitter space across which the universe grows by up to four e-folds, four factors of e, about 2.7 each. As the expansion passes one e-fold, all three readings fall, the first two fastest, so the chain test reads highest. The exponential rules show the opposite order, with the chain test lowest, so de Sitter expansion does not explain them. The flat row is a sprinkling across 0.01 e-fold, separate from those of Table 1.

**Table 7.** Sprinkled intervals of four-dimensional de Sitter space, 800 points (means of three sprinklings).

| Expansion across the interval | Ordering fraction | Midpoint | Longest chain |
|---|---:|---:|---:|
| None (flat) | 3.92 | 4.06 | 3.87 |
| 1 e-fold | 3.78 | 3.98 | 3.87 |
| 2 e-folds | 3.44 | 3.69 | 3.81 |
| 3 e-folds | 3.04 | 3.19 | 3.53 |
| 4 e-folds | 2.69 | 2.71 | 3.39 |

## 5. Why a rule can match four dimensions only at large scales

A link is a pair of related elements with nothing between them. In a rule's causal graph every link is a direct cause: if one update depends on another through a third, the third lies between them. An update consumes as many relations as its left side holds, each made by one earlier update, so it has at most that many direct causes. Each relation it makes is consumed at most once, so it has at most as many direct effects as its right side holds. No update therefore has more than five links in the two-into-three class, six in the two-into-four sample, or seven in the three-into-four sample, however large the graph grows. The long runs confirm the bound: no update had more than two direct causes or more than four direct effects.

A sprinkling behaves differently. In four-dimensional sprinkled intervals the mean number of links per element rises from 11.9 at 200 elements to 18.6 at 400, 31.8 at 800 and 47.6 at 1,600, and the largest number from 41 to 171 (Table 4). In a sprinkling of all of flat spacetime, every element has infinitely many links, reaching out along its light cone (Moore, 1988). Bombelli, Henson and Sorkin (2009) proved more. No graph with finitely many edges at each element can be attached to a sprinkling of flat spacetime in a way that respects Lorentz symmetry, the symmetry of special relativity. Any such graph picks out a preferred frame.

**Table 4.** Links per element: sprinkled four-dimensional intervals and two grids.

| Order | Elements | Mean links | Largest |
|---|---:|---:|---:|
| Sprinkled, d = 4 | 200 | 11.9 | 41 |
| Sprinkled, d = 4 | 400 | 18.6 | 72 |
| Sprinkled, d = 4 | 800 | 31.8 | 134 |
| Sprinkled, d = 4 | 1,600 | 47.6 | 171 |
| Light-cone lattice, 1+1 | 900 | 3.9 | 4 |
| Cellular automaton, 3+1 | 3,781 | 35.0 | 54 |
| Any rule, two into three | any | at most 5 | at most 5 |

So the causal graph of such a rule, update by update, is not a random four-dimensional causal set, and no search over these rules can find one that is. The light-cone lattice shows what a bounded number of links still allows: it reads two on all three tests with at most four links per element. What the bound costs is the frame. A lattice has a rest frame, and by the theorem of Bombelli, Henson and Sorkin, a rule of this kind whose updates sat in flat spacetime would pick one out too, at the scale of single updates.

The target is then a rule whose causal graph matches four-dimensional spacetime at large scales, after coarse-graining. One way to coarse-grain is to keep a random fraction of the updates, which leaves the first two tests nearly unchanged on average. The rule's small-scale frame must stay hidden at the scales and speeds that can be tested. That is a demanding requirement. Quantum corrections tend to carry a small-scale frame into low-energy physics unless something protects the symmetry (Collins, Perez, Sudarsky, Urrutia and Vucetich, 2004), and gamma-ray bursts bound such effects near the Planck scale (Abdo et al., 2009). A candidate rule must pass a frame test as well as the dimension tests.

The bound holds for rules that rewrite a fixed number of relations at a time. A rule escapes it if each update can draw on an unbounded part of its past. Classical sequential growth works this way: each new element of a causal set takes a whole set of earlier elements as its past, with probabilities fixed by a sequence of constants, in the simplest case one (Rideout and Sorkin, 2000). Such a rule can still be short. It gives up the locality of rewriting, in which each update reads only a few relations.

## 6. Limits of this pass

The classes are small. Most relations join two nodes, the left sides hold two or three relations, and the larger classes were only sampled. Wolfram's own examples often use relations of three nodes and larger rules (Wolfram, 2020).

One update order and one start were used. A rule whose causal graph does not depend on the update order is causally invariant (Wolfram, 2020; Gorard, 2020b). This pass did not test for it, and a different order can change the causal graph of a rule that lacks it.

The sizes are small. The longest runs reach 144,000 updates and intervals of about 59,000 elements, and the samples of the second pass ran only to 3,000 updates. In four dimensions the longest chain grows only as the fourth root of the number of elements, so an interval of 59,000 elements is about 33 steps tall, and one 200 steps tall needs about $10^8$. Dimension can depend on scale: several approaches to quantum gravity find the effective dimension, as a random walk measures it, falling toward two at the shortest distances (Ambjørn, Jurkiewicz and Loll, 2005; Carlip, 2017). If a rule's updates are events at the Planck scale, the intervals measured here are a few hundred Planck lengths across or less, and a rule that reads two or three there could read four at large scales. The test that counts is the trend with size: a rule whose readings climb toward four as its intervals grow is a candidate. No rule run long shows such a climb on all three tests; the chain readings of the near misses peak between 1.9 and 3.0 and fall back. A rule that stays three-dimensional at large scales, with two dimensions of space and one of time, would not describe our world. Gravity in three dimensions has no waves, and masses in it do not attract (Deser, Jackiw and 't Hooft, 1984).

Dimension is the first test only. A candidate must also show no preferred frame at large scales, curvature that obeys Einstein's equations, and matter. It needs a quantum version too: a unitary rule, one that keeps total probability fixed. Quantum cellular automata of this kind recover the Dirac equation of the electron (D'Ariano and Perinotti, 2014).

Slow rules can hit the 8-second limit, and then stop at slightly different points on another machine. In the full class 29 runs hit it, and in the sample 66. Twelve of them were single chains cut off below 300 updates; rerun with 30 seconds, all twelve reached 400 updates as single chains, and Table 2 counts them that way. Rerunning a random 60 rules of the full class reproduced 59 exactly; the other was one of the twelve. In the second pass 23 and 24 runs hit the limit. Eighteen of them stopped below 300 updates; rerun with 30 seconds, 9 passed 300 and 9 stopped again, each making at most one update per generation.

## 7. The next pass

The next pass should go up in length, shortest first: three relations into four and relations of three nodes in full, then larger rules. It should run several update orders and test causal invariance. It should follow each rule's readings as intervals grow, to $10^6$ or $10^8$ updates, which needs a compiled engine and tests that sample pairs of updates instead of storing every update's past. A rule whose readings climb toward four should count as a candidate. It should add a test for a preferred frame at large scales, and rules of the sequential-growth kind, whose updates draw on unbounded parts of their past. Each rule takes seconds, so the search divides easily among many computers. The code and data are in the supplements.

Enumeration has a reach. Each extra bit doubles the number of rules, so a blind search covers rules of tens of bits. A complete symbolic statement of the Standard Model with general relativity, its gauge groups, matter fields, action and 28 declared inputs, written out and compressed with a standard compressor, comes to about 11,300 bits, and a shorter encoding may exist. General relativity is about 600 bits of it, 5 to 6 percent: one term, (R - 2Λ)/(16πG), two constants, and the metric that every other term uses. The two theories meet only through the metric, and in a causal set the metric is the order of the elements plus their count (Surya, 2019). A rule whose causal graph is four-dimensional at large scales would supply that metric; gravity's own dynamics, the one term and two constants, would still have to come out of it. Most of the bits are not gravity: they are the gauge groups, the matter fields and their charges, and the inputs. If those bits must be written into the rule, no enumeration will reach it.

The search can then be guided from the other end: write the known laws as a discrete rule, remove their infinities, and shorten the result. Discrete geometry already removes one kind. Loop quantum cosmology replaces the big bang with a bounce at about 0.41 of the Planck density, and general relativity holds once the density falls below about one percent of the Planck density (Ashtekar, Pawlowski and Singh, 2006). Rovelli and Vidotto (2014) apply the same correction inside a collapsing star, which then bounces at the Planck density instead of reaching a singularity. The hard parts of this route are the ones this paper meets from the first end. The discrete rule must keep Lorentz symmetry at large scales, and it must carry the Standard Model's left-handed particles, which resist being placed on a grid (Nielsen and Ninomiya, 1981). The target is a rule reached from both ends at once: short, and recovering the known laws.

## 8. Conclusion

Among the shortest rewriting rules, under one update order, none grows a four-dimensional causal order. Of the 2,935 rules that rewrite two relations into three and samples of 3,500 from three larger classes, none reads four on all three of causal-set theory's dimension tests, in runs of up to 144,000 updates. Two-dimensional orders do appear, so the tests can find a dimension when a short rule has one. The near misses read four on two tests and 1.7 to 3.0 on the third: wide like four dimensions, with chains far longer than four dimensions allow, the opposite of an expanding universe. A rule that rewrites a few relations at a time bounds the links of every update, while a random four-dimensional causal set gives its elements ever more links. Such a rule can therefore match our spacetime only at large scales, where its small-scale frame must stay hidden. The search goes on with longer rules, following the readings as the intervals grow, and from the other end, with the known laws made discrete.

## References

Abdo, A. A., et al. (2009). A limit on the variation of the speed of light arising from quantum gravity effects. *Nature*, 462, 331-334. https://arxiv.org/abs/0908.1832

Ambjørn, J., Jurkiewicz, J., & Loll, R. (2005). Spectral dimension of the universe is scale dependent. *Physical Review Letters*, 95, 171301. https://arxiv.org/abs/hep-th/0505113

Ashtekar, A., Pawlowski, T., & Singh, P. (2006). Quantum nature of the big bang. *Physical Review Letters*, 96, 141301. https://arxiv.org/abs/gr-qc/0602086

Bolognesi, T. (2010). Causal sets from simple models of computation. *International Journal of Unconventional Computing*, 6(6), 489-524. https://arxiv.org/abs/1004.3128

Bombelli, L., Henson, J., & Sorkin, R. D. (2009). Discreteness without symmetry breaking: A theorem. *Modern Physics Letters A*, 24, 2579-2587. https://arxiv.org/abs/gr-qc/0605006

Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Space-time as a causal set. *Physical Review Letters*, 59, 521-524.

Brightwell, G., & Gregory, R. (1991). Structure of random discrete spacetime. *Physical Review Letters*, 66, 260-263.

Carlip, S. (2017). Dimension and dimensional reduction in quantum gravity. *Classical and Quantum Gravity*, 34, 193001. https://arxiv.org/abs/1705.05417

Collins, J., Perez, A., Sudarsky, D., Urrutia, L., & Vucetich, H. (2004). Lorentz invariance and quantum gravity: An additional fine-tuning problem? *Physical Review Letters*, 93, 191301. https://arxiv.org/abs/gr-qc/0403053

D'Ariano, G. M., & Perinotti, P. (2014). Derivation of the Dirac equation from principles of information processing. *Physical Review A*, 90, 062106. https://arxiv.org/abs/1306.1934

Deser, S., Jackiw, R., & 't Hooft, G. (1984). Three-dimensional Einstein gravity: Dynamics of flat space. *Annals of Physics*, 152, 220-235.

Gorard, J. (2020a). Algorithmic causal sets and the Wolfram model. arXiv:2011.12174. https://arxiv.org/abs/2011.12174

Gorard, J. (2020b). Some relativistic and gravitational properties of the Wolfram model. *Complex Systems*, 29(2), 599-654.

Meyer, D. A. (1988). *The dimension of causal sets* (Doctoral dissertation). Massachusetts Institute of Technology.

Moore, C. (1988). Comment on "Space-time as a causal set". *Physical Review Letters*, 60, 655.

Myrheim, J. (1978). *Statistical geometry* (CERN preprint TH-2538). CERN.

Nielsen, H. B., & Ninomiya, M. (1981). Absence of neutrinos on a lattice (I): Proof by homotopy theory. *Nuclear Physics B*, 185, 20-40.

Reid, D. D. (2003). Manifold dimension of a causal set: Tests in conformally flat spacetimes. *Physical Review D*, 67, 024034. https://arxiv.org/abs/gr-qc/0207103

Rideout, D. P., & Sorkin, R. D. (2000). Classical sequential growth dynamics for causal sets. *Physical Review D*, 61, 024002. https://arxiv.org/abs/gr-qc/9904062

Rovelli, C., & Vidotto, F. (2014). Planck stars. *International Journal of Modern Physics D*, 23, 1442026. https://arxiv.org/abs/1401.6562

Russell, B. (1903). *The Principles of Mathematics*. Cambridge University Press.

Solomonoff, R. J. (1964). A formal theory of inductive inference, Part I. *Information and Control*, 7(1), 1-22.

Surya, S. (2019). The causal set approach to quantum gravity. *Living Reviews in Relativity*, 22, 5. https://arxiv.org/abs/1903.11544

Wolfram, S. (2002). *A new kind of science*. Wolfram Media.

Wolfram, S. (2020). A class of models with the potential to represent fundamental physics. *Complex Systems*, 29(2), 107-536. https://arxiv.org/abs/2004.08210
