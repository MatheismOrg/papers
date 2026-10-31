# Searching for a Short Discrete Theory of Everything: Rewriting Rules and Supplied-Lattice Tests

## Abstract

If spacetime is the causal order of discrete updates, a theory of everything is an update rule, and the search can run shortest first, since each extra bit halves a rule's prior weight. This paper searches all 2,935 rules that rewrite two connected relations into three, and samples of 3,500 longer ones, reading each causal graph with three dimension tests from causal-set theory, fixed in advance with noise bands. Run oldest relations first, a 20-bit mirror pair grows two-dimensional orders; under random orders it reads near one, so a rule's update order belongs in its length. No rule reads four on all three tests. The closest read four on two, with chain-based dimension estimates 4.7 to 21 standard deviations below the four-dimensional reference. Part of the gap is a bound: an update has at most as many direct causes as its rule consumes relations. Sequential growth escapes the bound yet fails: transitive percolation reaches four on the chain at a million elements while its midpoint drifts to 5.8. From the top, one local unitary step reproduces the Dirac equation in curved spacetime, horizon crossing and pair creation in 2,488 bits. A proposed specification of the known laws on a supplied lattice, with declared dark matter and initial data, compresses to 16,272 bits. All four groups in a new 24-chain Regge audit fail the mixing screen; a complete quantum-gravity construction remains open. Code and data are supplied.

Keywords: causal sets; rewriting systems; Wolfram model; dimension estimators; Myrheim-Meyer dimension; discrete spacetime; Lorentz invariance; algorithmic probability; quantum walks; black holes; particle production; Schwinger model; semiclassical gravity; causal invariance; transitive percolation; sequential growth; preregistration; Ginsparg-Wilson relation; Regge calculus

## 1. Introduction

Two discrete approaches to spacetime point at the same object. Causal-set theory builds spacetime from events and their order: which events can influence which, and how many events there are, give the geometry (Bombelli, Lee, Meyer and Sorkin, 1987; Surya, 2019). Rewriting systems build the world from a rule applied again and again to its own output (Wolfram, 2002, 2020). Each application of the rule is an update, and each update uses relations that earlier updates made. The record of which update used which update's output is a causal graph, and its transitive closure, the relation "can influence through a chain of updates", is a causal order. So every rewriting rule grows a candidate spacetime. The question is which rule grows ours.

Both approaches take change to come in steps. This paper adopts that as its working hypothesis, motivated by an argument rather than a measurement: a cause produces its effect as the next state, and a continuum has no next state, since between any two points lies a third. The argument rests on the premise that change is produced, each state making the next; a law that only lists states at every time describes change without producing it. Russell's continuum, in which motion is only being at different places at different times (Russell, 1903), gives the emperor's new points: positions everywhere and nothing that moves them.

The rules can be listed in order of length. Write each rule as a string of bits, and let each extra bit halve its prior weight, as algorithmic probability does (Solomonoff, 1964). This weights individual short programs more strongly and motivates searching them first; it does not establish that the physical rule must be short. The discreteness argument supplies the premise of the search; the numerical tests below do not independently prove that premise. Others have looked at the causal orders of short rules. Wolfram (2020) surveyed the behavior of many short rules, including the growth of their causal graphs. Bolognesi (2010) built causal sets from Turing machines, string rewriting and other simple models and found growth dimensions up to about 3.5. Gorard (2020a, 2020b) connected Wolfram's models to causal sets and to relativity and related their dimension estimates to causal-set ones, with examples in two dimensions. None of these, as far as we can find, reports a short rule whose causal order reads four on the tests causal-set theory uses to recognize a four-dimensional spacetime.

This paper runs those tests on the smallest class of rules with room for more than a single line of causes, all 2,935 of its rules with at most two new nodes, and on samples of three larger classes. It calibrates three tests on random scatterings of points in flat spacetime and classifies 6,435 rules. It runs the most promising to 20,000 updates and ten of them past 100,000, compares them with what the equations give for an expanding universe, and supplies the code and data. It fixes a four-dimension criterion in advance, with noise bands from random point sets of matched size, tests every rule of the smallest class for causal invariance, reruns the long runs under four update orders, and reads sequential growth with the same tests. The result is negative for four dimensions and positive for two. Several rules grow two-dimensional causal orders, so the tests can find a dimension when a rule has one. No rule grows a four-dimensional order. The rules that come closest read four on two tests and between 1.7 and 3.0 on the third. The two-dimensional orders hold under the first passes' update order only, and sequential growth, the way around the obstacle below, does not read four on all three tests anywhere it was run.

The search also turned up an obstacle that holds for every rule that rewrites a fixed number of relations at a time. An update has at most as many direct causes as its rule consumes relations, so the causal graph of such a rule gives each update a bounded number of links, the related updates with nothing between it and them. A random four-dimensional causal set gives its elements ever more links as the region grows. No graph with a bounded number of edges per element can sit on one without picking out a preferred frame (Bombelli, Henson and Sorkin, 2009). A rule's causal graph can therefore match four-dimensional spacetime only at large scales. Section 5 states the bound and what it asks of the next pass.

The search therefore runs from both ends. From the bottom it lists rules shortest first, and enumeration reaches rules of a few tens of bits, while a complete symbolic statement of the Standard Model with general relativity comes to about 11,300 bits compressed (Section 7). From the top it takes the known laws, which treat their fields as continua, writes them as local unitary steps, checks the steps against the continuum where the spacing is small against every wavelength and curvature radius, and then shortens them while the match holds. Section 8 runs that route: one local update reproduces the Dirac equation in curved spacetime, gauge fields with Gauss's law exact, pair creation and transport across a black-hole horizon, and it exposes a lattice that passes every wave test yet doubles the source of gravity. Section 8.7 records a proposed specification of the known laws on a supplied four-dimensional lattice, with declared dark matter and initial data, at 16,272 bits. This description length is distinct from the length of a verified executable completion. A rule found from the top that recovers the known laws bounds from above the length the bottom search must reach, and each shortening lowers the bound. The target is the rule both routes reach. The paper's results are the census and its failures (Sections 3 to 6), the local-step benchmarks (Sections 8.1 to 8.6) and the Regge audit (Section 8.9); Section 8.7 is a proposal, a specification and its bit count, whose dynamics are untested as a whole.

## 2. Rules and their causal graphs

### 2.1 Rules

A state is a collection of relations, each an ordered pair of nodes; it is a directed graph that may contain loops and repeated edges. A rule has a left side, a pattern of relations written with variables, and a right side, relations written with the left side's variables and possibly new ones. An update finds relations that match the left side, removes them, and adds the right side, creating a fresh node for each new variable. Different variables may match the same node. For example, the rule

{{x,y},{x,z}} → {{x,y},{x,w},{z,x}}

takes two relations leaving the same node x. It keeps the first, adds a relation from x to a new node w, and reverses the second.

The search starts with rules that consume two relations and make three, for two reasons. A rule that consumes one relation gives each update one direct cause, the update that made that relation, so its causal order is a tree. A rule that consumes two and makes two keeps the number of relations fixed, so from a small start its updates can only follow one another in a line. Two into three is the smallest class with room for more.

The class searched in full has left sides of two relations, connected through a shared node, and right sides of three relations, connected, built from the left side's variables and at most two new nodes. Rules that differ only by renaming variables count once, and so do rules whose relations are listed in a different order. Under the fixed update order of Section 2.2 the listing is not neutral, since new relations are aged in the order the right side lists them: listed as {{x,y},{x,z}} → {{z,x},{x,w},{x,y}}, the two-dimensional rule of Section 4.2 grows a single chain. The census runs one listing of each rule, its canonical one. The class then holds 2,935 distinct rules. Reversing every relation turns a rule into its mirror image, which usually grows the same kind of causal graph; the count keeps both. A second set holds 1,500 distinct rules that rewrite two relations into four, drawn at random. The draw picks right sides at random before renamings are merged, so it is not uniform over distinct rules.

A second pass samples two larger classes, 1,000 distinct rules each, drawn the same way: rules that rewrite three binary relations into four, and rules that rewrite two relations of three nodes into three. A relation of three nodes, such as {x,y,z}, links three nodes at once, and a structure built from such links is a hypergraph, the kind many of Wolfram's examples use (Wolfram, 2020). The rule stays fixed; only the hypergraph it rewrites changes. These runs start from as many self-loops as the left side has relations. They use an engine for relations of any size and left sides of any length, which reproduces the first engine's causal graphs exactly on 277 rules of the first pass.

The rule-length score is the number of node symbols times the base-two logarithm of the number of distinct variables. It is a structural score, not a self-delimiting executable encoding; the scheduler, decoder and initial-state costs are additional. It is therefore not directly comparable to the compressed byte counts in Sections 7 and 8. The full class runs from 0 to 23.2 bits, the two-into-four sample to 27.9, and the second pass's samples to 36.2 and 42.1. By this count the classes overlap: 903 rules in the sample, at 19 bits or fewer, are shorter than most rules in the full class. A strict shortest-first search would interleave the classes; this pass takes whole classes, the smallest first.

### 2.2 Updating and the causal graph

The update order is fixed. Updates happen in generations. Each generation lists the relations present at its start, oldest first, and matches each one that is still free with its oldest compatible partner. The matches found in a generation share no relation, and all of them are applied at once. Section 4.6 varies this order. Every run of the first pass starts from two self-loops on one node, a minimal start of the kind Wolfram (2020) uses. It ends with the generation that passes 3,000 updates, or earlier if no match is left, 400 generations pass, 60,000 relations exist or 8 seconds elapse. The first pass's long runs go to 20,000 updates or 60 seconds; the second pass's go to about 120,000 updates or 600 seconds, without the limits on generations and relations.

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

A fixed window ignores how much each test scatters at a given size, and the scatter shrinks as intervals grow. Table 2 measures it on sprinkled intervals of 120 to 4,000 points, with 40 sprinklings at the smaller sizes and 6 at the largest. At 1,000 points a four-dimensional interval reads 3.96 ± 0.06 on the ordering fraction, 4.12 ± 0.13 on the midpoint and 4.02 ± 0.12 on the longest chain. From Section 4.4 on, an interval reads d-dimensional when each of its three readings lies within three standard deviations of the d-dimensional mean at its size. A rule is a four-dimensional candidate when its largest interval of at least 1,000 elements reads four-dimensional and so do the two largest intervals of its size series, so that its readings do not drift away from four as it grows. The size floor matters, since the bands are wide at 120 points.

This criterion, and the runs of Sections 4.6 and 5.1, were written down and fingerprinted with SHA-256 before they ran (04eca2085b16830e…; the full value and the protocol are in the supplements), so a near miss cannot be reread after the fact. Runs added after the protocol are marked exploratory.

**Table 2.** The scatter of the three tests on sprinkled flat intervals: mean ± standard deviation.

| Dimension | Points | Sprinklings | Ordering fraction | Midpoint | Longest chain |
|---:|---:|---:|---:|---:|---:|
| 4 | 120 | 40 | 3.66 ± 0.16 | 3.91 ± 0.24 | 3.98 ± 0.41 |
| 4 | 250 | 40 | 3.84 ± 0.14 | 4.09 ± 0.21 | 3.99 ± 0.23 |
| 4 | 500 | 40 | 3.92 ± 0.08 | 4.13 ± 0.16 | 4.02 ± 0.22 |
| 4 | 1,000 | 20 | 3.96 ± 0.06 | 4.12 ± 0.13 | 4.02 ± 0.12 |
| 4 | 2,000 | 12 | 3.99 ± 0.06 | 4.14 ± 0.12 | 4.08 ± 0.14 |
| 4 | 4,000 | 6 | 3.99 ± 0.04 | 4.07 ± 0.05 | 4.03 ± 0.10 |
| 3 | 1,000 | 20 | 2.99 ± 0.07 | 3.04 ± 0.12 | 2.96 ± 0.09 |
| 3 | 4,000 | 6 | 2.99 ± 0.02 | 3.02 ± 0.04 | 2.97 ± 0.06 |
| 2 | 1,000 | 20 | 1.99 ± 0.03 | 1.99 ± 0.06 | 2.02 ± 0.03 |
| 2 | 4,000 | 6 | 1.99 ± 0.01 | 1.99 ± 0.02 | 2.01 ± 0.01 |

## 4. Results

### 4.1 A census of the shortest rules

Table 3 sorts the rules by what they do in 3,000 updates. A third of the two-into-three rules stop or stay below 300 updates: the median one makes a single update and then finds no match. More than a quarter grow a single chain, one update per generation, which is one-dimensional. Another 28 percent grow exponentially, about 10 percent or more updates in each generation than in the one before, over the second half of the run; their median rate multiplies the updates by about 1.47 per generation. These rules spread so fast that pasts stay thin: the largest past in a run holds a median 1 percent of its updates, and only 42 of 825 make an interval of 120 elements. The remaining 374 grow polynomially. If the updates per generation grow as the generation number to a power p, the total grows like a volume of about p + 1 dimensions. Of these rules, 291 have p near 0, 75 near 1, 4 near 2 and 4 near 3.

**Table 3.** What the rules do in 3,000 updates.

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

Eighteen of the polynomially growing two-into-three rules read near two on the first two tests in their largest intervals, twelve of them in intervals of 500 elements or more, and five rules in the sample do. Six ran to 20,000 updates; Table 4 shows the two cleanest,

{{x,y},{x,z}} → {{x,y},{x,w},{z,x}}

and its mirror image {{x,y},{z,y}} → {{x,y},{y,z},{w,y}}, at 20 bits each. In about 250 generations they make intervals of about 1,000 to 9,400 elements, which read 1.86 to 2.09 on the first two tests and 1.86 to 1.93 on the longest chain; the check reads 1.83 to 2.04. Their updates per generation grow in proportion to the generation, as a two-dimensional volume does. These rules grow two-dimensional causal orders, flat at the sizes reached: spacetimes of one space dimension and time. Two other rules read two, with a bump to between 2.5 and 3.0 on the first two tests at middle sizes, and one rule in the sample reads between 1.8 and 2.5.

So the tests find a dimension when a short rule has one, and two dimensions are within reach of the shortest rules.

The noise bands of Table 2 sharpen this. The pair's orders are two-dimensional and are not two-dimensional sprinklings: their chains run 18 to 62 percent longer than a sprinkled interval's of the same size, 4.5 to 12 standard deviations below the band, and at their largest sizes the first two readings fall to 1.71 to 1.98. Their updates have few links (Section 5), so at the scale of single updates the order is closer to a lattice than to a sprinkling. Section 4.6 shows that the two-dimensional readings also depend on the update order.

**Table 4.** Long runs to 20,000 updates: the range of each test over a rule's size series.

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

The noise bands of Table 2 turn the miss into a distance. On the 12 intervals, the chain readings lie 4.7 to 18 standard deviations below the four-dimensional mean at their sizes. The runs past 100,000 updates add 9 intervals that read four on the first two tests; their chains lie 5.1 to 21 standard deviations below and run 1.7 to 11 times longer than four dimensions allow. No interval of the long runs, at any size, reads four-dimensional within the bands.

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

### 4.6 Other update orders, and causal invariance

A rule fixes its spacetime by itself only if its causal graph does not depend on the order of its updates; such a rule is causally invariant (Wolfram, 2020; Gorard, 2020b). The local test is the critical pair. Two updates compete for a relation; apply each in the smallest state that holds both, and ask whether further updates bring the two results to the same state, up to the names of new nodes. Every rule of the full class was tested with up to two further updates on each side. 677 of the 2,935 pass, 23 percent: 33 grow polynomially, 116 exponentially, 235 as single chains or trees, and 293 stop. With three further updates, 99 of the first 1,691 failing rules rechecked also pass, so the count is a floor. Of the 18 rules that read two, only the cleanest pair passes; 15 of the other 16 still fail with three further updates.

Table 8 reruns the 27 rules of the first pass's long runs to 20,000 updates under four orders: the first passes' order, which takes the oldest relations first and pairs each with its oldest partner; its reverse, newest first; and two random orders. No rule reads four on all three tests under any order. The highest chain reading under any order, 3.28 at 1,411 elements, belongs to {{x,y},{y,z}} → {{x,y},{y,x},{z,y},{w,y}} under newest first, still 5.9 standard deviations below four dimensions. Under newest first, 8 of the 23 candidates grow one-dimensional orders, one or two updates per generation, and the exponential near misses lose their readings near four.

The two-dimensional pair passes the critical-pair test, and its readings still change with the order. Under the two deterministic orders it reads 1.67 to 2.09 on all three tests. Under the two random orders the past of its latest update holds 95 to 96 percent of all updates, and its intervals read 1.08 to 1.28 on the first two tests and 1.97 to 2.43 on the chain: nearly one-dimensional on two tests. Listing its right side in another order makes the oldest-first engine grow a single chain (Section 2.1); under the random orders that listing grows the same nearly one-dimensional orders.

So a finite run shows the rule together with its order. Passing the critical-pair test is a local condition. When the rewriting runs without end, local resolution does not guarantee that every run grows the same causal graph, and here the runs do not: each grows only part of the graph, and different orders grow different parts. A rule's update order, and under a fixed order the listing of its right side, therefore belong in its description, and their bits count toward its length. The first passes' order costs few bits, so the pair's two-dimensional orders stand as the output of a short program, the rule run oldest first. A four-dimensional candidate must read four under every order, or carry its order in its length.

**Table 8.** The rules that read two and the near misses, run to 20,000 updates under four orders: the readings of the largest interval of each size series (ordering fraction, midpoint, longest chain), or the run's behavior when it grows a single chain. The second column gives the critical-pair test, with three further updates where two did not resolve.

| Rule | Critical pairs | Oldest first | Newest first | Random 1 | Random 2 |
|---|---|---:|---:|---:|---:|
| {{x,y},{x,z}} → {{x,y},{x,w},{z,x}} | passes | 1.86, 1.98, 1.86 | 1.67, 1.70, 1.84 | 1.09, 1.08, 1.99 | 1.09, 1.08, 1.97 |
| {{x,y},{z,y}} → {{x,y},{y,z},{w,y}} | passes | 1.87, 1.99, 1.89 | 1.96, 1.98, 1.85 | 1.09, 1.09, 1.99 | 1.09, 1.08, 1.99 |
| {{x,y},{x,z}} → {{y,x},{y,w},{z,x}} | fails | 2.04, 1.85, 1.87 | 1.01, 1.00, 1.29 | 2.92, 2.55, 2.27 | 3.18, 2.90, 2.37 |
| {{x,y},{z,y}} → {{y,x},{y,z},{w,x}} | fails | 1.94, 1.96, 1.85 | 2.65, 2.89, 2.06 | 2.75, 2.65, 2.29 | 3.15, 3.23, 2.37 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y}} | fails | 2.69, 3.04, 1.98 | 1.21, 1.19, 2.21 | 1.88, 1.93, 2.21 | 1.87, 1.92, 2.20 |
| {{x,y},{y,z}} → {{x,y},{y,w},{z,y},{w,y}} | fails | 3.00, 3.42, 1.87 | 1.53, 1.37, 2.18 | 3.16, 3.61, 2.63 | 3.13, 3.61, 2.52 |
| {{x,y},{y,z}} → {{x,y},{y,x},{z,y},{w,y}} | passes | 3.30, 3.94, 2.04 | 3.48, 3.94, 3.02 | 3.11, 3.54, 2.70 | 3.17, 3.67, 2.74 |
| {{x,y},{x,z}} → {{y,y},{z,w},{w,y}} | fails | 3.40, 4.01, 2.38 | 1.25, 1.23, 1.56 | 1.75, 1.78, 1.73 | 1.55, 1.39, 1.63 |
| {{x,y},{z,y}} → {{x,x},{x,w},{w,z}} | fails | 3.40, 4.01, 2.38 | 1.25, 1.23, 1.56 | 1.60, 1.54, 1.72 | 1.53, 1.37, 1.74 |
| {{x,y},{y,z}} → {{y,x},{y,w},{z,x}} | fails | 3.02, 3.23, 2.27 | single chain | 2.60, 2.68, 2.09 | 2.82, 2.74, 2.20 |
| {{x,y},{y,z}} → {{z,x},{z,w},{w,x}} | fails | 3.14, 3.74, 2.31 | single chain | 2.82, 3.27, 2.15 | 2.70, 3.09, 2.07 |
| {{x,y},{z,y}} → {{x,z},{x,w},{y,x},{w,y}} | fails | 3.31, 3.79, 2.61 | 2.44, 2.46, 2.19 | 3.08, 3.50, 2.55 | 3.23, 3.76, 2.58 |
| {{x,y},{x,z}} → {{x,y},{z,y},{z,w},{w,z}} | fails | 3.74, 4.27, 2.64 | 2.49, 2.73, 2.09 | 3.06, 3.29, 2.42 | 3.10, 3.48, 2.46 |

## 5. Why a rule can match four dimensions only at large scales

A link is a pair of related elements with nothing between them. In a rule's causal graph every link is a direct cause: if one update depends on another through a third, the third lies between them. An update consumes as many relations as its left side holds, each made by one earlier update, so it has at most that many direct causes. Each relation it makes is consumed at most once, so it has at most as many direct effects as its right side holds. No update therefore has more than five links in the two-into-three class, six in the two-into-four sample, or seven in the three-into-four sample, however large the graph grows. The long runs confirm the bound: no update had more than two direct causes or more than four direct effects.

A sprinkling behaves differently. In four-dimensional sprinkled intervals the mean number of links per element rises from 11.9 at 200 elements to 18.6 at 400, 31.8 at 800 and 47.6 at 1,600, and the largest number from 41 to 171 (Table 9). In a sprinkling of all of flat spacetime, every element has infinitely many links, reaching out along its light cone (Moore, 1988). Bombelli, Henson and Sorkin (2009) proved more. No graph with finitely many edges at each element can be attached to a sprinkling of flat spacetime in a way that respects Lorentz symmetry, the symmetry of special relativity. Any such graph picks out a preferred frame.

**Table 9.** Links per element: sprinkled intervals, two grids and the two-dimensional pair.

| Order | Elements | Mean links | Largest |
|---|---:|---:|---:|
| Sprinkled, d = 4 | 200 | 11.9 | 41 |
| Sprinkled, d = 4 | 400 | 18.6 | 72 |
| Sprinkled, d = 4 | 800 | 31.8 | 134 |
| Sprinkled, d = 4 | 1,600 | 47.6 | 171 |
| Sprinkled, d = 4 | 4,000 | 84.0 | 308 |
| Sprinkled, d = 2 | 1,000 | 11.0 | 24 |
| Sprinkled, d = 2 | 4,000 | 13.7 | 26 |
| Light-cone lattice, 1+1 | 900 | 3.9 | 4 |
| Cellular automaton, 3+1 | 3,781 | 35.0 | 54 |
| Any rule, two into three | any | at most 5 | at most 5 |
| Two-dimensional pair, four orders | 20,000 updates | 3.0 to 4.0 | 5 |

So the causal graph of such a rule, update by update, is not a random four-dimensional causal set, and no search over these rules can find one that is. The light-cone lattice shows what a bounded number of links still allows: it reads two on all three tests with at most four links per element. What the bound costs is the frame. A lattice has a rest frame, and by the theorem of Bombelli, Henson and Sorkin, a rule of this kind whose updates sat in flat spacetime would pick one out too, at the scale of single updates.

The two-dimensional pair shows the cost in numbers. Its updates have 3.0 to 4.0 links on average under the four orders, and at most five, while sprinkled two-dimensional intervals of 1,000 to 4,000 points give their elements 11.0 to 13.7 on average and up to 26 (Table 9).

The target is then a rule whose causal graph matches four-dimensional spacetime at large scales, after coarse-graining. One way to coarse-grain is to keep a random fraction of the updates, which leaves the first two tests nearly unchanged on average. The rule's small-scale frame must stay hidden at the scales and speeds that can be tested. That is a demanding requirement. Quantum corrections tend to carry a small-scale frame into low-energy physics unless something protects the symmetry (Collins, Perez, Sudarsky, Urrutia and Vucetich, 2004), and gamma-ray bursts bound such effects near the Planck scale (Abdo et al., 2009). A candidate rule must pass a frame test as well as the dimension tests.

The bound holds for rules that rewrite a fixed number of relations at a time. A rule escapes it if each update can draw on an unbounded part of its past. Classical sequential growth works this way: each new element of a causal set takes a whole set of earlier elements as its past, with probabilities fixed by a sequence of constants, in the simplest case one (Rideout and Sorkin, 2000). Such a rule can still be short. It gives up the locality of rewriting, in which each update reads only a few relations.

### 5.1 Sequential growth under the same tests

Sequential growth escapes the bound, so it gets the same tests. In transitive percolation, the simplest such dynamics, elements are born one at a time, and each earlier element becomes a direct cause of the newborn with probability p; the order is the transitive closure (Rideout and Sorkin, 2000). It has one constant, so it is short, and its elements can have as many links as they need. A root added below all elements serves as the bottom of the largest intervals.

Over the protocol's grid, from 1,000 to 8,000 elements and p from 0.0003 to 0.3, no causal set is a candidate. Dense ones, with pN of 30 or more, are nearly one long chain: their largest intervals read between 1.0 and 2.2. Sparse ones, with pN of 3 or less, fall apart into pieces too small to test, and at pN of 4 to 8 they read 3.4 to 6.5 on the first two tests and 2.2 to 3.4 on the chain. In between, at pN of 8 to 24, the intervals that read near four on the first two tests read 2.3 to 3.1 on the chain: wide and long at once, like the rewriting near misses.

Rideout and Sorkin (2001) found evidence that transitive percolation has a continuum limit when p falls as N grows so that the fraction of related pairs stays fixed, including at its four-dimensional value. Table 10 follows that trajectory as an exploratory run outside the protocol: p is tuned so that the largest interval reads four on the ordering fraction, and the causal set grows from 4,000 to 64,000 elements. The midpoint reads 4.26 to 4.97, above four, and the chain 2.90 to 3.46, below it; the chain lies 6 to 10 standard deviations below four dimensions. The chain reading climbs by about 0.1 each time the causal set doubles, while the midpoint drifts slowly up.

A compiled engine followed the trajectory further, to two million elements, reading the same largest interval: the chain exactly, the ordering fraction from 3,000 sampled ancestor counts, and the midpoint from the 600 elements that bottom-k sketches of their pasts and futures rank highest, each then counted exactly. At 16,000 and 128,000 elements this search finds the same midpoint as checking every element. The chain reading reaches 3.92 at one million elements and 4.09 at two million, through four, while the midpoint drifts on to 5.78 (Table 10). Sprinkled four-dimensional intervals of 8,000 to 128,000 points read 3.98 to 4.02 on the ordering fraction, 4.02 to 4.09 on the midpoint and 3.92 to 4.03 on the chain. So the trajectory reads four on the ordering fraction and on the chain at about a million elements, and misses the midpoint by more than one and a half: no element of its largest intervals has both a large past and a large future, as an element near the middle of a four-dimensional interval does. Transitive percolation is not a four-dimensional spacetime at any size reached, and the lead it offered closes on the midpoint.

Ahmed and Rideout (2010) fitted the chains and largest intervals of originary percolation, the variant in which every element descends from a first one, to four-dimensional de Sitter space for p from 0.0001 to 0.005. On the three tests it comes close to de Sitter space on two and misses on the third. At p = 0.001 and 32,000 elements its largest intervals read 1.91 on the ordering fraction and 1.79 on the midpoint, as a de Sitter interval of the same size does across about eight e-folds of expansion, 1.89 and 1.72, yet their longest chains run 86 steps where de Sitter's run 58 to 60, a chain reading of 2.72 against 3.02 to 3.05 in three sprinklings. Wherever originary percolation reads near four on the first two tests, its chain reads 2.3 to 3.1, while de Sitter space with the same ordering fraction reads above 3.7. Glaser and Surya (2013) reached the same verdict on transitive percolation by another route: it looks manifoldlike on coarse measures and fails a count of small intervals that sprinklings pass.

**Table 10.** Sequential growth: the medians of the largest intervals, means of two seeds for transitive percolation held at four on the ordering fraction up to 64,000 elements and of three seeds from 128,000 on (compiled engine), one seed for originary percolation, sprinkled flat intervals, and de Sitter space for comparison (means of three sprinklings).

| Model | Elements | p | Largest interval | Ordering fraction | Midpoint | Longest chain |
|---|---:|---:|---:|---:|---:|---:|
| Transitive, held at four | 4,000 | 0.00291 | 2,412 | 3.91 | 4.40 | 3.00 |
| Transitive, held at four | 8,000 | 0.00164 | 4,912 | 3.93 | 4.56 | 3.14 |
| Transitive, held at four | 16,000 | 0.00089 | 9,645 | 4.11 | 4.73 | 3.18 |
| Transitive, held at four | 32,000 | 0.00049 | 19,559 | 4.11 | 4.70 | 3.31 |
| Transitive, held at four | 64,000 | 0.00027 | 39,874 | 4.07 | 4.81 | 3.44 |
| Transitive, held at four | 128,000 | 0.000153 | 78,767 | 4.01 | 5.14 | 3.56 |
| Transitive, held at four | 256,000 | 0.000083 | 158,272 | 3.97 | 5.04 | 3.63 |
| Transitive, held at four | 512,000 | 0.000044 | 313,713 | 4.11 | 5.55 | 3.79 |
| Transitive, held at four | 1,024,000 | 0.000024 | 628,911 | 4.04 | 5.52 | 3.92 |
| Transitive, held at four | 2,048,000 | 0.000013 | 1,265,886 | 4.04 | 5.78 | 4.09 |
| Sprinkled, flat, d = 4 | 8,000 to 128,000 | | 8,002 to 128,002 | 3.98 to 4.02 | 4.02 to 4.09 | 3.92 to 4.03 |
| Originary | 32,000 | 0.001 | 26,717 | 1.91 | 1.79 | 2.72 |
| Originary | 32,000 | 0.0005 | 19,971 | 3.78 | 4.48 | 3.05 |
| de Sitter, 8 e-folds | 26,700 | | 26,702 | 1.89 | 1.72 | 3.04 |

## 6. Limits of this pass

The classes are small. Most relations join two nodes, the left sides hold two or three relations, and the larger classes were only sampled. Wolfram's own examples often use relations of three nodes and larger rules (Wolfram, 2020).

One start was used. The first passes used one update order and one listing of each rule; Section 4.6 adds three orders for the long runs and a critical-pair test, a local condition, for the full class. The census itself was not rerun under other orders.

The sizes are small. The longest runs reach 144,000 updates and intervals of about 59,000 elements, and the samples of the second pass ran only to 3,000 updates. In four dimensions the longest chain grows only as the fourth root of the number of elements, so an interval of 59,000 elements is about 33 steps tall, and one 200 steps tall needs about $10^8$. Dimension can depend on scale: several approaches to quantum gravity find the effective dimension, as a random walk measures it, falling toward two at the shortest distances (Ambjørn, Jurkiewicz and Loll, 2005; Carlip, 2017). If a rule's updates are events at the Planck scale, the intervals measured here are a few hundred Planck lengths across or less, and a rule that reads two or three there could read four at large scales. The test that counts is the trend with size: a rule whose readings climb toward four as its intervals grow is a candidate. No rule run long shows such a climb on all three tests; the chain readings of the near misses peak between 1.9 and 3.0 and fall back. Transitive percolation held at four on the first test climbs on the chain, through four at about a million elements, while its midpoint moves away (Section 5.1). A rule that stays three-dimensional at large scales, with two dimensions of space and one of time, would not describe our world. Gravity in three dimensions has no waves, and masses in it do not attract (Deser, Jackiw and 't Hooft, 1984).

Dimension is the first test only. A candidate must also show no preferred frame at large scales, curvature that obeys Einstein's equations, and matter. It needs a quantum version too: a unitary rule, one that keeps total probability fixed. Quantum cellular automata of this kind recover the Dirac equation of the electron (D'Ariano and Perinotti, 2014).

Slow rules can hit the 8-second limit, and then stop at slightly different points on another machine. In the full class 29 runs hit it, and in the sample 66. Twelve of them were single chains cut off below 300 updates; rerun with 30 seconds, all twelve reached 400 updates as single chains, and Table 3 counts them that way. Rerunning a random 60 rules of the full class reproduced 59 exactly; the other was one of the twelve. In the second pass 23 and 24 runs hit the limit. Eighteen of them stopped below 300 updates; rerun with 30 seconds, 9 passed 300 and 9 stopped again, each making at most one update per generation.

## 7. The next pass

The next pass should go up in length, shortest first: three relations into four and relations of three nodes in full, then larger rules. It should keep the critical-pair test as a first filter, read every candidate under several update orders, and count the order a rule needs in its length. It should follow each rule's readings as intervals grow, to $10^6$ or $10^8$ updates, which needs a compiled engine and tests that sample pairs of updates instead of storing every update's past. A rule whose readings climb toward four should count as a candidate. It should add a test for a preferred frame at large scales, and rules of the sequential-growth kind, whose updates draw on unbounded parts of their past. Transitive percolation, the lead this pass followed to two million elements, reads four on two tests and fails the midpoint. A growth law that favors manifold-like orders is the next candidate: weighting each causal set by the causal-set action of Benincasa and Dowker (2010) gives, in two dimensions, a phase whose orders look like a continuum (Surya, 2012). Each rule takes seconds, so the search divides easily among many computers. The code and data are in the supplements.

Enumeration has a reach. Each extra bit doubles the number of rules, so a blind search covers rules of tens of bits. A complete symbolic statement of the Standard Model with general relativity, its gauge groups, matter fields, action and 28 declared inputs, written out and compressed with a standard compressor, comes to about 11,300 bits, and a shorter encoding may exist. General relativity is about 600 bits of it, 5 to 6 percent: one term, (R - 2Λ)/(16πG), two constants, and the metric that every other term uses. The two theories meet only through the metric, and in a causal set the metric is the order of the elements plus their count (Surya, 2019). A rule whose causal graph is four-dimensional at large scales would supply that metric; gravity's own dynamics, the one term and two constants, would still have to come out of it. Most of the bits are not gravity: they are the gauge groups, the matter fields and their charges, and the inputs. If those bits must be written into the rule, no enumeration will reach it.

The search also runs from the other end. Every physical process is a sequence of steps, so general relativity and quantum field theory, which treat their fields as continua, describe a rule at the scales where they were measured. This route applies those continuum fields to rules and steps, checks that the outcome matches them at identifiable scales, where the spacing is small against every wavelength and curvature radius in the problem, and then shortens the rules while keeping the match. Section 8 runs the first tests. The hard parts of this route are the ones this paper meets from the first end. The discrete rule must keep Lorentz symmetry at large scales, and it must carry the Standard Model's left-handed particles, which resist being placed on a grid (Nielsen and Ninomiya, 1981). The target is a rule reached from both ends at once: short, and recovering the known laws. The two routes bound each other. A discrete rule found from the top that recovers the known laws caps the length of the shortest such rule, and each shortening lowers the cap; the enumeration from the bottom excludes candidates only within the classes, initial states, schedules and finite sizes it tests. The first cap is in hand for part of the known laws. Written as one self-contained program in Python, the local step of Section 8, which reproduces the Dirac equation on a curved static 1+1 background, prescribed electromagnetic phases and transport across a black-hole horizon, takes 451 source bytes, or 2,488 bits conditional on the Python interpreter, numerical dependencies and the specified compression decoder, after the raw LZMA2 compression that brings the symbolic statement of the Standard Model with general relativity to about 11,300 bits. It reproduces the steps behind those tests to 10^-15. Run unchanged on a background it was not built or checked for, an expanding 1+1 universe, with only its mass angle changing from step to step, the same program creates fermion pairs as the continuum does: 0.1274877 at 2,048 sites against 0.1274880, with errors falling 3.99 to 4.00 times per halving. The geometry and the gauge fields are supplied, not grown, so this caps a rule for matter on a given spacetime; growing the spacetime is what the search from the bottom must add, and two-dimensional orders already appear at 20 bits. Supplying a four-dimensional lattice as data costs few bits, and the statement's general relativity can bend it, as numerical relativity bends the metric's values on a grid; Section 8.6 runs a homogeneous version driven by quantum matter. What a fixed lattice leaves open is gravity at the scale of its own spacing, which a rule that grows its own spacetime would settle too. Section 8.7 writes this candidate out in full, with the lattice form of each field and the completion the Standard Model needs: 16,272 bits.

## 8. Known laws as local steps

One family of local steps already reproduces a broad range of known physics: the Dirac equation in curved spacetime, a quantum packet falling along a geodesic, electromagnetic phases, Maxwell waves, a quench that creates fermion pairs, and transport across a black-hole horizon. Each error falls about fourfold when the spacing halves, the mark of a second-order approximation. These tests supply the geometry, the fields and the starting state, and ask whether the steps reproduce the continuum answer.

### 8.1 The update

For a static metric $ds^2 = -N(x)^2dt^2 + dx^2$, the target is the Dirac Hamiltonian, with $c = \hbar = 1$:

$$H = -i\sigma_z\,[N\,\partial_x + (\partial_x N)/2] + mN\sigma_x.$$

Neighboring amplitudes are updated by real two-site rotations. With spacing $a$ and time step $\Delta t$, the bond angle is $\theta = \arcsin(\Delta t\,N/2a)$. Even and odd bonds are applied in half, full and half layers, and the parity that goes first alternates from step to step. The two chiral components move with opposite signs. Local mass rotations $\exp(-i\Delta t\,mN\sigma_x/2)$ bracket each transport step, and gauge links carry prescribed electromagnetic phases. The arcsine sets the long-wavelength speed exactly at finite $\Delta t/a$, and alternating the layer order cancels an error that depends on parity. One prescription serves every test, and nothing is fitted to a target outcome.

### 8.2 Results

**Table 11.** Local steps against continuum references.

| Situation | Result | Scope |
|---|---|---|
| Massive Dirac field, lapse $N = 1 + 0.25\cos x$ | State error $9.1\times10^{-4}$ at 1,024 sites, falling about fourfold per halving | Supplied 1+1 metric |
| Falling packet and the geodesic | At mass 64: circuit 1.69047, quantum reference 1.69057, geodesic 1.69257 | Finite-width packets |
| Prescribed electromagnetic field | Gauge covariance to $10^{-15}$; interference probability 0.206222 against 0.206107 | No dynamical photons |
| Fermion pairs from a mass quench | Eight modes, full Fock space of 256 states: state error $3.4\times10^{-5}$ at 2,048 sites; mean pairs 0.646932 against 0.646951 | Free fields, selected modes |
| Maxwell plane wave in a curved metric | Field error $7.1\times10^{-5}$ at 2,048 sites; conserved energy to $3\times10^{-16}$ | Plane-symmetric sector |
| Crossing a horizon | State error $5.4\times10^{-4}$ at 8,192 sites; 3.8 to 4.0 times smaller per halving | Radial, massless, fixed background |

The falling packet separates two effects. The circuit matches its quantum reference to $10^{-4}$, and the quantum packet approaches the geodesic as its mass grows and its width, $0.8/\sqrt{m}$, narrows: at masses 8, 16, 32 and 64 its center ends at 1.678, 1.685, 1.689 and 1.691, against the geodesic's 1.693.

### 8.3 Wrong rules fail

Agreement in one number is weak evidence, and wrong rules show why. Replacing the mass term $mN$ by $m$ keeps the light cone and the norm, but leaves massive-state errors of 0.014 to 0.16 that do not shrink, and the packet stays near its starting point instead of falling. Dropping the gauge links destroys the interference. In the quench, evolving with the old mass reproduces the pair count almost exactly while the full state is off by 1.18. Full amplitudes, gauge covariance and independent references catch what a single observable misses.

### 8.4 Black holes

A horizon is no barrier to local steps. In Painlevé-Gullstrand coordinates,

$$ds^2 = -dt^2 + (dr + \sqrt{r_s/r}\,dt)^2,$$

the radial characteristic speeds $v_\pm = \pm 1 - \sqrt{r_s/r}$ stay finite at the horizon $r = r_s$, and the same bond rule, with speed $v$ in place of $N$, carries a wave across it. With $r_s = 2$ on $0.25 \le r \le 6.25$, three packets were run at six resolutions from 256 to 8,192 sites and compared with independently integrated characteristics. An ingoing packet starting outside ends inside with probability 0.99999999992. Inside, even the outward branch falls inward, because $v_+ < 0$ there; outside, it escapes. Errors fall 3.8 to 4.0 times on the last halving, and the norm holds to $10^{-11}$.

The center needs the rule's own dynamics. Keep the classical geometry and stop at the first site, $r = a/2$, and the curvature there, $K = 12r_s^2/r^6$, grows 64-fold each time the spacing halves. A fixed physical spacing caps this sampled expression away from r = 0; it does not establish a regular center or bound the fields evolved by a discrete rule. What matter does at the center must come from the rule. Discrete quantum geometry shows what such dynamics can do. Loop quantum cosmology replaces the big bang with a bounce at about 0.41 of the Planck density, and general relativity holds once the density falls below about one percent of it (Ashtekar, Pawlowski and Singh, 2006). Rovelli and Vidotto (2014) apply the same correction inside a collapsing star, which then bounces at the Planck density instead of reaching a singularity. A bounce can survive the quantum matter it holds. Couple a published regular gravity (Bueno, Cano, Hennigar and Murcia, 2025) to the four-dimensional vacuum stress of free conformal fields, angular pressures included (Herzog and Huang, 2013), and a closed universe still contracts, bounces and expands. The response stays single-valued while the ratio of the quantum coupling to the gravity scale is below 8/(3√3), about 1.54. The Standard Model's free fields, 995.5 in anomaly units, give 0.88 at a gravity scale of one Planck area, and quantum matter then reaches 44 percent of the density; at 1.6 a contracting solution's acceleration diverges. The fields are free and massless. The coupling to curvature squared, which the quantum fields leave free, decides the outcome. Restored with the stable sign at up to 2 x 10^-3 of the gravity scale, it keeps the bounce, and peak curvature rises from 21 to 270 in those units; at 5 x 10^-3 and above, or with the opposite sign at any size tested, the contraction drives the curvature out of the theory's range before it can bounce. For the initial conditions and coupling values tested, the regular bounce survives only at the tested small positive couplings. The scan does not establish a universal necessity or a sharp threshold.

### 8.5 What a rule can settle beyond the cutoff

A specified lattice can regulate short-distance divergences. Discreteness alone does not guarantee bounded fields, a convergent path integral or a dynamically regular center. A rule can settle more, and black holes are where it can show it.

A black hole's entropy is one quarter of its horizon's area in Planck units. Counting the causal links that cross a horizon gives an entropy proportional to that area (Dou and Sorkin, 2003). A rule whose horizons give exactly one quarter would derive the Bekenstein-Hawking law by counting.

Hawking's derivation traces the outgoing radiation back to waves of unbounded frequency at the horizon, far beyond the Planck scale. A discrete rule has no such waves. Models with altered short-wavelength behavior still radiate a thermal spectrum (Jacobson, 1991; Unruh, 1995; Corley and Jacobson, 1996), so discreteness keeps Hawking's result while removing its weakest step. The horizon test above extends directly: give the update a short-wavelength dispersion and measure the converted modes.

A unitary rule loses no information, so an evaporating black hole must return what fell in, and a finite number of states per region makes that bookkeeping explicit. Beyond black holes, the number of elements in a causal set fluctuates by about its square root. That predicted a cosmological constant of about the observed size, $10^{-120}$ in Planck units, before it was measured, and makes it fluctuate in time (Sorkin, 1997; Ahmed, Dodelson, Greene and Sorkin, 2004).

### 8.6 Harder tests: pair creation, gauge fields and the source of gravity

Matching a known equation is the weakest test a rule can pass. Six harder ones show what a search must demand, and three of them catch rules that pass every easier check.

**Pair creation in an expanding universe.** In a 1+1 universe whose scale factor grows as 1.75 + 1.25 tanh(η/τ), the same gates, written in conformal time or in proper time, create fermion pairs as the continuum does: for τ = 0.15 and momentum 1 the occupation is 0.127487 in one clock and 0.127489 in the other, against 0.127488, with errors falling fourfold per halving and the two clocks converging on each other. Without expansion, or without mass, nothing is created. Fermion production in expanding 1+1 spacetimes is established (Fulgado-Claudio, Sánchez Velázquez and Bermudez, 2023); the test is whether local steps reproduce it, and with the wrong mass coefficient they do not, leaving an error of 2.49 that refinement keeps.

**A gauge field with Gauss's law.** In the lattice Schwinger model, fermions coupled to an electric field in one dimension of space (Lamm, Lawrence and Yamauchi, 2019; Shaw et al., 2020; Martinez et al., 2016), every gate preserves Gauss's law exactly across all 15,552 states, and each halving of the time step cuts the error fourfold, from 9.2 x 10^-7 to 5.7 x 10^-8 over two halvings. Yet an electric cutoff of one flux unit overstates pair creation by 30 percent, at full precision. A rule can converge, conserve and still be wrong. Updates that freeze the links, or wrap the flux around, keep the norm and the total charge and break Gauss's law.

**Wrong couplings pass consistency checks.** On a curved expanding background (Budd, Florio, Frenklakh and Mukherjee, 2026), the electric term must scale as the square of the conformal factor. Rules that scale it by the first power or not at all stay unitary, keep Gauss's law, agree with themselves across clocks and satisfy their own energy-work balance, and end in states 0.76 and 1.10 from the right one. Agreement under a change of clock is bookkeeping, not evidence of the right coupling to geometry.

**A continuum benchmark.** With 20 sites, every neutral state kept and the published lattice mass correction (Dempsey, Klebanov, Pufu and Zan, 2022), the interacting lattice reproduces the Schwinger boson's lowest gap within 0.26 to 0.62 percent, and the probability that the vacuum survives an expansion, 0.93264 against the continuum's 0.93280. The excitation energy is still 2 percent low, and the earlier six-site massive result changes with volume: convergence depends on what is measured.

**The source of gravity.** Three lattices share the long-wavelength wave equation and conserve energy to within 5 x 10^-10. One of them, the squared centered derivative, hides a second field at the edge of its momentum zone, and its quantum source of gravity, energy and pressure, comes out twice the single field's at twice the spacing, to within 2 x 10^-10: O_centered(ε) = 2 O_nearest(2ε). At the finest spacing its peak energy and pressure are 2.00004 and 2.00003 times the target. Lattice doublers are known (Buonanno, Cella and Curci, 1995); here a rule that passes every wave test would double the gravity its matter makes. A short rule must be tested on its source, not only on its waves.

**A step has its own vacuum.** A discrete time step has an invariant vacuum of its own. Start the field in the continuum's vacuum instead and the apparent energy after one step grows as the spacing shrinks, as 1/ε^2: 0.23, 0.91, 3.64 and 14.6 over three halvings, while the step's own vacuum stays empty, with occupations below 4 x 10^-30. With the matched vacuum, space and time refine together and the pressure error falls to 0.002 percent.

**Matter drives the expansion.** Coupled through one Hamiltonian to a homogeneous 3+1 universe (Husain and Singh, 2021), quantum matter in a vacuum, a squeezed state or with a vacuum energy drives the expansion with scale-factor errors near 10^-6 and the Friedmann constraint held within 3 x 10^-5, in proper and conformal time alike; the last halving cuts every error by 3.997 to 4.008. Freeze the feedback and the matter evolves perfectly while the Friedmann constraint fails by 85 percent and the two clocks disagree by 7.6 percent. The same model shows what remains open. At a fixed box the bare vacuum energy grows as 1/ε^4 and its pressure approaches a third of its energy, so canceling the energy with a constant leaves the pressure: the vacuum's contribution to gravity cannot be set by canceling one number.

**The end of evaporation.** In the solvable model of an evaporating black hole of Russo, Susskind and Thorlacius (1992), the traditional final mass is negative, close to minus one unit of the model's quantum scale, while the Hamiltonian charge of the complete semiclassical action (Barenboim, Frolov and Kunstatter, 2026) is the same amount positive. Pairing the complete-action charge with the traditional flux leaves an error of twice the quantum scale, 49 percent of the mass for a hole of four units and more for lighter holes, that no step size removes. No simple rule for the end, stopping at zero mass, switching branch or flooring the field, gives a regular ending: as at the center in Section 8.4, the endpoint needs the rule's own dynamics.

### 8.7 A proposed specification on a supplied lattice

Section 7 points to a candidate from the top: supply a four-dimensional lattice as data, run the known laws on it, and let general relativity bend it. This section writes the candidate out, completes it, tests the lattice forms it needs and counts its bits with the compressor of Section 7 (Table 12; the supplement's lattice, form-fermion and Regge scripts).

**Lattice form.** The lattice is hypercubic with the Planck length as its spacing, and it is physical, with no continuum limit to reach. Gauge fields become phases on its links (Wilson, 1974), and the metric becomes the lengths of the edges of a triangulation, so general relativity acts through the Regge action (Regge, 1961). Naive lattice fermions come with doublers (Nielsen and Ninomiya, 1981); the overlap operator (Neuberger, 1998) removes them and keeps an exact lattice chiral symmetry, the Ginsparg-Wilson relation (Ginsparg and Wilson, 1982; Lüscher, 1998). On a two-dimensional U(1) lattice the naive operator has 8 zero modes in the free field where the continuum has 2, and in 12 backgrounds of topological charge -2 to 3, smooth and rough (plaquette angles up to 3.1 radians), the overlap operator's index equals the charge every time while the relation holds to $4 \times 10^{-15}$. This verifies the overlap construction in the tested U(1) backgrounds; it does not test four-dimensional SU(3) dynamics. The weak force acts on left-handed fields only. Its lattice form exists for U(1) (Lüscher, 1999) and is proposed for the full Standard Model (Wang and Wen, 2020), against objections (Golterman and Shamir, 2023, 2026). In 1+1 dimensions the route has worked numerically: in the ungauged 3-4-5-0 model, whose chirality lies in its charge assignments, the mirror sector gaps without breaking the symmetry while the light sector stays gapless (Zeng, Zhu, Wang and You, 2022). Gauged versions have been proposed (Wang, 2022), but a numerical test with a fluctuating gauge field, and a selective mirror gap in 3+1 dimensions, are still open. Simplicity favors that route: a local lattice rule gives both handednesses in equal numbers, so a short rule keeps both and makes the mirror hand heavy at the lattice scale without breaking the symmetries, an echo of the mirror particles that Lee and Yang (1956) proposed. The anomalies must cancel, including the mod-16 anomaly of 5(B - L) - 4Y: the statement's 16 Weyl fermions per family cancel them all, where 15 would leave 13 (García-Etxebarria and Montero, 2019). Perturbative calculations can reproduce many measured weak-interaction observables. They do not complete the proposed nonperturbative chiral lattice theory; the count below measures its specification.

**Lorentz symmetry.** A lattice that treats its four directions alike forbids by hypercubic symmetry the terms that would give particles different top speeds, and its time step is the transfer matrix along one direction; a step that treats time apart, like that of Section 8.1, would need one tuned speed per particle type (Collins et al., 2004). The photon then has $E = p - a^2p^3/12$ along an axis and $E = p - a^2p^3/18$ along a body diagonal, so LHAASO's limit on quadratic dispersion, $E_{QG,2} > 6.9 \times 10^{11}$ GeV (Cao et al., 2024), allows spacings up to $7.0 \times 10^{-28}$ m, about 43 million Planck lengths.

**Completion.** The Standard Model with general relativity lacks dark matter and an initial state. The shortest completion declares them: a free scalar with a declared mass, and five numbers for the initial state. SMASH (Ballesteros, Redondo, Ringwald and Tamarit, 2017) derives them instead, from an axion, a heavy quark and Majorana neutrino masses, and predicts a tensor-to-scalar ratio of at least 0.004 and an axion of 50 to 200 μeV; it costs 928 more bits and would earn them only if those predictions are measured. Dark energy may need no entry: in the timescape model clocks run faster in emptier regions, so the averaged expansion only appears to accelerate (Wiltshire, 2007), which a recent supernova analysis mildly favors (Seifert, Lane, Galoppo, Ridden-Harper and Wiltshire, 2025) against arguments that such effects are too small (Green and Wald, 2014); then Λ is zero. Inputs may also become rules at the lattice scale: assuming asymptotically safe gravity drives the Higgs self-coupling close to zero there, Shaposhnikov and Wetterich (2010) predicted a Higgs mass near 126 GeV before its discovery, against 125.2 measured.

**What the lattice leaves open.** A physical lattice regulates short-distance modes without by itself establishing convergence of the gravitational integral. The parameter estimates also expose fine-tuning: the bare cosmological constant must cancel the vacuum energy of the lattice's fields to about 123 decimal places, and the bare Higgs mass parameter its shift to about 32. As numbers those digits would cost 408 and 106 bits; as the condition that the measured values result, as the statement writes them, they cost only that condition, so the renormalization condition has a short description. The target values, their precision and a procedure for determining the bare parameters must still be specified; compression does not explain those values. Every measurement lies far below the lattice scale, where general relativity as an effective quantum theory gives finite predictions (Donoghue, 1994), which motivates the low-energy target; an empirical fit of the complete supplied-lattice construction has not been performed here. Joining general relativity with quantum theory at black-hole centers and the first instants of expansion turns on gravity at the lattice's own scale, where the Euclidean action is unbounded below (Gibbons, Hawking and Perry, 1978). That test has begun. On a fixed 16⁴ Regge lattice Hamber and Williams (1995) found an attractive potential near a critical point, consistent with Newton's law within large errors; dynamical triangulations recover Newtonian binding within 20 to 30 percent, without matter acting back (Dai, Laiho, Schiffer and Unmuth-Yockey, 2021); and causal triangulations, which glue four-dimensional blocks causally, find a four-dimensional de Sitter universe (Ambjørn, Jurkiewicz and Loll, 2004; Ambjørn, Görlich, Jurkiewicz and Loll, 2008). On the supplied lattice itself, 4⁴ sites with λ = 1 and a uniform measure, the original supplement finds a long-lived negative-average-curvature regime after a flat start whose deficit angles are below $2 \times 10^{-15}$: the average curvature is negative and shrinks from -1.84 at k = 0 to -0.77 at k = 0.065, and above k ≈ 0.07 the lattice degenerates, near the critical value $k_c = 0.0639$ that Hamber (2015) finds on lattices up to 32⁴ with exponent $\nu = 0.334(4)$. On lattices this small the curvature does not reach zero first, so the run locates an apparent finite-run stability boundary, not a critical point. The independent execution in Section 8.9 tests its dependence on initial conditions and sampling history. Positive curvature-squared terms extend the smooth phase (Hamber, 1992), and Hamber (2019) proposes that the observed cosmological constant is set by the correlation length near the critical point, $\Lambda \approx 3/\xi^2$, which leaves that length to be fixed.

**Matter on causal triangulations.** Fermions there need frames on a changing geometry, and each change alters their determinant. Written as differential forms, Kähler-Dirac fields, they need no frames (Kähler, 1962; Becher and Joos, 1982; Catterall, Laiho and Unmuth-Yockey, 2018). On six sampled two-dimensional causal triangulations of a torus (Ambjørn and Loll, 1998) such a field has exactly 4 zero modes, the torus's Betti numbers, and low modes 0.92 to 1.10 times Weyl's law, in the measured low spectrum; those counts alone do not exclude every extra propagating species. On the supplied four-dimensional lattice, bent by random weights, its zero modes number 16. Each such field carries four copies of a Dirac fermion. Catterall (2021) proposed interactions that would leave 16 Weyl fermions if they gap one half and the continuum limit holds, tested in two dimensions, and four-dimensional gauge theories show symmetric gapped phases (Hasenfratz, 2022; Butt, Catterall and Hasenfratz, 2025), though none yet leaves the Standard Model's chiral spectrum. Free conformal fields in the Standard Model's spin counts give a tractable first proxy for quantum stress. Masses, curvature couplings, interactions and the quantum state remain necessary for a full matter calculation.

**Description length.** Written out in full (the statement, its lattice form and the declared completion), the candidate compresses to 16,272 bits, a description length in this encoding for the proposed specification with dark matter and initial data. A compressed description alone does not establish an executable physical model. The search from the bottom excludes candidates only within its tested classes and conditions (Sections 4 to 6), and no program can be proved the shortest for its output once its length exceeds a constant set by the axioms used (Chaitin, 1974). A demonstrated program gives an upper bound for its stated output and tests, while the present census excludes candidates within its tested classes. These are different bounds; neither establishes the minimum length of a complete physical theory.

**Table 12.** The proposed specification on a supplied four-dimensional lattice: each piece, its cost in compressed bits (raw LZMA2, as in Section 7) and its status.

| Piece | On the supplied lattice | Bits | Status |
|---|---|---:|---|
| Standard Model with general relativity, 28 inputs | the statement of Section 7 | 11,328 | complete as stated |
| Lattice form of each field | link phases, Regge edge lengths, overlap fermions, transfer matrix | +3,144 | 2D U(1) overlap checks; chiral completion unresolved |
| Dark matter and initial state, declared | free scalar; amplitude, tilt, two mass ratios, temperature | +1,800 | declared inputs; no new full-data fit |
| Total: specification | | 16,272 | description length; full dynamics unverified |
| Mechanisms instead (SMASH) | axion, seesaw, leptogenesis, inflation | +928 more | predicts r of at least 0.004, axion of 50 to 200 μeV |
| Fine-tuning of Λ and the Higgs mass | renormalization condition and declared targets | included in description (digit estimate: 514) | parameters not explained |
| Lorentz symmetry | hypercubic symmetry; corrections of order (aE)² | 0 | spacing up to 7.0 × 10⁻²⁸ m allowed |
| Gravity at the lattice scale | Regge path integral | 0 beyond the lattice form | open; finite-run regimes here; see Section 8.9; larger-lattice critical behavior in cited work |

### 8.8 Scope

These tests supply the geometry and the fields; the rule that grows them is what the first half of this paper searches for. Discrete does not mean finite. A short rule can add sites without end, and an infinite causal set can be locally finite (Surya, 2019); the arrays here are finite for computing, and their boundaries were checked. Earlier work reached the same kind of agreement for the Dirac equation (D'Ariano and Perinotti, 2014), for quantum walks outside and inside a horizon (Di Molfetta, Brachet and Debbasch, 2013) and for quantum electrodynamics in three dimensions of space (Eon, Di Molfetta, Magnifico and Arrighi, 2023). Section 8.6 runs the next tests listed in earlier passes: a time-dependent background against particle production, a dynamical gauge field that conserves Gauss's law, an interacting process, and matter acting back on the geometry. What remains is a geometry that propagates locally, a vacuum subtraction derived from the rule itself, and the endpoints that need the rule's own dynamics. The scripts, independent references and results are in the supplements.

### 8.9 Independent execution and equilibration audit

A new audit reran the supplied Regge C engine without changing its action, measure or update rules. The protocol was fixed before the main runs. At each of L = 3 and 4 and k = 0.05 and 0.065, three seeded chains began flat and three began from independently prepared rough configurations. Every rough preparation used the same prescribed 3,000-sweep evolution at k = 0.1; no configuration was selected for its observed outcome. Each of the 24 main chains used 3,000 equilibration sweeps followed by 12,000 measurement sweeps. All per-sweep data, starting preparations and ending configurations are retained. These are new executions and independently written diagnostics, not a second implementation of the gravitational dynamics.

The geometry self-tests reproduced, including the flat-start deficit angles and incremental/full recomputation check. The statistical audit used ordinary split R-hat across all six chains at each parameter pair, with a prespecified screening threshold of 1.05, and batch-mean, drift and start-dependence diagnostics. Passing this screen would not prove equilibrium; failing it prevents treating the pooled chains as a settled ensemble. Curvature means below are averages over each start class and are descriptive when the chains have not mixed.

| L | k | Flat-start mean R | Rough-start mean R | Split R-hat | Screen |
|---|---|---|---|---|---|
| 3 | 0.05 | -1.0470 | -0.6731 | 1.573 | fails |
| 3 | 0.065 | -0.6022 | 44.1443 | 13.032 | fails |
| 4 | 0.05 | -1.0655 | 31.2166 | 5.498 | fails |
| 4 | 0.065 | -0.6793 | 45.4855 | 17.716 | fails |

All four parameter groups fail the prespecified R-hat screen. The table and trajectories distinguish initial-state dependence from an established equilibrium phase. Negative mean curvature after a flat start is insufficient to demonstrate that the ensemble has equilibrated. These finite-run tests neither locate a critical point nor establish the nonexistence of an equilibrium limit; they restrict what can be claimed from the sampled histories.

All 24 chains satisfy the scaling identity in their measurement averages to within 0.09 percent, despite the failed mixing screen. This separates a successful overall-scale consistency check from equilibration of the geometry. It provides a concrete example of an implementation check that cannot establish the sampled phase.

The next tests are normalizability of the specified ensemble and convergence from both start classes at longer durations, followed by finite-size scaling and curvature correlations only where convergence is supported. A Newtonian potential and a Lorentzian matter-coupled theory require separate tests. The original run's histories are kept in its archive, and this audit supplies its own complete histories and trajectory figure in the audit archive.


## 9. Conclusion

Among the shortest rewriting rules, under the update orders tried, none grows a four-dimensional causal order. Of the 2,935 rules that rewrite two relations into three and samples of 3,500 from three larger classes, none reads four on all three of causal-set theory's dimension tests, in runs of up to 144,000 updates. Two-dimensional orders do appear under the first passes' order, so the tests can find a dimension when a short rule has one; under random orders the same rules grow orders that are nearly one-dimensional on two tests, so a rule's update order belongs in its length. The near misses read four on two tests and 1.7 to 3.0 on the third: wide like four dimensions, with chains far longer than four dimensions allow, the opposite of an expanding universe. A rule that rewrites a few relations at a time bounds the links of every update, while a random four-dimensional causal set gives its elements ever more links. Such a rule can therefore match our spacetime only at large scales, where its small-scale frame must stay hidden. A criterion fixed in advance with noise bands puts the near misses' chains 4.7 to 21 standard deviations from four dimensions. Sequential growth escapes the link bound and still misses: transitive percolation held at four on the first test, followed to two million elements, reaches four on the third while its midpoint drifts to 5.8, and originary percolation, read against de Sitter space, runs its chains too long. The search goes on with longer rules, following the readings as the intervals grow, and from the other end, with the known laws made discrete. From that end, one family of local unitary steps already reproduces quantum fields in curved spacetime and crossing a black-hole horizon, with errors falling fourfold as the spacing halves, its single-particle Dirac step fitting in a program of 2,488 bits, and the same steps, run as that program unchanged, create pairs in an expanding universe, while a separate 3+1 lattice model drives its expansion. The proposed specification on a supplied four-dimensional lattice, with dark matter and initial data, compresses to 16,272 bits. Component tests support parts of this construction; they do not establish a single executable model covering present data. Smooth large-scale geometry, a consistent chiral matter sector and the quantum-gravity dynamics remain separate requirements. The independent Regge audit fails the prespecified mixing screen in all four parameter groups, so the sampled regime is not established as an equilibrated smooth phase. Harder tests sharpen the search: a rule must be judged by the source of gravity it makes, since one lattice that matches every wave doubles it. What remains is the short rule that grows the geometry they need, and the dynamics it gives a black hole's center and a black hole's end.

## References

Abdo, A. A., et al. (2009). A limit on the variation of the speed of light arising from quantum gravity effects. *Nature*, 462, 331-334. https://arxiv.org/abs/0908.1832

Ahmed, M., & Rideout, D. (2010). Indications of de Sitter spacetime from classical sequential growth dynamics of causal sets. *Physical Review D*, 81, 083528. https://arxiv.org/abs/0909.4771

Ahmed, M., Dodelson, S., Greene, P. B., & Sorkin, R. (2004). Everpresent Λ. *Physical Review D*, 69, 103523. https://arxiv.org/abs/astro-ph/0209274

Ambjørn, J., & Loll, R. (1998). Non-perturbative Lorentzian quantum gravity, causality and topology change. *Nuclear Physics B*, 536, 407-434. https://arxiv.org/abs/hep-th/9805108

Ambjørn, J., Görlich, A., Jurkiewicz, J., & Loll, R. (2008). Planckian birth of a quantum de Sitter universe. *Physical Review Letters*, 100, 091304. https://arxiv.org/abs/0712.2485

Ambjørn, J., Jurkiewicz, J., & Loll, R. (2004). Emergence of a 4D world from causal quantum gravity. *Physical Review Letters*, 93, 131301. https://arxiv.org/abs/hep-th/0404156

Ambjørn, J., Jurkiewicz, J., & Loll, R. (2005). Spectral dimension of the universe is scale dependent. *Physical Review Letters*, 95, 171301. https://arxiv.org/abs/hep-th/0505113

Ashtekar, A., Pawlowski, T., & Singh, P. (2006). Quantum nature of the big bang. *Physical Review Letters*, 96, 141301. https://arxiv.org/abs/gr-qc/0602086

Ballesteros, G., Redondo, J., Ringwald, A., & Tamarit, C. (2017). Unifying inflation with the axion, dark matter, baryogenesis, and the seesaw mechanism. *Physical Review Letters*, 118, 071802. https://arxiv.org/abs/1608.05414

Barenboim, J., Frolov, A. V., & Kunstatter, G. (2026). Quantum-corrected Bondi mass for 2D Hawking radiation. *Physical Review D*, 114, 024067. https://arxiv.org/abs/2510.18298

Becher, P., & Joos, H. (1982). The Dirac-Kähler equation and fermions on the lattice. *Zeitschrift für Physik C*, 15, 343-365.

Benincasa, D. M. T., & Dowker, F. (2010). The scalar curvature of a causal set. *Physical Review Letters*, 104, 181301. https://arxiv.org/abs/1001.2725

Bolognesi, T. (2010). Causal sets from simple models of computation. *International Journal of Unconventional Computing*, 6(6), 489-524. https://arxiv.org/abs/1004.3128

Bombelli, L., Henson, J., & Sorkin, R. D. (2009). Discreteness without symmetry breaking: A theorem. *Modern Physics Letters A*, 24, 2579-2587. https://arxiv.org/abs/gr-qc/0605006

Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Space-time as a causal set. *Physical Review Letters*, 59, 521-524.

Brightwell, G., & Gregory, R. (1991). Structure of random discrete spacetime. *Physical Review Letters*, 66, 260-263.

Budd, E., Florio, A., Frenklakh, D., & Mukherjee, S. (2026). Quantum dynamics of cosmological particle production: Interacting quantum field theories with matrix product states. *Journal of High Energy Physics*, 2026(04), 183. https://arxiv.org/abs/2601.02331

Bueno, P., Cano, P. A., Hennigar, R. A., & Murcia, Á. J. (2025). Regular black hole formation in four-dimensional non-polynomial gravities. arXiv:2509.19016.

Buonanno, A., Cella, G., & Curci, G. (1995). Lattice energy-momentum tensor with Symanzik improved actions. *Physical Review D*, 51, 4494. https://arxiv.org/abs/hep-lat/9506008

Butt, N., Catterall, S., & Hasenfratz, A. (2025). Symmetric mass generation with four SU(2) doublet fermions. *Physical Review Letters*, 134, 031602. https://arxiv.org/abs/2409.02062

Cao, Z., et al. (LHAASO Collaboration). (2024). Stringent tests of Lorentz invariance violation from LHAASO observations of GRB 221009A. *Physical Review Letters*, 133, 071501. https://arxiv.org/abs/2402.06009

Carlip, S. (2017). Dimension and dimensional reduction in quantum gravity. *Classical and Quantum Gravity*, 34, 193001. https://arxiv.org/abs/1705.05417

Catterall, S. (2021). Chiral lattice fermions from staggered fields. *Physical Review D*, 104, 014503. https://arxiv.org/abs/2010.02290

Catterall, S., Laiho, J., & Unmuth-Yockey, J. (2018). Kähler-Dirac fermions on Euclidean dynamical triangulations. *Physical Review D*, 98, 114503. https://arxiv.org/abs/1810.10626

Chaitin, G. J. (1974). Information-theoretic limitations of formal systems. *Journal of the ACM*, 21, 403-424.

Collins, J., Perez, A., Sudarsky, D., Urrutia, L., & Vucetich, H. (2004). Lorentz invariance and quantum gravity: An additional fine-tuning problem? *Physical Review Letters*, 93, 191301. https://arxiv.org/abs/gr-qc/0403053

Corley, S., & Jacobson, T. (1996). Hawking spectrum and high frequency dispersion. *Physical Review D*, 54, 1568-1586. https://arxiv.org/abs/hep-th/9601073

D'Ariano, G. M., & Perinotti, P. (2014). Derivation of the Dirac equation from principles of information processing. *Physical Review A*, 90, 062106. https://arxiv.org/abs/1306.1934

Dai, M., Laiho, J., Schiffer, M., & Unmuth-Yockey, J. (2021). Newtonian binding from lattice quantum gravity. arXiv:2102.04492. https://arxiv.org/abs/2102.04492

Dempsey, R., Klebanov, I. R., Pufu, S. S., & Zan, B. (2022). Discrete chiral symmetry and mass shift in the lattice Hamiltonian approach to the Schwinger model. *Physical Review Research*, 4, 043133. https://arxiv.org/abs/2206.05308

Deser, S., Jackiw, R., & 't Hooft, G. (1984). Three-dimensional Einstein gravity: Dynamics of flat space. *Annals of Physics*, 152, 220-235.

Di Molfetta, G., Brachet, M., & Debbasch, F. (2013). Quantum walks as massless Dirac fermions in curved space-time. *Physical Review A*, 88, 042301. https://arxiv.org/abs/1212.5821

Donoghue, J. F. (1994). General relativity as an effective field theory: The leading quantum corrections. *Physical Review D*, 50, 3874-3888. https://arxiv.org/abs/gr-qc/9405057

Dou, D., & Sorkin, R. D. (2003). Black hole entropy as causal links. *Foundations of Physics*, 33, 279-296. https://arxiv.org/abs/gr-qc/0302009

Eon, N., Di Molfetta, G., Magnifico, G., & Arrighi, P. (2023). A relativistic discrete spacetime formulation of 3+1 QED. *Quantum*, 7, 1179. https://arxiv.org/abs/2205.03148

Fulgado-Claudio, C., Sánchez Velázquez, J. M., & Bermudez, A. (2023). Fermion production at the boundary of an expanding universe: A cold-atom gravitational analogue. *Quantum*, 7, 1042. https://arxiv.org/abs/2212.01355

García-Etxebarria, I., & Montero, M. (2019). Dai-Freed anomalies in particle physics. *Journal of High Energy Physics*, 2019(8), 003. https://arxiv.org/abs/1808.00009

Gibbons, G. W., Hawking, S. W., & Perry, M. J. (1978). Path integrals and the indefiniteness of the gravitational action. *Nuclear Physics B*, 138, 141-150.

Ginsparg, P. H., & Wilson, K. G. (1982). A remnant of chiral symmetry on the lattice. *Physical Review D*, 25, 2649-2657.

Glaser, L., & Surya, S. (2013). Towards a definition of locality in a manifoldlike causal set. *Physical Review D*, 88, 124026. https://arxiv.org/abs/1309.3403

Golterman, M., & Shamir, Y. (2023). Propagator zeros and lattice chiral gauge theories. arXiv:2311.12790. https://arxiv.org/abs/2311.12790

Golterman, M., & Shamir, Y. (2026). Symmetric mass generation and the Nielsen-Ninomiya theorem. arXiv:2603.15985. https://arxiv.org/abs/2603.15985

Gorard, J. (2020a). Algorithmic causal sets and the Wolfram model. arXiv:2011.12174. https://arxiv.org/abs/2011.12174

Gorard, J. (2020b). Some relativistic and gravitational properties of the Wolfram model. *Complex Systems*, 29(2), 599-654.

Green, S. R., & Wald, R. M. (2014). How well is our universe described by an FLRW model? *Classical and Quantum Gravity*, 31, 234003.

Hamber, H. W. (1992). Phases of four-dimensional simplicial quantum gravity. *Physical Review D*, 45, 507-512.

Hamber, H. W. (2015). Scaling exponents for lattice quantum gravity in four dimensions. *Physical Review D*, 92, 064017. https://arxiv.org/abs/1506.07795

Hamber, H. W. (2019). Vacuum condensate picture of quantum gravity. *Symmetry*, 11(1), 87. https://arxiv.org/abs/1707.08188

Hamber, H. W., & Williams, R. M. (1995). Newtonian potential in quantum Regge gravity. *Nuclear Physics B*, 435, 361-397. https://arxiv.org/abs/hep-th/9406163

Hasenfratz, A. (2022). Emergent strongly coupled ultraviolet fixed point in four dimensions with eight Kähler-Dirac fermions. *Physical Review D*, 106, 014513. https://arxiv.org/abs/2204.04801

Herzog, C. P., & Huang, K.-W. (2013). Stress tensors from trace anomalies in conformal field theories. *Physical Review D*, 87, 081901.

Husain, V., & Singh, S. (2021). Quantum backreaction on a classical universe. *Physical Review D*, 104, 124048. https://arxiv.org/abs/2109.12752

Jacobson, T. (1991). Black-hole evaporation and ultrashort distances. *Physical Review D*, 44, 1731-1739.

Kähler, E. (1962). Der innere Differentialkalkül. *Rendiconti di Matematica e delle sue Applicazioni*, 21, 425-523.

Lamm, H., Lawrence, S., & Yamauchi, Y. (2019). General methods for digital quantum simulation of gauge theories. *Physical Review D*, 100, 034518. https://arxiv.org/abs/1903.08807

Lee, T. D., & Yang, C. N. (1956). Question of parity conservation in weak interactions. *Physical Review*, 104, 254-258.

Lüscher, M. (1998). Exact chiral symmetry on the lattice and the Ginsparg-Wilson relation. *Physics Letters B*, 428, 342-345. https://arxiv.org/abs/hep-lat/9802011

Lüscher, M. (1999). Abelian chiral gauge theories on the lattice with exact gauge invariance. *Nuclear Physics B*, 549, 295-334. https://arxiv.org/abs/hep-lat/9811032

Martinez, E. A., et al. (2016). Real-time dynamics of lattice gauge theories with a few-qubit quantum computer. *Nature*, 534, 516-519. https://arxiv.org/abs/1605.04570

Meyer, D. A. (1988). *The dimension of causal sets* (Doctoral dissertation). Massachusetts Institute of Technology.

Moore, C. (1988). Comment on "Space-time as a causal set". *Physical Review Letters*, 60, 655.

Myrheim, J. (1978). *Statistical geometry* (CERN preprint TH-2538). CERN.

Neuberger, H. (1998). Exactly massless quarks on the lattice. *Physics Letters B*, 417, 141-144. https://arxiv.org/abs/hep-lat/9707022

Nielsen, H. B., & Ninomiya, M. (1981). Absence of neutrinos on a lattice (I): Proof by homotopy theory. *Nuclear Physics B*, 185, 20-40.

Regge, T. (1961). General relativity without coordinates. *Il Nuovo Cimento*, 19, 558-571.

Reid, D. D. (2003). Manifold dimension of a causal set: Tests in conformally flat spacetimes. *Physical Review D*, 67, 024034. https://arxiv.org/abs/gr-qc/0207103

Rideout, D. P., & Sorkin, R. D. (2000). Classical sequential growth dynamics for causal sets. *Physical Review D*, 61, 024002. https://arxiv.org/abs/gr-qc/9904062

Rideout, D. P., & Sorkin, R. D. (2001). Evidence for a continuum limit in causal set dynamics. *Physical Review D*, 63, 104011. https://arxiv.org/abs/gr-qc/0003117

Rovelli, C., & Vidotto, F. (2014). Planck stars. *International Journal of Modern Physics D*, 23, 1442026. https://arxiv.org/abs/1401.6562

Russell, B. (1903). *The Principles of Mathematics*. Cambridge University Press.

Russo, J. G., Susskind, L., & Thorlacius, L. (1992). End point of Hawking radiation. *Physical Review D*, 46, 3444. https://arxiv.org/abs/hep-th/9206070

Seifert, A., Lane, Z. G., Galoppo, M., Ridden-Harper, R., & Wiltshire, D. L. (2025). Supernovae evidence for foundational change to cosmological models. *Monthly Notices of the Royal Astronomical Society: Letters*, 537, L55-L60.

Shaposhnikov, M., & Wetterich, C. (2010). Asymptotic safety of gravity and the Higgs boson mass. *Physics Letters B*, 683, 196-200. https://arxiv.org/abs/0912.0208

Shaw, A. F., Lougovski, P., Stryker, J. R., & Wiebe, N. (2020). Quantum algorithms for simulating the lattice Schwinger model. *Quantum*, 4, 306. https://arxiv.org/abs/2002.11146

Solomonoff, R. J. (1964). A formal theory of inductive inference, Part I. *Information and Control*, 7(1), 1-22.

Sorkin, R. D. (1997). Forks in the road, on the way to quantum gravity. *International Journal of Theoretical Physics*, 36, 2759-2781. https://arxiv.org/abs/gr-qc/9706002

Surya, S. (2012). Evidence for the continuum in 2D causal set quantum gravity. *Classical and Quantum Gravity*, 29, 132001. https://arxiv.org/abs/1110.6244

Surya, S. (2019). The causal set approach to quantum gravity. *Living Reviews in Relativity*, 22, 5. https://arxiv.org/abs/1903.11544

Unruh, W. G. (1995). Sonic analogue of black holes and the effects of high frequencies on black hole evaporation. *Physical Review D*, 51, 2827-2838. https://arxiv.org/abs/gr-qc/9409008

Wang, J. (2022). CT or P problem and symmetric gapped fermion solution. *Physical Review D*, 106, 125007.

Wang, J., & Wen, X.-G. (2020). A non-perturbative definition of the Standard Models. *Physical Review Research*, 2, 023356. https://arxiv.org/abs/1809.11171

Wilson, K. G. (1974). Confinement of quarks. *Physical Review D*, 10, 2445-2459.

Wiltshire, D. L. (2007). Cosmic clocks, cosmic variance and cosmic averages. *New Journal of Physics*, 9, 377.

Wolfram, S. (2002). *A new kind of science*. Wolfram Media.

Wolfram, S. (2020). A class of models with the potential to represent fundamental physics. *Complex Systems*, 29(2), 107-536. https://arxiv.org/abs/2004.08210

Zeng, M., Zhu, Z., Wang, J., & You, Y.-Z. (2022). Symmetric mass generation in the 1+1 dimensional chiral fermion 3-4-5-0 model. *Physical Review Letters*, 128, 185301. https://arxiv.org/abs/2202.12355
