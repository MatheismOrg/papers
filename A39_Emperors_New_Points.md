# The Emperor's New Points: Zero-Length Points and Physical Space

*Continuum Models and Physical Reality*

Gary Abraham Bernstein

Independent researcher · ORCID: 0009-0009-1761-2867

## Abstract

A line has length although each of its points has none. How does the mathematics account for that, and what does it imply about physical space? Standard measure theory distinguishes taking a union of points from adding their individual lengths, so the zero-length argument is not a mathematical contradiction. The physical question remains: which structure supplies measured distances? This paper develops the case for discrete construction under an explicit generation premise and shows what it must recover from continuous models. Simple examples reveal both agreement in measured quantities and errors that can survive a finer grid.

## 1. The question that calculation alone does not answer

A line has length. Each individual point has zero length. To understand both statements, we must distinguish collecting points into a set from adding lengths. Resolving that mathematical question leaves a physical one: what structure accounts for the distances we measure?

The apparent contradiction can be stated sharply: if every constituent contributes zero width, how can the whole have nonzero width? Its force comes from a familiar rule. Split a rod into two pieces and add their lengths. Repeat finitely, and the rule still works. The mistake is to extend that rule to an arbitrary uncountable decomposition without checking what the mathematical operation means.

The appropriate response explains the operations. Merely saying “infinity” leaves the central question unanswered. Infinity does not turn a sum of zeros into one.

### 1.1 The position advanced here

My position is that physical continua cannot exist. Continuum mathematics can represent physical quantities, including exact collective quantities, without every coordinate in that representation being a physical location. I propose that physical extension rests on discrete relations and that smooth geometry describes their collective behaviour.

This is a rejection of a fundamental physical continuum, not a claim that the real interval is an undefined mathematical object. The two objections motivating the position are direct: a continuum has no next point, and zero singleton widths cannot add to a positive width. Their mathematical content must be stated accurately. Dense order permits the first property, while measure theory does not assert the disputed addition. These objections therefore motivate a demand for a physical account; they do not themselves establish physical impossibility.

The stronger conclusion requires an additional physical premise and its defence. The bearer proposal below asks what physically accounts for extension. Its force depends on whether it excludes a relational continuum for an independently justified reason. A demand for a physical account is shared by both candidates; a demand that the account must be discrete adds the disputed premise.

### 1.2 The bearer proposal and its scope

The bearer requirement is that measured extension have a physical account. In the usual construction of Lebesgue measure on the real line, interval lengths supply the starting size data. A singleton has measure zero because it fits inside intervals of arbitrarily small length. One can also define an uncountable nonnegative sum as the supremum of finite subsums; the sum of singleton measures is then zero. This operation is not the measure of their union.

The physical proposal is to account for distances through discrete relations and their dynamics. Calling those relations grains does not yet explain size. If a grain has positive extent, its size remains part of the account; if elementary sites have zero extent, their separations must come from relations. Discreteness alone supplies no unit of length.

A defender of continuous space can make the same relational move: points have no individual width, while a metric assigns distances between them. The interval and its metric structure are a candidate bearer of extension. They need not be conceived as a pile of miniature lengths. Thus “nothing bears extension” does not follow from “every singleton has zero length.” This is a substantive alternative, not merely permission to leave extension unexplained.

A stronger premise would require every extended physical region to decompose into elementary parts whose nonnegative individual extents sum to the region's extent. Together with a decomposition into zero extent points, that premise contradicts positive extension. The inference is valid conditional on the premise. But the premise imposes precisely the additive composition that the continuum model does not assert, and some discrete relational models do not assert it either. Its physical justification remains to be supplied.

I pursue the discrete account because it offers a concrete construction to investigate: specify elementary relations, derive measured geometry, and recover quantum and gravitational predictions from the same dynamics. The preference becomes an explanatory advantage when the construction earns those recoveries. The bearer requirement identifies the question; the model must provide the answer.

## 2. A union is not a sum of lengths

Work in the usual real line with Lebesgue measure, written λ. Then

$$[0,1]=\bigcup_{x\in[0,1]}\{x\},\qquad \lambda(\{x\})=0,\qquad \lambda([0,1])=1.$$

The first equality specifies membership: each member of the interval belongs to its singleton. It supplies no rule equating the measure of that union with an uncountable sum of singleton measures.

For any family of nonnegative numbers, the generalized sum is defined as the supremum of all finite subsums. Applying that definition gives

$$\sum_{x\in[0,1]}\lambda(\{x\})
=\sup_{F\subseteq[0,1],\ F\text{ finite}}\sum_{x\in F}0=0.$$

Consequently, the measure of this union is not the sum of those measures. Lebesgue measure is countably additive: for a sequence of disjoint measurable sets, the measure of the union equals the sum. The singletons of a real interval cannot be placed in such a sequence. This is the scope of the theorem, not an exception secretly introduced during a calculation. Hunter (2011) gives the standard measure theoretic framework.

Singleton values also fail to determine the measure. On the same interval, the measure defined by twice Lebesgue measure assigns zero to every singleton but assigns two to the whole interval. Specifying the points and their individual measures therefore leaves information unspecified. In ordinary geometry, the metric and the chosen normalization of length provide that information. A mathematical model includes these relations and conventions as well as its underlying set.

### 2.1 Ordered, but with no next point

An order does not require immediate neighbours. For any real positions a<b, the midpoint c=(a+b)/2 satisfies a<c<b. If b were the next point after a, c would contradict that claim. The conclusion is that no immediate successor exists in this order. The order itself remains defined.

Choosing two points selects endpoints; it does not create the order of the real line. That order is already specified by the relation <. Nor is the absence of a next point peculiar to an uncountable continuum. The rational numbers also have a rational midpoint between any two distinct members, despite being countable. Dense order and uncountability are different properties.

Rationals can be generated by a discrete process over integer ratios. For positive rationals, visit pairs of positive integers by increasing numerator plus denominator and keep only coprime pairs. Every positive rational then appears at a finite stage; zero and negative rationals can be included by interleaving. Each output has a next output in this process, although no rational has an immediate successor in numerical order. The two orders answer different questions. No sequence can enumerate all real numbers. If physical construction is required to consist of countably many successive outputs, this excludes enumerating every real position as a separate output. The physical requirement is an additional premise, not a consequence of numerical density.

The distinction matters physically. A procedure that advances by visiting the immediate successor of each real point cannot operate on a real interval, because those successors do not exist. But a continuous trajectory is not defined as that procedure. For example, x(t)=vt assigns a position to each time without requiring a next instant. This is a mathematical model, not a proof that physical time is continuous. Rejecting its physical realization requires a reason why nature must execute the unavailable successor procedure.

The same applies to a block account. A specified history is not required to come into existence by enumerating its positions. The proposed discrete ontology must be defended through its physical structure and consequences, not by attributing a point enumeration algorithm to every competing model.

### 2.2 Can motion cross a continuum?

In units where the speed is one, let x(t)=t for 0≤t≤1. The position is zero at time zero and one at time one. Every intermediate time has a corresponding position, although neither the times nor the positions have immediate successors. For every nonzero time increment within the domain, Δx/Δt=1. Taking its limit yields the same velocity; no division of zero by zero is performed.

One can partition the duration into successive intervals of lengths 1/2, 1/4, 1/8 and so on. The first n lengths sum to 1−2^(−n), and their series sums to one. These are positive durations, not zeros. The countably many intervals group uncountably many instants; their enumeration is not an enumeration of those instants.

This example establishes how motion is represented consistently in the continuum model. It does not establish that nature instantiates the model. A physical objection must explain why motion requires immediate successor steps or some other property absent from this representation. Infinite divisibility alone does not imply infinite elapsed time, and the absence of a last interval in this partition does not remove the endpoint t=1 from the model.

## 3. What an integral actually adds

Divide the unit interval into n cells of width 1/n. For the constant function f=1, every finite Riemann sum is

$$\sum_{i=1}^{n}1\cdot\frac1n=1.$$

The widths are positive at every finite stage. There is no last partition at which they become zero and the calculation continues by adding zeros. The limit concerns a sequence of whole sums. Replacing the summands by their individual limits while their number changes requires a justification that is absent here.

For a continuous function f, the right endpoint sums converge to its integral:

$$\lim_{n\to\infty}\frac1n\sum_{i=1}^{n}f(i/n)=\int_0^1f(x)\,dx.$$

The point i/n chooses a sampled value. The cell contributes the width. Confusing those jobs creates the impression that integration gives length to dimensionless samples. Introductory real analysis makes the limiting conditions explicit; see Lebl (2026), Chapter 5.

This construction presupposes an interval, its metric and a notion of length. It computes within that structure. It does not describe a physical process that manufactures space from point particles.

## 4. A finer grid still has zero length

Consider the nested grids

$$D_n=\{k/2^n:k=0,1,\ldots,2^n\},\qquad n=0,1,2,\ldots.$$

Each grid is finite and has Lebesgue measure zero. Their union D is countable, so it also has measure zero. Yet D is dense in the unit interval: every point of the interval is arbitrarily close to a point of D. Its closure is the entire interval, with measure one.

The geometric approximation can be quantified. The maximum distance from a point of [0,1] to the nearest point of D_n is half the grid spacing, or 2^(−n−1). Thus these compact grids converge to the interval in Hausdorff distance while their Lebesgue measures remain zero.

Nothing contradictory occurred. The increasing union is D, not its closure. Countable continuity of measure gives λ(D)=0. Hausdorff convergence and taking a closure do not generally preserve Lebesgue measure. The example identifies precisely which operation changes the answer.

A drawing of ever closer dots can suggest a solid segment. The drawing has finite resolution and physical marks with extent. The mathematical grid consists of singletons. Those are different objects, and neither picture licenses an uncountable addition of widths.

## 5. Discrete weights can approximate a continuous measure

A finite set of points can carry a different measure. Give each of n sample points weight 1/n:

$$\mu_n=\frac1n\sum_{i=1}^{n}\delta_{i/n},$$

where δ_a is the unit point mass at a. Every μ_n has total mass one. Its weights are assigned masses, not the Lebesgue lengths of its supporting singletons.

Suppose f satisfies |f(x)−f(y)|≤L|x−y| on [0,1]. Splitting the integral into the n cells gives

$$\left|\int f\,d\mu_n-\int_0^1f(x)\,dx\right|
\leq\sum_{i=1}^{n}\int_{(i-1)/n}^{i/n}|f(i/n)-f(x)|\,dx
\leq\frac{L}{2n}.$$

The last bound follows by integrating the distance from each point in a cell to its right endpoint. For a declared class of observables with a common Lipschitz bound L, choosing n sufficiently large makes the discrepancy uniformly small. This supplies a concrete comparison at stated resolution.

The restriction matters. Let D be the countable union of the supports of all these μ_n. Its indicator has integral one under every μ_n and zero under Lebesgue measure. Arbitrary measurable tests do not obey the preceding approximation guarantee. A physical application must justify which observations its apparatus permits, rather than selecting a convenient test class after seeing the result.

This example contains no zeros acquiring weight. It distinguishes three quantities: the length of a support, the mass assigned to that support, and the predictions obtained by integrating observables.

### 5.1 Exact agreement for an average

Representation is the broad term. Approximation describes agreement within a tolerance or through a limiting procedure. A representation can instead give an exact answer for a specified observable. For a>0, put equal probability at −a and +a. The mean position is exactly zero even in a model where only those two individual positions are allowed. The mean represents their weighted relation; it does not add zero as a third allowed position. This mathematical example does not establish which positions are available to a real particle.

Compare that discrete probability measure with the uniform probability measure on [−a,a]. Both have total mass one and mean zero. Their second moments differ:

$$\mathbb E_{\mathrm{discrete}}[X^2]=a^2,
\qquad \mathbb E_{\mathrm{uniform}}[X^2]=a^2/3.$$

Thus their agreement on the centroid is exact, while their agreement on the full distribution fails. Equal weights are an assumption; different weights change the mean. If the proposed particle geometry is centred on two positions, those positions determine a midpoint and separation under a specified metric. They do not uniquely determine a circle, surface or dynamical law.

This supplies a precise version of the suggested physical mechanism: discrete constituents can support exact collective quantities represented in a continuous description. A particle lattice proposal must specify its geometry, weights, observable map and dynamics before that agreement can be extended beyond the chosen average. An exact identity within a model and evidence that the model describes a real particle are separate achievements.

## 6. Specification, construction and physical realization

The constructive restriction can be stated precisely. Suppose every elementary spacetime event is introduced in a succession of stages indexed by the natural numbers, with finitely many new events at each stage. Label an event by its stage and its index within that stage. The generated events are then at most countable. If these events exhaust spacetime's physical constituents, an uncountable continuum of elementary events is excluded by the premise. This is a conditional countability argument. It neither derives a minimum spacing nor excludes continuous mathematical representations: the rational numbers already show that a countable set can be dense. The physical defence concerns why event generation exhausts the ontology.

“Construct the interval” can mean several things. In classical mathematics, a set can be specified through a definition within an accepted foundation. That does not require listing its members one at a time. A sequential enumeration cannot exhaust the reals, but inability to enumerate an object is not a contradiction in its definition.

Constructive mathematics introduces another meaning. It studies what existence and operations can be supported by explicit constructions and suitable approximation data. Identifying constructive mathematics with finitely many points would erase the constructive treatment of real numbers. The mathematical questions concern the permitted foundations and procedures; they are not settled by the physical limitations of a pen.

A circle specified by x²+y²=r² is exact within Euclidean geometry. A manufactured circular object must be assessed at some tolerance. The fact that manufacture does not establish infinite precision neither refutes the equation nor proves that the object inhabits discrete space. The same distinction applies to an ideal square and to the real numbers used in its description.

Physical realization adds empirical interpretation. A model needs an account of what its states, distances and observations correspond to. No physical assembly sequence follows merely from a set theoretic union. Conversely, a continuum need not have been assembled point by point in time in order to be proposed as physical. Imposing that requirement would add a premise about physical existence which itself needs defence.

## 7. The physical question survives the mathematical answer

The physical position advanced here rejects fundamental continua and proposes that continuum geometry represents very fine discrete relations. Some collective quantities may agree exactly; other predictions may agree only within a stated tolerance. On this account, smoothness belongs to an effective description, while the underlying structure supplies the interactions from which measured distances emerge. “Tiny slices” is one possible picture, but slices already assume an extended geometry; a more fundamental proposal must specify the relations that generate that geometry. This is a proposed physical direction, not a conclusion of the singleton calculation.

Successful continuum calculations establish their adequacy for particular observations. They do not, by that success alone, establish infinite physical divisibility. Nor does a finite measurement record prove discreteness. The task is to compare candidate structures, dynamics and predictions under declared conditions.

It is useful to separate three claims: a continuum is mathematically specified; a continuum model is empirically successful in a regime; a continuum is fundamental in nature. Each requires its own support. Calling the third claim an assumption can be justified where it has been adopted without a derivation. Calling it a logical contradiction requires an actual incompatible pair of premises. Singleton measure zero does not supply that pair.

Minimum measurable length and discrete spacetime are also different proposals. A limit on localization concerns operations available to physical observers; it does not by itself decide the underlying set of events. Hossenfelder (2013) surveys minimal length proposals and their implementation in quantum physics. Such work makes the physical question substantial without changing countable additivity.

### 7.1 Measure, geometry and the existing debate

The measure objection has a substantial history. Skyrms's “Zeno's Paradox of Measure” (1983) explicitly identifies Grünbaum's work on the standard resolution using measure theory and infinite sets. This paper does not claim that resolution as a discovery. Its question is what additional physical premises would justify moving from the mathematical resolution to an ontology of space. The bearer argument in Section 1.2 makes that extra commitment visible.

A relevant predecessor is Spinoza, [*Ethics* I.15, scholium](https://www.gutenberg.org/files/3800/3800-h/3800-h.htm). He rejects the attempt to compose a line from points, a surface from lines or a solid from surfaces. His conclusion is not discrete space: he defends extended substance as infinite, one and indivisible. This is not modern measure theory, but it already distinguishes extension from an additive assembly of dimensionless parts. It therefore cannot be cited as a proof of the discrete conclusion advanced here. The historical lesson is sharper: rejecting a sum of zero widths does not, by itself, choose between discrete relations and a continuous whole.

There is an independent challenge to discrete proposals. Weyl's tile argument asks whether counting elementary steps recovers the geometry measured at larger scales. Fritz's treatment (2013, corrected version) makes the difficulty precise for a class of periodic graphs: the attainable velocity set is a polytope and has distinguished directions. The scope matters. This constrains the specified graph dynamics; it is not a theorem excluding every discrete physical theory.

A simple control calculation shows why refining a grid is insufficient. On a square lattice of spacing h, allow only horizontal and vertical steps, each of length h. For a>0 with a/h an integer, the shortest path from (0,0) to (a,a) has length 2a. Its Euclidean distance is √2 a. Their ratio is √2 for every permitted h. Making the lattice arbitrarily fine does not remove that discrepancy. Assigning a Euclidean distance to the same sites changes the metric, but does not establish that a physical propagation rule measures that metric. A construction must connect its distances to its dynamics and observations.

Causal set theory pursues a different structure from a regular spatial lattice. In its continuum approximation, causal order is combined with counting related to volume. Poisson sprinkling avoids the regular lattice's built in preferred directions in a statistical Lorentz invariant construction. This is not a complete derivation of spacetime dynamics or proof that typical causal sets approximate manifolds; those remain substantive problems (Surya 2019). The proposal here must meet these recovery problems rather than cite causal sets as a universal answer to the tile argument.

Other discrete proposals make different commitments. Fredkin (2003) advocates digital mechanics based on computational rules. 't Hooft (2016) investigates a cellular automaton interpretation of quantum mechanics. Rovelli and Smolin (1995) derive discrete area and volume spectra in their quantum gravity framework. A discrete spectrum of geometric operators is not simply a lattice of physical locations. These approaches motivate different tests and cannot be treated as interchangeable proofs that continuous space is impossible.

### 7.2 What this paper contributes

The contribution is a structured comparison with three obligations. First, distinguish the measure of a set from sums of singleton measures. Second, state which observables an effective representation reproduces, exactly or within a bound. Third, test the geometry and dynamics claimed for a physical realization, including failures that persist as a mesh is refined. The centroid and second moment examples in Section 5 and the lattice calculation above make these obligations explicit.

The proposed advance is to make the physical commitment explicit and connect it to a constructive comparison. A discrete candidate earns its place by recovering the successful observables with fewer independent assumptions or by correctly predicting a departure. Exact agreement on selected observables can suffice for a representation without making its continuum coordinates fundamental constituents. The search in Section 8 turns that distinction into a programme for physics.

## 8. What would make discreteness a physical explanation?

A finite grid can approximate continuum measurements without proving that space is a grid. A physical alternative must name its events and relations, show how measured distance and direction arise from them, and recover the predictions for which continuous physics already works. It must also reproduce quantum interference. For two contributions of equal magnitude but opposite phase, their amplitudes cancel before their probability is calculated; merely counting positive branches would not capture that result.

The test is to publish the rules and initial conditions, derive the effective geometry and quantum behavior, and check for errors or preferred directions that experiments could detect. A rule that produces an attractive graph but fails these tests is not a replacement for continuum physics. This paper identifies the mathematical and physical gap; it does not select a unique discrete theory.

## 9. Implications for mathematical monism and teaching

For Matheism, the example sharpens what “structure” must include. A list of points together with their individual widths omits distances, measure and, in a physical application, dynamics and observation relations. Complete structure cannot mean a pile of isolated constituents. This clarifies the ontological proposal without proving it from a theorem about length. Both discrete and continuous structures are mathematical candidates.

A useful teaching exercise asks students to explain why four answers differ: the length of each singleton is zero; the unit interval has length one; all finite grids have length zero; the normalized atomic measures have mass one. The explanation must name the object, the measure and the operation in each case. An additional exercise asks why the real line and the rational numbers can both be ordered without having immediate successors. A final exercise asks which conclusion about physical space follows from those calculations. None selects a fundamental physical ontology without additional premises.

An educational study could compare procedural integration scores with performance on these explanation tasks, using a predefined rubric and a specified student population. Until such evidence exists, the defensible claim is that the distinction deserves explicit teaching. The mathematics resolves the apparent contradiction. The physical interpretation remains a question for a theory that earns its predictions.

### Observer recovery

A complete candidate must recover the organization needed for observers like us. If its boundary lacks that organization, a history containing those observers must develop it. A short rule can generate an elaborate observer without becoming a long rule; output size and runtime are distinct from description length. The search must account for the starting state and the observation map as well as the update law.

## Conclusion

Zero point-lengths do not make standard continuum mathematics contradictory: a union of uncountably many points is not a countable sum of their lengths. The physical explanation of extension remains a separate question.

This paper favors a discrete construction and states the extra generation premise that motivates it. The examples show why matching one average or shrinking a grid is not enough to recover all geometry. A successful physical account must derive the relevant distances and dynamics and demonstrate agreement with observations.

## References

Skyrms, B. 1983. “Zeno's Paradox of Measure.” In R. S. Cohen and L. Laudan, eds., *Physics, Philosophy and Psychoanalysis*, pp. 223 to 254. [Publisher abstract and bibliography](https://doi.org/10.1007/978-94-009-7055-7_12).

Fritz, T. 2013. “Velocity Polytopes of Periodic Graphs and a No-Go Theorem for Digital Physics.” *Discrete Mathematics* 313, 1289 to 1301. [Corrected arXiv version](https://arxiv.org/html/1109.1963v3). [Corrigendum](https://doi.org/10.1016/j.disc.2013.07.011). See also the author's [Weyl argument presentation](https://tobiasfritz.science/2012/weyl_argument.pdf).

Surya, S. 2019. “The causal set approach to quantum gravity.” *Living Reviews in Relativity* 22, 5. [Review](https://arxiv.org/abs/1903.11544).

Fredkin, E. 2003. “An Introduction to Digital Philosophy.” *International Journal of Theoretical Physics* 42, 189 to 247. [Article](https://doi.org/10.1023/A:1024443232206).

't Hooft, G. 2016. *The Cellular Automaton Interpretation of Quantum Mechanics*. [Author manuscript](https://arxiv.org/abs/1405.1548).

Rovelli, C., and L. Smolin. 1995. “Discreteness of area and volume in quantum gravity.” [Paper](https://arxiv.org/abs/gr-qc/9411005).

Gorard, J. 2020, revised 2021. “Some Relativistic and Gravitational Properties of the Wolfram Model.” *Complex Systems* 29(2), 599 to 654. [Paper](https://arxiv.org/abs/2004.14810).

Gorard, J., and J. Dannemann-Freitag. 2023. “Axiomatic Quantum Field Theory in Discrete Spacetime via Multiway Causal Structure: The Case of Entanglement Entropies.” [Paper](https://arxiv.org/abs/2301.12455).

Wolfram, S. 2020. “A Class of Models with the Potential to Represent Fundamental Physics.” [Paper](https://arxiv.org/abs/2004.08210). See also the project's proposed [quantum formalism](https://www.wolframphysics.org/technical-introduction/potential-relation-to-physics/quantum-formalism/).

Feynman, R. P., R. B. Leighton and M. Sands. *The Feynman Lectures on Physics*, Volume III, Chapter 16, “The Dependence of Amplitudes on Position.” [Online edition](https://www.feynmanlectures.caltech.edu/III_16.html).

Hossenfelder, S. 2013. “Minimal Length Scale Scenarios for Quantum Gravity.” *Living Reviews in Relativity* 16, 2. [Article](https://arxiv.org/abs/1203.6191).

Hunter, J. K. 2011. *Measure Theory*. University of California, Davis. [Lecture notes](https://www.math.ucdavis.edu/~hunter/measure_theory/measure_notes.pdf).

Lebl, J. 2026. *Basic Analysis I: Introduction to Real Analysis, Volume I*. Version 6.3, May 15. [Text](https://www.jirka.org/ra/realanal.pdf).

Sorkin, R. D. 2003. “Causal Sets: Discrete Gravity.” Notes for the Valdivia Summer School. [Paper](https://arxiv.org/abs/gr-qc/0309009).

## Generative AI assistance

ChatGPT assisted with drafting, source checking and elementary calculations. The proposed pedagogical study has not been conducted. This draft requires author review before submission. AI systems are not authors.
