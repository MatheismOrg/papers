# The Emperor's New Points: How a Line Gets Its Length

*How a line gets its length, and what that tells us about space*

Gary Abraham Bernstein

## Abstract

A line has length; its points have none. The apparent contradiction comes from treating a union of points as a sum of their lengths. This essay explains the distinction with small examples, then asks what it can and cannot tell us about physical space. A finer grid can approximate a continuum, but its geometry and physical interpretation still need to be supplied.

## A line and its points

A line has length. Each of its points has none. If we keep adding zeros, where does the length arrive?

“Infinity” sounds like an answer until we ask it to do the arithmetic. Ten zeros sum to zero. A million do too. An infinite series of zeros does not improve matters. Yet the unit interval has length one.

The surprise is real. The contradiction comes from asking a sum to do a union's job.

Put every singleton {x}, for x between zero and one, into a union. We obtain the unit interval. This tells us which points belong to the set. It does not tell us that the length of the union is the sum of the lengths of all its singleton subsets. That further rule would need justification.

Ordinary length, formalized as Lebesgue measure, is countably additive. For a sequence of disjoint measurable pieces, add their lengths to get the length of their union. The individual points of the real interval cannot be arranged into such a sequence. They form an uncountable set.

We can still define a nonnegative sum over that set: take the supremum of all its finite subsums. In this case every finite subsum is zero, and the result is zero. It is simply not the length of the union. No zero acquired a hidden positive contribution.

The distinction is old. Skyrms's discussion of Zeno's paradox of measure points to Grünbaum's treatment of its resolution through measure theory and infinite sets. What deserves renewed attention is how easily the distinction disappears when we move from calculating with a continuum to talking about what exists.

Spinoza had already rejected making a line from points in [*Ethics* I.15](https://www.gutenberg.org/files/3800/3800-h/3800-h.htm). He defended indivisible extension, not a world made of tiny grains. His answer was different from the modern one, but the warning survives: rejecting the sum does not settle what space is.

## What the calculation contains

A list of points and their individual lengths is incomplete information. Give the unit interval twice its ordinary measure. Every singleton still has measure zero, but the whole interval now has measure two. Singleton values alone cannot tell us which measure we chose.

In ordinary geometry, distances and a normalization supply the missing scale. The whole mathematical structure contains more than membership. Removing the relations and then asking where the length went is rather like removing the multiplication rule and asking what the numerals can calculate alone.

Order is not enough either. It tells us which position is before another. It does not establish how far apart they are. A continuum model includes the structure needed to answer that further question.

An integral does not evade this distinction. Divide the unit interval into n cells, each of width 1/n. Add their widths. At every finite stage the answer is one, and every width is positive. There is no final partition where the widths have become zero and the addition continues as before.

For a continuous function, a Riemann sum samples its value at a point in each cell and multiplies that value by the cell's width. The point chooses a value. The cell supplies a width. Taking the limit of the whole sums is not replacing every term with zero while keeping the total unchanged.

A student might learn to execute that calculation without being asked to explain these different jobs. That is a reason to ask the question explicitly. It is not evidence about how many students misunderstand it; that would require a study.

## More dots, still no length

Mark the endpoints of a unit interval. Then its midpoint. Then the quarter points, and continue by halving every gap. At each finite stage there are only finitely many marked points. Each finite set has length zero.

The union of all these stages is countable, so its length is also zero. Nevertheless, its points lie arbitrarily close to every point of the interval. Taking the closure produces the entire interval, with length one.

A union and a closure have different results. The new length did not emerge from summing the old singleton lengths. We changed the set through a different operation.

There is even a precise sense in which the finite grids approach the interval. The greatest distance from any position in the interval to the nearest grid point is half the grid spacing. That distance tends to zero as the grid is refined. Geometric closeness of this kind does not force the sets' measures to approach one another.

This makes a useful classroom question: how can the grids get arbitrarily close to a line while all having zero length? The answer names what is converging and what is being measured. “They become the line at infinity” obscures both.

## Exact answers from a different representation

Approximation is often the right word for an effective description. Representation is broader.

Suppose a model permits only two particle positions, minus one and plus one, with equal probability. The mean position is exactly zero. Zero need not be a third permitted position. It expresses a relation between the two alternatives.

Now compare a model that gives uniform probability to every position in the interval between minus one and plus one. Its mean is also zero. For this observable, the agreement is exact.

Square the position before taking the average and the models disagree. In the two position model the result is one. In the uniform model it is one third. Matching the mean does not match the distribution.

A smooth model can therefore give an exact collective answer without establishing that all its individual coordinates are physically available. Conversely, observing an exact average does not establish a discrete underlying world. We must specify which result matches and what the experiment measures.

## Smaller tiles do not settle geometry

Replacing a smooth line with grains brings another question: how do the grains determine distance?

Imagine a square grid. Permit motion only horizontally and vertically. From the origin to a point one unit right and one unit up, every shortest route travels two units. The Euclidean diagonal is √2 units long. Making the grid finer adds more shorter steps, but the total route length remains two. The relative discrepancy does not shrink.

This elementary example illustrates the kind of obstruction raised by Weyl's tile argument. Fritz gives a more general result for specified periodic graph dynamics, in which attainable velocities retain distinguished directions. The conditions of that result matter. It is not a proof against every discrete theory.

We could assign Euclidean distances to the grid points. But then we must explain how the allowed dynamics produces measurements of those distances. The mathematical assignment and the physical recovery are different tasks.

“Grainy” does not finish an explanation. A grain already has a size. A collection of sites still needs relations, a scale and dynamics. Continuous models need a physical interpretation of their metric structure as well. The demand applies to both.

## Conclusion: What remains physical

My position is that physical space has a discrete underlying account. The zero widths do not prove it. They explain why adding singleton lengths is the wrong construction and why a physical proposal must say more.

A continuum defender can place extension in metric relations among points. A discrete proposal can also place distances in relations. The comparison must ask what each structure explains, which assumptions it requires, and where their predictions differ.

A continuum also has no next point. Between any two distinct real positions lies another. This does not destroy their order, and it does not invalidate the mathematical trajectory x(t)=t. The trajectory assigns positions to times without enumerating immediate successors. Whether nature instantiates such a history is a further question.

The rationals make the distinction clear. A program can list integer ratios so that every rational eventually appears. There is a next output in the list, but no next rational in numerical order: another rational always fits between two distinct ones. A discrete generating process can therefore describe a dense ordering. All real numbers cannot be listed this way. Whether physical existence requires such a listing is the physical premise to defend.

We can investigate discrete models while keeping continuum mathematics as a working instrument. Some calculations may agree exactly with collective observables. Others may approximate them within a measured tolerance. A physical claim should identify which kind of agreement it has earned.

The useful lesson is precise. A union is not a sum. A limit of finite grids is not necessarily their union. An exact average is not an extra occupied position. And a finer lattice need not give the right geometry.

Calculus answers the questions posed to its structures. Physics must also account for the structure being measured. An exact answer for an average can come from a representation whose intermediate positions are never occupied. Understanding that possibility gives the search for a discrete account of space a clearer starting point.

## Further reading

John K. Hunter, *Measure Theory* (2011), Chapters 1 and 2. [Notes](https://www.math.ucdavis.edu/~hunter/measure_theory/measure_notes.pdf).

Jiří Lebl, *Basic Analysis I*, Chapter 5. [Text](https://www.jirka.org/ra/realanal.pdf).

Brian Skyrms, “Zeno's Paradox of Measure” (1983), pp. 223 to 254. [Publisher record and abstract](https://doi.org/10.1007/978-94-009-7055-7_12).

Tobias Fritz, “Velocity Polytopes of Periodic Graphs and a No-Go Theorem for Digital Physics” (2013). [Corrected text](https://arxiv.org/html/1109.1963v3).

## Author and assistance

Gary Abraham Bernstein is an independent researcher. ORCID: 0009-0009-1761-2867. ChatGPT assisted with drafting and source checking. The examples are standard mathematical constructions, not new experimental findings.
