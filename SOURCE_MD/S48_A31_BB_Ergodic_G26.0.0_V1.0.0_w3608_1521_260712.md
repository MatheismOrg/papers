# Boltzmann Brains Are Not Guaranteed by Infinite Time

## Abstract

Infinite time does not make every observer inevitable, because an endless sequence need not contain every possible event: the even numbers go on forever without producing an odd number. The same holds for Boltzmann brains, observers formed by rare fluctuations instead of biological development. A cosmology has to show that such a system can form, that it can work as an observer, and how much weight its observations receive. Infinite time settles none of these, and determinism alone does not exclude fluctuation observers either. Exact examples show where duration matters: a steady formation rate makes formation certain eventually, a rate that falls fast enough does not, a falling rate can still give infinitely many events in an expanding volume, and a steady rate can lose the race to vacuum decay. For observers assembled from matter that thins out as space expands, the expansion itself makes the rate fall fast enough: each doubling of distances makes each needed particle eight times rarer while the counted volume grows only eightfold, so in any region that expands with the universe, anything built from two or more particles has a finite expected count over endless time, as Zeno's halves, 1/2 + 1/4 + 1/8 and so on, sum to one. In an infinite space, a measure that compares the infinite totals region by region therefore gives such observers a finite ratio to ordinary ones. The paper separates the three requirements and shows that the Boltzmann-brain threat is a claim about a particular cosmology, to be checked inside it, and not a consequence of endless time.

## 1. Introduction

Infinite time alone cannot make a Boltzmann brain inevitable, because whether one forms, functions as an observer and dominates our evidence depends on the cosmology's laws, its criterion for observers and its weighting rule. A Boltzmann brain is a proposed observer assembled by a rare fluctuation rather than by biological development; the term recalls the physicist Ludwig Boltzmann's suggestion that our ordered world might itself be such a fluctuation (Boltzmann, 1895). Some cosmologies give these observers enough weight to threaten the reliability of the records we reason from, and infinite time is one input among those that decide the matter (Carroll, 2017).

The modern problem is sharpened by a positive cosmological constant, which leaves an expanding universe at late times with a fixed horizon and a finite entropy. Lisa Dyson, Matthew Kleban and Leonard Susskind (2002) examine what follows from treating such a de Sitter space (a universe expanding at a constant rate) as a finite-entropy system, including the difficulty of explaining our ordered observations through fluctuations. Andreas Albrecht and Lorenzo Sorbo (2004) compare fluctuated and inflationary origins. These arguments work within particular physical and statistical frameworks; infinite duration alone does not deliver them.

Existing responses include vacuum decay (empty space dropping to a lower-energy state), stationary vacuum accounts (in which empty space sits in a changeless quantum state), models with asymptotic freezing (late-time evolution slowing to a standstill), and different cosmological measures (rules for weighting observations across an infinite universe). Don Page (2008) derives a lower bound on the rate of vacuum decay that keeps the expected surviving four-volume (volume of space multiplied by time) of each comoving region (a region that grows with the expansion) finite, under his assumptions. Raphael Bousso, Ben Freivogel and I-Sheng Yang (2009) and Andrea De Simone and colleagues (2010) study how observers are weighted under the scale-factor measure, which counts events up to a cutoff in how much space has expanded, and related measures. Kimberly Boddy, Sean Carroll and Jason Pollack (2015) distinguish a projection (the part of a quantum state that matches a given configuration) associated with a putative observer from a process that forms and operates one. Roderich Tumulka (2018) studies freezing in a specified Bohmian field model. The present discussion builds on that work.

The useful question is narrower than whether every possible microstate (an exact arrangement of every part) appears. Does the specified cosmology realize histories with the organization attributed to a functioning observer, and with what weight?

This paper's contribution is to separate the three requirements a dominance claim needs, that the observer forms, that it functions as one and how much its observations weigh, and to show with exact cases where infinite time helps and where it does not: a steady formation rate makes formation certain eventually, while a rate that dies away fast enough leaves the total chance below one forever. Two further exact examples show that a falling local rate can still give infinitely many expected events in an expanding volume, and that a steady rate can fail to give any when the region can decay first. A dilution estimate then shows that for assembly from matter that thins out with the expansion, the local rate falls fast enough: the expected count over endless time is finite. Its decisive premise is that each requirement has to be met inside a specified cosmology. It rules out the slogan that infinite time makes every observer inevitable, and it rules out dismissing fluctuation observers by determinism alone.

## 2. Three requirements behind a dominance claim

A fluctuation observer argument can be organized as follows. Each step requires a model; none follows from the word *infinite*.

(P1) The cosmology supplies an appropriate duration and volume, together with a law for transitions or formation events.

(P2) Under that law and its boundary conditions (the conditions fixed at the start or edges), a specified class of candidate observer histories has positive occurrence weight or rate.

(P3) Those histories meet the adopted criterion for experience and records.

(C) Under a specified weighting and conditioning rule (a rule restricting probabilities to cases matching the evidence), their observations dominate those of observers with ordinary developmental histories.

The conclusion requires comparing observer weights, not merely establishing a nonzero rate. If both candidate populations are infinite, a regulator (a cutoff that makes infinite totals comparable) or other measure prescription is needed.

Under ergodicity, a typical history visits a region as often as its measure predicts, and that frequency says nothing about whether anything in the region functions as an observer. Stated precisely, a probability-preserving map T on (X, mu) is ergodic when every invariant measurable set (one the map carries onto itself) has measure zero or one. For a measurable target B, George Birkhoff's ergodic theorem then gives visit frequency mu(B) for almost every starting point (all except a set carrying zero weight) (Green, 2015).

If mu(B)>0, this establishes repeated visits for typical starts. It does not establish that every exact point is visited, or that the states in B constitute functioning observers.

## 3. Determinism does not settle reachability

### 3.1. Unitary evolution and branch structure

In Hugh Everett's many-worlds formulation (Everett, 1957), where every outcome of a measurement occurs in its own branch, the quantum state evolves unitarily: the evolution is reversible and conserves total probability. Which states it reaches depends on the Hamiltonian (the energy operator that drives the evolution) and the initial condition.

Unitary evolution can turn a zero coefficient into a nonzero one, so a configuration with no weight at the start can carry all of it later. In symbols, in units with hbar = 1, take H = omega sigma_x and the initial state |0>. Then |psi(t)> = cos(omega t)|0> - i sin(omega t)|1>. The probability of |1> begins at zero and reaches one at t = pi/(2 omega). A direction stays empty only if it lies in an invariant subspace (a family of states that the evolution maps among themselves) that started with zero weight. Starting at zero amplitude (the number whose squared size gives a probability) is not enough on its own.

A nonzero projection at an instant does not by itself give an operating observer. The physical claim is about the organization and the evolution attributed to the projected system.

### 3.2. Bohmian mechanics

Tumulka's concrete cosmological model uses a scalar field (a field with a single number at each point) in Bohmian mechanics, David Bohm's version of quantum theory in which definite configurations are guided by the wavefunction, on a de Sitter background, a universe emptied out and expanding forever under a cosmological constant. Its late evolution becomes too simple to support the functioning attributed to a brain. The argument concerns field dynamics and a specified class of states. It should not be summarized as a theorem that every configuration in every Bohmian cosmology freezes.

### 3.3. Deterministic examples

Determinism by itself settles nothing here: two deterministic rules on the same states differ in whether every state gets visited. One rule steps around a ring and reaches everything; the other leaves each state where it is and reaches nothing new. Written out, take a finite counter: X = {0,...,N-1}, T(x) = x+1 mod N, uniform measure. Every orbit (the sequence of states reached from one start) visits every state; the system is ergodic. On the same X, T(x) = x is not ergodic when N > 1. The rule decides the result.

A rule can visit every part of a space with exactly the right frequency and still miss almost every individual point. Stated precisely, an irrational rotation x -> x + alpha mod 1, stepping around a circle by an irrational fraction of a turn, is ergodic for uniform measure, although a single orbit is countable (its points can be listed one by one) and cannot contain every point of the circle. Being dense (coming arbitrarily close to every point) and visiting with the right frequencies is not the same as touching every point of the continuum.

## 4. Infinite does not mean everything

A process generating events one after another needs a starting state, since without one no first event is produced. That does not guarantee an infinite future. The process can stop, continue forever or repeat.

The phrase "given infinite time, every possible configuration will be realized" trades on an ambiguity. The even integers are an infinite subset of the natural numbers containing no odd integer. An infinite decimal expansion of 1/7 = 0.142857... repeats forever without producing "999."

A formation process can make an event almost certain eventually with no appeal to ergodicity, provided its rate holds up over the whole interval. Stated precisely, suppose a proposed formation process is Poisson with constant rate lambda > 0: events arrive independently, at a steady average rate. The probability of no event by time T is exp(-lambda T), tending to zero. That model gives almost sure eventual formation without global ergodicity. Two assumptions carry the result: the rate stays positive, and the process holds over the interval.

A time-dependent rate can give a different answer. If the rate dies away fast enough, the total opportunity stays finite however long the process runs. Written out, for lambda(t) = lambda_0 exp(-t/tau), the integrated rate is lambda_0 tau, and the probability of at least one event is 1 - exp(-lambda_0 tau), below one for finite lambda_0 tau. Infinite time did not make the integrated opportunity infinite.

## 5. Mathematical randomness and description length

Matheism treats regularity and randomness as mathematical. A probability measure is a mathematical structure even when it assigns probability one to a class of individually noncomputable infinite sequences (sequences no finite program can generate). The ontology, the account of what exists, does not select a cosmological transition law by itself.

### 5.1. Not all random distributions produce all configurations

A distribution's support, the set of values it can produce, matters. A process confined to [0,1] cannot output 2. Small weight is different from zero weight. A nondegenerate Gaussian (a bell curve with nonzero spread) has positive density throughout its domain, although distant regions receive little probability.

A probability attached to a state is not yet a claim that an observer forms and operates. For a projector P and state |psi(t)>, the Born probability (quantum theory's standard probability for an outcome) is <psi(t)|P|psi(t)>. Connecting that number to the formation and operation of an observer takes a model of histories and records. A stochastic (chance-driven) formation rate, a quantum projection probability and an algorithmic prior (a weighting of possibilities by description length) are distinct objects.

A simplicity comparison belongs at the weighting stage, and it has to treat both candidates alike. A developed observer may have a compact generative description; a fluctuation observer may also be generated by a compact cosmological model. It is invalid to charge the first a generator while forcing the second to list every constituent independently. An established description length difference Delta L gives a factor 2^(Delta L) only under the stated code prior, in which each bit is a 50/50 choice (Li and Vitanyi, 2019).

### 5.2. The scope of a sufficient reason argument

A strong Principle of Sufficient Reason, the thesis that everything has an explanation, may motivate the search for deeper explanations of outcomes. Accepting that principle proves nothing about nonergodicity, and neither does rejecting irreducible chance. The two deterministic maps in Section 3.3 already settle that logical point.

### 5.3. Declining rates and growing spacetime volume

A rate must be paired with the region in which it is counted. Let Gamma(t) be a rate per physical volume and V(t) the volume used by the chosen measure. In an illustrative inhomogeneous Poisson model, in which events occur independently at a rate that may change with time, the expected count through time T is

$$\Lambda(T)=\int_0^T\Gamma(t)V(t)\,dt,$$

and the probability of at least one event is 1 - exp(-Lambda(T)). Take Gamma(t) = Gamma_0 exp(-alpha t) and V(t) = V_0 exp(3Ht), with positive constants and H fixed. The infinite-time integral is finite when alpha > 3H and diverges when alpha is less than or equal to 3H. A decreasing rate per unit volume can therefore coexist with infinitely many expected events in an expanding global region.

This example does not choose a cosmological measure. It identifies why the choice matters. A fixed causal patch (the region one observer can ever exchange signals with), an expanding spatial slice and a regulated global ensemble need not give the same count. A finite-hazard result cannot pass from one region to another unless the rate and the counting domain are transformed together. Nor does an expected count by itself determine what a typical observer sees; ordinary observers and the reference class still have to enter the comparison.

### 5.4. A finite lifetime competes with fluctuations

Vacuum decay supplies a second distinction. Consider a fixed patch with a constant fluctuation rate beta and an independent terminal decay rate kappa, both positive. With exponential waiting times, the probability of at least one fluctuation before decay is

$$\Pr(B\ \mathrm{before\ decay})=\int_0^\infty\beta e^{-(\beta+\kappa)t}\,dt=\frac{\beta}{\beta+\kappa}.$$

If fluctuations do not end the patch, their expected number before decay is beta/kappa. These are ordinary competing-hazard calculations, and they show exactly where persistence enters the claim. A constant positive local rate does not make an event almost certain when the region has an independently finite lifetime. In the stable-patch limit, as kappa tends to zero, the probability tends to one and the expected count diverges. The familiar threat returns under the assumptions that support it.

Page's cosmological calculation adds the geometry of bubble nucleation, the birth of a region of lower-energy vacuum, to this distinction. Let A be a constant nucleation rate per spacetime four-volume, and let H be the late-time de Sitter expansion rate. In units with c = 1, A has units of inverse time to the fourth power, whereas H has units of inverse time. If a nucleated bubble destroys its causal future, a point survives with probability exp(-A V_4), where V_4 is its relevant past four-volume. At late times the effective survival hazard is 4 pi A/(3H^3). Finite expected surviving four-volume per comoving region requires this hazard to exceed the rate 3H at which spatial volume grows, giving Page's bound (Page, 2008):

$$A>A_{\min}=\frac{9H^4}{4\pi}.$$

The comparison with 3H is therefore between two rates per time; a nucleation rate per four-volume and a rate per time cannot be compared directly. The bound controls the integral of surviving volume; it does not directly compare decay with the Boltzmann-brain production rate. Nor is it a universal lifetime bound independent of the cosmology and weighting prescription. The fixed-patch calculation above illustrates the competition and leaves Page's geometric calculation in place.

The two examples close opposite shortcuts. Declining local activity need not suppress a global count if volume grows faster. A persistent local rate need not guarantee an event if the patch can end. A Boltzmann-brain argument must say which of these circumstances its cosmology and measure provide. These calculations clarify the premises; they do not refute a detailed model that has already established them.

### 5.5. Dilution: assembly from expanding matter

Section 5.3 left open how fast the local rate falls. For one class of fluctuation observers, the expansion itself fixes it. Consider an observer that must be assembled from N pre-existing particles, gathered into a region of fixed physical volume v, in a universe whose distances grow with the scale factor a(t). The particles thin out: their number density falls as n(t) = n_0 a(t)^(-3). In a dilute, uncorrelated gas, the chance that a given region of volume v holds the N required particles at once is about (nv)^N/N! when nv is small, so the formation rate per physical volume falls at least as fast as a^(-3N). The particles' own motions slow as space expands, and a cooling gas makes thermal assembly harder still; both only strengthen the bound. The physical volume of a fixed comoving region, one carried along by the expansion, grows as a^3. The expected count therefore accumulates at a rate

$$\Gamma(t)V(t)\propto a(t)^{-3(N-1)}.$$

With late-time exponential expansion, a(t) = exp(Ht), this is the case of Section 5.3 with alpha = 3NH, and alpha > 3H whenever N >= 2, so the infinite-time integral converges. Counted in doublings of the scale factor, each doubling multiplies the rate by 8^(-(N-1)), and the total over all future doublings is a geometric series: 1 + 1/8 + 1/64 + ... = 8/7 times the first term for a two-particle assembly, and incomparably less for a brain of some 10^26 atoms. Endless time leaves the total finite, as Zeno's halves, 1/2 + 1/4 + 1/8 + ..., sum to 1.

The estimate has a definite scope. It covers assembly from matter that dilutes with the expansion. Gravitationally bound matter does not dilute: stars, their remnants and black holes keep their densities while the space between them grows. Their lifetimes are finite on the usual accounts, through proton decay if it occurs, and otherwise through the slow collapse Freeman Dyson estimated for matter that cannot decay, followed by the evaporation of black holes (Dyson, 1979). A finite lifetime bounds their contribution, as in Section 5.4, though the bound can be vast. Nor does dilution touch fluctuations of a de Sitter horizon, whose temperature stays fixed as space expands. There the question returns to Sections 5.4 and 6: whether the late-time vacuum decays first, and whether a stationary vacuum state produces dynamical events at all. On the reading in which nothing happens without a cause, a stationary state produces no events, the conclusion Boddy, Carroll and Pollack (2015) reach under their assumptions about the state.

Where counts still diverge, divergence settles nothing by itself: two infinite populations have to be compared under a stated regulator. In an infinite space every epoch holds an infinite expected number of fluctuation events, so sums over epochs are infinite and a measure has to compare them. Counting per comoving region, the dilution result applies region by region: each region's expected total of assembled observers is the convergent sum above, a fixed multiple of its first term, and each region forms a finite number of ordinary observers from the same matter, so the ratio of the two populations is finite and is set by the earliest, densest epochs. Other cutoffs weigh epochs differently; under any of them a falling local rate lowers the fluctuation observers' count compared with a steady one, though it need not make their share small.

## 6. Comparison with existing approaches

Changing the measure can change the dominance criterion. Bousso, Freivogel and Yang (2009) give a local formulation of the scale-factor measure in expanding regions and compare its vacuum-stability conditions with those of the causal-diamond measure, which counts only what one observer's causal region contains. De Simone and colleagues (2010) derive conditions for avoiding Boltzmann-brain domination under the scale-factor cutoff. These results depend on the regulator, vacuum transitions, reference population and conditioning as well as on the volume counted, and no single Poisson integral with a new V(t) reproduces them all.

Each existing response holds under its own model's assumptions, and none of them supplies a general exclusion. Boddy, Carroll and Pollack's stationary vacuum proposal depends on assumptions about the state and its larger Hilbert space, the space of possible quantum states. Their earlier analysis distinguishes finite closed de Sitter systems where recurrences (returns arbitrarily close to earlier states) remain possible. Their conditional conclusions supply no blanket exclusion under unitary dynamics (Boddy, Carroll and Pollack, 2014, 2015).

Tumulka's conclusion depends on the field model and its late dynamics. His paper explicitly distinguishes configurations shaped like brains from configurations capable of functioning as brains. The present account identifies those dependencies and supplies elementary counterexamples to invalid general inferences.

## 7. Objections and replies

**Does allowing stochastic laws undermine the structural account?** No. Both deterministic maps and probability measures have mathematical structure. The present argument tests formation and weighting under whichever law the model specifies.

**Even if the exact configuration is not visited, configurations "close enough" will be.** That is the right way to pose it, and it has to be made precise. Specify a class of histories robust enough to implement the relevant organization. Ergodicity can imply repeated visits to a positive-measure target; nonergodicity alone cannot show a particular target is avoided.

**Does ordinary thermodynamic behavior require visiting every exact microstate?** No. Equilibration (settling toward steady values) of selected observables (measurable quantities) is different from exhaustive visitation.

**In Many-Worlds, every branch exists, so every nonzero-amplitude configuration is realized.** A branch is an account of correlated records and dynamics, not an arbitrary nonzero coefficient. Once a model establishes suitable observer histories, compare their weights using the same observation rule used for ordinary observers.

## 8. Implications

Three concrete tasks for a proposed cosmological account.

First, define the target as a class of functioning histories and calculate reachability or formation rates under the proposed state and dynamics.

Second, separate formation from dominance. Integrate rates over the relevant duration and volume, compare the rate at which ordinary observers are produced, and disclose any regulator needed for infinite populations.

Third, condition on the evidence consistently. David Wolpert, Carlo Rovelli and Jordan Scharnhorst (2025) distinguish the roles of assumptions about memories, temporal conditioning and entropy (roughly, a measure of disorder). A model should state which records it treats as data and which as reliable traces of a past.

## 9. Conclusion

Infinite time alone does not guarantee a fluctuation observer, because an endless process can keep avoiding a whole class of states, and determinism alone does not exclude one, because a deterministic rule can visit every state. The model must specify which histories are possible, which systems function as observers and how their observations are weighted. It must also count in a stated region: a falling rate can still produce infinitely many events in an expanding volume, and a steady one can lose the race to decay. For assembly from matter that thins out as space expands, the expansion itself makes the total in each comoving region finite, so a region-by-region measure gives such observers a finite ratio to ordinary ones even in an infinite space; bound matter and horizon fluctuations need their own accounts.

The even integers go on forever without producing an odd integer. A cyclic counter visits every state. Which analogy applies to a universe is a question for its dynamics.

## References

Albrecht, A., and Sorbo, L. (2004). Can the universe afford inflation? *Physical Review D*, 70, 063528. https://arxiv.org/abs/hep-th/0405270

Boddy, K. K., Carroll, S. M., and Pollack, J. (2014). De Sitter Space Without Dynamical Quantum Fluctuations. arXiv:1405.0298. https://arxiv.org/abs/1405.0298

Boddy, K. K., Carroll, S. M., and Pollack, J. (2015). Why Boltzmann Brains Don't Fluctuate Into Existence From the De Sitter Vacuum. arXiv:1505.02780. https://arxiv.org/abs/1505.02780

Boltzmann, L. (1895). On certain questions of the theory of gases. *Nature*, 51, 413-415.

Bousso, R., Freivogel, B., and Yang, I.-S. (2009). Properties of the scale factor measure. *Physical Review D*, 79, 063513. https://doi.org/10.1103/PhysRevD.79.063513

Carroll, S. M. (2017). "Why Boltzmann Brains Are Bad." arXiv:1702.00850.

De Simone, A., Guth, A. H., Linde, A., Noorbala, M., Salem, M. P., and Vilenkin, A. (2010). Boltzmann brains and the scale-factor cutoff measure of the multiverse. *Physical Review D*, 82, 063520. https://doi.org/10.1103/PhysRevD.82.063520

Dyson, F. J. (1979). Time without end: Physics and biology in an open universe. *Reviews of Modern Physics*, 51, 447-460.

Dyson, L., Kleban, M., and Susskind, L. (2002). Disturbing implications of a cosmological constant. *Journal of High Energy Physics*, 2002(10), 011. https://arxiv.org/abs/hep-th/0208013

Everett, H. (1957). "Relative State" formulation of quantum mechanics. *Reviews of Modern Physics*, 29, 454-462.

Green, B. (2015). Ergodic Theory. Lecture notes, section 3.4. https://people.maths.ox.ac.uk/greenbj/papers/ergodic-2015.pdf

Li, M. and Vitanyi, P. (2019). *An Introduction to Kolmogorov Complexity and Its Applications*. 4th ed. Springer.

Page, D. N. (2008). Is our universe likely to decay within 20 billion years? *Physical Review D*, 78, 063535. https://arxiv.org/abs/hep-th/0610079

Tumulka, R. (2018). The Problem of Boltzmann Brains and How Bohmian Mechanics Helps Solve It. arXiv:1812.01909. https://arxiv.org/abs/1812.01909

Wolpert, D., Rovelli, C., and Scharnhorst, J. (2025). Disentangling Boltzmann Brains, the Time-Asymmetry of Memory, and the Second Law. Entropy 27(12), 1227. https://arxiv.org/abs/2507.10959
