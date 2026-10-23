# Quantum Mechanics Without Collapse: The Simplicity Argument

## Abstract

Many worlds is the simpler account of quantum mechanics, because it is quantum mechanics with nothing added: one evolution rule covers particles, measuring devices and observers. Collapse theories add a second law that interrupts that evolution, and Bohmian mechanics adds an actual particle configuration and a law to guide it. Each added law is extra description, and each bit of description halves an account's weight, so the extras have to buy something the single rule cannot, and no observation so far requires either. Collapse also asks for an outcome that nothing selects, which this paper counts as an effect without a cause; indeterminists accept that cost, and the comparison of description length stands without it. Many worlds has no such event. The one rule produces every branch, and what looks random is which branch an observer is in, because no branch, observer or record holds the information of the whole wave. Given time and memory, one unitary rule can even run every computable program step by step, minds included. The branches need no rules of their own, because the one rule generates them all; starting conditions and the probabilities of recorded outcomes are counted for every rival alike. What many worlds owes is the Born rule's weights. The decision-theoretic derivations of David Deutsch and David Wallace and the self-locating derivation of Charles Sebens and Sean Carroll supply them from stated assumptions, and the comparison here counts those assumptions as part of the account's cost.

Keywords: Many-Worlds Interpretation, Simplicity Dominance Principle, Kolmogorov complexity, algorithmic probability, quantum foundations, Schrödinger equation, Born rule, measurement problem, mathematical monism

---

## 1. Introduction

This paper argues that a measurement needs no law of its own. One rule, unitary evolution (the deterministic, reversible change of a quantum state over time), already applies to the measured system, the apparatus and the observer. A collapse theory, in which the quantum state abruptly reduces to one outcome, changes that evolution partway through. A Bohmian theory keeps it and adds a selected configuration of particles with a rule guiding it. Each complete account then has to show what it explains for the information it requires.

A shorter complete description receives greater prior weight (its weight before any evidence is considered), because each of its bits is a 50/50 choice that no rule makes. Its posterior weight, the weight after evidence is considered, also depends on its observational likelihood, the probability it gives to what is observed. Using the same weight as a measure of what exists follows from the rule-random dichotomy (*Why These Simple Laws?*, Section 3.3).

## 2. Reality as mathematical structure

On the view taken here, all existence is either ruled or random, and both are mathematical. This motivates the monist position: causal powers and observable distinctions are structural, so the physical account can be identified with its structure rather than adding an unspecified bearer.

Isomorphism (a one-to-one match preserving all structure) alone does not imply numerical identity (being one and the same thing): distinct copies can be isomorphic. The identity is the interpretive claim, not a consequence of using the word isomorphism. The present comparison is conditional on explicit predictive models. It can be carried out whether or not a reader accepts the further ontological identity, an identity claim about what really exists.

## 3. The Simplicity Dominance Principle

Which structures observers should expect to inhabit is settled by a measure, and the measure comes from the same dichotomy that makes reality mathematical. A structure's weight halves with each bit of its complete description, because each bit of a shortest description is random, an unruled 50/50 choice. Stated precisely:

P(S) proportional to 2^(-K(S))

Here K, the length in bits of the shortest program that outputs S, is defined relative to a fixed reference machine (the computer chosen to run the programs) and an effective description of the complete model, not merely the length of its displayed equation. This weighting is the prior, to be combined with observational likelihoods.

Three limits belong with it. Changing a universal description language (one in which any program can be written) changes K by at most an additive constant. That constant does not guarantee that a ranking of finitely many models survives the change. Prior preference alone cannot rescue a model with poor likelihood. The extension from inductive inference (reasoning from observations to predictions) to physical realization rests on the dichotomy; the formal results about universal priors support it without proving it (Solomonoff, 1964; Li and Vitányi, 2008).

## 4. Application to quantum mechanics

### 4.1 Comparing interpretations

| Account | Complete specification includes | Explanatory task |
|---|---|---|
| Everett | State space, Hamiltonian, boundary, subsystem structure and record weights | Explain stable records and their probabilities within unitary dynamics |
| Operational quantum account | Preparations, operations and measurement probabilities | Predict observations with a specified scope of application |
| Bohmian account | Wave dynamics, configuration, guidance and preparation measure | Account for the selected configuration and observed statistics |
| Objective collapse | State, modified dynamics and collapse parameters | Account for outcomes through a distinct dynamical process |

Everett treats unitary evolution as complete dynamics and interprets emergent branches (effectively separate parts of the quantum state, each recording one outcome) as realized (Everett, 1957; Wallace, 2012). A measurement correlates the system with the apparatus and the environment; decoherence, the spreading of quantum correlations into the surroundings, can make interference between records effectively inaccessible, though it does not itself select a unique outcome or ontology (Zurek, 2003). Decoherence picks out record states only given a split into subsystems, a cost the table charges to Everett; collapse and Bohmian accounts settle the same question by writing a preferred quantity, position, into their laws.

The two main alternatives both carry structure that unitary evolution does without. Objective collapse proposals, such as the Ghirardi-Rimini-Weber (GRW) and continuous spontaneous localization (CSL) models, specify modified dynamics with extra parameters (Bassi et al., 2013). David Bohm's mechanics retains a wavefunction, the mathematical description of a quantum state, and adds a configuration and guidance law (Bohm, 1952). In quantum equilibrium, where particle positions are distributed according to the squared wavefunction, it reproduces the usual nonrelativistic statistics (Dürr, Goldstein and Zanghì, 1992). The parsimony objection asks: if unitary evolution already accounts for records and their probabilities, what requires these additions? Under the simplicity prior, a shorter complete account wins when likelihoods are equal.

One finite random measurement result has finite description length. A short stochastic law, a law that sets probabilities for outcomes, can assign a likelihood to a long record without listing that record as part of the law. Everett must likewise explain the weights of observed records. If we charge one account for the exact history it selects and let the other off, we are comparing different tasks.

### 4.2 Bell correlations

John Bell's inequalities (Bell, 1964), bounds on correlations between distant measurements, constrain combinations of assumptions including local factorization (each outcome depending only on its own setting and the shared past) and measurement independence (measurement settings statistically independent of the measured systems' hidden properties). Observed violations do not refute every deterministic model or independently confirm mathematical monism. An entangled state is a joint mathematical description; calling it structure does not answer Bell's questions about locality and independence. The correlations do not allow controllable signaling. That operational fact is distinct from whether a proposed ontology is local in Bell's sense.

### 4.3 The Born rule

Why do observed frequencies match |psi|^2? David Deutsch (1999) and David Wallace (2012) use decision-theoretic assumptions (premises about how rational agents choose) in Everettian quantum mechanics. Charles Sebens and Sean Carroll (2018) use self-locating uncertainty (uncertainty about which branch one is in) and an explicit separability principle (that an observer's probabilities depend only on its local situation). Both routes need substantive premises beyond the Schrödinger equation, the law governing how a quantum state evolves.

A short code assigned using known Born weights, the outcome probabilities quantum theory predicts, would encode those weights, not independently derive them. No general theorem identifies a branch's Kolmogorov complexity (its shortest description length) with the negative logarithm of its squared amplitude (its Born weight).

#### A refinement test for record weights

Counting branches and weighting them by squared norm (each branch's Born weight) come apart as soon as a branch is described more finely. Take two records an observer could tell apart, each correlated with its own state of the environment. With orthogonal (perfectly distinguishable) records |0> and |1> correlated with orthogonal environmental states:

$$|\Psi\rangle=\frac{\sqrt{3}}2|0,E_0\rangle+\frac12|1,E_1\rangle.$$

The squared norms are 3/4 and 1/4. Now split the environment that accompanies record 0 into several equally weighted substates, leaving the record itself alone. Refine only the environment correlated with record 0 by an isometry (a map preserving lengths and angles) sending E_0 to an equal superposition (a weighted sum of states) of m orthogonal substates E_0j:

$$|\Psi_m\rangle=\sum_{j=1}^{m}\frac{\sqrt{3}}{2\sqrt{m}}|0,E_{0j}\rangle+\frac12|1,E_1\rangle.$$

Tracing out the environment (summing over its states to describe the record alone) gives the same record state, diag(3/4, 1/4), for every m. Counting displayed terms instead assigns record 0 the fraction m/(m+1), changing with m while the observed distinction has not. Squared-norm weights satisfy refinement invariance; naive branch counting does not. That is the reason the comparison needs a weighting rule invariant under refinement. Wallace (2009) makes the required rationality assumptions explicit.

## 5. Objections

**"MWI adds infinite branches, that is complex."** Branches are consequences, not axioms. Complexity is measured by description length. A short program can generate infinite output.

**"This is untestable metaphysics."** Specified collapse models can make different predictions from unitary evolution, and those can be tested. A simplicity prior need not assign zero probability to complex models.

**"Why algorithmic probability?"** Because each bit of a shortest description is an unruled 50/50 choice, and universal induction also dominates every computable model class up to a constant factor. The choice of reference machine and observation rule still matters.

**"The ontological question is misguided."** Pragmatist and QBist interpretations (Fuchs, Mermin and Schack, 2014), which treat quantum states as agents' tools or beliefs, pose a different question about quantum states. The present proposal compares explicitly predictive models.

## 6. Complete comparison

### 6.1 Terms of comparison

A complete comparison charges every account for everything it needs. Use one explanatory task and encoding. The score includes a model description and the information cost of data under its likelihood. Joint description can exploit shared information.

### 6.2 What each account includes

Collapse: include the stochastic dynamics and parameters. A finite observed bit string is not an infinite specification of the law. Everett: include the dynamics, boundary data (the starting conditions) and the account assigning record weights. Saying that all outcomes occur does not remove the need to explain observed frequencies. Superdeterminism: include the common generator, preparation measure (the weighting of possible starting states) and response rules. In the simplest option, one rule generates both the source and the setting devices, reproducing Bell correlations while keeping local response and definite outcomes, at the cost of measurement dependence.

### 6.3 Record weights and the Born rule

An Everettian comparison has to pay for two things explicitly. The Born rule's assumptions carry a cost that cannot be declared zero because an argument is labeled a derivation. Single-world accounts pay a matching cost, a premise tying chance to rational belief (Saunders, 2010). Everett derives that link from premises about rational choice, such as indifference to branching, and collapse accounts add a law besides. The interactions and initial data underlying stable records may compress through reusable structure, but they cannot be omitted from one candidate and counted against another.

### 6.4 Causal completeness and boundary conditions

Collapse selects an outcome with nothing selecting it. Probability theory describes such a process consistently, as a list of chances, and that is all it does: it says how often, not why this one. On the view taken here, an outcome that nothing produces is an effect without a cause, and an account that allows it leaves each outcome unexplained; indeterministic interpretations reject that premise, and the simplicity comparison does not rest on it. Many worlds has no such event. Every branch is produced by the one rule. What looks random from inside is which branch an observer is in. Even an observer who knew the whole wave could not tell, before looking, which copy it is (Sebens and Carroll, 2018), and no observer knows the whole wave, since no part of it holds the information of the whole. A branch, an observer or a record carries only its share, while the rest stays in correlations with the environment, as a displayed sum of 2 hides which pair was added while the machine that computed it keeps both addends. Nothing then requires collapse to account for information an observer cannot reach. Deterministic accounts must still constrain their starting states, and boundary conditions, outcome probabilities and observational predictions remain part of the comparison.

On the view taken here, the universe starts simple and has to grow its observers, and growth takes causal steps: each state is produced from the one before, and since a continuum has no next state, the steps are discrete, an a priori proof in the logician's sense. A history grown one step at a time needs a first state, though no last one, and that is also why time runs forward: an endless history has no last state from which to run it backward, and a reversed history must also specify the point it turns back from, each bit of which halves its weight again. Many worlds fits this without strain, since its one rule produces each state from the last and leaves no event uncaused. The comparison here needs only that records form away from a low-entropy (highly ordered) start.

### 6.5 One rule for every mind

The one rule is enough for minds as well as measurements. A single quantum system with unbounded memory can run every computable program: a schedule gives each program its turn, reversible circuits keep the working information that each step would otherwise erase, and one unitary operator advances an internal clock and performs the next step. On the computational view of mind (that minds are computations) taken here, every computable mind is then present as a running computation, each state produced from the last, which is what separates a mind from a recording of one. Many worlds needs no new law for minds, any more than it needs one for measurements. Bohmian and collapse dynamics can run the same programs, so the point counts against accounts that give minds a law of their own.

## 7. Conclusion

Quantum mechanics without collapse is the simpler account, because it covers systems, apparatus and observers with one evolution rule. Many outcomes do not require a separate law for each outcome, and no outcome is left uncaused: the one rule produces every branch, and an observer who cannot predict an outcome is ignorant of which branch it is in, since no part of the wave holds the whole. The same kind of rule can run every computable mind. The simplicity argument becomes decisive only when complete descriptions and predicted probabilities are compared. Starting conditions count. The explanation of recorded outcomes counts. Under the declared weighting, where each bit is an unruled 50/50 choice that halves the weight, a shorter account with the same observational fit has the corresponding advantage. Everett's case rests there, and so does the standard its competitors can answer.

## References

Bassi, A., Lochan, K., Satin, S., Singh, T. P., and Ulbricht, H. (2013). Models of Wave-function Collapse, Underlying Theories, and Experimental Tests. Reviews of Modern Physics 85, 471-527. https://arxiv.org/abs/1204.4325

Bell, J. S. (1964). On the Einstein Podolsky Rosen paradox. Physics Physique Fizika, 1(3), 195-200.

Bohm, D. (1952). A Suggested Interpretation of the Quantum Theory in Terms of Hidden Variables. I. Physical Review 85, 166. https://doi.org/10.1103/PhysRev.85.166

Deutsch, D. (1999). Quantum theory of probability and decisions. Proceedings of the Royal Society A, 455, 3129-3137.

Dürr, D., Goldstein, S., and Zanghì, N. (1992). Quantum Equilibrium and the Origin of Absolute Uncertainty. Journal of Statistical Physics 67, 843-907. https://arxiv.org/abs/quant-ph/0308039

Everett, H. (1957). "Relative state" formulation of quantum mechanics. Reviews of Modern Physics, 29(3), 454-462.

Fuchs, C. A., Mermin, N. D., and Schack, R. (2014). An introduction to QBism with an application to the locality of quantum mechanics. American Journal of Physics, 82(8), 749-754.

Li, M., and Vitányi, P. (2008). An Introduction to Kolmogorov Complexity and Its Applications. Springer.

Saunders, S. (2010). Chance in the Everett interpretation. In S. Saunders, J. Barrett, A. Kent, and D. Wallace (Eds.), Many Worlds? Everett, Quantum Theory, and Reality. Oxford University Press.

Sebens, C. T., and Carroll, S. M. (2018). Self-locating uncertainty and the origin of probability in Everettian quantum mechanics. British Journal for the Philosophy of Science, 69(1), 25-74.

Solomonoff, R. J. (1964). A formal theory of inductive inference. Part I. Information and Control, 7(1), 1-22.

Wallace, D. (2009). A formal proof of the Born rule from decision-theoretic assumptions. https://arxiv.org/abs/0906.2718

Wallace, D. (2012). The Emergent Multiverse. Oxford University Press.

Zurek, W. H. (2003). Decoherence, einselection, and the quantum origins of the classical. Reviews of Modern Physics, 75(3), 715.
