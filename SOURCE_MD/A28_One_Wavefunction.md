# Can One Quantum Universe Run Every Computable Mind?

## Abstract

One mathematical quantum system can run every computable mind, because it executes every finite program through every stage it reaches, and a mind is its complete structure at work. The paper constructs the system. A schedule gives every program its turn, reversible circuits keep their working information, and one unitary operator advances an internal clock and performs the next stage; its inverse proves unitarity. Each step of every computation it runs is reproduced, so each computable mind is present in it as an ongoing computation, built step by step like any mind. The construction is many worlds at work: one unitary rule, no collapse, and every program's history held in one wave, with no branch holding the whole. Each mind in it is running, every state produced from the last, which is what separates a mind from a recording of one, and since a structure is one wherever it occurs, the bare form of being a subject is the same in all of them. Matheism's claim that every consistent structure exists is therefore no extravagance at the level of minds: one short quantum rule already runs every computable one, given unbounded time and memory. Whether our universe supplies those is a separate question, and computable structures are only part of mathematics.

Keywords: Many-Worlds Interpretation, Mathematical Universe Hypothesis, Tegmark Level IV, path integral, universal computation, Schmidhuber, Gödel encoding, Simplicity Dominance Principle, algorithmic probability, mathematical monism

## 1. The Setup

An ideal computer can run many programs by giving each a turn. With unbounded time and resources, a suitable schedule eventually reaches every finite step that each program executes. A first step does not require a beginningless past, and it does not guarantee an endless future. The question is whether a quantum system can realize the schedule, and which of the computations it runs are minds.

Matheism places minds and the physics implementing them inside one mathematical ontology, an account of what exists, so both fall under one comparison. Universal computation, one machine's ability to run any program, reaches only computable processes (those a finite program can carry out): it does not enumerate every mathematical structure, and it does not establish that every scheduled process feels.

A path integral adds up one contribution from every path to get a single transition amplitude (a number whose squared size gives the transition's probability); Richard Feynman (1948) showed that this matches ordinary quantum mechanics. Neither each integration path nor each Feynman diagram, a sketch standing for one term in the calculation, is automatically a literal particle trajectory or an Everettian world. Taking the amplitude structure as real is a proposed ontological reading; assigning every term its own traversing particle requires an additional argument.

The physicist Max Tegmark's Mathematical Universe Hypothesis (Tegmark, 1998, 2014) goes further: every consistent (contradiction-free) mathematical structure exists as a physical universe. In his four-level scheme that ensemble is Level IV, which includes structures with different constants, gauge groups (the symmetries that define a theory's forces), dimensions and mathematics. Level III is the branching of one quantum state in Hugh Everett's many-worlds interpretation, MWI (Everett, 1957; Wallace, 2012). The question this paper raises: is Level III already sufficient for everything we would want to call existence? Or does Level IV add something Level III cannot reach?

## 2. Why Level III Might Be Sufficient

Three considerations push toward Level III sufficiency.

**Code weights.** Simplicity can be priced once a description language is chosen: each bit is a 50/50 choice that no rule makes, so a shorter candidate gets more weight, by a factor of two for every bit it saves. Fix that language and a class of complete candidate models. Stated precisely, with assigned code weights w(H) proportional to 2^(-L(H)), a real difference of Delta L bits gives a prior ratio 2^(Delta L), the ratio of weights before any evidence. Adding up every program that produces an output, as Ray Solomonoff's algorithmic probability does (Solomonoff, 1964), gives a different number from weighting only its shortest program. Neither construction automatically defines a normalized measure (a weighting whose total is one) over every mathematical structure.

**Dovetailing.** A universal machine, Alan Turing's computer that can run any program it is given (Turing, 1936), can dovetail: run program 1 for a step, then programs 1 and 2 for a step, then 1, 2 and 3, so each program keeps getting turns. Every finite stage of each total computable history can then eventually be produced. Jürgen Schmidhuber built his proposal that one short program computes every computable universe on this observation (Schmidhuber, 1997, 2000). Dovetailing does not identify all total programs, finish every infinite history or show that our universe executes the enumeration.

**Arithmetic encoding.** Kurt Gödel's numbering (Gödel, 1931) assigns a number to each finite expression and proof, so claims about expressions become claims about numbers. Yuri Matiyasevich went further: every set a computer can list is the solution set of some whole-number polynomial equation (Matiyasevich, 1970), so plain arithmetic can describe any computation. It does not make arithmetic prove every encoded truth, implement every encoded system or represent every arbitrary infinite structure by a finite program. The forest can be described inside a tree without being physically grown there.

### 2.1 One quantum system executing every finite program

A single schedule can give every program as many instructions as anyone asks for, so long as the number asked for is finite. List all finite programs, including finite input data, as $p_1,p_2,\ldots$. Keep a separate working state for each. At stage s, advance programs 1 through s by one simulated instruction each. Leave a halted program idle. Each program that has not halted then picks up one more instruction at every later stage. Stated precisely, after stage S, every nonhalting program j with $j\leq S$ has received

$$N_j(S)=S-j+1$$

instructions. For any requested finite n, stage $S=j+n-1$ suffices.

**Make the stages reversible.** A step is reversible when running it backwards recovers exactly what went in. Retain the information that each simulated instruction would erase, together with the scheduler state, and every stage becomes reversible. Charles Bennett's reversible-computation construction supplies this general method (Bennett, 1973). The price is additional time and storage. The simulated computation itself is untouched.

Each stage can then be written as one circuit that runs backwards as readily as forwards. For each stage $t\geq0$, let $V_t$ be a finite reversible circuit, a fixed list of reversible gates (elementary logic operations), implementing stage $t+1$ and its record keeping. The circuits can be generated effectively (by a mechanical step-by-step procedure) from the scheduling rule. They act on a countable workspace (one whose cells can be numbered 1, 2, 3, and so on) with finite support at every finite stage: only finitely many cells are nonblank. Each circuit has an inverse obtained by reversing its gate order and inverting its gates. For $t<0$, let $V_t$ be the identity, the operation that changes nothing, so the operator defined below acts at every integer clock value.

**Make one autonomous unitary.** A unitary operator is a one-to-one map on states that preserves inner products (the overlaps between pairs of states), which is what quantum evolution requires. Carrying a clock label alongside the workspace lets one fixed operator do the work at every stage: it reads the clock, applies the circuit belonging to that stage, and advances the clock by one. In symbols, let C be the countable set of finite-support workspace configurations, include an integer clock register, take orthonormal basis states (perfectly distinguishable states that combine to give every state) $|t,c\rangle$ for $t\in\mathbb Z$ and $c\in C$, and define

$$U|t,c\rangle=|t+1,V_t c\rangle.$$

Every output has a unique predecessor, so the operator can be run backwards as well as forwards. Formally,

$$U^{-1}|t,c\rangle=|t-1,V_{t-1}^{-1}c\rangle.$$

U is a bijection (a perfect one-to-one pairing) of an orthonormal basis. Its linear extension (the same map applied to every combination of basis states) preserves inner products and is onto (reaches every state), so it is unitary on $\ell^2(\mathbb Z\times C)$. Applying it repeatedly runs the schedule forward one stage at a time. Written out, starting from $|0,c_0\rangle$ gives

$$U^S|0,c_0\rangle=|S,V_{S-1}\cdots V_0c_0\rangle.$$

This one fixed unitary and initial state execute every finite program through every finite stage it reaches. The clock labels describe the order of internal transitions; no outside clock is invoked.

**What follows for minds.** A mind is its complete structure, so a computation with a mind's organization is that mind, and the reversible simulation reproduces each of its steps while the extra bookkeeping leaves them untouched. Every computable mind with a finite program and input is therefore an ongoing computation inside this system. Interaction with an environment must be supplied as part of the computable program. Uncomputable input streams are outside this result.

The construction answers an ideal existence question. It does not show that our universe has the required initial state, resources or dynamics. The unitary has not been derived from a spatially local Hamiltonian, the energy operator that drives quantum evolution, or from gravity. It supplies no observer weighting and no reduction of noncomputable structures to computation.

## 3. What Would Survive Under Level III Alone

One wavefunction, the mathematical description of a quantum state, supplies the states reachable from its initial condition under its dynamics. A state can be consistent with the laws and still be unreachable from that condition. Combining the Standard Model, the established theory of particles and forces apart from gravity, with general relativity into a complete quantum cosmological dynamics is part of the problem being studied, not a supplied Hamiltonian in this paper.

Candidate phenomena include different measurement records, cosmological developments and physical implementations of minds. Each requires support in the specified evolution. Boltzmann fluctuations, brain-like structures thrown together by rare chance events, depend on the cosmology, its state and its dynamics. A small verbal addition to a description supplies no probability by itself. A false memory can represent an extended life without the represented events having occurred; a theory attributing experience to a brief fluctuation must specify the relevant organization and duration.

## 4. What Would Not Survive

Within a specified universal dynamics, every branch shares the governing theory. Branches can look different because symmetry-breaking fields (fields that fall into one of several equivalent lowest-energy states) settled differently, which changes the effective parameters, the constants that hold at accessible energies. Changing effective parameters is a smaller move than replacing the fundamental theory.

Structures requiring different fundamental dynamics are not automatically branches of the same state. These would require either different effective physics within inflationary regions (Level II, string landscape) or qualitatively different mathematical structures (Level IV).

The Schmidhuber/Gödel argument complicates this. If our universe's closure is computationally universal, it can run the dynamics of any other consistent computable structure. Does running it amount to ontological existence of the simulated structure, or does that structure also exist separately?

## 5. A Further Constraint: Decidability of Consistency

For effectively axiomatized theories (theories whose axioms a program can list) in suitable expressive settings, there is no algorithm deciding consistency in every case. Consistency is a property of a set of axioms relative to its logic. An actual mathematical structure satisfies its own truths. Any computational version of the question has to supply an effective input, a finite description a machine can read. It cannot assume that arbitrary infinite structures have one.

If a contradiction has a finite proof, proof search can find it. Failure to find one at a finite stage does not certify consistency. These are precise limitations on methods; they do not imply that every particular consistency question is undecidable.

Universal realization adds: mathematical existence already supplies the complete structure, and structural identity leaves no extra substance to make it real. A selector can distinguish parts of a larger structure; it cannot erase their mathematical existence. The commitment under dispute: existence needs no extra actuality condition (a further requirement for being real beyond mathematical existence). This paper holds that commitment, because a primitive actuality property (a basic, unanalyzable mark of being real) would either enter determinate relations with the structure, and so belong to it, or make no difference to anything.

A fictional sketch can be incomplete without being inconsistent. Describing a character as afraid does not build a system that feels fear. Completing that description would require an account of its relevant organization. No conclusion about the prevalence of conscious completions follows from the sketch alone.

## 6. Fictions as Experiences, Not Universes

A fiction is real as a description within mathematical reality. What it describes may have consistent completions, or may combine incompatible conditions. A separate question asks which completed structures contain the causal organization of an experiencing mind.

A fictional world need not be physically realized for someone to experience a representation of it. A dream is a process in an experiencing system; extending this to standalone minds requires both the universal realization premise and an implementation criterion for conscious dynamics. Under those additional premises, some otherwise unrealized scenarios could occur as represented experiences.

## 7. Caveats and Residuals

Universality alone does not license the conclusion that Level III swallows Level IV. Three distinctions block it.

*Ongoing computation.* A finite region cannot finish an infinite sequence. It can implement finite portions of an ongoing process. Nontermination (a program running forever) is not automatically zero probability; the probability construction must specify which outputs it counts.

*Noncomputable structures.* A prior over finite programs does not automatically cover arbitrary noncomputable structures. A computable stochastic law (a law that sets probabilities) can assign probability one to a class of individually noncomputable infinite outcomes. The zero probability of each singleton (each single outcome) does not give zero probability to the whole class.

*Implementation beyond representation.* A mathematical description can specify relations that a particular physical system does not implement. Even among computable candidates, resource limits and unreachable initial conditions can matter. These are substantive obstacles without introducing a nonstructural substance.

None of these distinctions licenses assigning zero total weight to the excluded class. A prior over finite outputs, a measure on infinite sequences and a physical observer measure have different domains.

## 8. Conclusion

One quantum system runs every computable mind, because one unitary evolution and one initial state execute every finite program to every finite stage it reaches, and a mind is its complete structure. Keeping the information each step would have erased makes every step reversible, and the explicit clock-and-workspace permutation proves unitarity.

The construction is many worlds at work. One unitary rule runs with no collapse, every program's history sits in one wave, and no branch holds the information of the whole, so what looks random from inside any one computation is which branch it is in, with no event left uncaused. Each mind in the wave is running, every state produced from the last, which is what separates a mind from a recording of one. A structure is one wherever it occurs, so the bare form of being a subject is the same in every one of them, with their memories and histories still distinct.

Whether our own cosmology realizes the construction, and how observations inside it are weighted, are questions about our universe, not about what one quantum rule can do.

## References

Bennett, C. H. (1973). Logical reversibility of computation. *IBM Journal of Research and Development* 17, 525-532. https://doi.org/10.1147/rd.176.0525

Everett, H. (1957). "Relative State" Formulation of Quantum Mechanics. *Reviews of Modern Physics* 29, 454.

Feynman, R. P. (1948). Space-Time Approach to Non-Relativistic Quantum Mechanics. *Reviews of Modern Physics* 20, 367.

Gödel, K. (1931). Über formal unentscheidbare Sätze der Principia Mathematica und verwandter Systeme I. *Monatshefte für Mathematik und Physik* 38, 173-198.

Matiyasevich, Y. (1970). Enumerable sets are Diophantine. *Doklady Akademii Nauk SSSR* 191, 279-282.

Schmidhuber, J. (1997). A Computer Scientist's View of Life, the Universe, and Everything. In *Foundations of Computer Science: Potential - Theory - Cognition*, ed. C. Freksa, Springer.

Schmidhuber, J. (2000). Algorithmic Theories of Everything. arXiv:quant-ph/0011122.

Solomonoff, R. (1964). A Formal Theory of Inductive Inference. *Information and Control* 7, 1-22 and 224-254.

Tegmark, M. (1998). Is "the theory of everything" merely the ultimate ensemble theory? *Annals of Physics* 270, 1.

Tegmark, M. (2014). *Our Mathematical Universe: My Quest for the Ultimate Nature of Reality.* Knopf.

Turing, A. M. (1936). On Computable Numbers, with an Application to the Entscheidungsproblem. *Proceedings of the London Mathematical Society* 42, 230-265.

Wallace, D. (2012). *The Emergent Multiverse: Quantum Theory According to the Everett Interpretation.* Oxford University Press.
