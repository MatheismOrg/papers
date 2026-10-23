# Which Mathematical Structure Does a Physical Theory Describe?

## Abstract

Every mathematical structure in a physical theory is real, the notation included: equations on a page are states of atoms, and the geometry they describe is structure in its own right. What a theory is about is the structure its predictions depend on, and four questions pick it out: what an element does, which changes preserve that role, what breaks when it is removed, and which observations tell the alternatives apart. Worked cases apply them. Paired real numbers can do the algebraic work of complex notation; a change of coordinates relabels a geometry without altering it; and a phase around a loop changes what electrons do even where the local field vanishes, so the loop phase is physical. One reality has many descriptions, and none of the mathematics in them is unreal.

## 1. Which structure?

A coordinate grid and the geometry it describes are both mathematical structures, and they answer different questions. Change the numbers assigned to events while preserving all distances and coincidences, and the physical predictions remain the same. Change the geometry itself, and clock comparisons or light propagation may change.

Both structures are real. The coordinate grid, whether printed on a page or stored in a computer, is a physical configuration of atoms or circuits, itself mathematical structure. The spacetime geometry is the structure whose relations determine what clocks and rulers measure. The grid does not become less real by being representational, and the geometry does not become more real by being represented. They are different mathematical structures playing different roles.

On the view that reality is mathematical structure (mathematical monism), a theory must still say which structure its predictions depend on. Fix a world model, meaning a structure together with its relevant relations and whatever history or boundary data (conditions fixed at a system's edges or starting time) the case needs. Then ask which of its mathematical relations the predictions depend on, and which choices express those relations in one notation rather than another.

Gauge freedom (the freedom to choose among descriptions giving identical predictions) and the dependence of results on representation are old lessons of physics [1, 5], and the recent debate over whether quantum theory needs complex numbers, which extend ordinary numbers with the square root of minus one, turns on the same question [2, 3, 4]. This paper's contribution is a method that works whether or not one accepts mathematical monism: four questions that separate the structure a theory is about from the structure it is written in. Its decisive premise is that two descriptions describe one world when they agree on every relation the predictions depend on, and physics decides which relations those are. It rules out reading ontology (what a theory says exists) off notation, and it rules out treating what is represented as something more than the structure its representation captures.

## 2. Closure, consistency and specification

Closure does not mean that removing any statement creates a contradiction. In formal logic, the deductive closure of a theory T contains the statements derivable from T. Deleting an axiom from a consistent classical theory cannot make the remaining axioms inconsistent. A theory with an independent statement P, one it can neither prove nor disprove, admits several possibilities before P is supplied. Adding P narrows those possibilities. Removing P restores alternatives; it does not destroy consistency.

A physical account is incomplete if it omits something its predictions depend on. Here, physical closure means supplying the structure, interactions and boundary information needed for the world under discussion. It is an explanatory requirement, distinct from logical closure. Two descriptions related by an accepted change of representation can still describe the same world.

Equivalent descriptions can be grouped together, but physics must first justify treating them as equivalent. Grouping them records the chosen equivalence; it does not establish that the choice preserves the relevant physical relations.

## 3. What an action principle supplies

An action principle, a law framed as a condition on a quantity called the action, by itself does not pick out one history. It constrains admissible histories through stationarity under specified variations. Boundary terms, allowed variations, constraints and initial data matter. Consider a free particle, one acted on by no forces, with action S[x] = integral of m(dx/dt)^2/2 over time. Its equation is d^2x/dt^2 = 0. Every trajectory x(t) = a + bt satisfies that equation. The action alone supplies neither a nor b. This blocks the inference from an internal action to a uniquely selected universe.

A mathematical world can include both its dynamics and its boundary data without an outside agent choosing either. A simplicity prior, a weighting that favors simpler descriptions before any evidence, may compare complete descriptions including their boundaries; it does not remove the need to identify what is being compared.

An action is not a prerequisite for mathematical structure or for every possible dynamics. A finite transition system, one with finitely many states, has a precise update rule whether or not anyone has found a useful variational formulation (a restatement of its rule as an action principle) for it. Where the action has a continuous symmetry (one holding for changes of any size, such as rotations), a conserved quantity can follow. Emmy Noether's first theorem connects suitable continuous variational symmetries with conservation laws under its hypotheses [1]. These results do not make every coordinate convention a physical field or guarantee a conserved global energy in every spacetime.

## 4. Four questions for identifying the referred structure

The following questions identify which mathematical structures a theory's predictions depend on. They are a method of analysis, not four independently established necessary and sufficient conditions.

Role. What relations or predictions does the element help specify? A boundary condition or fixed parameter can matter even when it is not varied in an action.

Equivalence. Which changes leave the predictions unchanged? Distinguish a relabeling from a change to the relations being labeled.

Removal. Does removing the element change predictions, leave the model underspecified, or permit an equivalent reformulation? An inconsistent reformulation is a fourth possibility; leaving an element out does not automatically produce one.

Discrimination. Do the alternatives differ logically, represent the same physics equivalently, or predict different observations? Different physical predictions invite empirical discrimination; identical predictions do not. They still leave room for an ontological argument, assessed a priori (by reasoning independent of observation) and through explanatory economy.

## 5. Worked cases

### 5.1 Complex amplitudes and real representations

Changing notation does not erase the role the notation performed. Complex arithmetic has a straightforward real representation: a + ib corresponds to the pair (a,b), and multiplication by i corresponds to J(a,b) = (-b,a). Applying J twice sends (a,b) to (-a,-b), which is the same as multiplying by -1, so J^2 = -I: two quarter-turns make a half-turn. Replacing complex scalars by paired real components can retain the relevant algebraic structure.

Experiments have tested one real-number formulation, and what they rule out depends on the assumptions built into it. Marc-Olivier Renou and collaborators [2] compared complex quantum theory with a real Hilbert space analogue (a version whose quantum state space uses only real numbers) in a network with independently prepared sources. Ming-Cheng Chen and collaborators [3] reported experimental correlations incompatible with the tested real formalism. The assumptions about composition (how separate systems are combined into one) and independence belong to the result; the experiments do not establish that every formulation using real numbers must fail. A 2026 preprint by Timothée Hoffreumon and Mischa Woods [4] explicitly distinguishes product state independence (requiring the sources' joint state to split into separate factors) from operational independence (no observable correlation between the sources) and claims an equivalence for finite protocols under operational independence.

For the method proposed here, the useful question is which transformations, correlations and composition rules an equivalent representation must preserve.

### 5.2 Coordinates and geometry

Relabeling is the clearest case where the structure the predictions depend on stays the same while the mathematical description changes. A coordinate change replaces a chart, one assignment of coordinates to points, and changes tensor components, the numbers representing geometric quantities in a chart, according to their transformation laws. The geometric object represented need not change. Moving the fields themselves is a different operation from relabeling the points under them. An active diffeomorphism, a smooth reversible remapping of spacetime points, acts on fields, and physical equivalence in general relativity concerns the appropriate transformation of the whole configuration, subject to boundary conditions.

The structure the predictions depend on is the geometry and relations of a spacetime model, not the numerical list g_mu_nu in one chart. Clock comparisons provide one way of probing that structure. Using a smooth manifold, a space that looks flat up close and supports calculus, specifies a mathematical model. Whether a physical continuum is fundamental or emerges from a discrete account remains a separate question.

### 5.3 Gauge structure and holonomy

Gauge theory shows that the boundary between the structures a theory's predictions depend on and the mathematical choices expressing them can fall in an unexpected place. In electromagnetism, adding the gradient of any function to the potential (the field from which electric and magnetic fields are derived), A_mu to A_mu + partial_mu lambda, leaves the field strength alone. In appropriate geometries, phase information (a wave's position within its cycle) around a loop can matter even when the local field strength vanishes along the accessible path, as in Yakir Aharonov and David Bohm's analysis [5]. For nonabelian gauge theories, those whose transformations give different results in different orders, the field strength transforms covariantly, changing in step with the transformation, rather than remaining componentwise invariant. There the loop quantity itself depends on the reference point chosen, while quantities built from it do not. A holonomy, the transformation accumulated around a loop, generally transforms by conjugation (being sandwiched between a transformation and its inverse), and suitable traces (sums of a matrix's diagonal entries) or conjugacy classes recover the gauge invariant information.

The method asks what survives an equivalent reformulation of the complete gauge structure.

## 6. Many structures, one world

Every mathematical structure involved in a physical theory is real, and they play different roles. The atoms forming ink on a page, the circuits storing coordinates in a computer, and the spacetime geometry those coordinates describe are all mathematical structures. They differ in what they contribute rather than in degree of reality: the theory's predictions depend on some of them and not on others.

A universe that contains pages, computers and brains performing calculations can take no more description than one without them: a short law and a simple starting state generate all of them, the way one rule generates every branch of a quantum history. Extra bits are charged only for what must be specified separately, and each such bit, a 50/50 choice that no rule makes, halves a description's weight. The representational structures are as real as anything they represent, and a world full of them is as probable as the short laws that produce it.

The further claim that all consistent mathematical structures exist rests on a separate argument, the rejection of a primitive actuality property, a basic unanalyzed feature making some structures actual and others not. Identifying which structure a particular theory's predictions depend on neither proves nor refutes that premise. It shows how the premise can coexist with the fact that a given world's physics depends on specific structures, while descriptions of those structures involve additional ones.

## 7. Limits and objections

Is the method circular? It is conditional on a model and its equivalences. That is appropriate for analyzing a theory already under discussion. Empirical comparison and philosophical argument retain separate work to do.

Does incompleteness, a theory's inability to settle every statement, prohibit closure? Deductive closure and completeness are different. A theory can contain all its consequences while leaving some sentences undecided. Kurt Gödel's results [6] concern appropriate consistent, effectively axiomatized theories, whose axioms a mechanical procedure can list. They do not show that every finite physical description is inadequate for every purpose.

Does the method derive physical realization? No. It identifies which structures a given account's predictions depend on. To derive a particular universe, one would still need admissibility conditions, dynamics, boundaries and a measure (a rule assigning probability weights) for comparing observer locations.

## 8. Conclusion

Every mathematical structure in a physical theory is real. Some of those structures are what the theory is about, and others are how the theory is written down. The four questions proposed here tell them apart: identify a role, test what preserves it, remove it and compare the consequences. A printed coordinate grid is atoms on a page; a spacetime geometry is what clocks and rulers measure. The theory's predictions depend on the geometry and not on which chart someone chose for it.

## References

[1] Noether, E. (1918). Invariante Variationsprobleme. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, 235-257.

[2] Renou, M. O., et al. (2021). Quantum theory based on real numbers can be experimentally falsified. Nature 600, 625-629. https://arxiv.org/abs/2101.10873

[3] Chen, M. C., et al. (2022). Ruling Out Real-Valued Standard Formalism of Quantum Theory. Physical Review Letters 128, 040403. https://arxiv.org/abs/2103.08123

[4] Hoffreumon, T., and Woods, M. P. (2026). Quantum theory based on real numbers cannot be experimentally falsified. Preprint, version 1. https://arxiv.org/abs/2603.19208

[5] Aharonov, Y., and Bohm, D. (1959). Significance of electromagnetic potentials in the quantum theory. Physical Review 115, 485-491.

[6] Gödel, K. (1931). Über formal unentscheidbare Sätze der Principia Mathematica und verwandter Systeme I. Monatshefte für Mathematik und Physik 38, 173-198.
