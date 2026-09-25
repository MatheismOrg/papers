# Reversible Physics, Incomplete Records: Why an Answer Can Hide Its Inputs

**Gary Abraham Bernstein**

Independent Researcher · ORCID: 0009-0009-1761-2867

## Abstract

The answer 2 does not tell us which numbers were added. Yet a machine can preserve its inputs even when its displayed answer omits them. This paper uses copying and memory examples to separate incomplete records from destruction of information in the complete system. It connects that distinction to constructor theory, which describes the tasks physical laws permit or forbid. A reusable device must retain its ability to work, not every microscopic detail of its previous state. These examples show why claims of irreversibility must identify the inputs, retained records and full system. They also clarify what an account of observer development must explain beyond merely displaying an answer.

**Keywords:** constructor theory, mathematical monism, reversible computation, records, physical possibility, observer development

## 1. The common question

What tasks can physical laws allow a machine to perform reliably? Constructor theory asks that question. A constructor is a system that can perform a task and retain the ability to do it again. Its microscopic details may change. This paper connects that approach to mathematical monism through simple examples of copying, recording and reconstructing information.

Deutsch and Marletto apply this framework to information, with informational properties grounded in physical possibilities rather than assumed as independent mathematical primitives [2]. A monist interpretation can identify those physical relations with mathematical structure. It should not misattribute that identification to constructor theory itself.

The useful meeting point is the task specification. A proposed transformation names an input, an output and constraints on the device and its environment. Treating the whole arrangement as structure makes those commitments explicit. It does not automatically determine whether the arrangement can be built under our laws.

## 2. Three distinct claims

**Consistency:** an input and output specification has a mathematical model under the stated definitions.

**Physical possibility:** the laws of the world under consideration permit a constructor for the task with the required accuracy and reliability.

**Realization:** the resources, history and actual arrangement include an implementation of the task.

A consistent map from an unknown quantum state to two copies is not thereby an allowed quantum operation. Conversely, a task permitted by the laws may never be attempted. These elementary distinctions prevent universal mathematical realism from making every task possible inside every universe.

A proposed world with different transformation laws can assign different possibilities. To compare it with our world, specify its laws and test its predictions. The statement that some consistent structure contains a transformation is weaker than the claim that our physics permits it.

## 3. The monist premise

Matheism argues that reality is mathematical structure. The rules and the objects they relate are parts of one account. Adding the word substance does not explain an additional causal power unless some role is specified for it.

The logical partition between pattern and absence of a chosen pattern is exhaustive, but it does not by itself establish an identity between mathematical description and described reality. Algorithmic randomness also presupposes an encoding and a definition of randomness. A patterned object can have high complexity under one description, and failure of one pattern is not a definition of maximal randomness. Constructor theory does not supply the missing ontological inference simply by classifying tasks.

I therefore use mathematical monism as the explicit interpretive premise here. If all physical distinctions are structural, a physically possible task is a relation supported by the structure of the relevant laws and substrates. This reading requires no additional agent to impose the laws from outside. It still requires the laws and substrates to be specified.

Universal realization is a further premise. The mathematical description of a filter that excludes no structures does not prove that every described structure is physically instantiated. Nor is existence itself an ordinary constructor task with physical input and output attributes. These arguments concern ontology; they cannot substitute for a constructor specification.

## 4. A reduced function and its implementation

Let f(a,b) = a + b over integers. Knowing the output 2 does not select one input pair. The inverse relation contains (0,2), (1,1), (2,0) and infinitely many others. The absence of a unique inverse is a fact about this function and its domain.

There are nevertheless computable selections from that relation: h(s) = (0,s) always returns a pair with sum s. This does not recover the historical pair that was actually supplied. The distinction is between producing a compatible input and recovering a particular past input. Saying there is nothing computable in the backward direction overstates the issue.

Now define G(a,b,c) = (a,b,c + a + b). Its inverse is G_inverse(a,b,d) = (a,b,d − a − b). With c initially zero, the third register records the sum. Projecting away the first two registers recovers the reduced output f(a,b). The full map is reversible while the reduced input to output map is not injective.

Both descriptions can be physically useful. Their mathematical reality does not make them identical descriptions with identical information. The reduced function's lost distinctions can remain in other degrees of freedom. That is why a local record or erasure does not prove that the full dynamics loses information.

## 5. Reusable constructors and records

Consider three bits: an unchanged device marker c, a substrate bit x and a memory bit m. Define U(c,x,m) = (c,x,m XOR x). U is its own inverse. For a prepared memory m = 0, its output is (c,x,x). The memory now records the substrate, the device marker is unchanged, and the map on all eight inputs is a permutation.

This finite example is not a model of a complete conscious observer. It refutes the logical claim that retaining a record and a reusable device state cannot coexist with reversible dynamics. The prepared memory is an explicit resource. Repeating the task on fresh memories requires further prepared substrates; hiding their preparation would hide part of the task.

Bennett's reversible computation construction retains enough information for inversion and can uncompute intermediate work while preserving the desired result and input [3]. It does not require every intermediate result to remain forever. Its relevance here is that temporary workspace can be restored without deleting all evidence of the computation.

Reversible learning is therefore not excluded by definition. A register can encode an updated estimate while the larger system retains the old estimate and evidence. Whether a particular update counts as useful learning depends on the model and its predictive performance. No requirement to erase all earlier alternatives follows from the ordinary concept of learning.

Landauer's principle concerns logically irreversible erasure under specified physical conditions [4]. The familiar k_B T ln 2 minimum concerns resetting an unbiased bit in the appropriate isothermal setting. Writing a correlation into prepared memory is not automatically that operation. Costs of preparation, storage, reset and environmental exchange must be accounted for in the full cycle.

## 6. What this establishes about time

A local task may discard distinctions accessible to an observer. Its repeatability may require a supply of ordered resources and disposal of entropy elsewhere. These are physical questions even when microscopic dynamics is reversible. Calling the local operation a coarse description does not make its practical constraints imaginary.

But many to one maps alone do not select which temporal boundary is low in entropy. An isolated equation f(x) = y supplies no cosmological probability measure or boundary condition. The direction in which ordinary records become reliable depends on the dynamics, the available resources and the statistical assumptions imposed on histories.

A short description alone does not establish low thermodynamic entropy. Nor does conditioning on observers by itself prove the conventional Past Hypothesis. The relation between simple boundaries and the relevant entropy must be demonstrated for the model being used.

An observer generated from a state lacking its organization necessarily acquires that organization along its developmental history. That implication is stronger than saying development is merely possible. It still does not require the laws themselves to grow more complex, or every measure of complexity to increase at every step. The claim concerns the organization needed for the observer.

Under global unitary dynamics, local record formation can occur without a collapse postulate. This compatibility does not establish that Many Worlds has lower total description length than every competitor. Such a comparison must include ontology, dynamics, boundaries, probability assignments and observational likelihoods, rather than count only one additional equation.

## 7. Generated and listed histories

Suppose a computable reversible dynamics F evolves a finite encoded state x_0 to x_t = F^t(x_0). Then K(x_t) ≤ K(F) + K(x_0) + K(t) + O(1). A complicated looking configuration can have a short generating description. When the inverse dynamics is computable, the same kind of description can recover earlier states from a complete later state and elapsed time.

A billiard illustration must specify whether it tracks an ideal conservative system or omits the table, friction and environment. In a reversible idealization, complete positions and momenta determine both directions. In a dissipative reduced model, the omitted degrees of freedom can carry the missing information. Neither case proves that the reverse description must list every particle independently.

The instruction to reverse a finite history also has a short description. An additional five bits would give a factor of 32 between specified description weights under a fixed binary prior. It is not a machine independent five bit penalty attached to the direction of time. One must show what is being encoded and include the likelihood and the full histories in the comparison.

A reverse dependency chain in which every operation requires a later output cannot begin if it has no final source and no independently supplied input. This is a constraint on initiated production. A constructor-theoretic application must state why the physical task requires that production chain: a jointly specified block is not generated merely by executing its reverse relations.

## 8. Quantum copying is a different obstruction

The quantum no cloning result concerns the impossibility of one physical operation exactly copying every arbitrary pure state [5]. For a unitary copier and blank |b>, preservation of inner products would require s = s², where s = <psi|phi>, for the desired mappings |psi>|b> to |psi>|psi> and |phi>|b> to |phi>|phi>. Distinct nonorthogonal states violate that condition.

This is not an inversion of addition or a general many to one obstruction. The proposed abstract map psi to (psi,psi) is injective; the issue is its compatibility with quantum transformations. A claim about allowed physics needs those additional constraints.

Copying distinguishable classical alternatives and cloning arbitrary quantum states are thus different tasks. Constructor theory organizes such distinctions. Mathematical monism can interpret the relevant constraints as structural facts without pretending that consistency alone derives the quantum state space or its transformation rules.

## 9. Conservation and biological development

Noether's theorem applies to suitable variational symmetries [6]. It does not show that any system supporting constructors must possess every familiar conservation law. A theory with explicit time dependence can lack the relevant time translation symmetry while still being a consistent theory of transformations. Conservation claims must specify the system and exchanges with its surroundings.

Biological systems can be studied through repeated transformations, resource use, maintenance and reproduction. Treating a cell as an approximate constructor can be useful for a specified task. A daughter cell is not an exact restored microscopic state of its parent, and evolution does not guarantee increasing complexity in every lineage.

The relevant developmental statement remains conditional and necessary: if a history begins without the organization required for a later modeler, that organization must develop before that modeler exists. A simple rule can generate it. This does not assign a probability to modelers without a measure, and it does not equate biological organization, thermodynamic entropy and Kolmogorov complexity.

## 10. Memory, prediction and observer orientation

A record is a physical correlation interpreted and used by a system. Memory of an event ordinarily includes an appropriate causal history; prediction is an inference about an event. Present correlation alone does not distinguish them. A clock can predict tomorrow's sunrise without recalling a tomorrow already experienced.

Simple initial states can contain extensive correlations. The finite string consisting of N identical bits has a short description for arbitrarily large N. Consequently, low description length does not entail an absence of initial correlations. A cosmological argument must specify which correlations its boundary assumptions permit.

A time reversed history of an entire reversible observer system cannot be ruled out merely by renaming its coordinate future our past. If the reversed history satisfies the theory's transformation and boundary requirements, claims about experience must track that observer's records and organization. This does not prove such observers are common or achievable in our environment. It prevents a coordinate choice from becoming an impossibility theorem.

## 11. What would advance the programme

Specify a substrate, admissible states, dynamics, resource states and the task to be repeated. Track the device, the output and the environment. Determine whether a sequence of implementations can achieve the required accuracy while retaining the task capability. Compare reverse and forward tasks under the same explicit resource accounting.

For observer development, supply an operational criterion of the relevant modelling organization and calculate when it appears. State separately the criterion for experience, if any. Self reference alone is insufficient as an established test: a simple recursive counter refers to its own state, but constructor theory does not identify it as conscious.

A proof that a task is possible or impossible under stated mathematical laws is a deductive result. It does not need a laboratory demonstration to become a theorem. Experimental work concerns whether those laws and resource assumptions apply to an actual system, or whether a proposed implementation achieves its stated behavior. The monist argument, in turn, concerns the nature of such laws and systems and must be evaluated on its own logical premises.

For simplicity comparisons, fix the description language and compare complete candidate accounts against the same observations. Under a specified prior, shorter descriptions can receive greater weight. This does not turn every compact mathematical construction into our physics, and it does not supply an observer measure automatically.

For computability, state an additional physical computability premise if the model needs one. Mathematical monism alone includes mathematical structures that are not computable. Discovering essential physical noncomputability would challenge a computable implementation thesis, not mathematical structure as such.

## 12. Conclusion

A displayed answer and the complete process producing it can retain different amounts of information. The worked examples show why reconstructing an input from one record can be impossible even when the underlying physical operation is reversible.

Constructor theory specifies which tasks laws permit. Mathematical monism identifies the systems and laws with structure. Their connection becomes useful when a model shows which resources and starting conditions allow records, reusable devices and observers to develop.

## References

[1] Deutsch, D. (2013). Constructor theory. Synthese 190, 4331 to 4359. https://arxiv.org/abs/1210.7439

[2] Deutsch, D., and Marletto, C. (2015). Constructor theory of information. Proceedings of the Royal Society A 471, 20140540. https://arxiv.org/abs/1405.5563

[3] Bennett, C. H. (1973). Logical reversibility of computation. IBM Journal of Research and Development 17, 525 to 532. https://www.cs.princeton.edu/courses/archive/fall04/cos576/papers/bennett73.html

[4] Landauer, R. (1961). Irreversibility and heat generation in the computing process. IBM Journal of Research and Development 5, 183 to 191.

[5] Wootters, W. K., and Zurek, W. H. (1982). A single quantum cannot be cloned. Nature 299, 802 to 803. https://doi.org/10.1038/299802a0

[6] Noether, E. (1918). Invariante Variationsprobleme. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, 235 to 257.
