# Reversible Physics, Incomplete Records: Why an Answer Can Hide Its Inputs

## Abstract

Reversible physics loses no information when an answer hides its inputs. A sum of 2 does not say which pair was added, yet a machine that keeps both addends beside the sum has thrown nothing away: the many-to-one function has no inverse on its outputs, while the reversible operation carrying it out still holds the input. Every verdict that a task is impossible must therefore cover the whole operation and what it retains, which is how constructor theory frames physical possibility, as the tasks the laws permit a reusable device to perform. Mathematical monism identifies the devices and laws with structure, and the distinction shows how records, copiers and observers can be built under reversible laws. Quantum measurement is the same case at the largest scale: no branch, observer or record holds the information of the whole wave, which stays in correlations, so an outcome that looks random from inside needs no uncaused collapse. Time's direction comes from production, not from incomplete records: a history grown step by step from a simple start needs a first state and has no last one to run back from, and a reversal must pay for the step it turns back from.

Keywords: constructor theory, mathematical monism, reversible computation, records, physical possibility, observer development

## 1. The common question

A machine can hand you an answer that does not determine its own inputs and still lose nothing: the inputs remain in the rest of its state. The map from inputs to the displayed answer can be many-to-one, sending different inputs to the same output, while the physical operation carrying it out stays reversible, able to run backward to recover its input. Constructor theory [1] supplies the setting for that distinction: it asks which tasks physical laws allow a machine to perform reliably. A constructor is a system that can perform a task and retain the ability to do it again. Its microscopic details may change. This paper connects that approach to mathematical monism, the view that reality is mathematical structure, through simple examples of copying, recording and reconstructing information.

David Deutsch and Chiara Marletto apply constructor theory to information: what a system can copy, store or distinguish follows from what the laws permit, rather than standing as an independent mathematical primitive, a basic notion defined from nothing else [2]. A monist interpretation can identify those physical relations with mathematical structure. It should not misattribute that identification to constructor theory itself.

The two approaches meet at the task specification. A proposed transformation has to state an input, an output and the constraints on the device and its environment. Treating the whole arrangement as structure keeps those commitments in view. The examples also bear on time. Its direction comes from generation: a history generated step by step has a first state and, being endless, no last one to run back from, and a reversal must pay in bits for its turning point, each bit halving its weight.

Rolf Landauer and Charles Bennett showed that computation can be made logically reversible [3, 4], and constructor theory recasts physics as a theory of which tasks are possible [1, 2]. This paper's contribution is to bring the first to bear on the second: a many-to-one record does not show that a process destroyed information, so a verdict of impossibility has to cover the whole operation and everything it keeps, as the worked examples of reversible addition and a reversible record make exact. The decisive premise is that information is lost only when no part of the complete system keeps it. It rules out reading destroyed information, or the direction of time, off an incomplete record.

## 2. Three distinct claims

Three claims about a task are easily run together and have to be kept apart.

Consistency: a specified input and output have a mathematical model under the stated definitions.

Physical possibility: the laws of the world under consideration permit a constructor for the task with the required accuracy and reliability.

Realization: the resources, the history and the actual arrangement contain something that performs the task.

A consistent map from an unknown quantum state to two copies is not thereby an allowed quantum operation. Conversely, a task permitted by the laws may go unattempted. A proposed world with different transformation laws can assign different possibilities; to compare it with our world, specify its laws and test its predictions.

## 3. The monist premise

Matheism argues that reality is mathematical structure. Its rule-random dichotomy covers both compressible and incompressible data, data with and without a shorter description. Its interaction argument requires any proposed extra causal ingredient to enter determinate relations with the physical system. Its identity thesis adds that the complete structure exhausts the relevant properties. A physical structure and a description of it still have to be told apart: when one chosen pattern fails, that does not establish algorithmic randomness (having no description shorter than the data itself), and classifying a task does not establish that anything implements it.

I use mathematical monism as the explicit interpretive premise here. Since all physical distinctions are structural, a physically possible task is a relation supported by the structure of the relevant laws and substrates, the physical systems that tasks act on. This reading requires no additional agent to impose the laws from outside, though it still requires the laws and substrates to be specified. Universal realization rests on a separate argument, the rejection of a primitive actuality property (a basic unanalyzed feature making some structures actual and others not), since describing a filter that excludes nothing does not by itself show that every structure it admits is instantiated. Existence itself is not an ordinary constructor task with physical input and output attributes.

## 4. A reduced function and its implementation

Addition is the plainest case of a function whose output does not say which inputs produced it. Written out, let f(a,b) = a + b over integers. Knowing the output 2 does not select one input pair; the inverse relation contains (0,2), (1,1), (2,0) and infinitely many others. There are nevertheless computable selections: h(s) = (0,s) returns a pair with sum s. This does not recover the historical pair that was actually supplied. Producing a compatible input and recovering the particular past input are two different achievements.

Keep the inputs alongside the answer and the addition can be reversed. In symbols, define G(a,b,c) = (a,b,c + a + b), whose inverse is G_inverse(a,b,d) = (a,b,d - a - b). With c initially zero, the third register, a slot that stores one value, records the sum. Projecting away the first two registers recovers the reduced output f(a,b). The full map is reversible while the reduced input-to-output map is not injective, meaning different inputs can share an output.

Both descriptions can be physically useful. The reduced function's lost distinctions can remain in other degrees of freedom, independent quantities that describe the system. A local record or erasure does not prove that the full dynamics loses information.

## 5. Reusable constructors and records

A reversible machine can take a record and still be ready for the next run. In symbols, with an unchanged device marker c, a substrate bit x and a memory bit m, define U(c,x,m) = (c,x,m XOR x). U is its own inverse. For a prepared memory m = 0, its output is (c,x,x). The memory now records the substrate, the device marker is unchanged, and the map on all eight inputs is a permutation, a reshuffling that sends each input to a distinct output. This finite example refutes the claim that a reversible machine cannot both keep a record and stay ready to run again.

Charles Bennett showed how to build a reversible computer that retains enough information for inversion and can uncompute intermediate work, running those steps backward to clear them, while preserving the desired result and the input [3]. Registers can encode updated estimates while the larger system retains earlier states. Whether a particular update counts as useful learning depends on the model and its predictive performance.

Rolf Landauer's principle, that erasing information releases a minimum amount of heat, concerns logically irreversible erasure under specified physical conditions [4]. The familiar k_B T ln 2 minimum concerns resetting an unbiased bit, one equally likely to be 0 or 1, in the appropriate isothermal setting, held at constant temperature. Writing a correlation into prepared memory is not automatically that operation. A full cycle must still account for the costs of preparation, storage, reset and environmental exchange.

## 6. What this establishes about time

Time runs forward for a different reason: a history generated step by step needs a first state, and an endless one has no last state from which to run it backward. A reversed history must also specify the step it turns back from, and each bit of that halves its weight again. Many-to-one maps alone do not select which temporal boundary is low in entropy, a measure of how many microscopic arrangements match the large-scale state. A local task may discard distinctions accessible to an observer. Running that task again may require a supply of ordered resources and somewhere to dump the entropy. These are physical questions even when microscopic dynamics is reversible. The direction in which ordinary records become reliable depends on the dynamics, the available resources and the statistical assumptions imposed on histories.

An observer generated from a state that lacks its organization necessarily acquires that organization as it develops. That implication is stronger than saying development is merely possible. It still does not require the laws themselves to grow more complex, or every measure of complexity to increase at every step.

## 7. Generated and listed histories

A short generator plus a start state is enough to place a ceiling on how much description a later state needs: describing the rule, the starting state and the elapsed time describes the later state as well. Stated precisely, suppose a computable reversible dynamics F evolves a finite encoded state x_0 to x_t = F^t(x_0). Then K(x_t) <= K(F) + K(x_0) + K(t) + O(1). A complicated looking configuration can have a short generating description. When the inverse dynamics is computable, the same kind of description can recover earlier states from a complete later state and elapsed time.

Using the forward generator to specify a finite reverse history avoids listing every coordinated particle separately. If the reversal operation adds k encoded bits to the unchanged forward program, its description receives 2^(-k) times the weight under the same binary measure, in which each bit is a 50/50 choice that no rule makes. Ten added bits give it one-1,024th the weight. Whether another program can absorb reversal without extra bits is the remaining question about shortest descriptions.

A process that produces each step from an earlier one needs at least one input supplied at the start. If every operation waits on an output none has yet produced, none is enabled. Anyone applying constructor theory here must say why the physical task requires such a production chain rather than relations specified jointly across a block, the whole of spacetime taken at once.

## 8. Quantum copying is a different obstruction

No single physical operation can exactly copy every arbitrary pure state, a quantum state described as completely as quantum theory permits. That is the quantum no-cloning result of William Wootters and Wojciech Zurek [5]. Written out, take a unitary copier, one acting by reversible quantum evolution, and blank |b>, with the desired mappings |psi>|b> to |psi>|psi> and |phi>|b> to |phi>|phi>. A unitary operation preserves inner products (measures of overlap between states), and preserving them here would require s = s squared, where s = <psi|phi>. Distinct nonorthogonal states, which partly overlap, violate that condition. The proposed abstract map psi to (psi,psi) is injective; the issue is whether quantum transformations allow it.

Constructor theory organizes such distinctions. Mathematical monism can interpret the relevant constraints as structural facts without pretending that consistency alone derives the quantum state space or its transformation rules.

Measurement under unitary evolution is the displayed sum at the largest scale. A measuring device that records an outcome becomes correlated with the system, and each branch of the joint state (a component corresponding to one outcome) carries one record. The record does not determine the whole superposition, the weighted combination of all the alternatives, and nothing has been destroyed: the information the record lacks stays in the correlations between the branches and their environment. An observer who sees one outcome holds one share of the whole wave, as the displayed 2 holds one share of the computation that produced it. What looks random to that observer is which share it holds. No outcome is selected by nothing, and no collapse, a sudden jump of the quantum state to one outcome, is needed to explain why an observer cannot reach information that the whole still keeps.

## 9. Conservation and biological development

A system that supports constructors need not possess every familiar conservation law. Emmy Noether's theorem, which links symmetries to conservation laws, delivers less than that: it applies to suitable variational symmetries [6]. A theory with explicit time dependence, laws that themselves change over time, can lack the relevant symmetry while still being a consistent theory of transformations.

Biological systems can be studied through repeated transformations, resource use, maintenance and reproduction. The claim about development is this: a history that begins without the organization a later modeler requires must develop that organization before the modeler exists. This does not assign a probability to modelers without a measure, and it does not equate biological organization, thermodynamic entropy and Kolmogorov complexity, the length of the shortest program producing a given output.

## 10. What would advance the program

Specify a substrate, admissible states, dynamics, resource states and the task to be repeated. Track the device, the output and the environment. For observer development, supply an operational criterion of the relevant modeling organization and calculate when it appears. State separately the criterion for experience, if any.

For simplicity comparisons, fix the description language and compare complete candidate accounts against the same observations. For computability, state an additional physical computability premise if the model needs one. Mathematical monism alone includes structures that are not computable.

## 11. Conclusion

The answer a machine displays and the process that produced it can hold different amounts of information, because the process can keep in its other parts what the displayed answer leaves out. The worked examples show why reconstructing an input from one record can be impossible even when the underlying physical operation is reversible.

An incomplete record is also how quantum measurement looks from inside a branch: the whole wave keeps what each record leaves out, so nothing needs to be destroyed or selected without a cause.

Time runs forward because a simple start has to grow its observers in causal steps: a history grown step by step needs a first state and, being endless, has no last state to run back from. The steps follow from an a priori proof in the logician's sense: a cause produces its effect as the next state, and a continuum has no next state. A reversed history must also specify the step it turns back from, and each bit of that, an unruled 50/50 choice, halves its weight. Many-to-one maps and incomplete records alone do not set that direction.

Constructor theory specifies which tasks laws permit. Mathematical monism identifies the systems and laws with structure. The two pay off together when a model shows which resources and starting conditions let records, reusable devices and observers develop.

## References

[1] Deutsch, D. (2013). Constructor theory. Synthese 190, 4331 to 4359. https://arxiv.org/abs/1210.7439

[2] Deutsch, D., and Marletto, C. (2015). Constructor theory of information. Proceedings of the Royal Society A 471, 20140540. https://arxiv.org/abs/1405.5563

[3] Bennett, C. H. (1973). Logical reversibility of computation. IBM Journal of Research and Development 17, 525 to 532. https://www.cs.princeton.edu/courses/archive/fall04/cos576/papers/bennett73.html

[4] Landauer, R. (1961). Irreversibility and heat generation in the computing process. IBM Journal of Research and Development 5, 183 to 191.

[5] Wootters, W. K., and Zurek, W. H. (1982). A single quantum cannot be cloned. Nature 299, 802 to 803. https://doi.org/10.1038/299802a0

[6] Noether, E. (1918). Invariante Variationsprobleme. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, 235 to 257.
