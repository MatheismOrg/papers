# Description length, entropy and physical prediction

Gary Abraham Bernstein

Academic consistency note · G47.5.0 · 22 September 2026

## What the simplicity proposal needs

Matheism identifies reality with mathematical structure. That ontology motivates a search for economical physical descriptions. A probability measure over those descriptions is an additional commitment. Its physical consequences require a specified encoding, candidate class, observational likelihood and rule for locating observers within the candidates. The ontology alone does not determine those choices.

Under a stated shortest description prior, w(T) is proportional to 2^(−K_U(T)). A program length L is a computable upper bound on K_U, not generally its measured value. Solomonoff's algorithmic probability instead sums weights over programs producing an output. These related constructions must not be treated as identical definitions. Infinite structures additionally require an effective representation and a measure on observations.

For two candidates compared using lengths in the same prefix code, a saving of ΔL bits gives prior odds 2^ΔL. Ten bits gives exactly 1,024; 100 bits gives approximately 1.27×10^30. Posterior odds multiply this by the likelihood ratio. The law, parameters, decoder, boundary data and observation map must all be charged consistently. A short rule producing a huge structure is not expensive merely because its output is huge. Runtime, memory and program length are different costs.

## Necessary observer development

Simple laws can remain unchanged while the structures they generate acquire the organization required for observers. Suppose the starting state lacks that organization and a later state contains observers who require it. The history must then develop the missing organization. This is a necessary consequence of those premises, not merely a possible outcome. A37 states it as an observer threshold crossing. The growing quantity is the observer relevant organization of the generated states; the law and the shortest description of the history can remain simple.

Observer conditioning supplies the second part of the argument. Histories that never contain the observers whose evidence is being explained have zero weight for that evidence. Within the class that begins below the required organization and later produces those observers, development is necessary. Comparing that class with histories that already contain observers at the chosen boundary is a separate selection calculation. Simple laws, simple total description and an initially unorganized state are three different premises.

For observers like us, the physical account also includes energy gradients, metabolism, maintained records and exchanges with an environment. Their organization can develop while entropy increases in the larger system. This connects the developmental claim to our thermodynamic setting without identifying organization with entropy or assuming that either grows at every step.

See [A37, Section 5](../../EDITORIAL/ACADEMIC_NOTES/A37_Arrow_of_Complexity.md), for the precise developmental result and its observer selection extension. A correction to an entropy or reversal argument must preserve this result wherever the developmental proposal is summarized.

## Finite reversal

Fix a finite effective encoding of a history h and a computable involution R that reverses it. A fixed wrapper runs a shortest program for h and then applies R. Hence K_U(Rh) ≤ K_U(h)+c_R. Apply the same wrapper again to obtain K_U(h) ≤ K_U(Rh)+c_R. Therefore |K_U(Rh)−K_U(h)| ≤ c_R. The constant depends on the language and encoding, not on the history's length.

This proves a bounded description overhead. It does not establish a numerically small prior ratio, an efficient reversal algorithm, or cheap laboratory control. An exact continuum microstate need not have a finite effective encoding; the result concerns the stated finite descriptions.

For computable evolution F and integer t, a program can compute the endpoint from F, x₀ and t. With a fixed tuple convention,

K_U(F^t(x₀)) ≤ K_U(F)+K_U(x₀)+K_U(t)+c.

A final state generated this way cannot be declared incompressible merely because it looks disordered. Restricting the reverse description to a literal list while allowing a generator for the forward description compares different tasks. Operational reversal costs can still be large; A13 specifies preparation models that make such a comparison meaningful.

## Low complexity and low entropy

Let N be even. Divide the N bit strings into macrostates by their number of ones. The Boltzmann entropy of the macrostate with m ones, in units with the entropy constant equal to one, is log binomial(N,m). The alternating string 0101… has K at most K(N)+c, yet belongs to the largest macroregion, m=N/2. The all zero string is also simple and belongs to a macroregion of size one.

Thus simple microstates can belong to both maximal and minimal entropy macroregions. A short description of a probability distribution is also not a short description of every draw from that distribution. An entropy argument must keep the distribution, macrostate and particular microstate separate.

The research target is P(low entropy boundary and reliable records | observations, model class), evaluated under explicit laws, macrostates and measure. A simplicity preference alone does not calculate that probability. Nor does preferring each short candidate establish that the class of short candidates dominates a larger class after summing its weights.

## Coarse graining and the arrow

A coarse map C sends a microstate x to a macrostate C(x). It may be many to one even when microscopic evolution F is invertible. The observed history is C(F^t(x)). The noninjectivity of C alone does not make its entropy increase with t.

For a direct counterexample to the information claim, a reset map sends every binary string to zero. For a uniform input its Shannon entropy falls from N bits to zero. A physical reset has an environment and resource costs; the mathematical map alone gives neither a total thermodynamic entropy change nor a temperature.

A Past Hypothesis, suitable dynamics and statistical assumptions can support a thermodynamic arrow. A simplicity measure is a proposed further explanation of boundary selection, not a universally necessary extra axiom for statistical mechanics. Alternatives to imposing a special initial boundary are also investigated; Lazarovici and Reichert (2020) discuss examples. Invertibility itself is not the same as invariance under a physical time reversal transformation.

Reliable records require a physical account of their creation and stability. Calling the observer's computation forward labels its orientation; it does not derive that orientation. Reversible computation and lossless encoding remain possible. A lossy record's failure to identify a unique input establishes an inference limitation, not the impossibility of observers in a reversed physical solution.

## Action, emergence and quantum rules

Noether's theorem links appropriate continuous variational symmetries to conservation laws within its hypotheses. It is not a theorem that any stable observer requires an action formulation. Coarse graining alone does not supply an action, complex quantum amplitudes or the Born rule. Where a stationary phase approximation applies, contributions near stationary action can dominate; multiple stationary paths and quantum corrections can remain.

A useful discrete search must identify its states, update rule, metric or causal structure, observation map and limiting regime. It must recover quantitative quantum and gravitational predictions, not just graphs that resemble branching histories or coordinate independence. A39 supplies concrete tests of representation and geometry. The shorter complete empirically adequate account earns the stipulated prior advantage.

## Companion work and references

[A0a](../ALTERNATIVES/A0a_Core_Matheism.md) states the ontology and distinguishes its physical extensions. [A13](../MANUSCRIPTS/A13_Arrow_of_Time_Description_Length_Asymmetry.md) develops finite reversal, records and preparation costs. [A32](../../EDITORIAL/ACADEMIC_NOTES/A32_Simplicity_Dominance_Three_Consequences.md) proposes applications of the simplicity measure. [A38](../MANUSCRIPTS/A38_Kolmogorov_Standard_Model.md) supplies a coding case study whose implementation counts are not exact Kolmogorov complexities. [A39](../MANUSCRIPTS/A39_Emperors_New_Points.md) distinguishes mathematical representation from physical realization. These are companion proposals and calculations, not independent validation of one another.

Bennett, C. H. 1973. Logical reversibility of computation. IBM Journal of Research and Development 17, 525 to 532. https://doi.org/10.1147/rd.176.0525

Noether, E. 1918. Invariant Variation Problems. English translation by M. A. Tavel. https://arxiv.org/abs/physics/0503066

Lazarovici, D., and P. Reichert. 2020. Arrow(s) of Time without a Past Hypothesis. https://arxiv.org/abs/1809.04646

Solomonoff, R. J. 1964. A formal theory of inductive inference, Parts I and II. Information and Control 7, 1 to 22 and 224 to 254.

The finite string examples above are elementary checks, not claims of new mathematical results.


## Observer implementation and cosmological reachability

[A30](../MANUSCRIPTS/S48_A30_Consciousness_G26.0.0_V1.0.0_w4197_1521_260712.md) separates the reality of an experience from the truth of its represented world. Its standalone mind argument requires both the realization premise and an adequate organizational criterion. Neither a narrative label nor one formally defined transition establishes consciousness. A12 supplies a proposal to investigate, not a completed sufficiency theorem.

[A31](../MANUSCRIPTS/S48_A31_BB_Ergodic_G26.0.0_V1.0.0_w3608_1521_260712.md) separates reachable histories, functioning observers and observational weight. Determinism neither excludes nor guarantees ergodicity. Unitary evolution does not preserve every zero coefficient in an arbitrary basis. Infinite duration alone does not determine integrated formation rates or the measure of our evidence.

The simplicity comparison must permit generative descriptions for both ordinary and fluctuation observers. A21 makes the additional observer sampling choice explicit. A28 distinguishes universal simulation capability from implementation in a particular cosmology. These dependencies prevent the companion papers from treating one another’s open questions as established results.


## Closure and constructor dependencies

A6 uses physical closure as a proposed complete account of internal relations, not the logical claim that deleting any axiom creates inconsistency. A partial specification can be consistent. An action may admit multiple histories and require boundary data. Complex quantum network results depend on model composition and source independence; they do not prohibit every real representation.

A7 separates mathematical consistency, task possibility under specified laws and actual realization. Its reversible copying example preserves both a record and a reusable device marker. No cloning is a quantum transformation constraint, not the absence of an inverse for addition. Necessary development of observer organization does not imply increasing law complexity, a universal complexity gradient, or a derivation of the Past Hypothesis. Infinite sequences without a last element do not prohibit all reverse oriented histories.


## Quantum model comparison and measurement independence

A8 and A9 compare complete predictive models, including preparation and observation rules. Outcome likelihood is not an automatic additional law cost. A finite random outcome has finite complexity. Shared determination need not produce statistical dependence, and a short measurement dependent rule need not recover quantum correlations. A9 includes a finite CHSH family to make the latter distinction explicit.

Specificity, coding complexity, causal faithfulness, robustness and typicality are different questions. Neither an attractor label nor a changed metric supplies a preparation measure by itself. A4's quantum section uses these distinctions. Claims that the cited relaxation proposal proves generic convergence, that collapse necessarily loses on infinite K, or that a block universe alone entails measurement dependence are superseded.


## Proof, physical identification and experiment

Do not replace established logical implications with merely empirical conjectures. Conditional observer development, the nonuniqueness of a lossy inverse, the finite history reversal bound and 2^delta description weight ratios have deductive status under their stated premises. A mathematical derivation of quantum or relativistic predictions is likewise possible without conducting an experiment.

Keep four questions separate: whether an inference is valid; whether its premises are justified, including by a priori argument; whether a particular model realizes those premises; and whether the model describes observed physics. Experiments address observational application and discrimination, not the validity of a correct deduction. Empirical equivalence does not prohibit a priori comparison, but a prior preference also does not become a uniquely established physical law without its missing argument.

The core Matheism identity argument is intended as an a priori argument. Calling it philosophical does not demote it to opinion, and calling it a proof does not establish unstated premises. Preserve its intended mode of justification while identifying the exact steps under assessment. Necessary development of observer organization remains necessary given its premises, not merely likely.
