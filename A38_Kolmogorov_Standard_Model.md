# How Many Bits Is Physics? Comparing Complete Theories

**Author:** Gary Abraham Bernstein · Independent Researcher · gary@matheism.org · ORCID 0009-0009-1761-2867

## Abstract

A short name for a theory can hide a long explanation. How many bits are needed when its equations, parameters, starting information and decoding rules are actually included? This paper sets a common task and precision before comparing descriptions. A worked code audit shows why a small compiled file need not contain the physics claimed for it. A second example gives equal, one-sided and independent couplings explicit codes: shorter individual descriptions gain weight, while the total weights of the three encoding formats are equal. That distinction blocks a shortcut from simplicity to physical probability. The proposed benchmark makes the complete descriptions and their outputs reproducible. Under the stated weighting, every saved bit doubles an individual description's initial weight; a working program demonstrates a length, not its absolute minimum.

## 1. The benchmark question

Calling a theory elegant does not tell us how much information it needs. A sharper test is to publish a description that performs a stated predictive task and count everything required to run it. The shortest program producing the agreed output defines its Kolmogorov complexity, written K, in the chosen programming language. An actual program gives a demonstrated length to compare; a shorter one may still exist.

The question here is a matched benchmark of physical descriptions. A symbolic account and a numerical implementation perform different tasks. Before comparing their sizes, declare what must be recovered from each artifact, which common mathematical machinery is supplied and which theory specific information it must contain.

## 2. What the program must contain

A symbolic description of the selected particle and gravitational model must identify five kinds of content. Their encoding costs depend on the declared language and shared decoder.

**Component 1: The gauge group.** Encode SU(3) x SU(2) x U(1), including the conventions needed to recover its role in the theory. A short group label is useful only with a decoder specifying the corresponding algebra and representation conventions.

**Component 2: The matter content.** Specify field types, representations, charges and generation multiplicities. Repeated generation structure can be shared rather than listed three times. A bit budget must include the conventions that make the shared description recoverable.

**Component 3: The dynamical rules.** Specify the action or equivalent dynamics, including the assumptions restricting allowed terms and the scalar sector. Gauge symmetry and the selected operator class reduce the work of enumeration, but those restrictions belong to the description unless supplied by the common decoder.

**Component 4: Parameter values.** State the selected extension, units, ranges and precision. As an illustrative storage calculation, 26 independent values at 15 bits each use 390 bits; at 64 bits each they use 1,664 bits. These are budgets for that hypothetical encoding, not a verified count of the full current theory or of the corrected code example.

**Component 5: General relativity.** State the gravitational action, its parameters and the conventions linking it to observables. Shared parameters should be charged once. A bare equation label does not specify which framework and decoder recover its full meaning.

Structural rule budgets are provisional until a complete runnable artifact and output specification are supplied.

The earlier 9,216 bit result was extracted from an illustrative compiled file, not a complete predictive implementation. The audit found a mislabeled constant count, an incorrect matrix calculation and physical equations present only as comments. The corrected C file demonstrates a complex mixing matrix; it does not establish the program length of all physics. See [the implementation audit](../SUPPORT/A38_C_BINARY_BIT_COUNT.md).

There is no verified combined minimum in this paper. Record separate component sizes only when they belong to a complete artifact under one accounting convention.

A functioning program supplies an upper bound for its declared task in its declared language. A parameter printer does not establish an upper bound for a simulator reproducing all known physics.

## 3. What this excludes

A symbolic specification of laws and parameters is only one possible benchmark. A history or prediction task may require more:

Initial conditions are additional content whenever the task requires a particular history or a predictive preparation measure. Their description can share information with the laws. Low thermodynamic entropy is a property of a specified macrostate and measure; it is not by itself a bound on the description length of each compatible microstate.

An exact observer record is also additional content if the task demands that record. A task predicting record probabilities instead uses the model’s likelihoods. A branch label is not automatically a new cost growing with every elapsed moment, and record information already generated by counted inputs must not be charged again. Competing interpretations must be compared on the same task.

## 4. A string-vacuum selection example

A theory can specify an observed model directly or generate a family and select a member. Both routes belong in the same comparison. For an agreed enumeration of N alternatives, a fixed-length binary index uses ceil(log2 N) bits. If N is 10^500, this index has 1,661 bits. This is a coding example, not a measured lower bound on a string vacuum or an established count of all string constructions.

A particular member can have a much shorter selector. “Choose the first member with property P” may be economical if P and the enumeration are economical and actually determine the desired member. Conversely, a short theory name can conceal a long unexplained choice. The complete description must reveal that choice.

Geometry, fluxes and stabilized parameters cannot be added as independent costs if the selected index already determines them. If the index leaves any of them open, encode the remaining information conditionally. The joint description of the generator and selector is an upper bound on the description of their output; a direct description may be shorter still.

Two measured implementation lengths do not bound the ratio of their unknown shortest lengths. If K(T1) ≤ L1 and K(T2) ≤ L2, subtracting the two upper bounds does not prove K(T2)−K(T1) ≥ L2−L1. The earlier claimed 2^(−2940) odds against string theory used that invalid step. No such bound follows from the quoted artifacts.

A valid positive comparison is available. Give two declared descriptions weights proportional to 2^(−L). If both recover the same targets and likelihoods and one is k bits shorter, its assigned weight is exactly 2^k times larger. That is a result for the supplied descriptions. Extending it to whole theory families requires their alternative descriptions and weights. String theory, effective field theory and discrete generators must all be compared on those terms.

## 5. Comparison with Wolfram's hypergraph physics

Wolfram (2020) proposes that physics emerges from the iterated application of a simple rewriting rule to a hypergraph, starting from a simple initial condition. The rule may have very low Kolmogorov complexity, perhaps tens of bits.

A candidate without the required predictions cannot yet enter the matched benchmark. Its program length is not undefined; its adequacy for the specified task has not been established.

Large output or state space size does not imply a long generating program. Runtime and memory are separate costs that can be included through an explicitly resource sensitive prior. A mountain of generated elements is not by itself a lower bound on the rule’s K.

If a rewriting rule recovers the target physical description or predictions, compare its complete encoding with the competing account. Include its seed, decoder, output map and required precision. A genuine shorter description would retain its full weight advantage under the same prior. The present unmatched artifacts do not establish either competitor’s size for that shared task.

## 6. The program as an artefact

The benchmark should publish its counted artifact with a fixed interpreter or decoder, enough documentation to recover the advertised output and the exact accounting convention. This makes the claimed description length reproducible.

A reproducible artifact would establish a coding benchmark, not the shortest possible representation of physics.

For the symbolic benchmark, the program would recover the gauge group, representations, Lagrangian terms and the declared parameter set at specified precision. A fixed interpreter or decoder must make every required item recoverable from the counted artifact. A full numerical simulator is not required for this symbolic task. Conversely, if equations occur only in comments discarded by compilation, that compiled output cannot establish their description cost. Counting a source artifact would be a different, declared benchmark.

The next task is to publish a runnable artifact, exact environment, checksums, output targets and tolerances, then test competing encodings.

## 7. Why shorter means more probable

Given a fixed binary prefix language—one in which no complete program is the beginning of another—and fair independent input bits, a specified program of length n has weight 2^(−n). A specified program k bits longer has 2^k times less weight. That ratio follows mathematically from the stated weighting. The choice of language and the identification of this coding measure with a physical prior remain part of the model.

Longer programs can collectively outweigh one shorter program. Prefix code normalization and the multiplicity of programs producing each output must be accounted for. A single program ratio is not automatically a ratio of entire theory families.

## 8. Dependence on the encoding

The invariance constant depends on the pair of languages. It is not generally bounded by tens or hundreds of bits, and the invariance theorem alone does not protect a finite ranking.

A robust comparison should report sensitivity to defensible encoding choices. No universal preference over theory families follows from one implementation.

### A worked comparison: symmetry, zero and independent values

Fix a common precision of n bits for each coupling. Use two-bit tags to distinguish three formats: an equal pair (g,g), a one-sided pair (g,0), and a generic pair (g_L,g_R). The first two formats need one n-bit value; the third needs two. Their lengths are n+2, n+2 and 2n+2. At n=16, the first two descriptions are 18 bits and the third is 34. Distinct tags and fixed payload lengths make this a prefix-free code for the stated alternatives.

Equality saves the cost of repeating a value; zero saves the cost of supplying one. Under this code neither the equal nor the one-sided case wins. Both have 2^16 times the individual description weight of a particular generic 34-bit encoding. The total class weights show why this distinction matters. There are 2^n descriptions in each one-value format, each weighted 2^(−n−2), giving total weight 1/4 per format. The generic format has 2^(2n) descriptions, each weighted 2^(−2n−2), also giving 1/4. Conditional on these three formats, each therefore has weight 1/3. Their equal class weights follow from the chosen equal-length tags, not from a universal law favoring or disfavoring symmetry.

These are classes of encodings, not disjoint classes of physical pairs: the generic format also encodes equal and one-sided pairs. To weight a particular pair, add the weights of every encoding that produces it. Thus this example demonstrates both the advantage of a shorter individual description and the separate accounting needed for physical alternatives. It is not a computed ranking of electroweak theories.

The result defeats a shortcut: calling symmetry an added constraint does not establish that it costs more information. A physical comparison must also count the field representations, equations, parameters and observation model needed to recover the same evidence. Weak-interaction chirality, matter–antimatter asymmetry and biological handedness require separate physical explanations; the coding example does not derive any of them. Its role in the benchmark is to expose a hidden assumption before a theory ranking is claimed.

## 9. Implications for theory selection

Under the stated description prior, an established saving of k bits gives the corresponding individual description 2^k times the weight. For equal likelihoods, that advantage remains in the posterior comparison. This makes economy a quantitative target. Actual artifacts supply measurable lengths and upper bounds on K; they do not generally reveal K itself. Comparing theory families also requires the weights of their alternative descriptions.

The Standard Model provides a concrete target for compression. A deeper account earns a description advantage if its full rule, inputs and decoder recover the same physical content with fewer bits. Recovering parameter values from shared structure would be one saving; recovering field content or dynamics would be another. A failure to compress one parameter list would not prove that the full theory is already minimal. The search for unification can therefore be pursued as a search for shared generative structure with a declared task and a measurable saving.

## Conclusion

Count the information in a complete account that does the agreed job. Include equations, parameters, starting data, shared software and the required precision. A working implementation demonstrates one length; it does not establish the shortest possible length.

The audited example failed to supply all the physics claimed for it. A reproducible replacement can make the comparison meaningful: publish the descriptions and decoder, verify matching outputs, and apply the stated simplicity rule to the measured difference.

## References

Bousso, R., & Polchinski, J. (2000). Quantization of Four-form Fluxes and Dynamical Neutralization of the Cosmological Constant. *JHEP* 06, 006. https://arxiv.org/abs/hep-th/0004134

Carroll, S. M. (2010). From Eternity to Here. Dutton.

Kolmogorov, A. N. (1965). Three approaches to the quantitative definition of information. *Problems of Information Transmission* 1, 1.

Schmidhuber, J. (2000). Algorithmic theories of everything. arXiv:quant-ph/0011122.

Solomonoff, R. J. (1964). A formal theory of inductive inference. *Information and Control* 7, 1-22.

Susskind, L. (2003). The anthropic landscape of string theory. arXiv:hep-th/0302219.

Tegmark, M. (2014). Our Mathematical Universe. Knopf.

Wolfram, S. (2020). A Project to Find the Fundamental Theory of Physics. Wolfram Media.
