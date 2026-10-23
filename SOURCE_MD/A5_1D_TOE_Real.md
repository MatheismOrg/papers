# Quantum Gravity in Fewer Effective Dimensions: What the Reductions Explain

## Abstract

Gravity and quantum theory have to fit one account, because reality is one mathematical structure, and several leading approaches to quantum gravity share a clue about how: at the smallest scales, spacetime behaves as if it had two dimensions. The reduction does not remove a coordinate. It comes from how the theory scales, measured by the spectral dimension, the dimension a diffusing particle would detect, and a worked example shows four coordinates giving spectral dimension two. Fewer dimensions is not automatically simpler. A lower-dimensional theory must carry the map back to the four dimensions we measure, as holography's boundary dictionary does, and the account is shorter only if theory and map together are. A 2026 theorem by Gianluca Calcagni and Fabio Briscese shows, under stated assumptions, that quantum probabilities stay consistent in a class of fractional field theories with this behavior. Seven approaches are compared by the physics they recover and the ingredients they still need.

Keywords: mathematical monism; spectral dimension; quantum gravity; unitarity; theory selection

## 1. One ontology, several dynamics

One ontology, a single account of what exists, does not supply one equation. The equations x-dot = x and x-dot = 2x are both mathematical, use the same kind of variable, and have different solutions. Reality is one mathematical structure, so gravity and quantum physics belong to the same account, and the common dynamics still has to be found.

Matheism's interaction argument asks how any additional substance could affect the physical world without definite relations. Together with the claim that complete structure exhausts physical properties, these arguments motivate one ontology. The metaphysical unification can therefore be argued a priori, by reasoning alone, while a derivation connecting general relativity and quantum theory still matters.

## 2. What is becoming lower dimensional?

The word dimension does several jobs here, and only some of those counts change at short distance. Steven Carlip's review (Carlip, 2017) separates the notions and collects the evidence for short-distance dimensional reduction across several quantum-gravity programs. Coordinate dimension, spectral dimension (the dimension a randomly wandering particle would detect), the dimension of a graph edge (a link in a network) and the dimension of a string worldsheet (the surface a string sweeps out over time) answer different questions.

### 2.1 A worked example

Change one scaling exponent and a theory with four coordinates measures as two-dimensional. Stated precisely, take a positive Euclidean operator with momentum dependence |k|^(2z), where z > 0, on d coordinate dimensions. Its heat-kernel trace per unit volume, read as a function of a diffusion scale s, has the form

$$P(s)=C\int_{\mathbb{R}^d}d^dk\,e^{-s|k|^{2z}}.$$

Rescale k = s^(-1/(2z))q. For s > 0 the radial integral converges and

$$P(s)=C' s^{-d/(2z)}.$$

Spectral dimension reads a dimension off the way that trace scales with s. Written out, it follows that

$$D_s=-2\frac{d\log P}{d\log s}=\frac{d}{z}.$$

With d = 4 and z = 2 it equals two although the integration still uses four coordinates. A kinetic operator (the part of the equations governing how a field propagates) that crosses from |k|^2 at low momentum to |k|^4 at high momentum can have different effective spectral dimensions at different diffusion scales. For the fourth-order operator, the off-diagonal heat kernel need not be nonnegative: P(s) is a spectral diagnostic, not automatically a classical return probability (Calcagni et al., 2016). Nor does this calculation establish Lorentzian unitarity (total probability conserved in physical spacetime).

### 2.2 Approaches to quantum gravity

| Approach | Relevant structure | What must be distinguished |
|---|---|---|
| String theory | One-dimensional strings sweeping two-dimensional worldsheets | String and worldsheet dimensions versus spacetime and spectral dimensions |
| Loop quantum gravity | Spin networks and geometric operators | Graph edges, area and volume spectra, and recovered spacetime geometry |
| Causal dynamical triangulations | A sum over triangulated geometries (Loll, 2020) | The lattice regulator and the proposed scaling limit |
| Asymptotic safety | Renormalization-group behavior | A fixed-point proposal and an effective dimension inferred from it |
| Horava-Lifshitz gravity | Anisotropic scaling of space and time | Scaling exponents versus the number of coordinates |
| Causal sets | A locally finite causal order | Discrete order, volume counting and manifold approximation |
| Fractional field theory | Fractional kinetic operators | Operator scaling, physical spectrum and amplitude prescription |

## 3. The recent theorem

Gianluca Calcagni and Fabio Briscese (2026) proved that specified fractional scalar, gauge and gravitational theories (whose equations use non-integer powers) keep unitarity, the rule that the probabilities of all outcomes add up to one, at every order of perturbation theory, the method of successive small corrections. The theorem covers that class of models and leaves the other programs untouched. Its results depend on the kinetic operator, exponent range and amplitude prescription. The hermitian construction uses the fakeon prescription, which removes selected modes from the physical asymptotic states (the particles that enter and leave interactions); the non-hermitian case has stated ranges compatible with the Feynman prescription, the standard rule in ordinary particle physics. The article is published as JHEP09(2026)206.

This establishes a consistency result within the analyzed class, not for every program in the table. The fakeon prescription gives up microcausality, the ban on faster-than-light influence, at short distances, as its own literature states (Anselmi, 2019). A program can now investigate that class without treating its previously unresolved unitarity issue as automatic grounds for dismissal. Ultraviolet behavior (how a theory acts at very high energies), physical limits and phenomenology (the link to experimental data) must still be answered for the chosen model and prescription.

## 4. Numbers specify particles; a parameter list is not a particle theory

The Standard Model, the established theory of the known particles and forces except gravity, is more than its numerical parameters. It specifies fields, transformation laws, interactions and a state space. An electron is not a scalar field, one with a single number at each point, and a gauge boson, a force-carrying particle such as the photon, does not become a scalar particle because its properties can be encoded by numbers. The physical content is the whole specification, and dropping the structures that give the parameters their meaning does not compress it to a parameter list.

## 5. What can be excluded a priori?

An explicit premise can exclude a candidate from a research program. A contradiction in the candidate's formalism can exclude that formulation. These are different exclusions.

Matheism holds that physical reality is discrete at its base, since physical succession generates each state from the one before and a continuum has no next state, an a priori proof in the logician's sense. A fundamental continuum ontology therefore falls outside its program: Russell's motion, only being at different places at different times (Russell, 1903), is the emperor's new points, a position for every instant with nothing generating the next. A lattice theory (one defined on a grid of points) that takes a continuum scaling limit (letting the grid spacing shrink to zero) may conflict with the premise at its proposed fundamental endpoint while still supplying useful effective calculations. As analysis, a description of many steps at once, the limit stays good mathematics. Whether a discrete alternative can recover the same successful limit becomes a constructive question.

The argument that zero-measure points cannot yield positive length does not itself prove the premise. Countable additivity (lengths of separate pieces adding up, even across an endless numbered list) does not extend to arbitrary uncountable unions (combining more pieces than can be numbered, such as a line's points). Using that requirement as a selection principle is stronger than attributing a contradiction to measure theory, the mathematics of assigning sizes to sets. Rejecting fundamental stochastic (chance-governed) outcomes rules out collapse laws (under which a measured quantum state jumps to one outcome) with those outcomes but does not select superdeterminism (the view that measurement settings correlate with what is measured) by logic alone. John Bell's statistical independence condition (Bell, 1964) concerns the joint distribution of hidden variables (unobserved properties underlying the outcomes) and settings, not whether either is produced by a deterministic algorithm.

Everett's unitary branching (the many-worlds picture, in which every outcome occurs) and a local single-outcome hidden-variable model relax different assumptions in familiar Bell arguments, proofs of limits on what local hidden-variable models can predict. They can be compared without treating their conjunction as forced by the block universe, the view that past, present and future equally exist. A short rule for digits of pi establishes computability, not normality (every digit string of a given length appearing equally often), cryptographic pseudorandomness (output no practical test can tell from true randomness) or a distribution matching all quantum experiments.

## 6. Parsimony with complete candidates

Fewer dimensions is not by itself a shorter description. The dimensions we measure carry a definite set of relations, and a theory set in fewer dimensions has to give them back: its description includes the map from its own variables to the four-dimensional world. Holography is the clearest case. In anti-de Sitter space (a model universe with saddle-like curvature), a quantum theory on the boundary, one dimension lower, reproduces the gravity inside (Maldacena, 1998), and the dictionary between the two is part of the account. The lower-dimensional account is simpler only if it and its map together are shorter than the four-dimensional one, and for our own universe no such boundary account is known.

A shorter description earns stated odds once both candidates cover the same observations, because each bit of a description is a 50/50 choice that no rule makes. Fix a prefix description language U, one in which no code word begins another, and specify candidates to the same observational precision. Stated precisely, under a shortest-description prior (a weighting set before any evidence) proportional to 2^(-K_U), an established difference of Delta bits gives an individual prior ratio of 2^Delta in favor of the shorter candidate.

A huge number of alternatives does not make our own expensive to specify. Counting 10^500 alternatives gives a fixed-length indexing cost of 1,661 bits. It does not prove that the particular alternative occupied by us needs such an index: a short rule could select it. Conversely, a short list of constants is not a complete competing account if interactions, geometry or boundary data have been left outside the charge. Mathematical monism, the view that reality is one mathematical structure, motivates this search without presupposing its winner.

## 7. Conclusion

A theory can behave as if it has fewer dimensions without losing coordinate directions, because its effective dimension comes from how its kernel scales, not from how many coordinates it uses. Fewer dimensions earns simplicity only when the reduced account and its map back to the measured world are shorter together. The probability-consistency theorem examined here establishes a further result for its specified class of models. To complete the search, a candidate must derive the relevant particles and interactions, recover the established limits and explain the observations. Mathematical monism gives that search one reality to explain.

## References

Anselmi, D. (2019). Fakeons, microcausality and the classical limit of quantum gravity. *Classical and Quantum Gravity*, 36, 065010. https://arxiv.org/abs/1809.05037

Bell, J. S. (1964). On the Einstein Podolsky Rosen paradox. *Physics*, 1, 195-200. https://doi.org/10.1103/PhysicsPhysiqueFizika.1.195

Calcagni, G., and Briscese, F. (2026). Perturbative unitarity of fractional field theories and gravity. *Journal of High Energy Physics*, 09, 206. https://arxiv.org/abs/2603.25709

Calcagni, G., Modesto, L., and Nardelli, G. (2016). Quantum spectral dimension in quantum field theory. *International Journal of Modern Physics D*, 25, 1650058. https://doi.org/10.1142/S0218271816500589

Carlip, S. (2017). Dimension and dimensional reduction in quantum gravity. *Classical and Quantum Gravity*, 34, 193001. https://arxiv.org/abs/1705.05417

Loll, R. (2020). Quantum gravity from causal dynamical triangulations: A review. *Classical and Quantum Gravity*, 37, 013002. https://arxiv.org/abs/1905.08669

Maldacena, J. (1998). The large N limit of superconformal field theories and supergravity. *Advances in Theoretical and Mathematical Physics*, 2, 231-252. https://arxiv.org/abs/hep-th/9711200

Russell, B. (1903). *The Principles of Mathematics*. Cambridge University Press.
