# Why These Simple Laws, Forward Time, Many-Worlds, and Superdeterminism

**Gary Abraham Bernstein**
Independent Researcher
ORCID: 0009-0009-1761-2867

## Abstract

If reality is mathematical structure (Bernstein, 2026c), all consistent structures exist. This paper derives the physical consequences. Algorithmic probability (Solomonoff, 1964) provides the natural measure over structure-space: simpler structures are exponentially more probable. Each additional bit of specification halves a structure's measure. From this single principle, four results follow. Simple laws: the laws we observe are compact because compact descriptions dominate the measure. Forward time: forward-described structures have lower complexity than backward-described ones. Many-Worlds: MWI adds zero specification cost per measurement while collapse adds unbounded cost. Superdeterminism: deterministic initial conditions have finite K while stochastic ones accumulate K per event. The fine-tuning problem dissolves: the constants are not arbitrary but among the simplest compatible with observers.

**Keywords:** algorithmic probability, mathematical structure, many-worlds, superdeterminism, arrow of time, fine-tuning, Kolmogorov complexity

---

## 1. Introduction

A companion paper (Bernstein, 2026c) argues that reality is identical to mathematical structure through the pattern-randomness dichotomy. This paper does not re-derive that conclusion. It takes it as given and asks: what follows for physics?

If all consistent mathematical structures exist, the question "why these laws?" transforms. It is no longer about what selected our physics from alternatives. It is about what measure weights one structure over another. The Level IV multiverse (Tegmark, 2014) follows from the null filter argument: any filter selecting among consistent structures is itself a consistent structure, so all filters exist, including the null filter that excludes nothing.

But a multiverse without a measure makes no predictions. The measure problem is the central challenge. This paper proposes algorithmic probability as the answer.

## 2. Algorithmic Probability as the Natural Measure

Under Solomonoff's (1964) universal prior, the probability of structure S is P(S) proportional to 2^(-K(S)), where K(S) is its Kolmogorov complexity. Simpler structures dominate exponentially. This is not one arbitrary choice among many. Solomonoff showed it dominates all computable alternatives. Li and Vitanyi (2019) proved it is the unique measure with this dominance property.

Since reality is structure, the measure is not imposed externally. It is discovered: algorithmic probability is the natural measure over structure-space. Occam's Razor becomes a theorem rather than a heuristic.

A structure requiring 100 bits of specification is 2^100 times less probable than one requiring none. This is why the laws of physics fit on a page: compact laws have maximal measure. Complicated laws exist in the multiverse but are exponentially rare among observer-containing structures.

### Limitations

The measure has known open questions. The choice of universal Turing machine introduces a machine-dependent constant. For structures whose complexity differs by hundreds or thousands of bits, this constant is negligible. For structures of similar complexity, it is not. The observation-measure problem, whether we weight by structure-count or observer-count, remains open (Bostrom, 2002). These limitations are real. They do not affect the qualitative conclusions.

## 3. Simple Laws

Why are the laws of physics simple? Under algorithmic probability, the answer is immediate: simple laws have low K and dominate the measure. A universe governed by 20 parameters has exponentially higher measure than one governed by 2,000. Among observer-permitting structures, the simplest ones overwhelm.

This predicts that outstanding puzzles will have simple solutions. Dark matter will be explained by lower-K descriptions over higher-K alternatives. The fundamental constants are not arbitrary but among the simplest values compatible with observers. String theory, with 3,000 to 4,700 bits of specification, is exponentially disfavored relative to the Standard Model's roughly 60 bits (Bernstein, 2026g).

### Fine-tuning resolved

The fine-tuning problem asks: why are the constants so precisely right for life? Under algorithmic probability, the question dissolves. All consistent structures exist. Observers find themselves in observer-permitting ones. Among those, algorithmic probability selects the simplest. The constants are not tuned by a designer. They are the simplest values that permit observers, and simplest values have the highest measure.

This is not the anthropic principle. The anthropic principle is a selection effect within a probability distribution. The present argument derives the distribution itself.

## 4. The Arrow of Time

Why does time have a direction? Under algorithmic probability, a universe described forward from initial conditions I_0 has total complexity K(L) + K(I_0), both low. Described backward from final state I_f, the total is K(L) + K(I_f), where K(I_f) may be arbitrarily large because specifying which of infinitely many inputs produced each output requires unbounded additional information.

Forward-described structures dominate the measure. The Past Hypothesis (Albert, 2000), that the universe began in a low-entropy state, follows as theorem rather than brute postulate.

A deeper structural argument grounds this: computation is many-to-one, the inverse function does not exist, and this irreversibility is the foundation algorithmic probability quantifies. Full development in Bernstein (2026h).

## 5. Many-Worlds Wins on Parsimony

The quantum measurement problem asks: what happens when a quantum system is measured? Three families of interpretation compete.

**Collapse interpretations** (Copenhagen, GRW, objective collapse) add a physical mechanism that reduces the wavefunction to one outcome. Each measurement event requires specifying which outcome occurs. For N measurements, the accumulated specification cost is K(O_1) + K(O_2) + ... + K(O_N). This grows without bound.

**Many-Worlds** (Everett, 1957) adds nothing. The wavefunction evolves unitarily. All branches exist. No outcome is selected because all outcomes occur. The specification cost per measurement is zero. K(MWI) = K(laws) + K(initial state). It does not grow with the number of measurements.

**Superdeterminism** adds nothing per measurement either. All outcomes are determined by initial conditions. K(SD) = K(laws) + K(initial conditions). Like MWI, it does not grow.

Under algorithmic probability, interpretations that accumulate specification cost are exponentially disfavored relative to those that do not. Collapse loses. MWI and superdeterminism survive. The choice between MWI and superdeterminism is not settled by parsimony alone: both have constant K.

### Bell's theorem as corroboration

Bell's theorem proves that quantum correlations cannot be explained by local hidden variables. Under mathematical monism, Bell violations are expected: entangled particles are one mathematical pattern, and patterns correlate perfectly across any spatial separation because space is structure within math. The correlations are not action at a distance. They are arithmetic.

### The Born Rule

MWI must recover the Born Rule: probabilities proportional to squared amplitudes. Under algorithmic probability, branches with simpler descriptions have higher measure. If the Born Rule is the unique probability assignment consistent with the structure of quantum mechanics (as Gleason's theorem establishes for Hilbert spaces of dimension greater than 2), then algorithmic probability and the Born Rule converge. The derivation is not complete, but the alignment is structural, not coincidental.

## 6. Superdeterminism Rehabilitated

The standard objection to superdeterminism is that reproducing quantum predictions requires "fine-tuned" initial conditions, specific correlations between hidden variables and measurement settings.

This objection conflates specificity with complexity. The correlations are specific (they must be exactly right), but specific does not mean complex. The digits of pi are specific to the last decimal but have Kolmogorov complexity of approximately zero: a short formula generates all of them. Under algorithmic probability, what matters is not whether the initial conditions are specific but whether they follow from a simple rule. If they do, K(initial conditions) is low and the theory is measure-favored regardless of how specific the outputs are.

The fine-tuning objection to superdeterminism is therefore a category error (Bernstein, 2026i). Full development in the companion paper.

## 7. Objections

**"This is untestable metaphysics."** The framework makes testable predictions: simpler explanations will be confirmed over complex ones. Dark matter will have a simple explanation. String theory will not be confirmed. These are falsifiable.

**"The measure is ad hoc."** Algorithmic probability is not ad hoc. It is the unique measure dominating all computable alternatives. It is not one arbitrary choice among many.

**"The measure problem is unsolvable."** The observation-measure problem (structure-count vs observer-count weighting) is genuinely open. But the qualitative prediction, simpler structures dominate, is robust across all reasonable measure choices.

**"The core premise is circular."** The argument does not assume mathematical structure to prove mathematical structure. It shows that non-mathematical existence is incoherent via the PRD (Bernstein, 2026c).

**"Occam's Razor is just a heuristic."** Under algorithmic probability, Occam's Razor is a theorem. Simpler descriptions have higher measure. The razor cuts because mathematics forces it to.

## 8. Conclusion

If reality is mathematical structure, algorithmic probability provides the natural measure over structure-space. From this single principle: simple laws dominate, the arrow of time points from simple to complex, Many-Worlds and superdeterminism survive while collapse does not, fine-tuning dissolves, and Occam's Razor becomes a theorem.

The framework is not complete. The observation-measure problem is open. The Born Rule derivation is incomplete. The exact K of the Standard Model's constants is unknown. These are research programs, not refutations.

What is established: the direction. Simpler structures dominate. The theory of everything must have lower total K than current physics, not higher. The theory that explains the most with the least specification wins.

---

## References

Albert, D. Z. (2000). *Time and Chance*. Harvard University Press.

Bernstein, G. A. (2026c). Reality is mathematical structure.

Bernstein, G. A. (2026g). String theory is complexity-disfavored: A description-length bound.

Bernstein, G. A. (2026h). The arrow of time is irreversible computation.

Bernstein, G. A. (2026i). The fine-tuning objection to superdeterminism conflates specificity with complexity.

Bostrom, N. (2002). *Anthropic Bias*. Routledge.

Everett, H. (1957). Relative state formulation of quantum mechanics. *Reviews of Modern Physics*, 29(3), 454-462.

Li, M., & Vitanyi, P. (2019). *An Introduction to Kolmogorov Complexity and Its Applications* (4th ed.). Springer.

Schmidhuber, J. (2000). Algorithmic theories of everything. arXiv:quant-ph/0011122.

Solomonoff, R. J. (1964). A formal theory of inductive inference. *Information and Control*, 7(1), 1-22.

Tegmark, M. (2014). *Our Mathematical Universe*. Knopf.

All companion papers available at https://matheism.org
