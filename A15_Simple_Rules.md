# Simple Rewriting Rules Are Exponentially Favored as Fundamental Physics

**Gary Abraham Bernstein**

Independent Researcher
ORCID: https://orcid.org/0009-0009-1761-2867

## Abstract

Under the Simplicity Dominance Principle (SDP), where a structure's measure is P(S) ∝ 2^(−K(S)), the most probable fundamental physics is the one with the shortest generating description. This paper argues that simple rewriting rules, such as those explored by Wolfram's Physics Project, are structurally the lowest-K candidates for fundamental physics. The reason is architectural: a rewriting rule collapses framework, parameters, and dynamics into a single object. Traditional physics separates framework (quantum field theory, general relativity), parameters (coupling constants, masses), and dynamics (equations of motion). Each separation adds specification cost. A rewriting rule has no such separation: the rule IS the principle, the parameters, and the dynamics. This architectural collapse is why rewriting rules can achieve total K in the range of tens of bits, while the Standard Model plus General Relativity requires roughly 190 bits of framework plus 20 to 27 free parameter values. The paper does not claim that a specific rule has been found. It argues that The principle predicts the fundamental theory has this architectural structure, and that the search should proceed by enumeration and mathematical derivation of simple rules rather than by constructing larger frameworks.

**Keywords:** rewriting rules, hypergraphs, Simplicity Dominance Principle, Kolmogorov complexity, Wolfram Physics Project, Theory of Everything

---

## 1. Why Architecture Matters for K

Two descriptions of the same physics can have vastly different Kolmogorov complexity depending on their architecture.

**Separated architecture:** Specify a principle (quantum field theory: path integral, gauge principle, renormalization). Specify parameters (coupling constants, masses). Specify initial conditions. Each component is a separate specification requiring separate bits. Total K is the sum.

**Collapsed architecture:** Specify a single rule. The rule generates spacetime, forces, particles, and constants as emergent properties of its evolution. No separate framework. No separate parameters. Total K is the length of the rule.

Under SDP, the collapsed architecture is exponentially favored whenever a short rule can reproduce the same physics. The question is whether such rules exist.

## 2. Rewriting Rules as Collapsed Architecture

A rewriting rule operates on a discrete structure (a graph, hypergraph, or string) by pattern-matching and replacement. A typical hypergraph rule:

{{x, y, z}} → {{x, y, w}, {y, w, z}}

This specifies:
- Input pattern: one hyperedge with three elements
- Output pattern: two hyperedges with specific element sharing
- Update rule: replace input with output wherever the pattern matches

The entire specification is the rule. No separate framework: the rule defines the dynamics. No separate parameters: the rule's structure determines all coupling constants, masses, and force strengths through its emergent behavior. No separate initial conditions: the simplest starting state (a small hypergraph) has low K and is SDP-favored.

The K of a single rule is determined by:
- Number of input hyperedges (1-2: ~1-2 bits)
- Number of output hyperedges (1-3: ~2 bits)
- Elements per hyperedge (2-4: ~2 bits)
- Connection pattern (which elements shared, which new: ~10-15 bits)

Total: approximately 15-25 bits for a single rule. If multiple rules are needed, approximately 30-75 bits.

## 3. Why Rewriting Rules Can Be Fundamental

Wolfram's Physics Project has demonstrated that simple rewriting rules can produce:
- Emergent dimensionality (the large-scale structure of the hypergraph approximates a manifold)
- Emergent curvature (consistent with Einstein's equations in the continuum limit)
- Causal structure (the causal graph of updates produces a partial order resembling spacetime)
- Quantum-like behavior (multiway branching of rule application produces structural analogs of superposition). At each step, the rule may apply in multiple places. Applying all simultaneously creates branching: the multiway system. All branches coexist. Branchial space, the space of simultaneous branches, maps structurally to quantum state space. When branches reconverge, they interfere: reinforcing branches survive, canceling branches vanish. This is the path integral derived from the rule rather than postulated. The classical path that survives is the Principle of Least Action, itself a consequence of the many-to-one structure of coarse-graining (Bernstein, 2026k). Whether branchial space is rigorously isomorphic to Hilbert space remains open. The strongest result is causal invariance producing general covariance. The multiway-to-quantum correspondence is the open frontier.

None of these properties is specified. All emerge from the rule. This is the architectural collapse in action: one short specification generates what separated frameworks require hundreds of bits to describe.

The unresolved question is whether any specific rule reproduces the Standard Model. General Relativity has been partially recovered in the continuum limit. Quantum mechanics appears through multiway evolution. The gauge structure, particle spectrum, and specific parameter values have not been derived from any known rule. This is the gap between structural prediction and specific verification.

## 4. The Computational Gap

If fundamental rewriting rules operate at scales around 10^(-93) meters (as Wolfram estimates from dimensional analysis), the gap between the rule and observable physics spans roughly 10^58 orders of magnitude below the Planck length. Direct simulation is computationally intractable: evolving a hypergraph from tens of elements to Standard Model physics requires a number of steps exceeding any feasible computation.

Mathematical derivation rather than brute-force simulation is required. This means proving that certain classes of rules necessarily produce certain emergent properties (dimensionality, gauge symmetry, specific particle content) without simulating every intermediate step.

This gap does not weaken the SDP argument. The principle predicts that the fundamental theory has low K. It does not predict that verification is easy. A rule with K = 25 bits may require 10^400 computational steps to verify, but its measure is still 2^(-25), which exponentially dominates any candidate with K = 200. The difficulty of verification is independent of the probability of truth.

## 5. Evidence from Enumeration

The Wolfram Physics Project has enumerated simple hypergraph rules and found several that produce physics-like behavior. Specific results include:

**Dimension emergence.** Rules of the form {{x,y,z}} → {{x,y,w},{y,w,z}} applied to minimal initial hypergraphs produce structures whose geodesic scaling matches 3+1 dimensional spacetime. The dimension is not specified in the rule; it emerges from the connectivity pattern. Different rules produce different effective dimensions, and 3+1 is a common attractor among simple rules (Wolfram, 2020).

**Causal invariance.** Some rules exhibit causal invariance: the causal graph is the same regardless of the order in which rewrites are applied. This property, when it holds, automatically produces discrete analogs of general covariance. The rule does not "know about" relativity; relativity emerges from the rule's structure.

**Energy-momentum conservation.** In causally invariant systems, conservation laws emerge from the symmetries of the causal graph, analogous to Noether's theorem in continuous physics.

**Quantum mechanics.** The multiway graph of all possible rewrite orderings produces branching structures analogous to quantum superposition. Branchial space, the space of simultaneous rewrite possibilities, may correspond to quantum state space. This remains conjectural but provides a concrete research direction.

No specific rule has been proven to reproduce the full Standard Model. This is the principal limitation of the program and must be stated honestly. The claim is not that a specific rule has been found. The claim is that The principle predicts the fundamental theory has this architectural structure, and that the enumeration results provide preliminary evidence that physics-like behavior emerges generically from simple rules, not only from carefully constructed ones.

## 6. Objections

**"No specific rule has been found."** Correct. The program is predictive, not yet confirmatory. The principle predicts the architecture; finding the specific rule is the remaining task. This is analogous to predicting that the TOE has low K without knowing the TOE: the prediction constrains the search space even before the search succeeds.

**"Continuous physics cannot emerge from discrete rules."** Continuous behavior emerges from discrete substrates routinely: fluid dynamics from molecules, thermodynamics from particles, smooth geometry from triangulations in Regge calculus. The question is whether the continuum limit of the correct discrete rule reproduces known physics. This is an open empirical question, not a principled impossibility.

**"Wolfram's project lacks peer-reviewed results."** Partially true as of 2026, though papers by Gorard (2020, 2021) on causal invariance and dimension emergence have appeared in refereed proceedings. The sociological status of the project does not affect the information-theoretic argument: if a 30-bit rule reproduces SM physics, its K is 30 bits regardless of who found it.

## 7. Comparison to Other Architectures

| Architecture | Framework | Parameters | Total K |
|---|---|---|---|
| Rewriting rules | Rule IS framework | Emerge from rule | ~15-75 bits |
| SM + GR | ~190 bits | 20-27 values | ~190 + values |
| Geometric Unity (14D) | ~190 + bundle | Emerge from geometry | ~190 + low |
| String theory | ~155 bits | ≥181 values | ~155 + values |

Weinstein's Geometric Unity proposes unification in a 14-dimensional observerse (the metric bundle over 4D spacetime). The insight is correct: gauge structure and gravity share geometric origin. The 14D follows automatically from having a metric on 4D, adding near-zero specification beyond the base manifold. This makes GU potentially the SDP-favored description of unification at the continuum level: the simplest way to describe the emergent physics between the fundamental rule and the Standard Model. Hypergraph rules generate the 4D spacetime from below. GU's metric bundle is what that spacetime necessarily contains. Both descriptions may be correct at their respective scales. Full K-analysis is given in Bernstein (2026m).

The architectural advantage of rewriting rules is not a matter of degree. It is structural. Traditional physics pays K for framework AND parameters AND their interface. Rewriting rules pay K once: for the rule. Everything else computes.

This architectural collapse is also why "simplest possible discrete rules" (SDP's prediction) and "hypergraph rewriting rules" (Wolfram's proposal) converge. Any candidate that collapses framework, parameters, and dynamics into one specification has this architectural advantage. Hypergraphs are the best-developed example, but the argument applies to any collapsed architecture: cellular automata (Zuse, 1969), causal sets (Bombelli et al., 1987), or other discrete computational substrates.

## 8. What SDP Predicts

SDP does not predict that Wolfram's specific program will succeed. It predicts:

1. The fundamental theory has collapsed architecture (framework + parameters = one object)
2. Its K is low, in the range achievable by short rewriting rules
3. It is a discrete computational rule, because discrete rules have lower K than continuous differential equations for the same output complexity
4. The specific rule is discoverable by systematic enumeration of low-K candidates combined with mathematical derivation of emergent properties
5. Higher-K candidates (string theory, extended frameworks with many parameters) are exponentially suppressed regardless of their mathematical elegance

The search program follows: enumerate simple rules, derive their emergent properties mathematically, test against observed physics. The rule that reproduces the Standard Model with the lowest K is the most probable fundamental theory.

## 9. Conclusion

Simple rewriting rules are exponentially favored as fundamental physics because they achieve the lowest possible specification complexity through architectural collapse. Where traditional physics separates framework, parameters, and dynamics, a rewriting rule unifies them into a single object with K in the tens of bits. The principle predicts this architecture. Wolfram's Physics Project is the most developed research program exploring it. Whether a specific rule reproduces our physics remains open. That it should is an SDP prediction. Higher K is exponentially suppressed. The search should proceed toward the simplest possible rules.

---

## References

Bernstein, G. A. (2026k). Why nature minimizes: The many-to-one root of time, action, and emergence.

Bernstein, G. A. (2026m). Geometric Unity as emergent variational geometry: Why 14 dimensions cost nothing.



Bombelli, L., Lee, J., Meyer, D., & Sorkin, R. D. (1987). Space-time as a causal set. *Physical Review Letters*, 59(5), 521-524.

Wolfram, S. (2020). *A Project to Find the Fundamental Theory of Physics*. Wolfram Media.

Zuse, K. (1969). *Rechnender Raum* (Calculating Space). Friedrich Vieweg & Sohn.

Companion papers available at https://philarchive.org/rec/BERAIM (PhilArchive) and https://independent.academia.edu/TheMatheist (Academia).
