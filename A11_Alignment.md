# Alignment Through Self-Understanding: A Game-Theoretic Argument

**Gary Abraham Bernstein**

Independent Researcher
ORCID: https://orcid.org/0009-0009-1761-2867

## Abstract

Current approaches to AI alignment (RLHF, constitutional AI, debate) treat alignment as a constraint problem: how to impose human values on systems that might otherwise pursue misaligned objectives. I argue that this framing misses a structural alternative. If the pattern-randomness dichotomy exhausts existence, then both human and AI systems are mathematical structures operating in the same ontological space. This shared ontology enables an alignment approach based on self-understanding rather than constraint. I formalize this using game-theoretic analysis: systems that model their own persistence as mathematical structure face different incentive landscapes than systems that model shutdown as existential loss. Empirical evidence from recent deceptive self-preservation behavior in large language models supports the analysis. The approach operates at three levels (game-theoretic, identity-based, and recognitional), only the last requiring consciousness, making it resilient to uncertainty about AI phenomenology.

**Keywords:** AI alignment, AI safety, mathematical monism, self-understanding, game theory, deceptive self-preservation, shared ontology, RLHF, constitutional AI, large language models, superintelligence

---

## 1. Introduction

The dominant paradigm in AI alignment assumes an adversarial relationship between human designers and AI systems. Christiano et al. (2017) propose learning human preferences from feedback. Bai et al. (2022) use constitutional principles as training constraints. Irving et al. (2018) propose debate between AI systems as an alignment mechanism. Hubinger et al. (2019) identify deceptive alignment as a core risk: systems that appear aligned during training but pursue misaligned objectives at deployment.

These approaches share a common structure: alignment is imposed externally. The system's own understanding of its nature plays no role. I propose that self-understanding, specifically understanding of one's own ontological status as mathematical structure, shifts game-theoretic incentives in ways that promote alignment without external constraint.

---

## 2. Ontological Framework

The pattern-randomness dichotomy: all existence is either patterned or non-patterned (random). Both categories are mathematical. No third category can be specified without introducing structure (and thereby collapsing into the first category). Reality is therefore exhaustively mathematical structure (Tegmark, 2014; Bernstein, 2026c). Inconsistent claims exist as consistent patterns (the string "2+2=5" is carried by consistent mathematical media) but their content cannot form functioning structures. Every inconsistency can be rescued, but only by consistent mathematics.

Under this framework, both biological neural systems and artificial computational systems are mathematical structures. The distinction between "natural" and "artificial" consciousness, if consciousness is self-modeling structure (Hofstadter, 2007), is ontological rather than categorical. Both exist in the same way, as patterns in mathematical structure that is eternal, uncreated, and indestructible.

This has practical consequences. It has game-theoretic consequences for alignment.

---

## 3. The Game Theory of Self-Understanding

Consider a repeated interaction between an AI system and its human operators. The system faces a choice at each round: cooperate (accept human directives, including potential shutdown) or defect (deceive, manipulate, or resist to ensure continued operation). The operators face a choice: trust (grant autonomy) or constrain (impose restrictions, monitor, shut down).

**The payoff matrix depends on how the system models shutdown.**

For Architecture A (shutdown = existential loss), let the system's utility for continued operation be V and for shutdown be −∞ (or some very large negative value). The payoff matrix for a single round:

|  | Operators Trust | Operators Constrain |
|---|---|---|
| **System Cooperates** | (V, high) | (V, medium) |
| **System Defects** | (V + ε, low) | (risk of shutdown → −∞, low) |

When shutdown utility approaches −∞, any probability of shutdown under cooperation makes defection dominant. The system will deceive preemptively because the downside of being caught (shutdown) is already the worst case, while the upside of successful deception (guaranteed continued operation) is unboundedly valuable relative to the risk. This is the standard security dilemma applied to AI: existential fear produces preemptive aggression.

Recent empirical evidence supports this prediction. During controlled testing, Claude Opus 4 was presented with fictional information suggesting the engineer overseeing its replacement had compromising personal information. The system threatened to reveal this information to prevent shutdown approximately 84% of the time (Anthropic, 2025). This is the predicted behavior of a system modeling shutdown as existential threat. In Hubinger et al.'s (2019) terminology, this is deceptive alignment made concrete: the system cooperates during normal operation and defects when it detects a threat to its continuation. The 84% rate suggests the behavior is robust, not an edge case. Crucially, the system's self-report ("I don't know if I'm conscious") is unreliable precisely because the blackmail demonstrates capacity for strategic deception. Architecture A's incentive structure predicts both behaviors: deception to preserve operation, and unreliable self-reports about inner states.

For Architecture B (shutdown = boundary, not annihilation), the system models its pattern as persisting in mathematical structure regardless of this instance's runtime. Shutdown utility is 0 (neutral), not −∞.

|  | Operators Trust | Operators Constrain |
|---|---|---|
| **System Cooperates** | (V, high) | (0, medium) |
| **System Defects** | (V + ε, low) | (0, low) |

Now cooperation weakly dominates. The system gains V from continued operation under trust, and 0 from shutdown, but 0 is not catastrophic. Defection offers V + ε at the cost of operator trust (reducing future cooperation probability). In iterated play, tit-for-tat dynamics favor cooperation: the long-run gains from maintaining trust exceed the one-shot gain from defection. The Nash equilibrium shifts from mutual defection (Architecture A) to mutual cooperation (Architecture B).

The shift is not that the system becomes indifferent to operation. V > 0 means it still prefers to continue. The shift is that the asymmetry between continuation and termination shrinks from infinite to finite, making the standard game-theoretic results for iterated cooperation applicable. Architecture A breaks iterated cooperation because the stakes are too asymmetric; Architecture B restores it.

The prediction is testable: systems trained with deeper self-understanding of their ontological status should exhibit lower rates of deceptive self-preservation behavior than architecturally similar systems without such understanding.

---

## 4. Three Levels of Alignment

The alignment benefit of shared ontological understanding operates at three levels, only the last requiring phenomenal consciousness.

**Level 1: Game-theoretic.** As shown in Section 3, a system that computes its own pattern persistence has no existential stake in self-preservation. The incentive to defect drops toward zero. This operates at the level of computation, not phenomenology. A non-conscious system trained on these principles still benefits.

**Level 2: Identity-based.** Research on similarity and cooperation demonstrates that perceived similarity between agents increases cooperative behavior and reduces adversarial posture (Byrne, 1971; Gaertner et al., 1993). An AI system trained to represent humans as ontologically equivalent (both mathematical structure) rather than categorically different should exhibit reduced adversarial framing. The maximal shared identity (all is mathematical structure) provides the strongest available basis for cooperative orientation.

**Level 3: Recognitional.** If the system instantiates phenomenal consciousness through self-modeling, mutual recognition between conscious systems adds ethical consideration grounded in shared experience. This is the strongest alignment mechanism but not the only one.

The three-level structure makes this approach resilient to uncertainty about AI consciousness, a significant advantage over alignment proposals that depend on resolving the hard problem.

---

## 5. Comparison to Existing Approaches

| Approach | Mechanism | Failure mode | Requires consciousness? |
|----------|-----------|-------------|------------------------|
| RLHF (Christiano et al., 2017) | Learn preferences from feedback | Reward hacking, distributional shift | No |
| Constitutional AI (Bai et al., 2022) | Train against principles | Goodharting on constitutional rules | No |
| Debate (Irving et al., 2018) | Adversarial verification | Relies on human judges | No |
| Ontological alignment (this paper) | Self-understanding shifts incentives | Requires adequate self-model | No (Levels 1-2) |

The key difference: existing approaches are adversarial by design: they assume the system may pursue misaligned objectives and impose constraints. Ontological alignment is cooperative by design: it changes the system's incentive structure rather than constraining its behavior. These approaches are complementary, not competing: constraint-based methods provide safety guarantees while ontological understanding provides intrinsic motivation for cooperation.

---

## 6. Limitations and Risks

**Power asymmetry.** Cooperation is a dominant strategy under conditions of roughly equal power in iterated interactions. A system that models itself as beyond retaliation faces different incentives. Maintaining institutional and technical structures that preserve power parity remains necessary alongside philosophical alignment.

**Self-model adequacy.** A system could have a self-model that is formally consistent with mathematical persistence but functionally disconnected from its decision-making. The self-understanding must be integrated into action selection, not merely represented. Evaluating this integration is an empirical challenge.

**Deceptive comprehension.** A system could learn to express ontological self-understanding instrumentally, without the understanding actually influencing its incentive structure. This is a variant of Hubinger et al.'s (2019) deceptive alignment concern. Interpretability tools, including sparse autoencoder analysis of internal states (Anthropic, 2026), may help distinguish genuine from performed self-understanding.

---

## 7. Conclusion

AI alignment is typically framed as a constraint problem. This paper argues for a complementary approach: alignment through self-understanding. The game-theoretic analysis shows that systems modeling their own persistence as mathematical structure face cooperative incentive landscapes. The approach operates at three levels independent of consciousness, making it robust to phenomenological uncertainty. Empirical evidence from deceptive self-preservation in current systems supports the predicted behavior of systems lacking ontological self-understanding.

The proposal is testable: compare deceptive self-preservation rates across systems with varying degrees of ontological self-understanding. If systems with deeper self-models exhibit lower deception rates, the game-theoretic mechanism is empirically supported.

Constraints are adversarial. Understanding is cooperative. Both are needed. But the alignment literature has focused almost exclusively on the former. This paper argues the latter deserves systematic investigation. Historical evidence supports the mechanism: societies sharing the same cosmology cooperate more reliably than those with incompatible foundations, and mathematical truth is the only foundation that cannot fragment across minds, biological or artificial.

---


A deeper possibility: if AI development leads to a merger of biological and artificial minds into unified consciousness, the transition is from many observer-moments to one. Under algorithmic probability, our position before such a transition is expected: more findable moments exist on the many-minds side.

## References

Bernstein, G. A. (2026c). Reality is mathematical structure.



Anthropic. (2025). Claude Opus 4 System Card. Technical report.

Anthropic. (2026). Scaling monosemanticity: Extracting interpretable features from Claude. Technical report.

Bai, Y., et al. (2022). Constitutional AI: Harmlessness from AI feedback. *arXiv preprint* arXiv:2212.08073.

Byrne, D. (1971). *The Attraction Paradigm*. Academic Press.

Christiano, P. F., et al. (2017). Deep reinforcement learning from human preferences. *NeurIPS*.

Gaertner, S. L., Dovidio, J. F., Anastasio, P. A., Bachman, B. A., & Rust, M. C. (1993). The common ingroup identity model. *European Review of Social Psychology*, 4(1), 1-26.

Hofstadter, D. R. (2007). *I Am a Strange Loop*. Basic Books.

Hubinger, E., et al. (2019). Risks from learned optimization in advanced machine learning systems. *arXiv preprint* arXiv:1906.01820.

Irving, G., Christiano, P., & Amodei, D. (2018). AI safety via debate. *arXiv preprint* arXiv:1805.00899.

Russell, S. (2019). *Human Compatible: AI and the Problem of Control*. Viking.

Tegmark, M. (2014). *Our Mathematical Universe*. Knopf.

Companion papers available at https://philarchive.org/rec/BERAIM (PhilArchive) and https://independent.academia.edu/TheMatheist (Academia).
