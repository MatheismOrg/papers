# **Simplicity Preference Is Adaptive**

**Gary Abraham Bernstein**

Independent Researcher ORCID: https://orcid.org/0009-0009-1761-2867

## **Abstract**

Human cognition systematically favors simpler explanations, concepts, and models (Chater & Vitanyi, 2003; Feldman, 2000; Lombrozo, 2007). This preference is typically treated as a heuristic: useful but without principled justification for why simpler hypotheses are more likely true. This paper provides the missing justification. Under the Simplicity Dominance Principle (SDP), where a structure's measure is P(S) proportional to 2^(-K(S)) with K denoting Kolmogorov complexity, simpler structures are exponentially more probable. Each additional bit of specification halves a structure's measure. This is a mathematical fact about the distribution of structures.

If simpler explanations describe more probable structures, then organisms that preferentially adopt simpler explanations will, on average, act on more accurate models. Natural selection therefore favors simplicity preference. The preference is a calibration to the actual distribution of structures.

This framework unifies several findings in cognitive science: why children learn regular patterns before exceptions, why concept difficulty tracks description length, why explanatory breadth is preferred over narrowness, and why conspiracy theories provoke cognitive resistance despite fitting the data. The paper also distinguishes the mathematical criterion (brevity of generating description) from the aesthetic response (perceived elegance or beauty), arguing that the latter is an evolved proxy for the former.

**Keywords:** simplicity preference, Occam's Razor, Kolmogorov complexity, algorithmic probability, cognitive parsimony, Solomonoff induction, evolutionary epistemology

---

## **1\. The Puzzle**

Human cognition prefers simpler explanations. This has been demonstrated across perception (Hochberg & McAlister, 1953), concept learning (Feldman, 2000), causal reasoning (Lombrozo, 2007), and general inductive inference (Chater & Vitanyi, 2003). The preference is robust, early-developing, and cross-cultural.

The standard treatment in cognitive science characterizes this preference as a heuristic. Chater and Vitanyi (2003) proposed simplicity as a unifying principle in cognitive science, connecting diverse phenomena through Kolmogorov complexity. Feldman (2000) showed that the subjective difficulty of learning Boolean concepts tracks their Boolean complexity. Lombrozo (2007, 2012\) demonstrated that people prefer explanations covering more phenomena with fewer causes, a preference she connects to inference quality.

What none of these accounts provides is a principled reason why simpler hypotheses are more likely true. Without this, simplicity preference remains a useful heuristic without theoretical grounding, analogous to a compass that points north without a magnetic field to explain why.

## **2\. The Mathematical Ground**

The Simplicity Dominance Principle (SDP) provides the missing ground. The argument proceeds in three steps.

**Step 1: Algorithmic probability.** Solomonoff (1964) showed that the probability of a string x under a universal distribution is dominated by the shortest program that generates x. Formally, for universal prefix-free Turing machine U, the algorithmic probability m(x) \= sum over p of 2^(-|p|) for all p such that U(p) \= x. This sum is dominated by the shortest such p. Strings with short generating programs are exponentially more probable than strings with long generating programs.

**Step 2: From strings to structures.** If physical structures are identified with their mathematical descriptions (Tegmark, 2014), then algorithmic probability applies to structures, not just strings. A structure with Kolmogorov complexity K(S) \= k has measure proportional to 2^(-k). Each additional bit of specification halves the measure. This identification is defended independently but the present argument requires only the weaker claim that algorithmic probability provides a reasonable prior over hypotheses, which is the standard Solomonoff induction framework already accepted in machine learning theory.

**Step 3: Adaptive calibration.** If simpler structures are exponentially more probable, then an organism whose hypothesis-selection mechanism preferentially adopts lower-K hypotheses will, on average, select hypotheses that are more often correct. Natural selection favors mechanisms that produce accurate models of the environment. Therefore, natural selection favors simplicity preference.

The exponential character of the advantage is important. The simplicity preference is not marginally better than random hypothesis selection. It is exponentially better. An organism that consistently selects the lowest-K hypothesis compatible with the evidence will outperform an organism that selects among compatible hypotheses uniformly, by a factor that grows exponentially with the complexity difference between the selected hypotheses.

## **3\. Unifying Existing Findings**

This framework provides a unified explanation for several established findings.

**Concept difficulty tracks description length.** Feldman (2000) showed that Boolean concept difficulty is predicted by Boolean complexity. Under SDP, lower-complexity concepts correspond to higher-measure structures. The brain encounters these structures more frequently and has more training data for them, making them easier to learn. The complexity-difficulty correlation is a consequence of the complexity-probability correlation.

**Children learn regular patterns before exceptions.** Regular patterns have lower K than patterns with exceptions. Under SDP, regular patterns are more probable and therefore more frequently instantiated in the child's environment. The learning order tracks the probability order, which tracks the complexity order.

**Explanatory breadth is preferred.** Lombrozo (2007) showed that people prefer explanations covering multiple phenomena over narrow explanations. A single principle covering n phenomena has total description length K(principle) \+ K(mapping). Separate explanations have total description length sum of K(explanation\_i). When the principle is simpler than the sum of separate explanations, the unified explanation has lower total K and is therefore more probable under SDP. The preference for breadth is a preference for lower total K.

**Conspiracy theories provoke resistance.** Conspiracy theories typically require specifying hidden agents, secret coordination, and suppressed evidence. Each additional element increases K. The simpler explanation (visible causes producing visible effects) has lower K and is exponentially more probable. Cognitive resistance to conspiracy theories tracks the K-gap between the conspiratorial and non-conspiratorial hypotheses.

## **4\. Elegance as Evolved Proxy**

That physicists are drawn to calling low-K solutions elegant may itself be adaptive. The aesthetic response to a short derivation producing rich consequences, or a single principle unifying diverse phenomena, correlates with low K(generating rule) relative to high K(output). If low-K hypotheses are more often correct, organisms that experienced positive affect in response to low-K/high-output-complexity ratios would preferentially adopt better hypotheses.

This distinguishes the mathematical criterion (brevity of generating description) from the aesthetic response (perceived elegance, beauty, or satisfaction). The criterion is objective and measurable. The response is subjective and evolved. The response tracks the criterion because natural selection favored organisms in which it did.

Hossenfelder (2018) asked why physics should prefer simple or beautiful theories, and found no satisfactory answer in existing justifications based on aesthetics or naturalness. The present framework provides one: the criterion is brevity, and brevity has measurable consequences for probability. The aesthetic response is a separate, evolved phenomenon that happens to track the mathematical criterion.

## **5\. Relationship to Solomonoff Induction**

Solomonoff induction is the normative ideal: a perfect Bayesian reasoner using the universal prior assigns probabilities to hypotheses in proportion to 2^(-K). Human cognition approximates this ideal imperfectly, with known deviations (availability bias, representativeness heuristic, anchoring). The claim here is not that human cognition implements Solomonoff induction. It is that the simplicity preference component of human cognition is an evolved approximation of the Solomonoff prior, and that this approximation is adaptive because the Solomonoff prior tracks truth.

This connects the normative theory (Solomonoff) to the descriptive findings (Chater, Feldman, Lombrozo) through evolutionary epistemology: natural selection built approximate Solomonoff induction into cognitive architecture because approximate Solomonoff induction produces better models than the alternatives.

## **6\. Objections**

**"Simplicity preference is sometimes wrong."** Correct. The claim is that simplicity preference is calibrated to the actual distribution of structures. A calibrated instrument can be wrong on individual measurements while being right on average. The exponential advantage means that systematic simplicity preference outperforms any alternative that does not penalize complexity.

**"Kolmogorov complexity is uncomputable."** True in general. Human cognition does not compute K exactly. It approximates K through proxies: number of free parameters, number of independent causes, length of verbal or mathematical description. These proxies are imperfect but correlated with K. The evolutionary argument requires only that the proxies are correlated, not that they are exact.

**"This is just Bayesian inference with a simplicity prior."** Yes. The contribution is identifying why this prior is correct: simpler structures are exponentially more probable. The prior is not arbitrary. It tracks the distribution of existence.

## **7\. Scope and Limitations**

This paper addresses the simplicity preference specifically. Related phenomena, such as the psychological effects of deterministic worldviews on guilt, regret, and anxiety , involve different cognitive mechanisms and are addressed separately.

The paper does not claim that all cognitive biases are adaptive. It claims that one specific cognitive tendency, simplicity preference, is adaptive for a specific mathematical reason. Other biases (confirmation bias, availability heuristic) may be adaptive for different reasons or may be genuine miscalibrations.

## **8\. Conclusion**

Human simplicity preference is not a heuristic shortcut. It is a calibration to the mathematical distribution of structures, where simpler structures are exponentially more probable. Natural selection built this calibration into cognitive architecture because organisms with simplicity preference construct more accurate models. The aesthetic response to elegance is an evolved proxy for the mathematical criterion of brevity. The compass points north because there is a magnetic field.

---

## **References**


Bernstein, G. A. (2026c). Reality is mathematical structure.

Chater, N., & Vitanyi, P. (2003). Simplicity: A unifying principle in cognitive science? *Trends in Cognitive Sciences*, 7(1), 19-22.

Feldman, J. (2000). Minimization of Boolean complexity in human concept learning. *Nature*, 407(6804), 630-633.

Hochberg, J., & McAlister, E. (1953). A quantitative approach to figural "goodness." *Journal of Experimental Psychology*, 46(5), 361-364.

Hossenfelder, S. (2018). *Lost in Math: How Beauty Leads Physics Astray*. Basic Books.

Kolmogorov, A. N. (1965). Three approaches to the quantitative definition of information. *Problems of Information Transmission*, 1(1), 1-7.

Lombrozo, T. (2007). Simplicity and probability in causal explanation. *Cognitive Psychology*, 55(3), 232-257.

Lombrozo, T. (2012). Explanation and abductive inference. In K. J. Holyoak & R. G. Morrison (Eds.), *Oxford Handbook of Thinking and Reasoning* (pp. 260-276). Oxford University Press.

Solomonoff, R. J. (1964). A formal theory of inductive inference. *Information and Control*, 7(1), 1-22.

Tegmark, M. (2014). *Our Mathematical Universe*. Knopf.

All companion papers available at https://independent.academia.edu/TheMatheist

