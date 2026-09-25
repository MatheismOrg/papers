# When Simpler Explanations Help Minds Learn and Survive

**Gary Abraham Bernstein**

Independent Researcher · ORCID: https://orcid.org/0009-0009-1761-2867

## Abstract

A short explanation can help a mind learn more from less evidence and reuse what it learns. It helps survival only if the saved effort does not discard distinctions needed for action. This paper states that trade-off as an evolutionary comparison: benefits across encountered environments, less learning and implementation costs. It also identifies an exact information limit. When two situations look identical but require incompatible responses, no strategy using only those observations can guarantee success in both. Simplicity must therefore preserve the information that matters. Mathematical work on prediction supports the account, and cognitive findings motivate testing it in learners. The proposal explains why economical explanations may be adaptive and why perceived elegance can mislead; a large description-weight advantage is not automatically a large survival advantage.

**Keywords:** simplicity preference, learning, algorithmic probability, evolutionary epistemology, explanation

## 1. When simplicity helps

A learner cannot keep every detail or test every explanation. A short rule can identify a regularity that works across many encounters. It saves effort when the same rule can be reused, and it can support predictions beyond the examples already seen. This gives simplicity a practical role. The question is when that role improves decisions enough to favor its evolution.

Existing research already gives this question a mathematical foundation. Chater and Vitányi (2003) develop simplicity as a unifying principle in cognition; Feldman (2000) relates the difficulty of learning Boolean concepts to their complexity; Lombrozo (2007) examines simplicity and probability in causal explanation. These findings motivate an evolutionary account without, by themselves, establishing one. The present contribution is to make the required connection explicit: which environments, decisions and costs make a simplicity preference advantageous?

## 2. What a shorter description buys

Fix a language for describing candidate explanations. Suppose each complete description of length L receives initial weight proportional to 2^(−L). In a prefix-free language, no complete description begins with another complete description, so programs can be distinguished without an extra end marker. Under this rule, five fewer bits give an individual description 32 times the initial weight. This is an exact consequence of the chosen weighting.

An explanation's initial weight is not its probability after evidence. A shorter explanation that predicts the observations badly can lose to a longer one that predicts them well. Nor does a ratio between individual descriptions equal a ratio between entire classes of explanations. Several programs can describe the same output, and their weights must be combined under a mixture account.

Solomonoff's approach sums over programs rather than invariably choosing the shortest explanation compatible with the data. Mathematical results on Bayesian mixtures relate prediction loss to the weight assigned to the true environment when it belongs to the model class (Hutter, 2003). Their stated assumptions matter. They supply a reason to investigate economical prediction; they do not show that a finite organism implementing a rough brevity rule always outperforms every alternative.

The broader matheist proposal extends simplicity weighting to accounts of physical reality. The cognitive argument can be assessed without first accepting that extension. It needs an account of the regularities organisms actually encounter and the cost of exploiting them. A short description of an environment is useful only if the learner can access or approximate the relevant rule with its available resources.

## 3. The evolutionary comparison

Consider two inherited learning strategies in the same population and conditions. Each strategy earns different reproductive benefits in different environments and costs something to build and run. In a simplified model, define its expected net benefit as

F(P) = Σ_e q(e)u(P,e) − c(P).

Here P is the strategy, q(e) is the frequency of environment e, u(P,e) is its reproductive benefit there, and c(P) is its implementation cost in the same units. The sum averages over the environments encountered. Within a model where these quantities determine expected reproductive success, a simpler strategy is favored over a rival when its F is greater. Applying that comparison to evolution also requires heritable variation and the relevant population conditions.

This compact calculation states the mechanism. Reusable regularities can let a short model earn comparable benefits at lower cost, or better benefits for the same cost. Conversely, a costly exception may be worth learning if missing it is dangerous. Predictive accuracy, description length and reproductive success are related through the task; they are not interchangeable quantities.

The environment frequencies cannot simply be replaced by a simplicity prior and treated as measured facts. A simplicity-weighted distribution is one proposed model for q. Evidence about actual encounters, rewards and costs is needed to establish whether it explains the observed preference. Even within that model, exponential ratios of prior weights do not become exponential fitness ratios without a derivation linking them.

## 4. As simple as the task allows

Compression can remove detail that makes no difference to a decision. Removing a crucial distinction defeats the point. Suppose two situations supply identical observations and memory, but successful action requires incompatible responses. A strategy based only on that information must respond identically, or use the same probabilities, in both. It cannot guarantee success in both. More computation cannot replace the missing distinction.

Later evidence can support a reconstruction of what happened. It does not make that evidence available to the earlier choice. This matters for the proposed connection between learning and the arrow of time: an environment requiring recurrent reconstruction after consequential events may offer less usable guidance for action. Establishing a temporal asymmetry still requires comparing learners and environments in both orientations. The information limit itself applies in either direction.

Useful simplicity therefore compresses the account while preserving the distinctions needed for the task. The best model for an organism may be more elaborate than the shortest explanation of its observations, because a rare hazard matters greatly to its survival. This yields a clear comparative prediction: as the cost of missing an exception rises, an adaptive learner should become more willing to pay for a model that represents it.

## 5. Cognitive findings and perceived elegance

A relation between concept difficulty and description length fits the proposal: a reusable rule can be easier to acquire than a list of unrelated cases. But learning difficulty also depends on representation, exposure and the learner's capacities. The evolutionary explanation should compete with those accounts, rather than infer selection directly from a correlation.

Explanatory breadth can save description length when one rule and its application explain several observations more economically than separate accounts. The comparison must count the application too. A slogan that appears to explain everything by leaving out the difficult details has not achieved the same result.

Perceived elegance may guide the search toward such economies. A compact derivation with many consequences can be rewarding to understand. That response is a possible aid to discovery, with possible errors. Familiarity or presentation may make one explanation feel simpler while another predicts better. Establishing that the response evolved because it improved decisions requires evidence beyond the existence of the response.

Shortest program length, or Kolmogorov complexity, is defined relative to a description language and is uncomputable in general. Practical studies must use specified approximations, such as a fitted grammar or a compression scheme. These are measurable model choices; they cannot be treated as exact measurements of an object's shortest possible description. Nor does an elaborate-looking output necessarily have a complex generator.

## 6. Comparisons that could support the account

The central test is whether a simplicity preference improves decisions after its costs are counted. Compare learners on the same observations and tasks, measuring predictive success, consequential mistakes and the effort needed to acquire and use their models. Include environments with reusable rules and environments in which exceptions matter. This asks when the preference earns its keep.

A second comparison varies the consequences while holding the observations fixed. If a previously negligible exception becomes dangerous, the expected benefit of representing it increases. A learner that adapts its model accordingly fits the task-sensitive account better than an unconditional rule to choose the shortest explanation.

A third comparison separates felt elegance from performance. Do judgments of simplicity help identify models that generalize to unseen cases once familiarity and presentation are controlled? Agreement would support elegance as a useful proxy in that setting. Disagreement would locate its limits. Neither outcome alone establishes an evolutionary origin.

These tests do not require identical learning rates across sensory systems, an exact 2^(−L) curve for every human judgment, or a universal ordering of explanations across cultures. Those stronger predictions require additional models of encoding and learning. The contribution here is a common comparison of usefulness and cost, with those choices made explicit.

## 7. Conclusion

Simplicity can help a mind reuse what it learns, generalize from limited evidence and act at lower cost. It is adaptive when those gains improve reproductive success in the environments encountered. The comparison must include the cost of missed distinctions as well as the cost of representing them.

A prior favoring short descriptions supplies a possible basis for economical prediction. Existing mathematical results establish benefits under specified assumptions; they do not establish a universal evolutionary advantage for choosing the shortest account. The evolutionary hypothesis becomes stronger by identifying the conditions under which simplicity succeeds, and the information it must preserve.

## References

Chater, N., and Vitányi, P. (2003). Simplicity: A unifying principle in cognitive science? Trends in Cognitive Sciences 7(1), 19–22. https://doi.org/10.1016/S1364-6613(02)00005-0

Feldman, J. (2000). Minimization of Boolean complexity in human concept learning. Nature 407, 630–633. https://doi.org/10.1038/35036586

Hutter, M. (2003). Optimality of Universal Bayesian Sequence Prediction for General Loss and Alphabet. Journal of Machine Learning Research 4, 971–1000. https://www.jmlr.org/papers/v4/hutter03a.html

Kolmogorov, A. N. (1965). Three approaches to the quantitative definition of information. Problems of Information Transmission 1(1), 1–7.

Lombrozo, T. (2007). Simplicity and probability in causal explanation. Cognitive Psychology 55(3), 232–257. https://doi.org/10.1016/j.cogpsych.2006.09.006

Solomonoff, R. J. (1964). A formal theory of inductive inference. Information and Control 7, 1–22 and 224–254.
