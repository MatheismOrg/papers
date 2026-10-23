# When Simpler Explanations Help Minds Learn and Survive

## Abstract

Minds prefer simple explanations because the world is grown from simple laws. A world built in causal steps from a short rule is compressible, so a learner that bets on the shortest rule consistent with what it has seen predicts new cases better on average, learns from less evidence and reuses what it learns. Selection rewards that bet wherever its savings in learning and action exceed its costs. Simplicity has one limit: compression must keep every distinction that action needs, since two situations that look identical but demand incompatible responses defeat any strategy that cannot tell them apart. The account explains why a taste for elegance is adaptive and why elegance can still mislead, and it states the comparisons in real learners that would test it.

Keywords: simplicity preference, learning, algorithmic probability, evolutionary epistemology, explanation

## 1. When simplicity helps

A short rule earns its keep on cases the learner has not seen. Compression stops helping when it drops a distinction that successful action needs. A learner has limited memory and time, so the advantage depends on which details a reusable explanation preserves.

Nick Chater and Paul Vitanyi (2003) develop simplicity as a principle in cognition, Jacob Feldman (2000) relates concept-learning difficulty to complexity, and Tania Lombrozo (2007) examines simplicity in causal explanation. These findings motivate the evolutionary question: in which environments do the predictive benefits exceed the learning and implementation costs?

## 2. What a shorter description buys

A shorter description starts with more weight, before any evidence is brought in, because each of its bits is a 50/50 choice that no rule makes: each bit saved doubles the starting weight. Stated precisely, fix a language for describing candidate explanations and suppose each complete description of length L receives initial weight proportional to 2^(-L). In a prefix-free language, no complete description begins with another, so programs can be distinguished without an extra end marker. Five fewer bits give an individual description 32 times the initial weight. The weighting follows from treating each unruled bit as 50/50 (*Why These Simple Laws?*, Section 3.3).

An explanation's initial weight is not its probability after evidence. A shorter explanation that predicts the observations badly can lose to a longer one that predicts them well. Nor does a ratio between individual descriptions equal a ratio between entire classes of explanations. Several programs can describe the same output, and their weights must be combined under a mixture account.

Ray Solomonoff's approach (Solomonoff, 1964) sums over programs rather than choosing the shortest explanation compatible with the data. Mathematical results on Bayesian mixtures, weighted blends of candidate models updated by evidence, relate prediction loss to the weight assigned to the true environment when it belongs to the model class (Hutter, 2003). Their stated assumptions matter: they supply a reason to investigate economical prediction, not a proof that a finite organism implementing a rough brevity rule outperforms every alternative. A short description of an environment is useful only if the learner can access or approximate the relevant rule with its available resources.

## 3. The evolutionary comparison

Whether a simpler strategy wins depends on what it earns across the environments it meets, minus what it costs to run. Consider two inherited learning strategies in the same population and conditions. Written out, expected net benefit is

F(P) = Sigma_e q(e)u(P,e) - c(P).

Here P is the strategy, q(e) is the frequency of environment e, u(P,e) is its reproductive benefit there, and c(P) is its implementation cost in the same units. A simpler strategy is favored when its F is greater. Applying that comparison to evolution also requires heritable variation, differences passed from parents to offspring, and the relevant population conditions.

Reusable regularities can let a short model earn comparable benefits at lower cost, or better benefits for the same cost. Conversely, a costly exception may be worth learning if missing it is dangerous. Predictive accuracy, description length and reproductive success are related through the task; they are not interchangeable quantities. A simplicity-weighted distribution is one proposed model for q. Only evidence about actual encounters, rewards and costs can establish whether it explains the observed preference.

## 4. As simple as the task allows

Compression can remove detail that makes no difference to a decision. Removing a crucial distinction defeats the point. Suppose two situations supply identical observations and memory, but successful action requires incompatible responses. A strategy based only on that information must respond identically in both. It cannot guarantee success in both. More computation cannot replace the missing distinction.

Useful simplicity therefore compresses the account while preserving the distinctions needed for the task. The best model for an organism may be more elaborate than the shortest explanation of its observations, because a rare hazard matters greatly to survival. This yields a clear comparative prediction: as the cost of missing an exception rises, an adaptive learner should become more willing to pay for a model that represents it.

## 5. Cognitive findings and perceived elegance

The cognitive findings fit the proposal without establishing it. A relation between concept difficulty and description length is consistent with it: a reusable rule can be easier to acquire than a list of unrelated cases. The evolutionary explanation should compete with accounts based on representation, exposure and the learner's capacities. A correlation by itself does not license an inference to selection.

Explanatory breadth can save description length when one rule explains several observations more economically than separate accounts. The comparison must count the cost of applying the rule too. A slogan that appears to explain everything by leaving out the difficult details has not achieved the same economy. Perceived elegance may guide the search toward such economies, but familiarity or presentation may make one explanation feel simpler while another predicts better. That the response exists does not show it evolved because it improved decisions. Showing that takes further evidence.

Any empirical test of the proposal has to work with a stand-in for shortest description length. Shortest program length is defined relative to a description language and is uncomputable in general (Kolmogorov, 1965). Practical studies must use specified approximations. These are measurable model choices; they cannot be treated as exact measurements of an object's shortest possible description.

## 6. Comparisons that could support the account

The central test is whether a simplicity preference improves decisions after its costs are counted. Compare learners on the same observations and tasks, measuring predictive success, consequential mistakes and effort. Include environments with reusable rules and environments in which exceptions matter. This asks when the preference earns its keep.

A second comparison varies the consequences while holding the observations fixed. If a previously negligible exception becomes dangerous, the expected benefit of representing it increases. A learner that adapts its model accordingly fits the task-sensitive account better than an unconditional rule to choose the shortest explanation. A third comparison separates felt elegance from performance: do judgments of simplicity help identify models that generalize to unseen cases once familiarity and presentation are controlled?

## 7. Conclusion

A short model lets a mind reuse what it learns, generalize from limited evidence and act at lower cost. It is adaptive when those gains improve reproductive success in the environments encountered. The comparison must include the cost of missed distinctions as well as the cost of representing them. Existing mathematical results establish benefits under specified assumptions. To strengthen the evolutionary hypothesis, specify the conditions under which simplicity succeeds and the information it must preserve.

## References

Chater, N., and Vitanyi, P. (2003). Simplicity: A unifying principle in cognitive science? Trends in Cognitive Sciences 7(1), 19-22. https://doi.org/10.1016/S1364-6613(02)00005-0

Feldman, J. (2000). Minimization of Boolean complexity in human concept learning. Nature 407, 630-633. https://doi.org/10.1038/35036586

Hutter, M. (2003). Optimality of Universal Bayesian Sequence Prediction for General Loss and Alphabet. Journal of Machine Learning Research 4, 971-1000. https://www.jmlr.org/papers/v4/hutter03a.html

Kolmogorov, A. N. (1965). Three approaches to the quantitative definition of information. Problems of Information Transmission 1(1), 1-7.

Lombrozo, T. (2007). Simplicity and probability in causal explanation. Cognitive Psychology 55(3), 232-257. https://doi.org/10.1016/j.cogpsych.2006.09.006

Solomonoff, R. J. (1964). A formal theory of inductive inference. Information and Control 7, 1-22 and 224-254.
