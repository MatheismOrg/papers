# Superdeterminism in a Block Universe: The Case for a Common Generator

## Abstract

Superdeterminism asks nothing beyond determinism except a shared cause: the devices that choose measurement settings are physical systems too, so they can share a past with the particles they measure. One short generating rule for the particle source and the setting devices then gives a local model with definite outcomes that reproduces Bell correlations, the way a short program fixes a billion digits of pi without listing them. Bell's measurement-independence condition forbids exactly that shared cause, so a local model reproducing Bell violations has to relax it, and a finite deterministic example shows how: the two setting devices stay statistically independent of each other while each stays correlated with the hidden information that decides the outcomes. A worked derivation of Michael Hall's model for a pair of spins recovers the quantum probabilities for every pair of measurement directions, and a source-first factorization turns that conditional density into one joint preparation of source and settings. What a full theory needs is a physical rule that produces that joint preparation for every experiment.

Keywords: superdeterminism, measurement independence, Bell correlations, description length, model comparison, fine tuning

## 1. Specificity is not complexity

A billion specified digits can come from a short generating rule. The digits of pi illustrate the point: a finite prefix, an initial run of the digits, can be described by an algorithm together with the requested precision. Specificity is not complexity.

This answers a specific objection to superdeterminism, the view that measurement settings correlate with the hidden states being measured. Correlations between settings and measured systems need not be separately entered in a cosmic lookup table. A shared generating rule can coordinate them.

Fine tuning has several meanings: description length, sensitivity to perturbations, small measure of admissible preparations, or a causal model's failure of faithfulness, the requirement that every observed independence follow from causal structure. Refuting an inference from specificity to long description does not answer all four. Christopher Wood and Robert Spekkens [1] formulate a problem involving observed independence and causal parameters. Their result is not the claim that each initial particle coordinate must be listed separately.

Matheism puts the source, apparatus, setting devices and observer inside the same mathematical world. This motivates a common generator account in which experiments are physical processes throughout. The generator must explain both what is measured and how the measurement conditions arise.

## 2. The conditional case for superdeterminism

Keep local causes and definite outcomes inside a hidden-variable account (a theory in which unobserved variables underlie measurement outcomes), require that account to match the observed violations, and measurement dependence, a statistical link between hidden states and chosen settings, follows. That conditional claim is what this section sets out.

Bell's comparison depends on two different independence claims: how the source is related to the chosen settings, and how each outcome depends on the distant apparatus. Measurement independence says the hidden state carries no statistical trace of how the dials were set. Local factorization says each outcome depends on its own dial and on the hidden state, and on nothing else. Written out, with a and b the settings, A and B the outcomes and lambda the model's relevant underlying state, measurement independence is P(lambda | a,b) = P(lambda) and local factorization is P(A,B | a,b,lambda) = P(A | a,lambda)P(B | b,lambda).

This paper starts from two commitments: causes act locally and deterministically, and each measurement has one definite outcome. Add correlations that violate a Bell inequality (a limit on correlations that local models with independent settings obey) and measurement dependence becomes necessary inside that hidden-variable model. A common past generator is a natural first candidate: the source and setting devices inherit their coordination from physical conditions in the same world. The preference is conditional, and a theorem stands behind the dependence it requires. The theorem does not select a unique generator.

Michael Hall [2] supplies explicit local deterministic models using measurement dependence, including a singlet correlation model (for two spins whose total spin is zero). This is established prior work. It shows why the task is to specify and assess the dependence rather than reject or accept it by name.

### The Bell bound and the role of independence

The familiar bound comes from two assumptions working together, and writing the derivation out shows the one step where independence does its work. The factorized form stated above, which covers stochastic (probabilistic) local responses, restricts attainable correlations once measurement independence is added to it. The deterministic special case makes the CHSH bound, the Bell inequality used throughout this paper, transparent. On any single hidden state all four responses are already fixed, so one of the two combinations of Bob's responses drops out. In symbols, for settings a,b in {0,1} and outcomes in {-1,+1}, write A_0(lambda), A_1(lambda), B_0(lambda), B_1(lambda). For each lambda,

$$C(\lambda)=A_0(B_0+B_1)+A_1(B_0-B_1).$$

Exactly one parenthesis vanishes and the other is +/-2. Averaging all four correlations with the same normalized rho(lambda), a weighting of hidden states that sums to one, gives

$$|E_{00}+E_{01}+E_{10}-E_{11}|\leq 2.$$

Measurement independence enters at exactly one step: the same distribution rho(lambda) averages all four setting pairs. Write rho(lambda|a,b) instead and the four correlations need not be averages over the same ensemble (population of hidden states), so the bound no longer follows.

Retain definite outcomes, local response functions and empirical agreement with Bell violations. Independence cannot remain alongside them. That is the affirmative argument for investigating measurement dependence.

## 3. Which variables must depend on which?

Two quantities can come from the same source and still be independent. In symbols, let U be uniform on {0,1,2,3} and define a = U mod 2 and lambda = floor(U/2). Each pair (a,lambda) occurs once, so P(a,lambda) = 1/4 = P(a)P(lambda): both are deterministic functions of the same U, and they are independent. Set lambda = a instead, with a balanced: then P(a=0,lambda=0) = 1/2, not 1/4. Both constructions have short descriptions, so common origin and shortness select neither independence nor dependence.

Alice's setting a and Bob's setting b can be mutually independent even when the pair is correlated with lambda. Bell's measurement independence concerns P(lambda | a,b), not P(a,b). A local deterministic model retaining measurement independence still obeys the relevant Bell bound.

A single fixed history is not itself an ensemble distribution. The distribution might describe repeated trials, uncertainty about initial data or another explicitly defined measure. Spread probability uniformly over a finite product space (every combination of values of several variables) and its coordinate projections (the separate variables) come out independent; a different mapping, or a conditioning, can produce dependence.

Superdeterminism is not simply another word for universal determination. It concerns dependence between settings and underlying states in a specified statistical account. Everettian (many-worlds) and superdeterministic approaches make different commitments about outcomes and their explanation; they are not two conclusions automatically delivered by one equation.

Take out fundamental randomness and determinism follows only if the specified laws and state fix the evolution. Standard Bohmian mechanics, a deterministic quantum theory with particles guided by a wave, supplies a counterexample: deterministic configurations follow a nonlocal guidance law, in which each particle's motion depends instantly on distant particles' positions, while settings independence can hold [4, 6]. The sharper implication follows under the assumptions of Bell's theorem: for a hidden-variable account with definite outcomes, retain local factorization and the other conditions; if that account reproduces a violation, measurement independence must fail [7].

The preference argued here is for one line of investigation: whether one local generator can do the same predictive work as a nonlocal alternative at lower cost. A shorter complete description with equal likelihood (the same probability for the observed data) has the prior advantage (the edge held before any data are weighed) specified in Section 5.

### Einstein, Bell and the block account

Albert Einstein, Boris Podolsky and Nathan Rosen supply the historical motivation: a demand for completeness tied to separated systems. A block perspective, which treats all of spacetime as one existing whole, includes the setting apparatus within the same world as the particles. Bell identifies conditions a local completion (a local theory supplementing quantum mechanics) must satisfy. Together these motivate treating the experimental settings as physical events inside the theory; they do not establish that Einstein unknowingly proved measurement dependence.

Bell's theorem should not be explained away as a failure to understand relativity. Its mathematical content survives the independence challenge: relaxing an assumption is one of the possibilities the theorem makes explicit.

## 4. A deterministic construction and its correlation measure

The two setting devices can be independent of each other while the hidden state depends on both. Let the deterministic trial counter n cycle through 0 to 7. Set a = floor(n/4), b = floor(n/2) mod 2 and r = n mod 2. All eight triples occur once, so a and b are balanced and mutually independent over a cycle.

Each side then answers using only its own dial and the state the source sent it. In symbols, the source supplies lambda = (a,b,r), with A(a,lambda) = r and B(b,lambda) = r XOR (lambda_a AND b). On the supported trials A XOR B = ab. Bell's measurement independence fails because lambda contains the settings, despite the independence of the setting devices from each other.

Scoring that construction shows how far a short deterministic rule can overshoot quantum theory. Write E_ab for the mean of (-1)^(A XOR B). Then E_00 = E_01 = E_10 = 1 and E_11 = -1. The CHSH expression is 4, exceeding the quantum maximum 2 sqrt(2). A very short, measurement-dependent rule can be nonsignalling, meaning that neither side's outcome statistics reveal the other side's setting, and still be empirically wrong for quantum experiments. Compactness and nonsignalling together are insufficient to recover quantum theory.

One extra deterministic bit tunes the construction down to the quantum value, and leaves the choice of that value unexplained. Extend lambda with a bit e, independent of a,b,r, with P(e=1) = q. Replace B by B XOR e. Now E_ab = (1-2q)(-1)^ab and CHSH = 4(1-2q). Choosing q = (2-sqrt(2))/4 reproduces 2 sqrt(2). For a computable q, assign successive eight-trial blocks the bit e_k = floor((k+1)q) - floor(kq). After N blocks the CHSH value is 4(1 - 2 floor(Nq)/N), converging to 4(1-2q). A wholly deterministic schedule supplies the stated limiting correlations. Choosing the quantum value of q still requires a physical explanation.

This match addresses only the specified four setting pairs. It does not recover full angle dependence, arbitrary preparations, sequential measurements or field theory.

### 4.1 A worked singlet model for every pair of directions

Hall's model [2] supplies the full singlet probabilities. Each side reports whether the hidden vector lies in the hemisphere around its own measurement direction, with Bob's report flipped in sign. In symbols, let $a,b$ be unit measurement directions, $\theta=\arccos(a\cdot b)$, and $\lambda$ a unit vector, and define local responses

$$A(a,\lambda)=\operatorname{sgn}(a\cdot\lambda),\qquad
B(b,\lambda)=-\operatorname{sgn}(b\cdot\lambda).$$

The hidden vector is not drawn from a fixed distribution: its density depends on the pair of directions being measured, and it weights the region where the two signs agree differently from the region where they disagree. Stated precisely, for $0<\theta<\pi$, Hall's density with respect to solid angle (area on the unit sphere of directions) is

$$\rho(\lambda\mid a,b)=\begin{cases}
\dfrac{1+a\cdot b}{8(\pi-\theta)},&\operatorname{sgn}(a\cdot\lambda)=\operatorname{sgn}(b\cdot\lambda),\\[6pt]
\dfrac{1-a\cdot b}{8\theta},&\operatorname{sgn}(a\cdot\lambda)\ne\operatorname{sgn}(b\cdot\lambda).
\end{cases}$$

The sizes of the two regions are what make this density normalize and return the quantum correlation. The two equal-sign regions have area $4(\pi-\theta)$; the unequal-sign regions have area $4\theta$. Normalization and correlation follow directly:

$$\int\rho\,d\Omega=\frac{1+\cos\theta}{2}+\frac{1-\cos\theta}{2}=1,$$
$$E(a,b)=-\frac{1+\cos\theta}{2}+\frac{1-\cos\theta}{2}=-\cos\theta.$$

Consequently, for $\alpha,\beta\in\{-1,+1\}$,

$$P(A=\alpha,B=\beta\mid a,b)=\frac{1-\alpha\beta\,a\cdot b}{4}.$$

At $\theta=0$ or $\pi$, uniform density $1/(4\pi)$ gives perfect anticorrelation or correlation. This recovers every singlet direction pair and nonsignalling marginals, the statistics of each side taken alone. Local deterministic responses with measurement dependence are therefore adequate to these data; the density itself still does not follow from any law of cosmology.

### 4.2 Turning the conditional model into a common preparation

A density conditioned on the settings can be turned around into one joint preparation: draw the hidden state first, then the settings. For a finite menu of setting pairs $j$, let $p_j>0$ be their preparation weights with $\sum_jp_j=1$. Write $\rho_j(\lambda)$ for the preceding density for pair j. Written out, the joint distribution is

$$P(j,\lambda)=p_j\rho_j(\lambda),\qquad
\bar\rho(\lambda)=\sum_jp_j\rho_j(\lambda),$$
$$P(j\mid\lambda)=\frac{p_j\rho_j(\lambda)}{\bar\rho(\lambda)}\quad\text{where }\bar\rho(\lambda)>0.$$

Generate a common preparation by selecting $\lambda$ with density $\bar\rho$ and selecting the setting pair using the displayed conditional probabilities. Sum and condition, and $P(j)=p_j$ and $P(\lambda\mid j)=\rho_j(\lambda)$ come back exactly. A physical account must explain the preparation and what remains true when the menu or apparatus changes.

## 5. Compare predictive models with predictive models

A comparison is between complete accounts, so everything each one needs has to be written into it before anything is counted. Let M include the state space (the set of possible states), dynamics, preparation distribution, observation map and experimental protocol. Fix an effective encoding convention U, a mechanical scheme for writing accounts as bit strings. The shorter a complete account is in that convention, the more prior weight it receives, because each bit of its description is a 50/50 choice that no rule makes: the prior is w(M) proportional to 2^(-K_U(M)) [3]. Two accounts are then compared by multiplying the ratio of those prior weights by the ratio of how well each predicted the data. Written out, the posterior comparison is

P(M_1 | D) / P(M_2 | D) = [w(M_1) / w(M_2)] [P(D | M_1) / P(D | M_2)].

A delta-bit advantage gives the factor 2^delta. When the models give identical likelihoods on the compared evidence, the prior advantage remains.

Shared information should be counted once. A joint description of the laws, preparation and prediction rules can reuse information that separate descriptions would repeat. The shortest program length is generally uncomputable (no algorithm can calculate it in every case), so an exhibited program supplies an upper bound. Changing the reference machine, the standard computer used to measure program length, can also change finite rankings. Input precision and the preparation measure must be included before a short program can count as a short physical explanation.

## 6. A stochastic law is not a list of its outcomes

This paper asks why one particular event occurs rather than another. A probability distribution answers a different question: how outcomes are weighted under specified conditions. If the selection is declared irreducible, a deterministic explanation has been set aside. This is a demand for causal completeness, not a contradiction in the probability axioms.

A fair coin law has a short specification. For N tosses it assigns each string probability 2^(-N). Describing one particular string exactly may require about N bits, but that is a data cost. An all-zero string also has probability 2^(-N) and a short exact description, so it is false that every randomly produced finite record is incompressible, impossible to describe in fewer bits than it contains.

A finite binary measurement outcome does not need infinitely many bits merely because it was stochastic. Deterministic continuum models can also place unlimited information in their initial data. The comparison must treat both sides at the same observational resolution.

An Everettian model must account for the weights of expected observations. A deterministic model with unknown initial conditions also needs a distribution over possibilities. For accounts assigning the same probabilities to D, their likelihood terms agree. Giving one account a data cost while giving its rival none would compare different explanatory tasks.

## 7. Existing proposals and their limits

Ciepielewski, Okon and Sudarsky [4] construct a local pilot-wave account (in which a guiding wave steers the particles) with specially homogeneous initial data. This establishes that a compact constraint can express a highly coordinated preparation. The proposed relaxation dynamics, a process by which the system settles toward a target distribution, leaves competition between terms unresolved, so generic convergence cannot yet be used as the explanation.

Tim Palmer [5] argues that the geometry of the admissible states changes how a counterfactual perturbation should be judged, the question of what would have followed had a dial been set differently. A metric defines distances; a probability measure assigns weights. Changing the first does not by itself derive the second.

An attractor, a state or set of states that evolution tends toward, does not erase every initial-condition cost. A basin (the set of starting states drawn to the attractor), convergence rate, tolerance and observation time must be specified. These distinctions apply to any proposed relaxation mechanism.

## 8. The stronger fine-tuning objection

Wood and Spekkens [1] connect Bell correlations with faithfulness in classical causal models: observed independence can require restrictions on causal parameters beyond those encoded by the graph. A short formula may impose such a restriction and still be fine tuned in that technical sense.

A useful answer must explain why the permitted physical variations preserve the observed independence and quantum statistics. One option is an independently motivated law restricting the parameter family. Another is a dynamical mechanism with a specified preparation measure.

The finite example makes this concrete. A balanced common bit makes both local outcomes unbiased. Bias that bit and Bob's marginal at b=1 generally depends on a. Whether that perturbation is admissible belongs to the model. Coding cost is one diagnostic alongside causal explanation, robustness and empirical fit.

## 9. Scientific inference and experimental variations

Science does not require settings to be chosen by an uncaused will. It does require a justified route from actual samples and interventions to the conclusions drawn from them. A measurement-dependent model must explain that route and identify which experimental variations preserve its predictions.

Specify the setting generator, source preparation and data selection procedure. Predict what changes when the setting protocol is altered. The empirical target includes repeated preparations, nonsignalling marginals, the angle dependence of correlations and consistency with other physical observations.

## 10. Conclusion

One physical process can lie behind the source, the setting devices and their correlations, because a short rule can produce a highly specific history without listing its events one by one. Apparent randomness does not establish the absence of an underlying generator.

The finite example separates independence of the settings from their joint dependence on hidden source information. The Hall construction recovers the complete singlet joint probabilities, with local deterministic responses and unbiased marginals. The source-first factorization supplies a common preparation for a declared setting menu. A complete physical theory must additionally explain this preparation and its stability across experimental changes.

## References

[1] Wood, C. J., and Spekkens, R. W. (2015). The lesson of causal discovery algorithms for quantum correlations: causal explanations of Bell inequality violations require fine tuning. New Journal of Physics 17, 033002. https://arxiv.org/abs/1208.4119

[2] Hall, M. J. W. (2010). Local deterministic model of singlet state correlations based on relaxing measurement independence. Physical Review Letters 105, 250404. https://arxiv.org/abs/1007.5518

[3] Hutter, M., Legg, S., and Vitanyi, P. M. B. (2007). Algorithmic probability. Scholarpedia 2(8), 2572. https://www.hutter1.net/ai/algprob.pdf

[4] Ciepielewski, G. S., Okon, E., and Sudarsky, D. (2023). On Superdeterministic Rejections of Settings Independence. The British Journal for the Philosophy of Science 74(2). Preprint: https://arxiv.org/abs/2008.00631

[5] Palmer, T. N. (2024). Superdeterminism without Conspiracy. Universe 10(1), 47. https://arxiv.org/abs/2308.11262

[6] Bohm, D. (1952). A Suggested Interpretation of the Quantum Theory in Terms of Hidden Variables. I. Physical Review 85, 166. https://doi.org/10.1103/PhysRev.85.166

[7] Hossenfelder, S., and Palmer, T. (2020). Rethinking Superdeterminism. Frontiers in Physics 8, 139. https://doi.org/10.3389/fphy.2020.00139
