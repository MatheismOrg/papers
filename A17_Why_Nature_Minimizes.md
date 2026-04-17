# Why Forward Time, Least Action, and Emergence Are the Same Theorem

**Gary Abraham Bernstein**

Independent Researcher
ORCID: https://orcid.org/0009-0009-1761-2867

## Abstract

The Principle of Least Action, the arrow of time, and the emergence of smooth physics from discrete rules are treated as separate results requiring separate explanations. This paper argues they are three instances of one structural fact: many-to-one mappings. Computation is many-to-one: 1+1=2 is deterministic, 2=x+y is not. The inverse does not exist. This is the arrow of computation (Bernstein, 2026h). Feynman's path integral is many-to-one: all quantum paths exist, but paths far from the extremal cancel under phase interference. What survives is the classical trajectory. The Principle of Least Action is what remains when everything else cancels. Coarse-graining (averaging over fine details) is many-to-one: many microstates map to one macrostate. When discrete rules are averaged at large scales, only extremal macro behavior survives, for the same structural reason. Variational emergent physics is forced by coarse-graining regardless of whether the fundamental rule is variational.

A forcing argument completes the picture. Observers require stable structures. Stable structures require conservation laws. Noether's theorem proves conservation laws require symmetries of a variational formulation. The chain from observers to variational physics runs one direction with no exit. Under algorithmic probability (Bernstein, 2026d), the simplest observer-supporting frameworks dominate. If all such frameworks produce variational physics through the Noether chain, algorithmic probability derives the Principle of Least Action for any universe we could inhabit.

**Keywords:** Principle of Least Action, many-to-one, arrow of time, path integral, coarse-graining, Noether's theorem, variational physics, algorithmic probability, mathematical monism

---

## 1. One Structural Fact, Three Physical Consequences

The Principle of Least Action, the arrow of time, and the emergence of variational physics from simple rules are the same mathematical pattern wearing different clothes. All three reduce to one fact: many-to-one mappings have a deterministic direction in which information survives, and an underdetermined direction in which it does not. The forward direction produces a unique output. The backward direction has infinitely many candidates and no procedure to choose among them.

This paper traces the pattern across three domains, shows why they converge, and argues that variational physics is forced for any observer-supporting universe. In five words: everything happens, most cancels out.

## 2. Many-to-One: The Structural Fact

Addition: 1+1=2 is deterministic. 2=x+y is not. The inverse function does not exist. Multiplication: 3×4=12 is deterministic. 12=x×y is not. Hashing: hash("hello")=0x5d41 is deterministic. 0x5d41=hash(?) is not. Neural computation: sensory inputs produce a percept. The percept does not uniquely determine which inputs produced it.

Nearly all functions are many-to-one: multiple inputs produce the same output. The forward direction is determined. The backward direction is underdetermined. This is a property of mathematics, true in every possible world. Under mathematical monism, the function is the territory. This structural asymmetry is as fundamental as anything in physics.

## 3. Instance One: The Arrow of Computation

Full development is given in Bernstein (2026h). The summary:

A system running backward cannot compute, because every operation has infinitely many possible inputs. It cannot predict, because prediction requires computing consequences from causes, and causes are underdetermined from consequences. It cannot model, because modeling is compression, compression is many-to-one, and reversible computation is bijective and cannot compress.

In a block universe where all states exist, backward states are present but nothing can compute through them. This is why time has a direction. Forward descriptions from simple initial conditions have low Kolmogorov complexity. The same universe described backward from its current state has high complexity, because many-to-one functions ran forward and made the end state complex. Algorithmic probability quantifies this asymmetry. The computational arrow is the structural foundation. Algorithmic probability measures its consequence.

## 4. Instance Two: The Principle of Least Action

Variational physics means specifying all of physics through one principle: systems follow the path that extremizes the action. One rule generates every trajectory. The alternative, listing each permitted trajectory separately, requires infinite specification. That is what variational means: one principle instead of an infinite list.

A thrown ball does not know to take the optimal path. In quantum mechanics, it takes all paths simultaneously. Feynman's path integral assigns every possible trajectory a phase factor proportional to e^(iS/ℏ), where S is the action along that path. Paths near the extremal have similar phases. They reinforce. Paths far from the extremal have rapidly varying phases. They cancel.

What survives the cancellation is the classical trajectory. The Principle of Least Action is what remains when everything else cancels.

This is many-to-one. Many paths contribute. One classical trajectory survives. The surviving path is determined. The inverse question, given the classical trajectory, which quantum paths contributed, has infinitely many answers. The structure is identical to 1+1=2: deterministic forward, underdetermined backward.

The analogy is not metaphorical. Both are instances of the same mathematical operation: a map from a large space to a small space, where the forward direction has a unique image and the backward direction does not.

## 5. Instance Three: Coarse-Graining and Variational Emergence

A simple discrete rule, like a cellular automaton or a hypergraph rewriting rule, operates step by step at the micro level. Each step is local: update this cell, rewrite this edge. The rule may have very low Kolmogorov complexity, perhaps just a few bits.

Zoom out. Average over many steps and many cells. This is coarse-graining: replacing a detailed microstate with a summarized macrostate. Many microstates map to the same macrostate. Coarse-graining is many-to-one by definition.

What survives coarse-graining? The same thing that survives phase cancellation in the path integral: extremal behavior. Non-extremal microstates vary in all directions. Their contributions average out, the way randomly oriented vectors sum to zero. Extremal microstates are aligned. Their contributions reinforce.

The emergent large-scale physics is variational even when the fundamental rule is not. The underlying cellular automaton has no action functional (the total kinetic-minus-potential along a path), no Lagrangian, no variational principle. The emergent physics does, because coarse-graining killed everything except the extremal survivor.

This supports the hypothesis that our physics emerges from extremely simple discrete rules (Wolfram, 2020). The simplest possible specification, a rewriting rule expressible in a few bits, may require a variational emergent physics not because anyone designed it that way, but because coarse-graining forces it. The fundamental rule compresses the description. The many-to-one structure of coarse-graining produces the variational form at the scale observers inhabit.

The connection to quantum mechanics is direct. Wolfram's multiway system applies the rewriting rule in all possible locations simultaneously. All branches coexist. When branches reconverge, they interfere: reinforcing branches survive, canceling branches vanish. This is coarse-graining over the space of branches. The multiway system IS the path integral, derived from the rule rather than postulated. The classical trajectory that survives IS what remains after many-to-one cancellation. The strongest result here is causal invariance producing general covariance (the analog of general relativity). The full correspondence between branchial space and Hilbert space remains open. The structural case is suggestive and narrowing.

## 6. The Noether Chain: Why Observers Require Variational Physics

The three instances converge on a forcing argument.

Observers require stable structures. A mind that models its environment needs atoms that hold their bonds, molecules that persist, chemistry that works. Stability means things are conserved: energy stays in orbitals, momentum carries through collisions, charge is neither created nor destroyed.

Conservation laws require symmetries. This is Noether's theorem (1918): every continuous symmetry of the action yields a conservation law. Conservation of energy comes from time-translation symmetry. Conservation of momentum comes from spatial-translation symmetry. Conservation of angular momentum comes from rotational symmetry. No symmetry, no conservation. No conservation, no stability. No stability, no chemistry. No chemistry, no observers.

Symmetries of the action require a variational formulation. Noether's theorem does not say "every symmetry yields a conservation law." It says "every continuous symmetry of the action functional yields a conservation law." The action functional is the variational object. Without it, the theorem has no input.

The chain runs one direction: observers → stability → conservation → symmetry → variational formulation. There is no exit. Any observer-supporting universe must have variational emergent physics, regardless of whether its fundamental rule is variational, discrete, or otherwise.

Can an observer-supporting universe lack variational physics? Two candidates fail:

Purely discrete physics with no continuum limit. Low K, but atoms need smooth orbitals to form stable bonds. Without a continuum limit, electron clouds cannot form. No orbitals, no chemistry, no molecular stability, no observers.

Fundamentally nonlocal dynamics with no local action. Things affect nearby things in our universe. Your hand pushes a ball, the ball hits a wall. Cause travels through space step by step. Nonlocal means every point affects every other point instantly with no connection through the space between. If effects do not travel through space, spatial structure collapses. "Here" and "there" stop meaning anything. No spatial structure, no molecules, no chemistry, no observers.

Both fail at observer-support. Variational physics appears forced for any universe containing minds.

## 7. Unification and Relationship to algorithmic probability

Algorithmic probability says structures with lower Kolmogorov complexity dominate exponentially. It selects for the simplest observer-supporting frameworks. If the Noether chain is correct, that all observer-supporting frameworks produce variational emergent physics, then algorithmic probability derives the Principle of Least Action. The derivation:

1. Algorithmic probability selects for lowest-K observer-supporting structures.
2. Observer-supporting structures require stability, conservation, and symmetry (Noether chain).
3. Conservation and symmetry require a variational formulation.
4. Therefore observers inhabit structures with variational physics.
5. Among variational structures, the simplest dominate.
6. The simplest wave equation (Schrödinger) produces the path integral.
7. The path integral produces PoLA as the coherent classical limit.

The claim that algorithmic probability "derives" PoLA rests on step 2-3: that the Noether chain has no exit. Whether there exists a non-variational observer-supporting framework remains an open question. The structural arguments against it (no stable chemistry without conservation, no conservation without Noether, no Noether without action) are strong. A proof of impossibility would close the gap. Until then, the derivation is conditional: algorithmic probability derives PoLA given that all observer-supporting physics is variational.

The three instances are connected through algorithmic probability as follows:

The computational arrow (instance one) is the structural foundation. Many-to-one functions create the K-asymmetry that algorithmic probability quantifies as the arrow of time.

The Principle of Least Action (instance two) is the coherent survivor. Many-to-one phase cancellation in the path integral produces variational classical physics from simple quantum mechanics.

Coarse-graining emergence (instance three) is the bridge. Many-to-one averaging produces variational macro physics from simple discrete rules, connecting the fundamental simplicity algorithmic probability selects for to the variational physics observers inhabit.

One mathematical structure. Three physical faces. One forcing argument from observers.

## 8. Implications

The unification has several consequences.

First, it connects papers that were previously separate: the computational arrow (Bernstein, 2026h), the algorithmic probability overview (Bernstein, 2026a), the flagship derivation (Bernstein, 2026d) . Many-to-one is the thread.

Second, it supports the hypothesis that our physics emerges from simple discrete rewriting rules. If coarse-graining forces variational emergence, the simplest possible rule (lowest K) will produce the physics we observe at our scale. The rule does not need to encode the Standard Model. It needs to be simple. The rest is coarse-graining.

Third, it strengthens the case that algorithmic probability is a single principle with convergent consequences rather than an ad hoc hypothesis applied separately to different domains. The arrow of time, the Principle of Least Action, and variational emergence are not three separate applications of algorithmic probability. They are one application seen at three scales.

Fourth, it clarifies what would count as evidence against the account. If a stable observer-supporting physics were found that is provably non-variational, with no action functional and no Noether symmetries yet with conservation laws arising from some other mechanism, that would challenge the Noether chain and weaken the derivation of PoLA from algorithmic probability.

---

## References

Bernstein, G. A. (2026a). Why these simple laws, forward time, many-worlds, and superdeterminism, in brief.

Bernstein, G. A. (2026d). Why these simple laws, forward time, many-worlds, and superdeterminism.

Bernstein, G. A. (2026h). The arrow of time is irreversible computation.




Bennett, C. H. (1973). Logical reversibility of computation. *IBM Journal of Research and Development*, 17(6), 525-532.

Feynman, R. P. (1948). Space-time approach to non-relativistic quantum mechanics. *Reviews of Modern Physics*, 20(2), 367-387.

Noether, E. (1918). Invariante Variationsprobleme. *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse*, 1918, 235–257.

Wolfram, S. (2020). *A Project to Find the Fundamental Theory of Physics*. Wolfram Media.

All companion papers available at https://matheism.org
