# 14 to Infinite Dimensions from a 4D Metric: Can Geometric Harmonics Explain Particles?

**Gary Abraham Bernstein**

Independent Researcher · ORCID: 0009-0009-1761-2867

## Abstract

Can one geometric rule explain many particles, as one instrument produces many tones? A four-dimensional spacetime metric has ten independent components, giving fourteen coordinates when positions and possible metric values are counted together. Adding derivatives of every order gives an infinite jet tower: a hierarchy recording how the metric changes. A short reusable rule specifies that hierarchy. The paper derives its coordinate counts and uses an oscillator to show how equations of motion, rather than derivative order alone, select physical modes. A worked quadratic-gravity example then connects orbital motion with light deflection. Its extra forces fade at large distance, exposing a limitation for that proposed dark-matter mechanism. The aim is to replace separate physical ingredients with consequences of one construction; the infinite mathematical hierarchy is established, while a derived particle spectrum remains the target.

## 1. One construction, many possible modes

The particle zoo invites a musical analogy. A small set of rules might produce many distinct excitations, just as one instrument produces many tones. The appeal is explanatory: replace a catalogue with a construction. A successful construction would account for both the existence of the modes and the pattern of their interactions.

An instrument, however, does not acquire its frequencies merely because its displacements can be differentiated. Its material, dynamics and boundary conditions determine its modes. A geometric theory has the same obligation. The metric bundle identifies possible configurations; an action and its solutions determine which disturbances propagate and how they interact.

The proposed term geometric harmony refers to that explanatory ambition. It does not identify each derivative of the metric with another overtone. Nor does the existence of an infinite sequence of derivatives show that infinitely many physical excitations are present with decreasing amplitudes.

The attraction of geometric unification is that apparently independent ingredients may follow from one construction. A successful account would derive interactions and relations among parameters that otherwise have to be entered separately. Counting the consequences as separate axioms would miss that economy.

Weinstein’s Geometric Unity motivates the fourteen-coordinate starting point. The count describes the metric bundle, a space containing possible metrics over spacetime. It is not fourteen observed spacetime directions. From that starting point, this paper follows one chain: metric components, their derivatives, physical modes selected by dynamics, and a calculable gravitational response.

## 2. Why the metric bundle has fourteen dimensions

Let M be a smooth n-dimensional manifold. At x in M, a metric is a nondegenerate symmetric bilinear form on the tangent space. In a coordinate chart, the symmetric matrix has n diagonal and n(n−1)/2 independent off-diagonal entries. Restricting to a fixed nondegenerate signature selects an open subset of that space and leaves its local dimension unchanged. Thus

$$\dim E=n+\frac{n(n+1)}2.$$

For n=4, there are four base coordinates and ten fibre coordinates: fourteen in total. A metric field g on M is a section of E. Assuming such a global field exists is a separate global geometric condition; the local count does not establish it.

Choosing one section does not make all fibre coordinates independent physical excitations. Constraints, gauge equivalence and the equations of motion determine physical degrees of freedom. Higher jets likewise organize derivatives without assigning a new particle to each coordinate. The following sections develop that distinction.

## 3. Jet bundles: organizing the derivatives

At a point in a four-dimensional manifold, a symmetric metric tensor has ten independent components. The bundle E of metrics of fixed Lorentzian signature therefore has fourteen local coordinates, including the four base coordinates. A spacetime metric is a section of E.

The n-th jet bundle records the metric's derivatives through order n. The number of coordinates is

$$\dim J^nE=4+10\binom{n+4}{4}.$$

The resulting dimensions begin 14, 54 and 154. The increase reflects the number of partial derivatives, including their index symmetries. It is not an increase in the dimension of the spacetime base.

For a fixed physical metric, its derivatives are derived quantities. A jet description makes them available as coordinates for local differential equations and variational constructions. It does not declare them to be independently excitable physical fields. Constraints and gauge relations must be accounted for before counting propagating degrees of freedom.

The construction rule is reusable. For an effective encoding, specifying order n and applying a fixed jet construction gives an upper bound of the form K of E plus K of n plus a constant. It does not follow that every finite extension adds literally zero bits, or that an action involving higher derivatives has the same description length as Einstein gravity.

### From fourteen to an infinite tower

The counts 14, 54, 154 and so on have no finite upper bound as the derivative order increases. Keeping compatible derivative data at every finite order defines the infinite jet tower. Each projection to a lower order simply forgets the higher derivatives. Its base remains four-dimensional spacetime; the additional coordinates describe field values and derivatives, rather than new spacetime directions.

This gives the infinity claim its precise meaning. A finite rule specifies the hierarchy without listing infinitely many coordinates separately. It does not specify an arbitrary value for every coordinate, nor give every coordinate an independent physical excitation. Infinite jet bundles are an established geometric construction; their proposed use in deriving particle content is the question here (De Vecchi, 2017).

The harmonic-series analogy expresses the ambition that one structure supports many related physical modes. Here it is not the numerical series 1 + 1/2 + 1/3 + …, and the jet orders are not automatically the integer-frequency overtones of an instrument. A physical harmonic spectrum would have to follow from the action and its admissible solutions. That distinction keeps both the infinite hierarchy and the unifying physical proposal visible.

### When derivatives determine a field

An infinite derivative list at a point does not generally determine a smooth function away from that point. For example, define f(x) to be exp(-1/x squared) for nonzero x and zero at x = 0. This function is smooth and all its derivatives at zero vanish, although it is nonzero away from zero. The zero function has the same infinite jet at the origin and is a different function.

For an analytic function, a convergent Taylor expansion determines the function in a neighborhood of the point. Global continuation needs additional domain and consistency conditions. The distinction matters when a proposal moves from smooth metrics to claims that a single point encodes an entire geometry.

Neither local Taylor reconstruction nor a code for a geometric object is equivalent to a physical holographic duality. Such a duality must relate observables and dynamics, not merely represent information in another notation. The present programme does not derive holography from the existence of derivatives.

## 4. The description cost of generated geometry

Use a fixed universal prefix machine and a finite effective description D of the input geometry. This may describe a particular computable atlas or a finite symbolic presentation; the argument does not assume every smooth geometry has a finite exact code. Let F be a fixed effective procedure producing the chosen representation of the metric bundle from D. The program “apply F to the supplied input” gives

$$K(F(D)\mid D)\leq c_F.$$

The constant depends on the procedure, reference machine and representation, not on the particular input D within the procedure’s domain. Similarly,

$$K(D,F(D))\leq K(D)+c'_F.$$

This is an upper bound, not a computed minimum. It establishes a useful positive result: the generated construction need not be encoded component by component. It also identifies the limit of the phrase “free dimensions.” A bounded reusable cost need not be zero, and a chosen section, action or boundary condition is not automatically an output of F.

For a finite jet construction F_n, a fixed algorithm can instead read both D and an encoded order n, giving K(F_n(D) | D,n) bounded by a constant. Increasing the order adds the cost of specifying n when it is not already supplied. The growing coordinate list need not imply comparable growth in the description of the rule that generates it.

## 5. The same geometry can support different laws

An explicit example shows what the bundle construction leaves open. On a four-dimensional Lorentzian base, consider the two actions

$$S_0[g]=\frac{1}{16\pi G}\int d^4x\sqrt{-g}\,R,$$

$$S_\Lambda[g]=\frac{1}{16\pi G}\int d^4x\sqrt{-g}\,(R-2\Lambda).$$

For compactly supported metric variations, the corresponding vacuum equations are G_ab=0 and G_ab+Lambda g_ab=0. The two theories use the same metric configuration bundle. A flat metric solves the first vacuum equation, while a nonzero Lambda excludes that flat vacuum in the second. Four plus ten therefore does not select even this elementary difference in dynamics.

The example is a mathematical underdetermination result. No experiment is needed to establish that these actions differ. Observations become relevant when identifying which action, or deeper construction generating it, describes a physical system. A priori mathematical monism can motivate the search for that construction without supplying the missing selection merely by naming the ontology.

## 6. Which disturbances become physical modes?

Take one oscillatory function with a fixed frequency. Differentiating it changes its amplitude and phase in a prescribed way; it does not produce a new independent oscillator at a new frequency for every derivative. This elementary example shows why a jet tower cannot be read directly as a particle spectrum.

Physical modes come from a kinetic operator and its constraints about a chosen background. Boundary conditions determine which solutions are admissible. Interactions determine how excitations exchange energy and whether they are stable. A higher-derivative action may introduce extra poles, alter constraints or produce pathological modes; it need not add one healthy mode per derivative order.

The instrument analogy becomes useful at this dynamical stage. The action plays the role of the instrument's construction. Its linearized solutions supply candidate tones. Interaction and preparation determine amplitudes. The geometry alone does not select any of those quantities.

### One oscillator, infinitely many derivatives

Let q(t) = A cos(omega t). Its n-th derivative is A omega^n cos(omega t + n pi/2). Every derivative has the same frequency. Their phases and amplitudes are fixed by A and omega. Treating q, dot q, double-dot q and so on as a list creates infinitely many mathematical entries without creating infinitely many freely excitable oscillators. This is the simplest counterexample to identifying derivative order with an overtone number.

The distinction remains when a Lagrangian contains higher derivatives. Consider the illustrative mechanical action

$$S[q]=\frac12\int dt\left[\dot q^2-\omega^2q^2-\eta\ddot q^2\right],\qquad \eta>0.$$

The higher-derivative Euler-Lagrange equation is

$$\eta q^{(4)}+\ddot q+\omega^2q=0.$$

Substitution of q proportional to exp(i nu t) gives eta nu^4 − nu^2 + omega^2 = 0. With 0 < 4 eta omega^2 < 1, the two squared frequencies are

$$\nu_\pm^2=\frac{1\pm\sqrt{1-4\eta\omega^2}}{2\eta}.$$

The extra frequency is supplied by eta and the action, not by the existence of q's derivatives. As eta tends to zero, the lower frequency tends to omega and the higher diverges. Moreover, the inverse quadratic operator has opposite residues at its two simple frequency poles. Real frequencies therefore do not by themselves establish a positive-energy physical spectrum. This toy model is a diagnostic, not a proposed viable gravitational theory.

For a nondegenerate higher-derivative action, the additional initial data and stability issue must be assessed. Constraints, degeneracies and quantization prescriptions can change the outcome, but they must be displayed. Infinite smoothness is compatible with a second-order equation having only a small number of physical modes. Conversely, one higher-derivative term can introduce an unwanted mode. The instrument is the dynamics, not the derivative list.

## 7. The action fixes the particle content

An effective field theory can organize operators by their relevance at a specified scale. A schematic contribution from an operator O of mass dimension d in four dimensions has a coefficient c divided by a scale Lambda raised to d minus four. For d greater than four, effects may be suppressed at energies well below Lambda under the assumptions of a controlled expansion.

That suppression requires a scale, coefficients and a regime of validity. Symmetries, light degrees of freedom, thresholds and possible nonlocal effects all matter. An arbitrary derivative tower does not supply the scale or guarantee modest coefficients.

Moreover, operator mass dimension and jet order are not identical classifications. Products of fields and derivatives can yield many operators at the same order, subject to identities and redundancies. A proposed correspondence must construct an operator basis and show how the geometric action fixes its coefficients.

Geometric harmony could become explanatory if a specific action determines a stable spectrum and relates apparently independent low-energy parameters. It has not yet derived effective field theory or the Standard Model merely by naming the derivative hierarchy a harmonic series.

### A geometric operator expansion with explicit costs

A possible four-dimensional effective action begins with

$$S=\int d^4x\sqrt{-g}\left[\frac{M_P^2}{2}R-\Lambda_{\mathrm{vac}}+\alpha R^2+\beta C^2+\frac{\gamma}{M^2}R\Box R+\cdots\right]+S_m.$$

The notation exposes the missing ingredients. M is a suppression scale; alpha, beta and gamma are coefficients; S_m determines matter and its coupling. Different contractions and identities affect the independent operator basis. The number of jet coordinates neither selects these coefficients nor sets Lambda_vac. Even when a deeper theory derives them from one parameter, that derivation must be supplied before the apparent economy can be credited.

The effective interpretation also matters. If a higher-derivative term is used only as a perturbative low-energy correction, treating a high-frequency pole of a truncated expression as an exact new particle may exceed the approximation's validity. If the action is proposed as fundamental, its full propagator and stability need assessment. These are different uses of the same written term.

This gives geometric harmony a definite task: derive a restricted operator pattern with a controlled low-energy regime, rather than call every possible coefficient an amplitude available to fit afterward. A successful result could predict ratios or forbid terms that a generic expansion allows. Such restrictions would make the musical analogy explanatory.

### Stability and the physical interpretation

Anselmi's fakeon prescription provides one approach to the unitarity problem by excluding selected modes from asymptotic physical states. It also involves violations of microcausality; preserving microcausality is not one of the conclusions of the cited paper (Anselmi, 2019).

The absence of an asymptotic particle is not equivalent to an independently clustering dark-matter fluid. An internal contribution to a gravitational response needs to be translated into field equations and observables. The choice of prescription, its classical limit and its scales are physical content that must be assessed.

The term 'wave forms' can describe candidate disturbances without presupposing a WIMP interpretation. It cannot establish that they are massive, stable, electromagnetically dark or able to cluster. Those properties must follow from the action and its solutions. The present paper does not identify every higher-derivative mode with a fakeon or every fakeon with dark matter.

## 8. Dark matter: a worked gravitational response

Consider the standard static weak-field point-source potential in the usual quadratic-gravity benchmark, expressed with positive inverse lengths m0 and m2 for the scalar and spin-2 contributions:

$$\Phi(r)=-\frac{GM}{r}\left[1+\frac13e^{-m_0r}-\frac43e^{-m_2r}\right].$$

The scalar contribution has coefficient +1/3; the spin-2 contribution has coefficient −4/3. This benchmark is not a derivation of an identical macroscopic potential for every fakeon prescription, background or nonlocal action. Each proposed theory needs its own response calculation.

For a circular orbit in this weak-field point-source setting, the squared circular speed is r times the radial derivative of the potential. Differentiating gives

$$v_c^2(r)=\frac{GM}{r}\left[1+\frac13(1+m_0r)e^{-m_0r}-\frac43(1+m_2r)e^{-m_2r}\right].$$

For positive fixed m0 and m2, both exponential terms vanish as r grows. Therefore the large-radius behavior is v_c squared approaching GM/r, and v_c approaches zero. This potential does not generate an asymptotically flat rotation curve for a finite point-source mass merely by containing higher derivatives.

The result does not exclude a useful fit over a finite radial interval. It also does not establish that one exists. A galactic calculation must integrate the appropriate response over an extended mass distribution, determine boundary conditions, compare the radial range actually observed and account for uncertainties in the baryonic model.

Lensing adds another requirement. In a metric theory, light deflection and nonrelativistic motions can depend on different combinations of the weak-field metric potentials. Matching a circular-speed profile is insufficient. One must derive both potentials and their relation to the matter source.

This elementary check is valuable because it turns a broad geometric suggestion into a constraint on a candidate mechanism. A different long-range behavior would require additional dynamics or a different regime, and those changes must be charged and tested explicitly.

### Extended sources and the lensing potential

The point-source expression can be extended in the standard linear static benchmark by superposition. The two-potential decomposition follows the classical, zero-logarithmic-correction limit discussed by Burzillà and colleagues (2021); the source integration and effective-density checks below are carried out explicitly here. Use the weak-field convention

$$ds^2=-(1+2\Phi)dt^2+(1-2\Psi)d\mathbf{x}^2,$$

with c = 1, an asymptotically flat background and a localized, slowly moving source of density rho. For s = |x − x-prime|, the two potentials in this benchmark are

$$\Phi(\mathbf{x})=-G\int d^3x'\frac{\rho(\mathbf{x}')}{s}\left[1+\frac13 e^{-m_0s}-\frac43e^{-m_2s}\right],$$

$$\Psi(\mathbf{x})=-G\int d^3x'\frac{\rho(\mathbf{x}')}{s}\left[1-\frac13 e^{-m_0s}-\frac23e^{-m_2s}\right].$$

The sum controlling leading weak lensing is consequently

$$\Phi+\Psi=-2G\int d^3x'\frac{\rho(\mathbf{x}')}{s}\left[1-e^{-m_2s}\right].$$

The scalar Yukawa contribution cancels in this sum although it remains in the nonrelativistic potential Phi. This is an immediate reason not to infer a lensing fit from a rotation-curve fit. The pair of equations belongs to the specified linear quadratic-gravity response; it is not asserted for arbitrary nonlocal theories or for every finalized fakeon classical equation. The physical prescription has to reproduce this response before these formulas can be used for its data analysis.

For an extended source, the Yukawa response does not obey the Newtonian shell theorem in the same way as the 1/s kernel. Integrating over the measured density is therefore necessary; replacing an entire galaxy by its total mass can erase the radial information the proposed test needs. Outside a source with bounded support and finite mass, the finite-range terms nevertheless decay at large distance. The ordinary far-field mass term survives. Thus the extended-source version does not rescue an asymptotically constant circular speed from positive fixed inverse ranges alone.

### The source equation: a calculable effective density

For the same benchmark use the distributional identities

$$\nabla^2\frac{1}{s}=-4\pi\delta^{(3)}(\mathbf{s}),\qquad
\nabla^2\frac{e^{-ms}}{s}=m^2\frac{e^{-ms}}{s}-4\pi\delta^{(3)}(\mathbf{s}).$$

Since 1 + 1/3 − 4/3 = 0, the contact terms cancel in the Laplacian of Phi. Defining an effective dynamical density by nabla^2 Phi = 4 pi G rho_eff gives

$$\rho_{\mathrm{eff}}(\mathbf{x})=\frac{1}{4\pi}\int d^3x'\frac{\rho(\mathbf{x}')}{s}\left[\frac43m_2^2e^{-m_2s}-\frac13m_0^2e^{-m_0s}\right].$$

This is a convolution, not a universal rule that a halo follows the local gradient of visible density. Its sign and spatial profile depend on both inverse ranges and the source. Moreover, the effective density is defined from one potential; it is not automatically a material fluid with its own conservation law, positive density and collisionless merger dynamics.

For positive m, the integral of m^2 exp(−ms)/(4 pi s) over all space is one. Consequently the total integral of rho_eff for a localized source is M(4/3 − 1/3) = M. This benchmark redistributes the inferred dynamical density but does not create an extra total far-field mass. The effective excess rho_eff − rho integrates to zero. It may matter over a finite radial interval; it cannot be described as a derived positive halo of arbitrarily larger total mass.

That result is a useful failure condition rather than a verdict on all geometric alternatives. To obtain a different total or asymptotic behavior one must change a stated ingredient: the kernel, boundary conditions, background, source content or regime. The change then belongs in the model and in its specification cost.

### A reproducible first observational test

For a disk galaxy, take a finite baryonic density model constrained by photometry and gas data. Evaluate Phi in the disk plane using the source integral, differentiate it to obtain v_c^2(R) = R partial_R Phi, and compare at the actual measured radii. Use consistent units and state whether m_0 and m_2 are shared across all galaxies. A negative v_c^2 where a circular orbit is assumed is a model failure for that proposed orbit, not a parameter value to silently omit.

Then evaluate the line-of-sight derivative of Phi + Psi for a system with lensing data. Use the same action and common parameters; do not fit an independent halo to the lensing calculation after claiming that geometry removed it. Stellar mass-to-light ratios, distances, inclination and gas modeling are nuisance quantities whose uncertainties must be carried in the likelihood. A successful fit to one profile is not the joint test.

The included calculation script checks algebraic and limiting claims, not galaxy data. No SPARC likelihood has been run in this revision. The advance is a specified response and a reproducible failure test, which make the next data calculation possible without pretending it has already happened.

## 9. Compare complete geometric and particle accounts

Fix the scope of the comparison: observables, precision, preparation conditions and the treatment of boundary data. Count the geometry, fields, equations of motion, starting or boundary conditions, and rule connecting predictions to observations together. Their joint description length is the length of a program specifying the whole account. It can be shorter than separate descriptions of each ingredient because one ingredient may determine another.

This supplies a concrete research target. If one construction derives relations among couplings that a rival must specify separately, those relations can reduce the joint cost. Count the derivation once and do not then charge its consequences again. Conversely, an unspecified action cannot receive a zero cost simply because a rich configuration space has been defined.

String theory and a metric-bundle proposal must face this same comparison. A generic fixed-length index into N alternatives uses ceil(log2 N) bits, but a particular member may have a short selector. Neither a landscape count nor the size of a software implementation determines the minimum description of a complete candidate.

Under a stipulated description prior and equal observational likelihoods, a demonstrated delta-bit saving gives a factor 2^delta in prior and posterior odds. Establishing the saving is the outstanding comparison, rather than a conclusion available from the coordinate counts alone.

A WIMP mass, interaction strength and production history are not necessarily independent knobs; a particular particle model can relate them. An axion's mass and couplings can likewise be related by its underlying construction. Assigning each a generic fifty or one hundred bits while treating a new gravitational kernel as free does not compare complete explanations. The proposed metric theory must pay for m_0, m_2, its action and its prescription at whatever precision is required for the same observations.

Conversely, a geometric theory should receive credit if a single set of coefficients explains observations for which a comparison model needs independently fitted freedom. The relevant quantity is the full predictive code or evidence, not a prejudice for particles or against them. The fixed jet-construction algorithm can be shared across geometric candidates without being confused with the physical content of every candidate.

For a practical two-part code, define L(T,theta) for the laws and specified parameters and add minus log_2 P(D|T,theta) for the residual data under an explicitly normalized observation model. This is one approximation to model comparison, not a computation of universal K. If nuisance parameters are integrated out instead, state their priors. A numerical superiority claim needs actual candidate codes and likelihoods; this paper supplies neither a measured bit advantage nor posterior odds.

## 10. From the construction to a physical theory

The worked response identifies what the next model must supply. Specify an action and background, derive the independent propagating modes, and check whether disturbances grow without bound or carry an unacceptable sign of energy. Then determine their coupling to matter and which parameters the construction fixes.

The first observational target is a joint account of orbital motion and lensing using the same coefficients, as set out in Section 8. The SPARC database provides resolved rotation curves and photometry for the orbital test (Lelli, McGaugh and Schombert, 2016). A successful result would connect observations through shared dynamics; fitting an unrelated halo afterward would remove that claimed economy.

A galaxy fit would not by itself establish agreement with cluster mergers, cosmology or gravitational waves. Those extensions require solutions in their respective regimes. They are subsequent tests, not additional results claimed here. The present contribution is the geometric hierarchy and a worked gravitational response that exposes a concrete limitation of the finite-range example.

## 11. Conclusion

The fourteen-dimensional metric bundle and its finite jets follow from specified geometric constructions. Their growing coordinate lists can share a short generating rule. An action, constraints and a background then determine the actual physical modes.

The quadratic-gravity calculation supplies a definite test. With the stated positive inverse ranges, its additional forces decay at large distance, its effective total mass equals the source mass, and its lensing potential differs from its orbital potential. A finite-range fit may still be investigated, but the construction alone is not a derived dark halo.

Geometric unification becomes explanatory when one action supplies relations among particles, interactions and observations that otherwise require independent inputs. The coordinate counts, mode calculation and source response specify the next derivation and comparison.

## References

Weinstein, E. (2021). *Geometric Unity: Author’s Working Draft, v1.0*. https://geometricunity.nyc3.digitaloceanspaces.com/Geometric_Unity-Draft-April-1st-2021.pdf

Gorard, J. (2020). Some relativistic and gravitational properties of the Wolfram model. https://arxiv.org/abs/2004.14810

Saunders, D. J. (1989). The Geometry of Jet Bundles. Cambridge University Press.

Alvarez-Gaume, L., Kehagias, A., Kounnas, C., Lüst, D., and Riotto, A. (2015). Aspects of Quadratic Gravity. https://arxiv.org/abs/1505.07657

Solomonoff, R. J. (1964). A formal theory of inductive inference. Information and Control 7, 1-22 and 224-254.

Burzillà, N., Giacchini, B. L., de Paula Netto, T., and Modesto, L. (2021). Newtonian potential in higher-derivative quantum gravity. https://arxiv.org/abs/2012.06254

Anselmi, D. (2019). Fakeons, Microcausality and the Classical Limit of Quantum Gravity. Classical and Quantum Gravity 36, 065010. https://arxiv.org/abs/1809.05037

Lelli, F., McGaugh, S. S., and Schombert, J. M. (2016). SPARC: Mass models for 175 disk galaxies with Spitzer photometry and accurate rotation curves. Astronomical Journal 152, 157.

Planck Collaboration (2020). Planck 2018 results VI. Cosmological parameters. Astronomy and Astrophysics 641, A6.

Stelle, K. S. (1977). Renormalization of higher-derivative quantum gravity. Physical Review D 16, 953.

Wiltshire, D. L. (2007). Cosmic clocks, cosmic variance and cosmic averages. New Journal of Physics 9, 377.

Wiltshire, D. L. (2009). Average observational quantities in the timescape cosmology. Physical Review D 80, 123512.

De Vecchi, F. C. (2017). Finite dimensional solutions to SPDEs and the geometry of infinite jet bundles. arXiv:1712.08490. https://arxiv.org/abs/1712.08490
