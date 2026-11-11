"""Make explicit local revisions to T18 and insert the independently generated audit."""
from pathlib import Path
import json
D=Path(__file__).resolve().parent;root=D.parent
s=(root/'T18_ORIGINAL/PAPER/A42_Shortest_Rules_Search.md').read_text();changes=[]
def replace(old,new):
 global s
 assert s.count(old)==1,(old[:100],s.count(old))
 s=s.replace(old,new);changes.append(dict(old=old,new=new))
replace('# The Theory of Everything Must Be Short and Discrete: A Search of Rewriting Rules for Four-Dimensional Spacetime','# Searching for a Short Discrete Theory of Everything: Rewriting Rules and Supplied-Lattice Tests\n\nGary Abraham Bernstein  \nORCID: 0009-0009-1761-2867  \nT19 research revision, October 31, 2026')
replace('with chains 4.7 to 21 standard deviations short of four dimensions','with chain-based dimension estimates 4.7 to 21 standard deviations below the four-dimensional reference')
replace("The full known laws, set on a supplied four-dimensional lattice and completed with dark matter and an initial state, compress to 16,272 bits; their lattice quarks keep an exact chiral symmetry, and the Planck spacing lies 43 million times inside the LHAASO Lorentz bound. Its open test is smooth spacetime from the lattice's own quantum gravity.","A proposed specification of the known laws on a supplied lattice, with declared dark matter and initial data, compresses to 16,272 bits. Component tests and a new 24-chain Regge audit distinguish description length from demonstrated dynamics; a complete quantum-gravity construction remains open.")
replace('Then the short rules carry most of the weight, and a search should start with them.','This weights individual short programs more strongly and motivates searching them first; it does not establish that the physical rule must be short. The discreteness argument supplies the premise of the search; the numerical tests below do not independently prove that premise.')
replace('Section 8.7 sets the known laws in full on a supplied four-dimensional lattice, completed with dark matter and an initial state, at 16,272 bits.','Section 8.7 records a proposed specification of the known laws on a supplied four-dimensional lattice, with declared dark matter and initial data, at 16,272 bits. This description length is distinct from the length of a verified executable completion.')
replace('the enumeration from the bottom raises the floor with each class it rules out under its tests.','the enumeration from the bottom excludes candidates only within the classes, initial states, schedules and finite sizes it tests.')
replace("A fixed physical spacing caps it, so in a discrete world the curvature stays finite, and what matter does at the center must come from the rule.","A fixed physical spacing caps this sampled expression away from r = 0; it does not establish a regular center or bound the fields evolved by a discrete rule. What matter does at the center must come from the rule.")
replace('A regular bounce needs that coupling small and positive.','For the initial conditions and coupling values tested, the regular bounce survives only at the tested small positive couplings. The scan does not establish a universal necessity or a sharp threshold.')
replace('Discreteness alone removes the infinities of point particles and the ultraviolet divergences of quantum fields.','A specified lattice can regulate short-distance divergences. Discreteness alone does not guarantee bounded fields, a convergent path integral or a dynamically regular center.')
replace('That settles the strong force.','This verifies the overlap construction in the tested U(1) backgrounds; it does not test four-dimensional SU(3) dynamics.')
replace('The weak force is weak, so a program can compute its effects in perturbation theory to any measured precision, and the cap below holds however this ends.','Perturbative calculations can reproduce many measured weak-interaction observables. They do not complete the proposed nonperturbative chiral lattice theory; the count below measures its specification.')
replace("so fine-tuning, like the strong CP problem's θ = 0, costs no bits.","so the renormalization condition has a short description. The target values, their precision and a procedure for determining the bare parameters must still be specified; compression does not explain those values.")
replace('so no doubled modes; on the supplied four-dimensional lattice','in the measured low spectrum; those counts alone do not exclude every extra propagating species. On the supplied four-dimensional lattice')
replace('Since the Standard Model\'s masses lie far below the lattice scale, its fields act back on the geometry mainly through their numbers of each spin, so a first test, a proxy for the full matter sector, needs only free fields in those numbers.','Free conformal fields in the Standard Model\'s spin counts give a tractable first proxy for quantum stress. Masses, curvature couplings, interactions and the quantum state remain necessary for a full matter calculation.')
replace('an upper bound in this encoding on a complete statement of the known laws with dark matter and an initial state.','a description length in this encoding for the proposed specification with dark matter and initial data.')
replace('A theory of everything can therefore be capped from the top, bounded from below within classes and tested against data; it cannot be proved minimal.','A demonstrated program gives an upper bound for its stated output and tests, while the present census excludes candidates within its tested classes. These are different bounds; neither establishes the minimum length of a complete physical theory.')
replace('| Total: the cap | | 16,272 | upper bound in this encoding |','| Total: specification | | 16,272 | description length; full dynamics unverified |')
replace('| Fine-tuning of Λ and the Higgs mass | the condition that the measured values result | 0 (as digits, 514) | puzzles of chance |','| Fine-tuning of Λ and the Higgs mass | renormalization condition and declared targets | included in description (digit estimate: 514) | parameters not explained |')
replace('Written out in full on a supplied four-dimensional lattice, with dark matter and an initial state, the known laws compress to 16,272 bits, enough for present data within its known tensions; they join general relativity with quantum theory if the lattice\'s own quantum gravity yields smooth four-dimensional spacetime, the test that remains.','The proposed specification on a supplied four-dimensional lattice, with dark matter and initial data, compresses to 16,272 bits. Component tests support parts of this construction; they do not establish a single executable model covering present data. Smooth large-scale geometry, a consistent chiral matter sector and the quantum-gravity dynamics remain separate requirements.')
replace('A physical lattice makes the path integral finite at the price of fine-tuning:', 'A physical lattice regulates short-distance modes without by itself establishing convergence of the gravitational integral. The parameter estimates also expose fine-tuning:')
replace('so the statement covers present data.', 'which motivates the low-energy target; an empirical fit of the complete supplied-lattice construction has not been performed here.')
replace('the supplement finds a smooth phase after a flat start', 'the original supplement finds a long-lived negative-average-curvature regime after a flat start')
replace('so the run locates the edge of the smooth phase, not the critical point.', 'so the run locates an apparent finite-run stability boundary, not a critical point. The independent execution in Section 8.9 tests its dependence on initial conditions and sampling history.')
replace('**The cap.**', '**Description length.**')
replace('An upper bound is what this route can prove. The search from the bottom proves lower bounds only within the classes it runs (Sections 4 to 6),', 'A compressed description alone does not establish an executable physical model. The search from the bottom excludes candidates only within its tested classes and conditions (Sections 4 to 6),')
replace('**Table 12.** The known laws on a supplied four-dimensional lattice:', '**Table 12.** The proposed specification on a supplied four-dimensional lattice:')
replace('strong force exact; weak force proposed, disputed', '2D U(1) overlap checks; chiral completion unresolved')
replace('matches present data, up to known tensions', 'declared inputs; no new full-data fit')
replace('open; smooth phase here up to k ≈ 0.07, critical point on 32⁴ lattices', 'open; finite-run regimes here; see Section 8.9; larger-lattice critical behavior in cited work')
replace('takes 451 bytes, or 2,488 bits after the raw LZMA2 compression', 'takes 451 source bytes, or 2,488 bits conditional on the Python interpreter, numerical dependencies and the specified compression decoder, after the raw LZMA2 compression')

replace("A rule's length in bits is the number of symbols it uses, one for each node of each relation, times the bits each symbol needs, the base-two logarithm of the number of distinct variables.",'The rule-length score is the number of node symbols times the base-two logarithm of the number of distinct variables. It is a structural score, not a self-delimiting executable encoding; the scheduler, decoder and initial-state costs are additional. It is therefore not directly comparable to the compressed byte counts in Sections 7 and 8.')
result=json.loads((D/'results.json').read_text())
nfail=sum(not g['passes_rhat_screen'] for g in result['groups'])
fail_label='All four' if nfail==4 else str(nfail)+' of four'
replace('Component tests and a new 24-chain Regge audit distinguish description length from demonstrated dynamics;', f'{fail_label} groups in a new 24-chain Regge audit fail the mixing screen;')
replace('Harder tests sharpen the search:', f'The independent Regge audit fails the prespecified mixing screen in {nfail} of four parameter groups, so the sampled regime is not established as an equilibrated smooth phase. Harder tests sharpen the search:')

# Audit text is generated from the new results; no numerical outcomes are preset here.
audit=(D/'paper_audit_section.md').read_text()
replace('## 9. Conclusion',audit+'\n\n## 9. Conclusion')
(root/'PAPER').mkdir(exist_ok=True)
(root/'PAPER/A42_T19.md').write_text(s)
(root/'PAPER/CHANGES.json').write_text(json.dumps(changes,indent=2,ensure_ascii=False)+'\n')
(root/'PAPER/CHANGES.md').write_text('# T19 changes to T18\n\nEvery original file is retained in T18_ORIGINAL.\n\n'+'\n\n'.join(f'## Change {i+1}\n\nPrevious:\n\n{c["old"]}\n\nRevised:\n\n**{c["new"]}**' for i,c in enumerate(changes)))
print('Applied',len(changes),'checked replacements')
