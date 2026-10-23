# A38 implementation audit

The earlier comparison did not establish that all known physics takes 9,216 bits or that string theory requires an irreducible extra 1,661 bits. Original sources and claims remain in the archive.

## What the files contain

The C example contains 27 illustrative constants despite its original label of 26. Its CKM routine read the wrong array positions and omitted the imaginary components of the matrix. The corrected example uses complex entries and prints the result so it can be checked. Its historical inputs are not advertised as a current precision fit.

GR and quantum equations in that file are comments. A compiler removes comments. The program does not solve those equations or recover a Standard Model prediction from arbitrary physical inputs. Extracting only .text and .data also omits other executable sections, runtime information and libraries. The old procedure did not count a complete self contained decoder.

The Python file produces symbolic descriptions requiring an external interpreter. Its minimal Standard Model entry now explicitly uses massless neutrinos and nineteen physical inputs including the strong CP angle. This is a different model scope from the C data list.

The earlier string C file added arrays without deriving a consistent vacuum or its low energy predictions. It is retained as an archived historical artifact. Its active entrypoint now refuses compilation so it cannot silently regenerate a misleading benchmark.

## What a package size establishes

The earlier CYTools and JAXVacua sizes were reported implementation measurements for software performing different tasks. They were not independently reproduced in this pass and are not lower bounds on string theory. No 616 fold physical complexity ratio is retained.

For a catalogue of 10^500 entries, a generic fixed length address needs 1,661 bits given the catalogue. A particular entry may have a shorter description. Using an address can replace a parameter list; charging both can double count the same information.

## Next benchmark

Specify a shared predictive task, inputs, output precision and decoder. Include libraries and data, test outputs against independent reference results and state what is excluded. Compare the resulting complete implementations as upper bounds, not as exact Kolmogorov complexities. Report sensitivity to languages and tolerances before assigning a prior ratio.

Verification in this release checks the corrected CKM matrix for finite entries and unitarity, not the correctness of a theory of everything. See G47_15_0_VERIFICATION.json in REVIEW.

Parameterization source: Particle Data Group, CKM Quark Mixing Matrix, https://pdg.lbl.gov/2025/reviews/rpp2025-rev-ckm-matrix.pdf.
