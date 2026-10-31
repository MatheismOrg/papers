Globular-cluster structural control
==================================

Run: python build_control.py
Dependencies: Python 3, numpy. Network is needed only if the two frozen inputs
are absent. Saved inputs and their SHA-256 hashes permit exact reproduction.

Outputs
- globular_cluster_controls.csv: all 171 rows in the primary author table;
  descriptive_control_selected marks the 151 rows used in the summary.
- globular_summary.json: selection, descriptive statistics, provenance, hashes.

Primary author sources, downloaded 27 October 2026:
https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/parameter.html
https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/combined_table.txt

This is the updated author database, not the original 112-row 2018 publication
table. The page's introductory number (169) is stale relative to the 171 data
rows actually downloaded. Methods originate in Baumgardt & Hilker (2018),
MNRAS 478, 1520-1557, DOI 10.1093/mnras/sty1057, arXiv:1804.08359. The author
site also credits later distance, photometry, mass-function and kinematic
updates; its saved HTML contains individual source tags.

The descriptive subset requires a blue fitted-system object link, a reported
positive mass uncertainty, and at least 20 radial-velocity measurements or 20
proper-motion measurements. The source HTML marks luminosity-estimated masses
with different object-link styling; these are not treated as independent
kinematic determinations. The script reads that flag from the name link, not
from italic aliases or estimated mass-function slopes.

Selected objects have median V-band dynamical-model M/L 1.81 solar units;
the 16th-84th percentile interval is 1.41-2.42. Their projected half-light
radii range from 1.05 to 29.04 pc. 133/151 have M/L between 1 and 3.
These statistics show the different scale and modest typical M/L of this
population, but are not a quantitative dark-matter exclusion.

NO RAR FIT OR TRANSFER ERROR IS CALCULATED FOR THESE ROWS.
The table's sigma0 is central, mass-weighted one-dimensional velocity
dispersion, not the luminosity-weighted global dispersion required by a Wolf
half-light estimator. In addition, half-mass and half-light radii differ due
to mass segregation, and N-body masses already rely on a stellar dynamical
model. Reusing those masses as independent baryonic and dynamical inputs
would be circular. A defensible future RAR test needs the full kinematic and
surface-brightness profiles plus independently constrained stellar/remnant
mass distributions, binaries, tidal effects, and external gravitational field.

Do not include the 151 controls in the number of independently fitted
acceleration systems or in a statement that one formula fits all dark matter.
