A43 additional fixed-radius cluster transfer

Run from the complete A43 supplement:
  python dm_unified_stress/run.py
Dependencies: NumPy, SciPy. Uses sibling dm_unified_limits and dm_fixed_radius.

All ten models are frozen at their all-data, five-population coefficients.
No X-COP target trains any parameter here. The seven clusters/21 profile points
were inspected in earlier work, so this is an exploratory transfer stress test,
not a new independent discovery sample. See PROTOCOL.txt and the original
fixed-radius README for source-model assumptions and the nine sensitivity cases.
The inputs are correlated hydrostatic reconstructions; star contributions are
approximate. No universal-force or formal model-exclusion inference follows.
