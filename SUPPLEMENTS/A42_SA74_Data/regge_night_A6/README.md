# Plan A6: the 24⁴ and 32⁴ chains of A42, Section 8.10 (paper version T37)

Per-sweep histories of the five chains of plan A6, as they ran on an Apple M4 with the multi-core engine (November 9 to 11, 2026):

- `A6_L20_k0.244_r2.dat`: 20⁴, k = 0.244, flat start, 2,500 + 3,000 sweeps
- `A6_L24_k0.244.dat`: 24⁴, k = 0.244, flat start, 2,500 + 3,000 sweeps
- `A6_L24_k0.248.dat`, `A6_L24_k0.252.dat`: 24⁴, each started from the previous chain's final configuration, 1,500 + 3,000 sweeps
- `A6_L32_k0.244.dat`: 32⁴, k = 0.244, flat start, 2,500 + 3,000 sweeps; `A6_L32_k0.244.log` is its run log

`queue3.txt` is the plan file; plan A6 was appended to it before the runs. The columns of the histories are named in each file's header.
The analysis is `analyze_a1.py` in ../../A42_Reproduce (sa74/regge/r2/night), run with `MHFILE=multihist_T30.json OUTSUFFIX=_T37`.
