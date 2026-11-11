# A42 Supplement: Known Laws as Local Steps

These scripts back Section 8 of "Searching for a Short Discrete Theory of Everything: Rewriting Rules and Supplied-Lattice Tests". Each test supplies a geometry, fields and a starting state, runs one family of local unitary steps, and compares the result with an independent continuum calculation. Every script writes its results as JSON beside itself and checks its own assertions.

## Requirements

Python 3 with NumPy and SciPy; matplotlib only for the one plot. Each script runs in a few seconds. Run each from its own folder, for example `cd A42_Supplied_Physics/horizon && python3 horizon_test.py`.

## The tests

- `massive/massive_curved.py`: the massive Dirac field on the static metric $ds^2 = -N(x)^2dt^2 + dx^2$ with $N = 1 + 0.25\cos x$, against an independent Fourier reference. For mass 2 and packet center 0.9 the state error is $9.1\times10^{-4}$ at 1,024 sites, about four times smaller per halving of the spacing. A flat-space case is included.
- `geodesic/circuit_geodesic.py` and `geodesic/geodesic_packet.py`: a falling packet against its quantum reference and the classical geodesic. At mass 64 and 8,192 sites the circuit center is 1.690468, the quantum reference 1.690566 and the geodesic 1.692566.
- `gauge/gauge_ring.py`: prescribed electromagnetic phases on a ring; gauge covariance to about $10^{-15}$ and two-arm interference probability 0.206222 against 0.206107.
- `qft/fock_mass_quench.py`: a free-fermion mass quench from 0.4 to 3 with eight modes and a full Fock space of 256 states; many-fermion state error $3.4\times10^{-5}$ at 2,048 sites and mean pair count 0.646932 against the analytic 0.646951.
- `maxwell/maxwell_gravity.py`: transverse Maxwell plane waves on a curved plane-symmetric metric; relative field error about $7\times10^{-5}$ at 2,048 sites, conserved energy to about $3\times10^{-16}$.
- `horizon/horizon_test.py`: radial massless transport across a Schwarzschild horizon in Painlevé-Gullstrand coordinates, $r_s = 2$, three packets at six resolutions from 256 to 8,192 sites; errors fall 3.8 to 4.0 times on the last halving, and the norm holds to $10^{-11}$. It also computes the curvature at the first site of a grid that stops short of the center, which grows 64-fold per halving.
- `audit/`: independent checks of the mass coupling, the geodesic comparison and the quench overlaps, including the wrong rules that fail: replacing the mass term $mN$ by $m$, dropping the gauge links, and evolving with the old mass after the quench.

## Scope

The tests supply the metric, the field content and the quantum state space. They test whether local steps reproduce known equations at identifiable scales; they do not grow the geometry, and they include no interacting fields or matter acting back on the geometry. Earlier work on the same kind of agreement: D'Ariano and Perinotti (2014) for the Dirac equation, Di Molfetta, Brachet and Debbasch (2013) for quantum walks at a horizon, and Eon, Di Molfetta, Magnifico and Arrighi (2023) for quantum electrodynamics.
