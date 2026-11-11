A42 Pass 36 with the curvature-squared coupling restored (30 October 2026)

Model: closed FLRW, square-root NPQT gravity (Bueno, Cano, Hennigar and Murcia 2025) with the
four-dimensional conformal vacuum stress of free fields (Astra Pass 36, lambda = 0.880, dust constant
3.2237), plus (1/2 kappa) beta R^2. Units alpha = 1; b = beta/alpha. b > 0 is the stable sign.

derive2.py   derives the R^2 contribution to the Friedmann constraint from the minisuperspace action
             (sympy) and checks it against Einstein-Hilbert and flat de Sitter; writes QP.pkl.
r2_bounce.py evolves the fourth-order system: the time derivative of the constraint divided by H,
             regular through the bounce; the constraint is monitored, not imposed (held to about 1e-11
             until the action's branch edge). Initial data: Pass 36's midpoint state, with the third
             derivative fixed by the constraint.
robust.py    repeats the b = 0.002 and 0.005 cases with two integrators.
r2_scan.json results for b from -0.1 to 1.

Result: the bounce survives for 0 < b <= 0.002, with peak curvature K alpha^2 rising from 20.9 to 270;
at b >= 0.005, and for every negative b tested (-1e-6 to -0.1), the contraction reaches the edge of
the square-root action (alpha Phi -> 1) before it can bounce. The R^2 coupling is a free renormalized
parameter; this constrains it within the model. Free massless fields only; homogeneous sector only.
