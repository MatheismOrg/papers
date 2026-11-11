"""A42 bounce test with the curvature-squared coupling restored.
Closed FLRW, square-root NPQT gravity + free conformal 4D quantum matter (Astra Pass 36) + b R^2.
Units alpha=1; b = beta/alpha, where the action gains (1/2kappa) beta R^2. b>0 is the stable (Starobinsky) sign.
The R^2 terms come from the minisuperspace action (derive2.py); the evolution is the time derivative of the
constraint divided by H, which is regular through the bounce; the constraint is monitored, not imposed."""
import sys, os, json, math, pickle, numpy as np, sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.integrate import solve_ivp
from pass36_model import F, Fp, turning_points, y_from_a, rhs as rhs2
Qs, Ps = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QP.pkl'), 'rb'))
a0s, a1s, a2s, a3s, a4s, ks, bs = sp.symbols('a0 a1 a2 a3 a4 k b')
Qexpr = sp.sympify(Qs).subs(ks, 1); Pexpr = sp.sympify(Ps).subs(ks, 1)
Qf = sp.lambdify((a0s, a1s, a2s, a3s, bs), Qexpr, 'numpy')
P0f = sp.lambdify((a0s, a1s, a2s, a3s, bs), Pexpr.subs(a4s, 0), 'numpy')
cQ3 = sp.lambdify((a0s, a1s, bs), sp.diff(Qexpr, a3s), 'numpy')        # coefficient of a3 in Q
Q_no3 = sp.lambdify((a0s, a1s, a2s, bs), Qexpr.subs(a3s, 0), 'numpy')
LAM = 0.88021525471101; C = 3.22370979547
def rhs4(t, z, b, lam=LAM, c=C):
    a, a1, a2, a3 = z; y = (1 + a1 * a1) / a ** 2
    if not 0 <= y < 1: raise ValueError('left the real branch, y>=1')
    fp = Fp(y, lam)
    a4 = -a / (12 * b) * (2 * fp * (a2 / a - y) + 3 * c / a ** 3 - P0f(a, a1, a2, a3, b))
    return [a1, a2, a3, a4]
def constraint(z, b, lam=LAM, c=C):
    a, a1, a2, a3 = z; y = (1 + a1 * a1) / a ** 2
    return F(y, lam) - Qf(a, a1, a2, a3, b) - c / a ** 3
def run(b, lam=LAM, c=C):
    amin, amax = turning_points(c, lam); a = (amin + amax) / 2
    y = y_from_a(a, c, lam); a1 = -math.sqrt(a * a * y - 1)
    a2 = a * (y - 3 * c / a ** 3 / (2 * Fp(y, lam)))
    out = dict(b=b, amin_b0=amin, amax_b0=amax)
    if b == 0:
        def bb(t, z): return z[1]
        bb.direction = 1
        def back0(t, z): return z[0] - a
        back0.direction = 1; back0.terminal = True
        sol = solve_ivp(lambda t, z: rhs2(t, z, c, lam), (0, 100 * (amax + 1)), [a, a1], method='DOP853', rtol=1e-11, atol=1e-12, events=[bb, back0], dense_output=True)
        T = np.linspace(0, sol.t[-1], 20001); A, V = sol.sol(T); Y = (1 + V * V) / A ** 2; D = Y - 3 * c / A ** 3 / (2 * Fp(Y, lam))
        out.update(bounced=len(sol.t_events[0]) == 1, returned=len(sol.t_events[1]) == 1, a_min=float(A.min()), K_max=float((12 * (Y * Y + D * D)).max()), y_max=float(Y.max()), t_end=float(sol.t[-1]))
        return out
    a3 = -(Q_no3(a, a1, a2, b) - (F(y, lam) - c / a ** 3)) / cQ3(a, a1, b)   # constraint fixes a3
    z0 = [a, a1, a2, a3]
    def bounce(t, z): return z[1]
    bounce.direction = 1
    def back(t, z): return z[0] - a
    back.direction = 1; back.terminal = True
    def edge(t, z): return 1 - 1e-9 - (1 + z[1] ** 2) / z[0] ** 2
    edge.terminal = True
    def tiny(t, z): return z[0] - 1e-3
    tiny.terminal = True
    try:
        sol = solve_ivp(lambda t, z: rhs4(t, z, b, lam, c), (0, 60 * (amax + 1)), z0, method='Radau', rtol=1e-10, atol=1e-12,
                        events=[bounce, back, edge, tiny], dense_output=False)
        status = sol.message
    except Exception as e:
        return dict(out, bounced=False, failure=str(e))
    A, A1, A2, A3 = sol.y; Y = (1 + A1 ** 2) / A ** 2; D = A2 / A
    res = np.array([constraint(sol.y[:, i], b, lam, c) for i in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 400))])
    out.update(bounced=len(sol.t_events[0]) >= 1, returned=len(sol.t_events[1]) == 1, hit_branch_edge=len(sol.t_events[2]) > 0,
               collapsed=len(sol.t_events[3]) > 0, n_bounces=len(sol.t_events[0]), a_min=float(A.min()), K_max=float((12 * (Y * Y + D * D)).max()),
               y_max=float(Y.max()), max_constraint_residual=float(np.max(np.abs(res))), t_end=float(sol.t[-1]), solver=status)
    return out
if __name__ == '__main__':
    bs_ = [float(x) for x in sys.argv[1:]] or [0, 1e-3, 1e-2, 0.1, 0.5, 1, 3, -1e-3, -1e-2, -0.1, -1]
    R = [run(b) for b in bs_]
    for r in R: print(json.dumps(r))
    json.dump(R, open('r2_scan.json', 'w'), indent=1)
