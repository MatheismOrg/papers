"""Transfer test: the frozen 451-byte A42_Step_Rule.py, unchanged, on a background it was not built or checked for:
a 1+1 universe that expands in conformal time, a(eta) = 1.75 + 1.25 tanh(eta/tau), mass 0.7, momentum 1, tau = 0.15.
The program is given a time-dependent mass angle each step; nothing else changes. A negative-energy mode at the start
is evolved on a ring of n sites, and the probability that it ends with positive energy (a created pair) is compared
with the continuum, integrated independently with scipy."""
import sys, os, json, math, time
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A42_Step_Rule import step
AI, AF, MASS, BOUND = 0.5, 3.0, 0.7, 10.0
X = np.array([[0, 1], [1, 0]], complex); Z = np.diag([1, -1]).astype(complex)
def scale(eta, tau): return (AI + AF) / 2 + (AF - AI) / 2 * math.tanh(eta / tau)

def spinors(s, k, A):
    w, v = np.linalg.eigh(s * k * Z + MASS * A * X)
    return v[:, 0], v[:, 1]          # negative energy, positive energy

def continuum(s, k, tau):
    def rhs(e, y): return -1j * (s * k * Z + MASS * scale(e, tau) * X) @ y
    neg, _ = spinors(s, k, scale(-BOUND * tau, tau))
    sol = solve_ivp(rhs, (-BOUND * tau, BOUND * tau), neg.astype(complex), method='DOP853', rtol=1e-12, atol=1e-14)
    _, pos = spinors(s, k, scale(BOUND * tau, tau))
    return abs(np.vdot(pos, sol.y[:, -1])) ** 2

def lattice_sign(n=64, k=1):
    """The sign s for which the program's free evolution of a right-mover is exp(-i s k t)."""
    dx = 2 * math.pi / n; dt = 0.25 * dx; x = np.arange(n) * dx
    F = [list(np.exp(1j * k * x)), [0j] * n]
    for j in range(2): step(F, [[1.0] * n, [-1.0] * n], None, None, dt, dx, j)
    ph = np.angle(np.vdot(np.exp(1j * k * x), np.array(F[0])))
    return 1 if ph < 0 else -1

def lattice(n, k, tau, s):
    dx = 2 * math.pi / n; start, end = -BOUND * tau, BOUND * tau
    steps = int(math.ceil((end - start) / (0.5 * dx))); steps += steps % 2; dt = (end - start) / steps
    x = np.arange(n) * dx; wave = np.exp(1j * k * x) / math.sqrt(n)
    neg, _ = spinors(s, k, scale(start, tau))
    F = [list(neg[0] * wave), list(neg[1] * wave)]
    V = [[1.0] * n, [-1.0] * n]
    for j in range(steps):
        A = scale(start + (j + 0.5) * dt, tau)
        step(F, V, None, [0.5 * dt * MASS * A] * n, dt, dx, j)
    _, pos = spinors(s, k, scale(end, tau))
    u, w = np.array(F[0]), np.array(F[1])
    amp = np.vdot(pos[0] * wave, u) + np.vdot(pos[1] * wave, w)
    norm = np.vdot(u, u).real + np.vdot(w, w).real
    return abs(amp) ** 2, norm, steps

if __name__ == '__main__':
    tau, k = 0.15, 1
    s = lattice_sign()
    ref = continuum(s, k, tau)
    out = {'program_sha256': __import__('hashlib').sha256(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'A42_Step_Rule.py'), 'rb').read()).hexdigest(),
           'tau': tau, 'k': k, 'mass': MASS, 'sign': s, 'continuum_pair_probability': ref, 'runs': []}
    sizes = [int(a) for a in sys.argv[1:]] or [64, 128, 256, 512, 1024, 2048]
    for n in sizes:
        t0 = time.time(); P, norm, steps = lattice(n, k, tau, s)
        row = {'sites': n, 'steps': steps, 'pair_probability': P, 'error': P - ref, 'norm_error': norm - 1, 'sec': round(time.time() - t0, 1)}
        out['runs'].append(row); print(json.dumps(row), flush=True)
        json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'transfer_expanding_results.json'), 'w'), indent=1)
    errs = [abs(r['error']) for r in out['runs']]
    out['error_ratios'] = [errs[i] / errs[i + 1] for i in range(len(errs) - 1) if errs[i + 1] > 0]
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'transfer_expanding_results.json'), 'w'), indent=1)
    print('continuum', ref, 'ratios', out['error_ratios'])
