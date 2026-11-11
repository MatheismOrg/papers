"""SA73 T46: self-interacting dark matter (SIDM) halos with baryons, the isothermal Jeans model of Kaplinghat, Tulin and
Yu (2016, PRL 116, 041302) and Ren et al. (2019, PRX 9, 031020): inside r1, where a particle has scattered about once in
t_age [(sigma/m) <v> rho_NFW(r1) t_age = 1, <v> = (4/sqrt(pi)) sigma0], the dark matter is isothermal in the total
potential, sigma0^2 dln(rho)/dr = -G (M_dm + M_b)/r^2; outside it is the NFW halo; density and enclosed mass match at r1.
Units: kpc, Msun, km/s. Baryons spherical (as in Ren et al.'s MCMC fits). Fixed-step RK4 in ln r, 3-D root find."""
import numpy as np
from scipy.optimize import root
G = 4.30091e-6                    # kpc (km/s)^2 / Msun
CM2G = 2.0893e-10                 # 1 cm^2/g in kpc^2/Msun
KMS_KPC_GYR = 1.02271             # 1 km/s/kpc in 1/Gyr
NS = 240; T = np.linspace(np.log(1e-4), 0., 2 * NS + 1)   # grid in ln(r/r1), with half steps
def ff(x): return np.log1p(x) - x / (1 + x)
class NFW:
    def __init__(s, mdm200, r200, c):
        s.rs = r200 / c; s.rhos = mdm200 / (4 * np.pi * s.rs ** 3 * ff(c))
    def rho(s, r): x = r / s.rs; return s.rhos / (x * (1 + x) ** 2)
    def mass(s, r): return 4 * np.pi * s.rhos * s.rs ** 3 * ff(r / s.rs)
def xs_vdep(vmean, s0=2.0, w=344.0):
    """sigma/m in cm^2/g: s0 at low speed, falling as 1/(1 + (v/w)^2); with s0 = 2, w = 344 km/s it is 0.1 at <v> = 1500 km/s
    (Kaplinghat, Tulin and Yu 2016: 1.9 cm^2/g for dwarfs and LSB galaxies, 0.10 for clusters)."""
    return s0 / (1 + (vmean / w) ** 2)
def _integ(lrho0, lsig0, r1, mbg):
    """RK4 from r = 1e-4 r1 to r1 in x = ln r. Returns (h(r1), m(r1)) and the grid solution for m."""
    rho0 = np.exp(lrho0); s2 = np.exp(2 * lsig0); rg = r1 * np.exp(T); dx = T[2] - T[0]
    k = 4 * np.pi * rho0
    h = 0.; m = k / 3 * rg[0] ** 3; ms = np.empty(NS + 1); ms[0] = m
    for i in range(NS):
        r0, rm, r2 = rg[2 * i], rg[2 * i + 1], rg[2 * i + 2]; b0, bm, b2 = mbg[2 * i], mbg[2 * i + 1], mbg[2 * i + 2]
        a1 = -G * (m + b0) / (s2 * r0); c1 = k * r0 ** 3 * np.exp(h)
        hh = h + .5 * dx * a1; mm = m + .5 * dx * c1
        a2 = -G * (mm + bm) / (s2 * rm); c2 = k * rm ** 3 * np.exp(hh)
        hh = h + .5 * dx * a2; mm = m + .5 * dx * c2
        a3 = -G * (mm + bm) / (s2 * rm); c3 = k * rm ** 3 * np.exp(hh)
        hh = h + dx * a3; mm = m + dx * c3
        a4 = -G * (mm + b2) / (s2 * r2); c4 = k * r2 ** 3 * np.exp(hh)
        h += dx / 6 * (a1 + 2 * a2 + 2 * a3 + a4); m += dx / 6 * (c1 + 2 * c2 + 2 * c3 + c4); ms[i + 1] = m
        if not (np.isfinite(h) and np.isfinite(m)) or h < -700: return np.nan, np.nan, ms
    return h, m, ms
def solve(mdm200, r200, c, mb, xs, t_age, rmin_frac=1e-3):
    """mb(r): baryon mass inside r (vectorized). xs(vmean): cross-section in cm^2/g. Returns dict with r1, rho0, sigma0 and
    mdm(r) (dark mass inside r); the NFW halo alone when r1 < rmin_frac * r_s or no solution (flag)."""
    h = NFW(mdm200, r200, c)
    def rate_log(lr1, lsig0):
        sig0 = np.exp(lsig0); vm = 4 / np.sqrt(np.pi) * sig0; r1 = np.exp(lr1)
        return np.log(h.rho(r1) * xs(vm) * CM2G * vm * KMS_KPC_GYR * t_age)
    def resid(x):
        lrho0, lsig0, lr1 = x; r1 = np.exp(lr1)
        if not (np.isfinite(r1) and 1e-6 < r1 < 100 * r200): return [50., 50., 50.]
        hh, m, _ = _integ(lrho0, lsig0, r1, mb(r1 * np.exp(T)))
        if not np.isfinite(hh) or m <= 0: return [50., 50., 50.]
        return [lrho0 + hh - np.log(h.rho(r1)), np.log(m / h.mass(r1)), rate_log(lr1, lsig0)]
    sig = np.sqrt(G * (h.mass(h.rs) + mb(np.array(h.rs))) / h.rs / 2)
    # r1 for this dispersion (NFW density is monotonic)
    lo, hi = 1e-6 * h.rs, 50 * r200
    if rate_log(np.log(lo), np.log(sig)) < 0: return dict(r1=0., flag='nfw', mdm=h.mass)
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if rate_log(np.log(mid), np.log(sig)) > 0: lo = mid
        else: hi = mid
    r1 = np.sqrt(lo * hi)
    if r1 < rmin_frac * h.rs: return dict(r1=float(r1), flag='nfw', mdm=h.mass)
    best = None
    for fac in (2., 5., 1.2, 20., 100.):
        sol = root(resid, [np.log(fac * h.rho(r1)), np.log(sig), np.log(r1)], method='hybr', options={'xtol': 1e-10})
        if sol.success and np.max(np.abs(sol.fun)) < 1e-7 and (best is None or sol.x[0] < best[0]): best = sol.x
        if best is not None and fac == 1.2: break
    if best is None: return dict(r1=float(r1), flag='fail', mdm=h.mass)
    lrho0, lsig0, lr1 = best; r1 = float(np.exp(lr1)); rho0 = float(np.exp(lrho0))
    hh, m, ms = _integ(lrho0, lsig0, r1, mb(r1 * np.exp(T)))
    lr_grid = np.log(r1) + T[::2]; lm_grid = np.log(ms)
    def mdm(r):
        r = np.asarray(r, float); out = np.array(h.mass(r), float); inner = r < r1
        if np.any(inner):
            lr = np.log(r[inner]); v = np.interp(lr, lr_grid, lm_grid)
            v = np.where(lr < lr_grid[0], np.log(4 / 3 * np.pi * rho0) + 3 * lr, v); out[inner] = np.exp(v)
        return out
    return dict(r1=r1, rho0=rho0, sigma0=float(np.exp(lsig0)), flag='ok', mdm=mdm, nfw=h)

if __name__ == '__main__':
    import time
    def hern(M, a): return lambda r: M * r ** 2 / (r + a) ** 2
    def plum(M, a): return lambda r: M * r ** 3 / (r * r + a * a) ** 1.5
    def r200(m, hz=0.0671): return (G * m / (100 * hz ** 2)) ** (1 / 3)   # H in km/s/kpc
    fb = 0.1543
    for name, m, c, mb, ta in [('MW-like', 1.5e12, 8., hern(6e10, 3.), 10.), ('dSph', 3e9, 15., plum(1e6, .3), 10.),
                               ('cluster', 1e15, 4., hern(1.5e12, 20.), 5.)]:
        for s0 in (0.0001, 0.5, 2., 10.):
            t0 = time.time(); out = solve((1 - fb) * m, r200(m), c, mb, lambda v: xs_vdep(v, s0), ta)
            rr = np.array([.3, 1., 3., 10., 30., 100.])
            print(f"{name:8s} s0={s0:6.4f} r1={out['r1']:8.3f} flag={out['flag']:4s} sigma0={out.get('sigma0', 0):7.1f} "
                  f"M_dm/M_nfw at {rr.tolist()}: {np.round(out['mdm'](rr) / NFW((1 - fb) * m, r200(m), c).mass(rr), 3).tolist()}  {time.time() - t0:.2f}s", flush=True)
