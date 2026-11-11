"""SA73 T64: self-interacting dark matter on every rotation-curve point. The isothermal Jeans model of T46 (dmfit/sidm46.py:
Kaplinghat, Tulin and Yu 2016; Ren et al. 2019; baryons included, spherical) in the halos of the formula without its core
(T48; A, B, C at the all-disk values 11.7481, 0.8178, -0.7119, M_t 11.339, as T46 did, so the held-out disk scores are
slightly optimistic for every cross-section alike). Cross-section sigma/m = s0 / (1 + (v/344 km/s)^2) cm^2/g, the
published fall with speed (s0 = 2 gives about 2 in dwarfs and 0.1 in clusters); t_age 10 Gyr for galaxies, 5 Gyr for
cluster-scale halos. Disks: all 175 SPARC galaxies of harness_T61.ALL, every point (with its fractional velocity error);
dwarfs, early types, MaNGA and CLASH: the benchmark's points. Writes res64_sidm_points.json (per host and s0, the errors).
Usage: PYTHONPATH=../dmfit:../t59 python3 sidm64.py <processes> [<max hosts per population, for timing>] [<t59 folder, to add Iorio's eight>]"""
import sys, json, time
import numpy as np
from multiprocessing import Pool
import harness_T64 as M
import sidm46 as SD
H = M.H; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
S0S = [0.0, 0.03, 0.1, 0.3, 1.0, 2.0, 5.0, 10.0]
cl = H.cluster_rows()
POPS = {'disk': M.ALL, 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
def iorio_rows(t59):
    import os, model59 as L
    gal = json.load(open(os.path.join(t59, 'lt59_inputs.json'))); KL, KR = 1.172, 1.424; out = []
    for nm, g in gal.items():
        if g['sparc'] is not None or g['f36'] is None or g['inc'] < 30: continue
        R = np.array(g['R']); Vc = np.array(g['Vc']); eV = np.array(g['eVc'])
        Ms = 0.5 * L.l36(g['f36'][0], g['D']) / KL; Re = 1.678 * g['read']['Rstar'] * g['D'] / g['read']['D'] / KR
        vg2, mhi, mcut = L.gas(R, g['Sd'], 'exp'); vb2 = vg2 + L.stars_v2(R, Ms, Re / 1.678); k = (R > 0) & (Vc > 0) & (vb2 > 0); Mp = Ms + 1.33 * mhi
        out += [dict(name=nm, kind='disk', r=x, gb=y, go=w, ms=Ms, mp=Mp, re=Re, z=0.) for x, y, w in zip(R[k], vb2[k] / R[k] * H.ACC, Vc[k] ** 2 / R[k] * H.ACC)]
    return H.arr(out)
if len(sys.argv) > 3 or (len(sys.argv) == 3 and not sys.argv[2].isdigit()): POPS['iorio'] = iorio_rows(sys.argv[-1])
def mb_from_points(r, mb, a_hern=None):
    o = np.argsort(r); r, mb = r[o], np.maximum.accumulate(mb[o]); lr, lm = np.log(r), np.log(mb)
    def f(x):
        x = np.asarray(x, float); out = np.exp(np.interp(np.log(np.maximum(x, 1e-12)), lr, lm)); lo = x < r[0]
        if np.any(lo):
            if a_hern is None: out = np.where(lo, mb[0] * (x / r[0]) ** 2, out)
            else: out = np.where(lo, mb[0] * (x / (x + a_hern)) ** 2 / (r[0] / (r[0] + a_hern)) ** 2, out)
        return out
    return f
def baryons(pop, a):
    r = a['r']; mb = a['gb'] / H.ACC * r ** 2 / H.G
    if pop == 'dwarfs':
        ap = a['re'][0]; norm = mb[0] / (r[0] ** 3 / (r[0] ** 2 + ap ** 2) ** 1.5); return lambda x: norm * np.asarray(x, float) ** 3 / (np.asarray(x, float) ** 2 + ap ** 2) ** 1.5
    if pop == 'MaNGA':
        ah = a['re'][0] / 1.8153; norm = mb[0] / (r[0] / (r[0] + ah)) ** 2; return lambda x: norm * (np.asarray(x, float) / (np.asarray(x, float) + ah)) ** 2
    if pop in ('ETG', 'CLASH'): return mb_from_points(r, mb, a_hern=a['re'][0] / 1.8153)
    return mb_from_points(r, mb)
def jobs(nmax=None):
    out = []
    for pop, a in POPS.items():
        for n in list(np.unique(a['name']))[:nmax]:
            k = a['name'] == n; out.append((pop, str(n), {key: np.asarray(v)[k] for key, v in a.items()}))
    return out
def run(job):
    pop, n, a = job; mb = baryons(pop, a); m, c, rv = (float(np.atleast_1d(v)[0]) for v in M.halo(TH, a)); t_age = 5. if pop in ('MaNGA', 'CLASH') else 10.; out = {}
    for s0 in S0S:
        if s0 == 0: sol = dict(r1=0., flag='nfw', mdm=SD.NFW((1 - H.FB) * m, rv, c).mass)
        else: sol = SD.solve((1 - H.FB) * m, rv, c, mb, lambda v, s0=s0: SD.xs_vdep(v, s0), t_age)
        g = a['gb'] + H.G * sol['mdm'](a['r']) / a['r'] ** 2 * H.ACC
        out[str(s0)] = dict(er=np.log10(a['go'] / g).tolist(), r1=sol['r1'], flag=sol['flag'], sigma0=sol.get('sigma0'))
    return pop, n, dict(by_s0=out, kind=[str(x) for x in a['kind']], fe=a['fe'].tolist() if 'fe' in a else None, r=a['r'].tolist(), logm=float(np.log10(m)), c=c)
if __name__ == '__main__':
    t0 = time.time(); timing = len(sys.argv) > 2 and sys.argv[2].isdigit(); js = jobs(int(sys.argv[2]) if timing else None)
    with Pool(int(sys.argv[1])) as pool: res = pool.map(run, js, chunksize=2)
    if timing:
        print('timing run:', len(js), 'hosts in', round(time.time() - t0, 1), 's; flags', {f: sum(o['by_s0'][s]['flag'] == f for p, n, o in res for s in o['by_s0']) for f in ('ok', 'nfw', 'fail')})
    else:
        json.dump([[p, n, o] for p, n, o in res], open('res64_sidm_points.json', 'w'))
        print('hosts', len(res), 'sec', round(time.time() - t0), 'flags', {f: sum(o['by_s0'][s]['flag'] == f for p, n, o in res for s in o['by_s0']) for f in ('ok', 'nfw', 'fail')}, flush=True)
