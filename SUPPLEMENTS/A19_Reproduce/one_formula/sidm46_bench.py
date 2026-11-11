"""SA73 T46 exploration (all data, no protocol yet): SIDM inner halos (sidm46.py) on the benchmark, for two halo-mass
models: 'std' = standard halos (Moster expected halo with each sample's selection, Dutton-Maccio c, no subhalo factor)
and 'formula' = the final formula's halos (T45, all-disk A, B, C, c with the Moline factor). Cross-section
sigma/m = s0 / (1 + (<v>/344 km/s)^2) cm^2/g (s0 = 2 gives Kaplinghat, Tulin and Yu's 1.9 for galaxies and 0.1 for
clusters); t_age 10 Gyr for galaxies, 5 Gyr for cluster-scale halos (MaNGA BCGs, CLASH), as in that paper.
Baryons: spherical, from the measured g_bar where a system has several radii (disks, CLASH, early types), otherwise
a Plummer (dwarfs) or Hernquist (BCGs) profile through the measured point. Writes res46_sidm.json (per-host mse)."""
import sys, json, time
import numpy as np
from multiprocessing import Pool
import harness_T45 as M
import sidm46 as SD
H = M.H; S = M.S; K = H.K; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
S0S = [0.0, 0.3, 1.0, 2.0, 5.0]
gs = H.disks(); cl = H.cluster_rows()
P = {'disk': H.arr(H.disk_rows(gs)), 'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()),
     'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']), 'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']),
     'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
CLASH = H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])
def mb_from_points(r, mb, a_hern=None, inner='r2'):
    o = np.argsort(r); r, mb = r[o], np.maximum.accumulate(mb[o]); lr, lm = np.log(r), np.log(mb)
    def f(x):
        x = np.asarray(x, float); out = np.exp(np.interp(np.log(np.maximum(x, 1e-12)), lr, lm))
        lo = x < r[0]
        if np.any(lo):
            if a_hern is None: out = np.where(lo, mb[0] * (x / r[0]) ** 2, out)
            else: out = np.where(lo, mb[0] * (x / (x + a_hern)) ** 2 / (r[0] / (r[0] + a_hern)) ** 2, out)
        return out
    return f
def host_jobs():
    jobs = []
    for pop in ['disk', 'dwarfs', 'ETG', 'MaNGA', 'CLASH']:
        a = CLASH if pop == 'CLASH' else P[pop]
        for n in np.unique(a['name']):
            k = a['name'] == n; jobs.append((pop, n, {key: np.asarray(v)[k] for key, v in a.items()}))
    return jobs
def baryons(pop, a):
    r = a['r']; mb = a['gb'] / H.ACC * r ** 2 / H.G
    if pop == 'dwarfs':
        ap = a['re'][0]; norm = mb[0] / (r[0] ** 3 / (r[0] ** 2 + ap ** 2) ** 1.5)
        return lambda x: norm * np.asarray(x, float) ** 3 / (np.asarray(x, float) ** 2 + ap ** 2) ** 1.5
    if pop == 'MaNGA':
        ah = a['re'][0] / 1.8153; norm = mb[0] / (r[0] / (r[0] + ah)) ** 2
        return lambda x: norm * (np.asarray(x, float) / (np.asarray(x, float) + ah)) ** 2
    if pop in ('ETG', 'CLASH'): return mb_from_points(r, mb, a_hern=a['re'][0] / 1.8153)
    return mb_from_points(r, mb)
def halos(pop, a):
    thr = S.thresholds(a)[:1] if np.ndim(S.thresholds(a)) else np.array([S.thresholds(a)])
    ms, mp, re, z = a['ms'][:1], a['mp'][:1], a['re'][:1], a['z'][:1]
    m_std = S.m_exp(ms, z, thr); m_std = (m_std ** 8 + (10. ** S.LFLOOR) ** 8) ** (1 / 8)
    m_f = M.P.G.mass(TH, ms, mp, re, z, thr)
    fac = a['fac'][0] if 'fac' in a else 1.
    return {'std': (float(m_std[0]), float(K.c_dm14(m_std, z)[0])), 'formula': (float(m_f[0]), float(K.c_dm14(m_f, z)[0] * fac))}
def run(job):
    pop, n, a = job; mb = baryons(pop, a); z = a['z'][0]; hz = H.hz(np.array([z]))[0]
    t_age = 5. if pop in ('MaNGA', 'CLASH') else 10.
    out = {}
    for model, (m, c) in halos(pop, a).items():
        rv = (H.G * m / (100 * hz ** 2)) ** (1 / 3)
        for s0 in S0S:
            if s0 == 0: sol = dict(r1=0., flag='nfw', mdm=SD.NFW((1 - H.FB) * m, rv, c).mass)
            else: sol = SD.solve((1 - H.FB) * m, rv, c, mb, lambda v, s0=s0: SD.xs_vdep(v, s0), t_age)
            g = a['gb'] + H.G * sol['mdm'](a['r']) / a['r'] ** 2 * H.ACC
            er = np.log10(a['go'] / g)
            out[f'{model}|{s0}'] = dict(er=er.tolist(), kind=list(a['kind']), r1=sol['r1'], flag=sol['flag'])
    return pop, n, out
if __name__ == '__main__':
    t0 = time.time(); jobs = host_jobs()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 2) as pool: res = pool.map(run, jobs, chunksize=4)
    json.dump([[p, n, o] for p, n, o in res], open('res46_sidm.json', 'w'))
    # all-data scores
    for model in ('std', 'formula'):
        for s0 in S0S:
            key = f'{model}|{s0}'; acc = {}
            for p, n, o in res:
                er = np.array(o[key]['er']); kinds = o[key]['kind']
                if p == 'CLASH':
                    for kk in ('CLASH_centre', 'CLASH_outer'):
                        sel = np.array([x == kk for x in kinds]); acc.setdefault(kk, []).append(np.mean(er[sel] ** 2))
                else: acc.setdefault(p, []).append(np.mean(er ** 2))
            fl = [o[key]['flag'] for p, n, o in res]
            print(key, {k: round(float(np.sqrt(np.mean(v))), 4) for k, v in acc.items()}, 'fails', fl.count('fail'), flush=True)
    print('sec', round(time.time() - t0))
