"""Independent re-computation of protocol 9's primary numbers (written separately from A19_KiDS_v9.py).
NFW lensing from colossus (dark mass (1 - f_b) M_ref), Simpson integration over each annulus on a linear grid,
coefficients typed from the protocol, data parsed from the raw release files and the compact covariances."""
import json, sys, numpy as np
from scipy.integrate import simpson
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw
cosmology.addCosmology('A19p', {'flat': True, 'H0': 67.1, 'Om0': 0.3175, 'Ob0': 0.049, 'sigma8': 0.8344, 'ns': 0.9624})
cosmology.setCosmology('A19p')
D = sys.argv[1]
h = 0.671; fb = 0.0490 / 0.3175; A, B, C = 11.7481, 0.8178, -0.7119
G_SI, MSUN, MPC = 6.67430e-11, 1.98841e30, 3.0856775814913673e22
edges = 10 ** np.linspace(-15, np.log10(5e-12), 16)
vdw_x, vdw_d = [9.25, 9.75, 10.25, 10.75], [-0.22, -0.33, -0.28, -0.16]
def mref(lm, early):
    ms = 10 ** lm; re = 10 ** (-1.6721 + 0.2245 * lm)
    if early: re *= 10 ** np.interp(lm, vdw_x, vdw_d); mp = ms
    else: mp = ms + 1.33 * 10 ** (4.1156 + 0.5396 * lm)
    return mp / fb + 10 ** A * (ms / 1e10) ** B * (re / 3) ** C, mp
cache = {}
def lens(lm, early):
    key = (lm, early)
    if key in cache: return cache[key]
    m, mp = mref(lm, early); c = 10 ** 0.905 * (m * h / 1e12) ** -0.101
    p = profile_nfw.NFWProfile(M=m * h, c=c, z=0.25, mdef='200c')
    mg = 10 ** lm * (1 + 10 ** (-0.69 * lm + 6.63))
    out = []
    for k in range(15):
        r1 = np.sqrt(G_SI * mg * MSUN / edges[k + 1]) / MPC; r2 = np.sqrt(G_SI * mg * MSUN / edges[k]) / MPC
        R = np.linspace(r1, r2, 801)
        ds = (1 - fb) * p.deltaSigma(R * 1e3 * h) * h / 1e6 + mp / (np.pi * (R * 1e6) ** 2)   # suns / pc^2
        out.append(simpson(ds * R, x=R) / ((r2 ** 2 - r1 ** 2) / 2))
    cache[key] = (np.array(out), mg); return cache[key]
f7 = json.load(open(D + '/brouwer2021_fig7_digitized.json'))['bins']
def stack(items):
    num = np.zeros(15); den = 0; ns = 0; mgs = 0
    for lm, n, early in items:
        e, mg = lens(round(lm, 3), early); num += n * mg * e; den += n * mg; ns += n; mgs += n * mg
    return num / den, np.log10(mgs / ns)
mid = lambda b: round((b['lo'] + b['hi']) / 2, 3)
matched = [(mid(b), b['matched']) for b in f7 if b['matched'] > 0]
late, lmg8 = stack([(x, n, False) for x, n in matched]); early, _ = stack([(x, n, True) for x, n in matched])
lims = [8.5, 10.3, 10.6, 10.8, 11.0]
models = {'blue': late, 'disc': late, 'red': early, 'bulge': early}
for i in range(4):
    models['fig9_bin%d' % (i + 1)], _ = stack([(mid(b), b['blue'], False) for b in f7 if lims[i] <= mid(b) < lims[i + 1] and b['blue'] > 0] +
                                             [(mid(b), b['red'], True) for b in f7 if lims[i] <= mid(b) < lims[i + 1] and b['red'] > 0])
models['full'], lmgf = stack([(mid(b), b['blue'], False) for b in f7 if 8.5 <= mid(b) < 11 and b['blue'] > 0] + [(mid(b), b['red'], True) for b in f7 if 8.5 <= mid(b) < 11 and b['red'] > 0])
def inner(lmg): return edges[:-1] >= G_SI * 10 ** lmg * MSUN / (0.3 * MPC) ** 2
def prof(f):
    d = np.loadtxt(D + '/' + f); return d[:, 1] / d[:, 4]
def cov(f):
    J = json.load(open(D + '/' + f)); nb = J['nb']; no = len(J['obs_min']); n = nb * no; M = np.zeros((n, n)); it = iter(J['v'])
    for a in range(no):
        for b in ([a] if J['blocks_only'] else range(a, no)):
            for i in range(nb):
                for j in range(nb):
                    if a == b and j < i: continue
                    M[a * nb + i, b * nb + j] = M[b * nb + j, a * nb + i] = next(it)
    return M
Cc, Cs, C9, Cf = cov('Fig-8_Colorbins_covmatrix_compact.json'), cov('Fig-8_Sersicbins_covmatrix_compact.json'), cov('Fig-9_Massbins_covmatrix_compact.json'), cov('Fig-4-5-C1_covmatrix_compact.json')
data = {'blue': (prof('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'), Cc[:15, :15]), 'red': (prof('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'), Cc[15:, 15:]),
        'disc': (prof('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'), Cs[:15, :15]), 'bulge': (prof('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt'), Cs[15:, 15:]),
        'full': (prof('Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt'), Cf)}
for i in range(4): data['fig9_bin%d' % (i + 1)] = (prof('Fig-9_RAR-KiDS-isolated_Massbin-%d.txt' % (i + 1)), C9[15 * i:15 * i + 15, 15 * i:15 * i + 15])
mg9 = [10.14, 10.57, 10.78, 10.96]
out = {}
for s, (d, Cv) in data.items():
    k = inner(mg9[int(s[-1]) - 1] if s.startswith('fig9') else (lmgf if s == 'full' else lmg8))
    m = models[s]; W = np.linalg.inv(Cv[np.ix_(k, k)]); a = (m[k] @ W @ d[k]) / (m[k] @ W @ m[k])
    out[s] = dict(log10_a=round(float(np.log10(a)), 4), n_inner=int(k.sum()))
k8 = inner(lmg8)
S = lambda e, l, k: float(np.log10(np.mean(e[k] / l[k])))
out['split_model'] = dict(inner=round(S(early, late, k8), 4), all=round(S(early, late, np.ones(15, bool)), 4))
out['split_data'] = {lab: dict(inner=round(S(data[e][0], data[l][0], k8), 4), all=round(S(data[e][0], data[l][0], np.ones(15, bool)), 4)) for lab, e, l in (('colour', 'red', 'blue'), ('sersic', 'bulge', 'disc'))}
out['log_mgal_matched'] = round(float(lmg8), 4); out['log_mgal_full'] = round(float(lmgf), 4)
print(json.dumps(out, indent=1))
