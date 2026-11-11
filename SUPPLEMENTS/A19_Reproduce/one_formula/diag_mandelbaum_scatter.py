"""T31 (PROTOCOL_T31.txt): T30 with the Moster term inverted with its scatter and the Tinker et al. (2008) halo mass
function, as in the paper's sixth protocol, for galaxies selected by stellar mass."""
import sys, csv, json
import numpy as np
import harness_T26 as M
F = M.F; H = M.H; K = H.K; B = __import__('harness_T15b')
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8, A19_Clusters_v5 as V5
from colossus.lss import mass_function
B.FROZEN['abc'] = abc = np.array([11.7481, .8178, -.7119]); th = np.array([11.1868])
(p1, p0), (q1, q0) = H.pop_coeff(H.disks()); HH = 0.671; GRID = np.arange(9.0, 16.5 + 1e-9, 0.005); SIG = np.hypot(.15, .10)
t30 = {(d['color'], d['logMs']): d for d in json.load(open('res30_mandelbaum.json'))['bins']}
rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
out = []
for r in rows:
    lm, z, obs, ep, em = float(r['logMs_eff_VAGC']), float(r['z_eff']), float(r['logM200m_hinv']), float(r['err_plus']), float(r['err_minus'])
    red = r['color'] == 'red'; ms = 10 ** lm
    prior = mass_function.massFunction(10 ** GRID * HH, z, mdef='200c', model='tinker08', q_out='dndlnM')
    like = np.exp(-.5 * ((lm - np.log10(V5.moster_mstar(10 ** GRID, z))) / SIG) ** 2)
    w = prior * like; w /= w.sum(); k = w > 1e-7 * w.max()
    mscat = float(np.sum(w[k] * 10 ** GRID[k]) / w[k].sum())                      # posterior mean M200c, suns
    m200m = [10 ** L8.to_m200m(10 ** g, K.c_dm14(10 ** g, z), z) for g in GRID[k]]
    moster_scatter = float(np.log10(np.sum(w[k] * np.array(m200m)) / w[k].sum()))
    re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
    q = float(F.Z.switch('MOG', th, np.array([ms]), np.array([mp]), np.array([re]))[0])
    mref = mp / H.FB + 10 ** abc[0] * (ms / 1e10) ** abc[1] * (re / 3) ** abc[2] + q * mscat
    pred = L8.to_m200m(mref, K.c_dm14(mref, z), z); dev = pred - obs; tol = max(.2, 2 * (ep if dev > 0 else em))
    out.append(dict(color=r['color'], logMs=lm, obs=obs, formula_T30=t30[(r['color'], lm)]['formula'], moster_mean=t30[(r['color'], lm)]['moster'],
                    moster_scatter=round(moster_scatter, 3), formula_scatter=round(pred, 3), dev=round(dev, 3), tol=round(tol, 3), ok=bool(abs(dev) <= tol), q=round(q, 3)))
print('color logM*  obs   | T30 formula  Moster mean | Moster scatter  formula scatter  dev    ok')
for d in out: print('%-5s %5.2f %5.2f | %7.2f %11.2f | %10.2f %13.2f %+7.2f  %s' % (d['color'], d['logMs'], d['obs'], d['formula_T30'], d['moster_mean'], d['moster_scatter'], d['formula_scatter'], d['dev'], d['ok']))
top = [d for d in out if d['color'] == 'red' and d['logMs'] > 11.2]
red = [d for d in out if d['color'] == 'red']; blue = [d for d in out if d['color'] == 'blue' and d['logMs'] < 11.4]
rms = lambda L, k: round(float(np.sqrt(np.mean([(d[k] - d['obs']) ** 2 for d in L]))), 3)
res = dict(selection_reading_supported=all(d['ok'] for d in top), red_rms=dict(T30=rms(red, 'formula_T30'), scatter=rms(red, 'formula_scatter')),
           blue_rms=dict(T30=rms(blue, 'formula_T30'), scatter=rms(blue, 'formula_scatter')), top_red_dev=[d['dev'] for d in top])
print(json.dumps(res)); json.dump(dict(bins=out, summary=res), open('res31_mandelbaum_scatter.json', 'w'), indent=1)
