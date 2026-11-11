"""T38 (PROTOCOL_T38.txt): the frozen T37 formula against the unfitted checks: Mandelbaum et al. 2016 central halos,
the KiDS early/late split, and the Milky Way halo. Nothing fitted."""
import sys, csv, json
import numpy as np
import harness_T37 as S
H = S.H; K = H.K; B = S.B; F = S.F
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
import kids_split as KS
B.FROZEN['abc'] = abc = np.array([11.7481, .8178, -.7119]); th = np.array([11.1768])
gs = H.disks(); (p1, p0), (q1, q0) = H.pop_coeff(gs)
t31 = {(d['color'], d['logMs']): d for d in json.load(open('res31_mandelbaum_scatter.json'))['bins']}
out = []
for r in csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')):
    lm, z, obs, ep, em = float(r['logMs_eff_VAGC']), float(r['z_eff']), float(r['logM200m_hinv']), float(r['err_plus']), float(r['err_minus'])
    red = r['color'] == 'red'; ms = 10 ** lm
    re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
    m = float(S.mass(th, np.array([ms]), np.array([mp]), np.array([re]), np.array([z]), 0.)[0])
    pred = L8.to_m200m(m, K.c_dm14(m, z), z); dev = pred - obs; tol = max(.2, 2 * (ep if dev > 0 else em))
    d31 = t31[(r['color'], lm)]
    out.append(dict(color=r['color'], logMs=lm, obs=obs, T26=d31['formula_T30'], T31=d31['formula_scatter'], T37=round(float(pred), 3),
                    dev=round(float(dev), 3), tol=round(tol, 3), ok=bool(abs(dev) <= tol)))
red = [d for d in out if d['color'] == 'red']; blue = [d for d in out if d['color'] == 'blue' and d['logMs'] < 11.4]
rms = lambda L, k: round(float(np.sqrt(np.mean([(d[k] - d['obs']) ** 2 for d in L]))), 3)
mand = dict(red_rms={k: rms(red, k) for k in ('T26', 'T31', 'T37')}, blue_rms={k: rms(blue, k) for k in ('T26', 'T31', 'T37')},
            bins_ok=sum(d['ok'] for d in out), bins=len(out), failed=[(d['color'], d['logMs'], d['dev']) for d in out if not d['ok']])
for d in out: print('%-5s %5.2f obs %5.2f | T26 %6.2f T31 %6.2f T37 %6.2f dev %+5.2f tol %.2f %s' % (d['color'], d['logMs'], d['obs'], d['T26'], d['T31'], d['T37'], d['dev'], d['tol'], d['ok']))
print(json.dumps(mand))
SS = KS.setup(gs); kid = dict(T26=KS.check('MOGF0IE', np.array([11.1868]), SS), T37=KS.check('MOGF0IES', th, SS)); print(json.dumps(kid))
mstar, mp = 5.43e10, None
a_ = np.linspace(1e-3, 60, 6000); discs = {'HI': (53.1, 7., 4., .085), 'H2': (2180., 1.5, 12., .045)}
mgas = sum(2 * np.pi * np.sum(a_ * p[0] * 1e6 * np.exp(-p[2] / a_ - a_ / p[1])) * (a_[1] - a_[0]) for p in discs.values()); mp = mstar + mgas
args = (np.array([mstar]), np.array([mp]), np.array([5.75]))
mw = dict(q=float(F.Z.switch('MOG', th, *args)[0]), M_halo_T37=float(S.mass(th, *args, np.array([0.]), 0.)[0]),
          M_halo_T26=float(F.mass('MOGF0IE', np.array([11.1868]), *args, np.array([0.]))[0]))
print(json.dumps(mw))
ok = (mand['red_rms']['T37'] <= mand['red_rms']['T26'] + .02 and mand['blue_rms']['T37'] <= mand['blue_rms']['T26'] + .02
      and kid['T37']['S_model'] >= kid['T26']['S_model'] - .01 and abs(mw['M_halo_T37'] / mw['M_halo_T26'] - 1) < .01)
json.dump(dict(mandelbaum=dict(summary=mand, bins=out), kids_split=kid, milky_way=mw, passes=bool(ok)), open('res38_frozen_T37_checks.json', 'w'), indent=1)
print('T38 passes:', ok)
