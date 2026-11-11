"""T30 (PROTOCOL_T30.txt): the frozen T26 formula against Mandelbaum et al. (2016) lensing halo masses of SDSS central
galaxies, red and blue. Reuses the paper's eighth-protocol code (A19_Lensing_v8.py) for the mass-definition conversion
and reads its saved predictions for the paper's rule and the Moster relation, for comparison."""
import sys, csv, json
import numpy as np
import harness_T26 as M
F = M.F; H = M.H; K = H.K; B = __import__('harness_T15b')
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); th = np.array([11.1868])
(p1, p0), (q1, q0) = H.pop_coeff(H.disks())
rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
old = {(r['color'], float(r['logMs'])): r for r in csv.DictReader(open(H.PKG / 'blind_test/A19_Lensing_v8_Bins.csv')) if r['stellar_masses'] == 'VAGC'}
def formula(lm, z, red, size_factor):
    ms = 10 ** lm; re = 10 ** (p0 + p1 * lm) * size_factor; mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
    m = float(F.mass('MOGF0IE', th, np.array([ms]), np.array([mp]), np.array([re]), np.array([z]))[0])
    q = float(F.Z.switch('MOG', th, np.array([ms]), np.array([mp]), np.array([re]))[0])
    return L8.to_m200m(m, K.c_dm14(m, z), z), q
out = []
for r in rows:
    lm, z, obs, ep, em = float(r['logMs_eff_VAGC']), float(r['z_eff']), float(r['logM200m_hinv']), float(r['err_plus']), float(r['err_minus'])
    red = r['color'] == 'red'; pred, q = formula(lm, z, red, .5 if red else 1.)
    o = old[(r['color'], lm)]; d = dict(color=r['color'], logMs=lm, z=z, obs=obs, ep=ep, em=em, formula=round(pred, 3), switch_q=round(q, 3),
             rule=round(float(o['rule_half_size'] if red else o['rule']), 3), moster=round(float(o['moster']), 3))
    if red: d['formula_disk_size'] = round(formula(lm, z, True, 1.)[0], 3)
    dev = pred - obs; d['dev'] = round(dev, 3); d['tol'] = round(max(.2, 2 * (ep if dev > 0 else em)), 3); d['ok'] = abs(dev) <= d['tol']
    out.append(d)
print('color logM*  obs    formula  dev    tol   ok  | rule   moster | q')
for d in out:
    print('%-5s %5.2f %6.2f %7.2f %+6.2f %5.2f %-5s | %6.2f %6.2f | %.2f' % (d['color'], d['logMs'], d['obs'], d['formula'], d['dev'], d['tol'], d['ok'], d['rule'], d['moster'], d['switch_q']))
red = [d for d in out if d['color'] == 'red']; blue = [d for d in out if d['color'] == 'blue' and d['logMs'] < 11.4]
rms = lambda L, k: float(np.sqrt(np.mean([(d[k] - d['obs']) ** 2 for d in L])))
res = dict(R_red_all_within=all(d['ok'] for d in red), B_blue_below_11p4_within=all(d['ok'] for d in blue),
           rms_red=dict(formula=round(rms(red, 'formula'), 3), rule=round(rms(red, 'rule'), 3), moster=round(rms(red, 'moster'), 3)),
           rms_blue_below_11p4=dict(formula=round(rms(blue, 'formula'), 3), rule=round(rms(blue, 'rule'), 3), moster=round(rms(blue, 'moster'), 3)),
           red_disk_size_rms=round(float(np.sqrt(np.mean([(d['formula_disk_size'] - d['obs']) ** 2 for d in red]))), 3))
split = []
for d in red:
    b = [x for x in out if x['color'] == 'blue' and x['logMs'] == d['logMs']] or [x for x in out if x['color'] == 'blue' and abs(x['logMs'] - d['logMs']) < .05]
    if b and 10.7 <= d['logMs'] <= 11.2: split.append(dict(logMs=d['logMs'], measured=round(d['obs'] - b[0]['obs'], 3), formula=round(d['formula'] - b[0]['formula'], 3), rule=round(d['rule'] - b[0]['rule'], 3)))
res['S_split'] = split
print(json.dumps(res, indent=1)); json.dump(dict(bins=out, summary=res), open('res30_mandelbaum.json', 'w'), indent=1)
