"""T41 step 2: Mandelbaum red-central rms (as T38) and HeCS numbers (as T40, selected branch) for the mass function in A19_HMF,
at the all-data threshold given on the command line. Usage: A19_HMF=<model> python3 check41.py <log10 M_t> <hecs_t2.csv> <hecs_t3.csv>"""
import sys, csv, json
import numpy as np
import harness_T41 as W
S = __import__('harness_T37'); H = S.H; K = H.K; B = S.B
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
B.FROZEN['abc'] = ABC = np.array([11.7481, .8178, -.7119]); th = np.array([float(sys.argv[1])])
gs = H.disks(); (p1, p0), (q1, q0) = H.pop_coeff(gs); dev = []
for r in csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')):
    if r['color'] != 'red': continue
    lm, z, obs = float(r['logMs_eff_VAGC']), float(r['z_eff']), float(r['logM200m_hinv']); ms = 10 ** lm
    re = 10 ** (p0 + p1 * lm) * .5; m = float(S.mass(th, np.array([ms]), np.array([ms]), np.array([re]), np.array([z]), 0.)[0])
    dev.append(L8.to_m200m(m, K.c_dm14(m, z), z) - obs)
t2 = {r['ID']: r for r in csv.DictReader(open(sys.argv[2]))}; t3 = {r['ID']: r for r in csv.DictReader(open(sys.argv[3]))}
man = {(r['ms'], r['re']) for r in H.cluster_rows() if r['kind'] == 'MaNGA_BCG'}; P = np.polyfit(np.log10([a for a, _ in man]), np.log10([b for _, b in man]), 1)
rows = [(float(r['Mass']), float(t2[k]['zspec']), float(t2[k]['M200']) * 1e14) for k, r in t3.items() if k in t2 and float(r['Mass']) > 0 and float(t2[k]['M200']) > 0]
lm = np.array([x[0] for x in rows]); z = np.array([x[1] for x in rows]); obs = np.log10([x[2] for x in rows]); ms = 10 ** lm
d = np.log10(S.mass(th, ms, ms, 10 ** np.polyval(P, lm), z, np.full(len(ms), 1e14))) - obs
print(json.dumps(dict(hmf=W.MOD, logMt=th[0], mandelbaum_red_rms=round(float(np.sqrt(np.mean(np.square(dev)))), 3), mandelbaum_red_dev=[round(float(x), 2) for x in dev],
                      hecs_N=len(d), hecs_mean=round(float(d.mean()), 3), hecs_median=round(float(np.median(d)), 3), hecs_rms=round(float(np.sqrt(np.mean(d ** 2))), 3))))
