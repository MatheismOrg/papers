"""SA73 T57 all-data look (no protocol): cold gas making halos bigger (mu) and looser (s), each alone and together, fitted on
all disks, against the final formula. Usage: PYTHONPATH=<dmfit> python3 describe57.py"""
import sys, csv, json, numpy as np
import harness_T57 as M7
import harness_T48 as M8
P = M7.P; H = M7.H; B = M7.B; S = M7.S; K = H.K
sys.path.insert(0, str(H.PKG / 'code')); import A19_Lensing_v8 as L8
gs = H.disks(); cl = H.cluster_rows(); d = H.arr(H.disk_rows(gs)); LENS = {False: H.Lens(gs, False), True: H.Lens(gs, True)}
ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
       'CLASH_centre': H.arr([r for r in cl if r['kind'] == 'CLASH_centre']), 'CLASH_outer': H.arr([r for r in cl if r['kind'] == 'CLASH_outer'])}
(p1, p0), (q1, q0) = H.pop_coeff(gs); rows = list(csv.DictReader(open(H.PKG / 'data/Lensing/mandelbaum2016_central_halos.csv')))
un, gi = np.unique(d['name'], return_inverse=True); ct = np.bincount(gi); lx = np.log10(d['r'] / d['re']); gmass = np.array([d['ms'][gi == i][0] for i in range(len(un))])
def score(model, th, a): return round(float(np.sqrt(np.mean(list(H.losses(model, th, a).values())))), 4)
def look(lab, free):
    M7.U['free'] = free; f = H.fit(M7.MODEL, .1, d, ext, LENS[False], [0, 1, 2]); th = np.array(f['th'])
    r = dict(abc=np.round(f['abc'], 3).tolist(), mu=round(f['mu'], 2), s=round(f['s'], 3), disk=score(M7.MODEL, th, d), **{k: score(M7.MODEL, th, a) for k, a in ext.items()})
    for drop in (False, True):
        c2 = [round(float(np.sum(LENS[drop].white(M7.MODEL, th, [j]) ** 2)), 1) for j in range(4)]; r['lens_' + ('strict' if drop else 'primary')] = round(sum(c2), 1)
    dw = ext['dwarfs']; e = np.log10(dw['go'] / H.gpred(M7.MODEL, th, dw)); X = np.c_[np.ones(len(e)), np.log10(dw['ms']), np.log10(dw['re'])]
    r['dwarf_slope'] = round(float(np.linalg.lstsq(X, e, rcond=None)[0][2]), 2)
    e = np.log10(d['go'] / H.gpred(M7.MODEL, th, d)); gm = np.bincount(gi, e) / ct; w = e - gm[gi]; r['within_disk'] = round(float(np.sqrt(np.mean(np.bincount(gi, w ** 2) / ct))), 4)
    sl = np.array([np.polyfit(lx[gi == i], w[gi == i], 1)[0] if (ct[i] >= 5 and np.ptp(lx[gi == i]) > .3) else np.nan for i in range(len(un))]); ok = np.isfinite(sl)
    r['slope_vs_mass'] = round(float(np.corrcoef(np.log10(gmass)[ok], sl[ok])[0, 1]), 2); r['inner_stack'] = round(float(np.mean(w[lx < np.log10(.6)])), 3)
    c2 = {'red': 0., 'blue': 0.}
    for w_ in rows:
        lm, z, obs, ep, em = float(w_['logMs_eff_VAGC']), float(w_['z_eff']), float(w_['logM200m_hinv']), float(w_['err_plus']), float(w_['err_minus']); red = w_['color'] == 'red'
        ms = 10 ** lm; re = 10 ** (p0 + p1 * lm) * (.5 if red else 1.); mp = ms if red else ms + 1.33 * 10 ** (q0 + q1 * lm)
        m = float(M7.mass(th, np.array([ms]), np.array([mp]), np.array([re]), np.array([z]), 0.)[0]); dev = L8.to_m200m(m, K.c_dm14(m, z), z) - obs; c2[w_['color']] += (dev / (ep if dev > 0 else em)) ** 2
    r['red_chi2'] = round(c2['red'], 1); r['blue_chi2'] = round(c2['blue'], 1)
    print(f'{lab:22s}', json.dumps(r), flush=True); return r
out = {lab: look(lab, free) for lab, free in (('final formula', (False, False)), ('bigger (mu)', (True, False)), ('looser (s)', (False, True)), ('bigger and looser', (True, True)))}
json.dump(out, open('res57_describe.json', 'w'), indent=1)
