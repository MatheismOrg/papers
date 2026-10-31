"""A19 supplement: learning curve. Hold out a random 5 to 90 percent of the 149 SPARC disks, refit every
rule on the rest (each galaxy weighted equally), and score the held-out disks, the original 44 dwarf
spheroidals and the 27 blind-test dwarfs. 50 random draws per training fraction, fixed seeds.
Usage: python3 A19_Learning_Curve.py SPARC_Lelli2016c.mrt Rotmod_LTG dwarfs_benchmark.csv A19_Blind_Dwarfs.csv out_prefix"""
import sys, os, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from A19_Disk_Fit_SPARC import table1, load, fit_halo, fit_a0, resid, g_nfw, rar, FB, G, ACC
rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
rmsg = lambda lst: float(np.sqrt(np.mean([np.mean(v ** 2) for v in lst])))
def dwarf_set(path, kind):
    r = list(csv.DictReader(open(path))); f = lambda k: np.array([float(x[k]) for x in r])
    if kind == 'mcconnachie':
        R = f('radius_kpc'); gb = f('gbar_m_s2'); lo = np.log10(f('gobs_m_s2')); Re = f('projected_Re_circularized_kpc'); Ms = f('stellar_mass_Msun')
    else:
        Re = f('R_e_circ_kpc'); Ms = 2 * f('L_V'); R = 4 / 3 * Re; gb = G * (Ms / 2) / R ** 2 * ACC; lo = f('log_g_obs')
    return dict(R=R, gb=gb, lo=lo, Re=Re, Ms=Ms, n=len(r))
def score_dwarfs(d, star, bar, a0):
    H = lambda p: d['gb'] + g_nfw(d['R'], d['Ms'] / FB + 10 ** p[0] * (d['Ms'] / 1e10) ** p[1] * (d['Re'] / 3) ** p[2])
    return {'stellar': rms(d['lo'] - np.log10(H(star))), 'baryonic': rms(d['lo'] - np.log10(H(bar))), 'rar': rms(d['lo'] - np.log10(rar(d['gb'], a0)))}
def main(t1, rot, mcc, blind, out):
    gal = load(table1(t1), rot); n = len(gal)
    D = {'original_44': dwarf_set(mcc, 'mcconnachie'), 'blind_27': dwarf_set(blind, 'lvdb')}
    fracs = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]; draws = 50
    res = {'galaxies': n, 'draws_per_fraction': draws, 'fractions': {}}
    for f in fracs:
        ntr = int(round(f * n)); rows = []
        for d in range(draws):
            rng = np.random.default_rng(1000 * int(f * 100) + d); idx = rng.permutation(n)
            tr = [gal[i] for i in idx[:ntr]]; te = [gal[i] for i in idx[ntr:]]
            star = fit_halo(tr, 'star'); bar = fit_halo(tr, 'baryon'); a0 = fit_a0(tr)
            row = {'coef_star': list(map(float, star)), 'a0': float(a0),
                   'disk_heldout': {'stellar': rmsg(resid('star', te, star)), 'baryonic': rmsg(resid('baryon', te, bar)), 'rar': rmsg(resid('rar', te, a0))}}
            for k, ds in D.items(): row[k] = score_dwarfs(ds, star, bar, a0)
            rows.append(row)
        q = lambda vals: [round(float(np.percentile(vals, p)), 4) for p in (16, 50, 84)]
        summ = {'training_galaxies': ntr, 'heldout_galaxies': n - ntr}
        for key in ['disk_heldout'] + list(D):
            summ[key] = {m: q([r[key][m] for r in rows]) for m in ('stellar', 'baryonic', 'rar')}
            if key != 'disk_heldout':
                summ[key]['draws_stellar_beats_rar'] = int(sum(r[key]['stellar'] < r[key]['rar'] for r in rows))
                summ[key]['draws_stellar_beats_baryonic'] = int(sum(r[key]['stellar'] < r[key]['baryonic'] for r in rows))
        summ['coef_star_median'] = [round(float(np.median([r['coef_star'][j] for r in rows])), 3) for j in range(3)]
        summ['coef_star_spread_16_84'] = [round(float(np.percentile([r['coef_star'][j] for r in rows], 84) - np.percentile([r['coef_star'][j] for r in rows], 16)), 3) for j in range(3)]
        res['fractions'][f'{f:.2f}'] = summ
        print(f, json.dumps({k: summ[k] for k in ['disk_heldout', 'original_44', 'blind_27', 'coef_star_median', 'coef_star_spread_16_84']}), flush=True)
    json.dump(res, open(out + '_Results.json', 'w'), indent=1)
if __name__ == '__main__': main(*sys.argv[1:6])
