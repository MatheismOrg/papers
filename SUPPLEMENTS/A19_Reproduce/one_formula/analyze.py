"""Summaries, pass marks and paired host bootstraps for the T15 cross-scale run (written before the run finished)."""
import sys, json
import numpy as np
SO73 = {  # SO73 published held-out scores on the same benchmark and partitions (primary mask unless noted)
 'SO73 original rule (T54)': dict(disk=.153298, dwarfs=.406391, ETG=.203097, MaNGA=.321920, CLASH_outer=.774209, lens=78.967700, lens_strict=54.093669),
 'SO73 Moster':              dict(disk=.171114, dwarfs=.439850, ETG=.206083, MaNGA=.181830, CLASH_outer=.203556, lens=149.991505, lens_strict=120.443462),
 'SO73 T54 mixture, 14 coef': dict(disk=.151259, dwarfs=.405366, ETG=.198573, MaNGA=.220365, CLASH_outer=.233694, lens=49.976224, lens_strict=30.735803),
 'SO73 T58 bridge, 18 coef': dict(disk=.151401, dwarfs=.379705, ETG=.198549, MaNGA=.205282, CLASH_outer=.194427, lens=50.127820, lens_strict=30.837878),
 'SO73 T62 constrained, 18 coef': dict(disk=.151825, dwarfs=.376910, ETG=.198755, MaNGA=.205270, CLASH_outer=.191533, lens=49.091270, lens_strict=30.173139),
}
TARGET = dict(disk=.151401, dwarfs=.376910, ETG=.198549, MaNGA=.205270, CLASH_outer=.191533, lens=49.091270)
NCOEF = {'ORIG': 3, 'KNEE': 5, 'KNEES': 6, 'KNEEC': 7, 'KNEESC': 8, 'KNEEE': 7}
MAIN = ['disk', 'dwarfs', 'ETG', 'MaNGA', 'CLASH_outer']
TOL = dict(disk=.01, dwarfs=.02, ETG=.05, MaNGA=.05, CLASH_outer=.05, lens=.05)

def table(res):
    sm = res['summary']; rows = []
    for cfg, d in sm.items():
        r = {k: d['scores'][k] for k in MAIN}; r['lens'] = d['lensing_chi2']; r['CLASH_centre'] = d['scores']['CLASH_centre']
        r['bound'] = d['bound_hit_partitions']; r['failed'] = d['failed_attempts']; rows.append((cfg, r))
    return dict(rows)
def passes(t, base):
    out = {}
    for cfg, r in t.items():
        better = sum(r[k] < base[k] for k in MAIN + ['lens'])
        worse_tol = [k for k in MAIN + ['lens'] if r[k] > base[k] * (1 + TOL[k])]
        match = all(r[k] <= TARGET[k] * 1.03 for k in MAIN + ['lens'])
        beats = sum(r[k] < TARGET[k] for k in MAIN + ['lens'])
        out[cfg] = dict(better_than_ORIG=better, worse_beyond_tolerance=worse_tol, improvement=better >= 4 and not worse_tol,
                        matches_SO73=match, beats_SO73_count=beats, beats_SO73=beats >= 4 and NCOEF[cfg.split('@')[0]] <= 8)
    return out
def host_losses(res, cfg, pop):
    pooled = {}
    for r in res['records']:
        if r['config'] != cfg: continue
        for host, mse in r['losses'][pop].items(): pooled.setdefault(host, []).append(mse)
    return {h: float(np.mean(v)) for h, v in pooled.items()}
def boot(res, cfg, base='ORIG@disk', n=10000, seed=15):
    rng = np.random.default_rng(seed); out = {}
    for pop in MAIN + ['CLASH_centre']:
        a = host_losses(res, cfg, pop); b = host_losses(res, base, pop); hs = sorted(a); A = np.array([a[h] for h in hs]); B = np.array([b[h] for h in hs])
        ii = rng.integers(0, len(hs), (n, len(hs))); d = np.sqrt(A[ii].mean(1)) - np.sqrt(B[ii].mean(1))
        out[pop] = dict(delta=round(float(np.sqrt(A.mean()) - np.sqrt(B.mean())), 6), ci95=[round(float(x), 6) for x in np.percentile(d, [2.5, 97.5])])
    return out
if __name__ == '__main__':
    for f in sys.argv[1:]:
        res = json.load(open(f)); t = table(res); base = t['ORIG@disk']
        print('==', res['mask'])
        for cfg, r in t.items():
            print('%-11s' % cfg, ' '.join('%s %.4f' % (k, r[k]) for k in MAIN), 'lens %.2f' % r['lens'], 'centre %.4f' % r['CLASH_centre'], 'bound', r['bound'], 'failed', r['failed'])
        for cfg, p in passes(t, base).items(): print('%-11s' % cfg, p)
