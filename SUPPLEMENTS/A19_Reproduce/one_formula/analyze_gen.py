"""General comparison: the paper's rule, T17b MOSS, T19 MOGC/MOG and any further result files; progress against a baseline.
Usage: python3 analyze_gen.py <mask> <baseline cfg> <file pattern with %s for the mask> ..."""
import sys, json, glob
import numpy as np
import analyze as A
A.NCOEF.update({k: 5 for k in ['SUBH', 'MOSS', 'CAUD', 'MOG', 'MOGC', 'MOSC', 'MOGa', 'MOGb', 'MOGd', 'MOGe', 'MOGC2', 'MOGF', 'MOGF85', 'MOGF0', 'MOGH5', 'MOGHf', 'MOGFH5', 'MOGFI', 'MOGFK']}); A.NCOEF.update({k: 4 for k in ['MOGF0I', 'MOGF0IE', 'MOGF0IEX', 'MOGF0IEY', 'MOGF0IES', 'MOGF0IESY', 'MOGF0IESG', 'MOGF0IESR']}); A.NCOEF.update({'MOGF0IESM': 3, 'MOGF0IESB': 3})
SIX = A.MAIN + ['lens']
mask, basecfg = sys.argv[1], sys.argv[2]
recs, summ, new = [], {}, []
for f, keep in [('res_%s.json', ['ORIG@disk']), ('res17b_%s.json', ['MOSS@0.1']), ('res19_%s.json', ['MOG@0.1', 'MOGC@0.1'])] + [(p, None) for p in sys.argv[3:]]:
    r = json.load(open(f % mask)); kk = keep or r['configs']
    if keep is None: new += [c for c in kk if c not in new]
    recs += [x for x in r['records'] if x['config'] in kk]; summ.update({c: r['summary'][c] for c in kk})
res = dict(mask=mask, summary=summ, records=recs); t = A.table(res); base = t['ORIG@disk']; bl = t[basecfg]
print('==', mask)
for cfg, r in t.items():
    print('%-11s' % cfg, ' '.join('%s %.4f' % (k, r[k]) for k in A.MAIN), 'lens %.2f' % r['lens'], 'centre %.4f' % r['CLASH_centre'], 'bound', r['bound'], 'by group', summ[cfg]['lensing_by_group'])
for cfg in new:
    r = t[cfg]; better = [k for k in SIX if r[k] < bl[k]]; worse = [k for k in SIX if r[k] > bl[k] * (1 + A.TOL[k])]
    th = np.array([x['fit']['th'] for x in recs if x['config'] == cfg])
    print('%-11s vs %s: better %s; worse beyond tol %s; progress %s | vs paper %s | coef median %s min %s max %s' % (cfg, basecfg, better, worse, len(better) >= 2 and not worse,
          A.passes({cfg: r}, base)[cfg]['improvement'], np.round(np.median(th, 0), 3), np.round(th.min(0), 3), np.round(th.max(0), 3)))
if mask == 'primary':
    for cfg in new:
        b = A.boot(res, cfg, base=basecfg); print('  %-9s bootstrap vs %s:' % (cfg, basecfg), ' '.join('%s %+.4f [%+.4f, %+.4f]' % (k, v['delta'], *v['ci95']) for k, v in b.items()))
