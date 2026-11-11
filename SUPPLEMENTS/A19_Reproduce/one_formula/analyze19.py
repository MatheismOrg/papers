"""T19 analysis: the paper's rule (T15), SUBH@0.1 (T16), MOSS@0.1 (T17b), CAUD@0.1 (T18) and the T19 runs; same partitions.
Usage: python3 analyze19.py <mask>"""
import sys, json
import numpy as np
import analyze as A
A.NCOEF.update({'SUBH': 6, 'MOSS': 4, 'CAUD': 5, 'MOG': 4, 'MOGC': 5, 'MOSC': 5})
SIX = A.MAIN + ['lens']
mask = sys.argv[1]
recs, summ = [], {}
for f, keep in [('res_%s.json', ['ORIG@disk']), ('res16_%s.json', ['SUBH@0.1']), ('res17b_%s.json', ['MOSS@0.1']),
                ('res18_%s.json', ['CAUD@0.1']), ('res19_%s.json', None)]:
    r = json.load(open(f % mask)); kk = keep or r['configs']
    recs += [x for x in r['records'] if x['config'] in kk]; summ.update({c: r['summary'][c] for c in kk})
res = dict(mask=mask, summary=summ, records=recs)
t = A.table(res); base = t['ORIG@disk']; moss = t['MOSS@0.1']
print('==', mask)
for cfg, r in t.items():
    print('%-10s' % cfg, ' '.join('%s %.4f' % (k, r[k]) for k in A.MAIN), 'lens %.2f' % r['lens'], 'centre %.4f' % r['CLASH_centre'],
          'bound', r['bound'], 'by group', summ[cfg]['lensing_by_group'])
new = [c for c in t if c.split('@')[0] in ('MOG', 'MOGC', 'MOSC')]
for cfg in new:
    r = t[cfg]; better = [k for k in SIX if r[k] < moss[k]]; worse = [k for k in SIX if r[k] > moss[k] * (1 + A.TOL[k])]
    print('%-10s vs MOSS: better on %s; worse beyond tolerance on %s; progress %s' % (cfg, better, worse, len(better) >= 2 and not worse))
    print('%-10s vs paper rule:' % cfg, A.passes({cfg: r}, base)[cfg])
    th = np.array([x['fit']['th'] for x in recs if x['config'] == cfg])
    print('%-10s coefficients median %s min %s max %s' % (cfg, np.round(np.median(th, 0), 3), np.round(th.min(0), 3), np.round(th.max(0), 3)))
if mask == 'primary':
    for cfg in new: print(cfg, 'bootstrap vs MOSS@0.1:', json.dumps(A.boot(res, cfg, base='MOSS@0.1')))
