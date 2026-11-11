"""T18 analysis: the paper's rule (T15 run), SUBH@0.1 (T16), MOSS@0.1 (T17b) and the T18 runs; same partitions.
Usage: python3 analyze18.py <mask>"""
import sys, json
import numpy as np
import analyze as A
A.NCOEF.update({'SUBH': 6, 'MOSS': 4, 'CAUD': 5, 'CAU5': 4, 'CAU1': 4})
SIX = A.MAIN + ['lens']
mask = sys.argv[1]
recs, summ = [], {}
for f, keep in [('res_%s.json', ['ORIG@disk']), ('res16_%s.json', ['SUBH@0.1']), ('res17b_%s.json', ['MOSS@0.1']), ('res18_%s.json', None)]:
    r = json.load(open(f % mask)); kk = keep or r['configs']
    recs += [x for x in r['records'] if x['config'] in kk]; summ.update({c: r['summary'][c] for c in kk})
res = dict(mask=mask, summary=summ, records=recs)
t = A.table(res); base = t['ORIG@disk']; moss = t['MOSS@0.1']
print('==', mask)
for cfg, r in t.items():
    print('%-11s' % cfg, ' '.join('%s %.4f' % (k, r[k]) for k in A.MAIN), 'lens %.2f' % r['lens'], 'centre %.4f' % r['CLASH_centre'],
          'bound', r['bound'], 'by group', summ[cfg]['lensing_by_group'])
for cfg, r in t.items():
    if not cfg.startswith('CAU'): continue
    better = [k for k in SIX if r[k] < moss[k]]; worse = [k for k in SIX if r[k] > moss[k] * (1 + A.TOL[k])]
    print('%-11s vs MOSS: better on %s; worse beyond tolerance on %s; progress %s' % (cfg, better, worse, len(better) >= 2 and not worse))
    print('%-11s vs paper rule:' % cfg, A.passes({cfg: r}, base)[cfg])
    ff = [x['fit'] for x in recs if x['config'] == cfg]
    s = np.array([f['disk_fit']['s'] for f in ff]); rt = np.array([f['th'][0] for f in ff])
    print('%-11s s median %.3f range %.3f-%.3f, s at bound %d; log R_t median %.3f range %.3f-%.3f' % (cfg, np.median(s), s.min(), s.max(),
          sum(f['disk_fit']['s_at_bound'] for f in ff), np.median(rt), rt.min(), rt.max()))
if mask == 'primary':
    for cfg in t:
        if cfg.startswith('CAU'): print(cfg, 'bootstrap vs MOSS@0.1:', json.dumps(A.boot(res, cfg, base='MOSS@0.1')))
