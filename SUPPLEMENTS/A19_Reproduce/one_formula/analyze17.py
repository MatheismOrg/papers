"""T17 analysis: merges the paper's rule (ORIG@disk, T15 run) and SUBH@0.1 (T16 run) with the T17 runs; same partitions.
Usage: python3 analyze17.py <mask> <res17 file> [more result files]"""
import sys, json
import numpy as np
import analyze as A
A.NCOEF.update({'SUBH': 6, 'MZ': 6, 'GASC': 7, 'MZGASC': 7, 'ZEVO': 7, 'BEST1': 7, 'BEST1J': 7, 'ORIGZ': 3, 'KNEESFZ': 6, 'MOSS': 4, 'MOSH': 5})
SIX = A.MAIN + ['lens']
mask = sys.argv[1]
recs, summ = [], {}
for f in ['res_%s.json' % mask, 'res16_%s.json' % mask] + sys.argv[2:]:
    r = json.load(open(f)); keep = [c for c in r['configs'] if f.startswith('res17') or f.startswith('res18') or c in ('ORIG@disk', 'SUBH@0.1')]
    recs += [x for x in r['records'] if x['config'] in keep]; summ.update({c: r['summary'][c] for c in keep})
res = dict(mask=mask, summary=summ, records=recs)
t = A.table(res); base = t['ORIG@disk']; subh = t['SUBH@0.1']
print('==', mask)
for cfg, r in t.items():
    print('%-12s' % cfg, ' '.join('%s %.4f' % (k, r[k]) for k in A.MAIN), 'lens %.2f' % r['lens'], 'centre %.4f' % r['CLASH_centre'], 'bound', r['bound'], 'failed', r['failed'])
for cfg, r in t.items():
    if cfg in ('ORIG@disk', 'SUBH@0.1'): continue
    better = [k for k in SIX if r[k] < subh[k]]; worse = [k for k in SIX if r[k] > subh[k] * (1 + A.TOL[k])]
    print('%-12s vs SUBH: better on %s; worse beyond tolerance on %s; progress %s' % (cfg, better, worse, len(better) >= 2 and not worse))
    print('%-12s vs paper rule:' % cfg, A.passes({cfg: r}, base)[cfg])
if mask == 'primary':
    for cfg in t:
        if cfg in ('ORIG@disk', 'SUBH@0.1'): continue
        print(cfg, 'bootstrap vs SUBH@0.1:', json.dumps(A.boot(res, cfg, base='SUBH@0.1')))
