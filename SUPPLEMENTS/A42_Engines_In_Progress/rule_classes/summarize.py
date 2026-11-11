#!/usr/bin/env python3
"""Summary of the class runs in out/: blocks present, rules, growth classes, readings near four, the closest rules."""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import toehome as T
for cls in sys.argv[1:] or T.ORDER[:4]:
    want = [T.uid(u) for u in T.unit_list(cls)]; got = {}
    for path in sorted(glob.glob(os.path.join(HERE, 'out', cls + '.w*.jsonl'))):
        for line in open(path):
            try: r = json.loads(line)
            except ValueError: continue
            got[r['uid']] = r
    recs = [rec for u in want if u in got for rec in got[u]['rules']]
    meas = [r for r in recs if 'root' in r or 'mid' in r]
    near = [(r['rule'], T.readings(r)) for r in recs if T.near_four(r)]
    zs = sorted((min(x[3] for x in T.readings(r)), r['rule'], min(T.readings(r), key=lambda x: x[3])) for r in recs if T.readings(r))
    two = sum(1 for r in recs if any(all(3.5 <= t <= 4.5 for t in x[2][:2]) for x in T.readings(r)))
    print('%s: blocks %d of %d, rules %d, cut by the step budget %d, with an interval to read %d, core seconds %.0f' % (cls, len(got), len(want), len(recs), sum(r['cut'] for r in recs), len(meas), sum(g['sec'] for g in got.values())))
    print('   stops or stays small %d; trees %d; all three tests between 3.5 and 4.5: %d; first two tests between 3.5 and 4.5: %d' % (sum(r['ev'] < T.MIN_EVENTS for r in recs), sum(bool(r.get('tree')) for r in recs), len(near), two))
    for z, rule, x in zs[:5]: print('   closest: %.1f sd  %s  %s n=%d  %.2f / %.2f / %.2f' % (z, rule, x[0], x[1], x[2][0], x[2][1], x[2][2]))
    for rule, rd in near[:20]: print('   NEAR FOUR:', rule, rd)
