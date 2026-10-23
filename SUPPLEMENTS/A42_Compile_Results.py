#!/usr/bin/env python3
"""A42 supplement: compile the run files of A42_Rule_Search.py into A42_Search_Results.json (validation, census of both
classes, long runs) and A42_All_Rules.csv (one row per rule). Run it in the folder holding the run files:
validation.json, res_all.jsonl and rules_22_32.json (the full class), s42_all.jsonl and rules_22_42_sample.json (the sample),
deep2.jsonl and deep2_rules.json (the 27 long runs, each rule with its group), and tl_*.jsonl (reruns of slow single chains
with 30 seconds). The second pass adds s3242.jsonl and s2333.jsonl (samples of three relations into four and two
three-node relations into three, with reruns in tl2_*.jsonl), long.jsonl and long_rules.json (ten rules to 120,000
updates) and desitter.json, when present. usage: python3 -I A42_Compile_Results.py OUT_DIR"""
import json, csv, sys, glob
from collections import Counter
OUT = sys.argv[1]
def load(p): return [json.loads(l) for l in open(p)]
def cls(r):
    if r['status'] != 'measured': return 'stops or stays small'
    if r.get('tree'): return 'tree'
    return 'exponential' if (r.get('fit') or {}).get('exp_rate', 0) > 0.1 else 'polynomial'
def band(p):
    if p is None: return 'no fit'
    return 'near 0 (volume about 1D)' if p < 0.5 else 'near 1 (about 2D)' if p < 1.5 else 'near 2 (about 3D)' if p < 2.5 else 'near 3 (about 4D)' if p < 3.5 else '3.5 or more'
near = lambda d, a, b: d is not None and a <= d <= b
census = {}; allrows = []
import os
classes = [('two relations into three, all', 'res_all.jsonl', 'tl_*.jsonl'),
           ('two relations into four, random sample', 's42_all.jsonl', 'tl_*.jsonl')]
if os.path.exists('s3242.jsonl'): classes.append(('three relations into four, random sample', 's3242.jsonl', 'tl2_*.jsonl'))
if os.path.exists('s2333.jsonl'): classes.append(('two three-node relations into three, random sample', 's2333.jsonl', 'tl2_*.jsonl'))
for label, path, rerun_glob in classes:
    rows = load(path); seen = set(); uniq = []
    # slow single chains cut off by the 8-second limit below 300 updates were rerun with 30 seconds; count them by the rerun
    rerun = {r['rule']: r for f in sorted(glob.glob(rerun_glob)) for r in load(f)}
    rows = [dict(rerun[r['rule']], rerun_seconds=30) if r['rule'] in rerun else r for r in rows]
    for r in rows:
        if r['rule'] not in seen: seen.add(r['rule']); uniq.append(r)
    rows = uniq
    c = Counter(cls(r) for r in rows)
    poly = [r for r in rows if cls(r) == 'polynomial']
    census[label] = {
        'rules': len(rows), 'classes': dict(c),
        'polynomial_growth_power': dict(Counter(band((r.get('fit') or {}).get('power')) for r in poly)),
        'exponential_with_interval_120_or_more': sum(1 for r in rows if cls(r) == 'exponential' and r.get('root')),
        'trees_that_are_single_chains': sum(1 for r in rows if cls(r) == 'tree' and max(r['per_gen']) <= 1),
        'read_four_on_mm_and_midpoint': [{'rule': r['rule'], 'interval': lab, **r[lab]} for r in rows for lab in ('root', 'mid')
                                         if r.get(lab) and near(r[lab]['mm'], 3.5, 4.5) and near(r[lab]['mid'], 3.5, 4.5)],
        'read_two_on_mm_and_midpoint_in_largest_intervals': [{'rule': r['rule'], **r['root']} for r in rows if cls(r) == 'polynomial'
                                         and r.get('root') and near(r['root']['mm'], 1.75, 2.25) and near(r['root']['mid'], 1.75, 2.25)],
        'bits_range': [min(r['bits'] for r in rows), max(r['bits'] for r in rows)],
        'hit_8_second_limit': sum(1 for r in rows if r['sec'] >= 7.9 and 'rerun_seconds' not in r),
        'rerun_with_30_seconds': sum(1 for r in rows if 'rerun_seconds' in r)}
    for r in rows:
        f = r.get('fit') or {}; ro = r.get('root') or {}; mi = r.get('mid') or {}
        allrows.append([label, r['rule'], r['bits'], cls(r), r['events'], r['gens'], f.get('exp_rate'), f.get('power'),
                        ro.get('n'), ro.get('mm'), ro.get('mid'), mi.get('n'), mi.get('mm'), mi.get('mid')])
groups = {r['rule']: r['group'] for r in json.load(open('deep2_rules.json'))}
g2 = {r['rule']: r.get('group') for r in json.load(open('long_rules.json'))} if os.path.exists('long_rules.json') else {}
deep = []
for r in load('deep2.jsonl'):
    r = dict(r); r['group'] = groups[r['rule']]; r.pop('per_gen', None); deep.append(r)
res = {
    'about': 'A42 search, first and second passes (SO70 T213 and T214, October 11, 2026). Produced by A42_Rule_Search.py and compiled by A42_Compile_Results.py; see A42_Search_README.md.',
    'classes': {'stops or stays small': 'fewer than 300 updates before no match is left or the run ends',
                'tree': 'every update has one direct cause, so the causal order is a tree (dimension 1)',
                'exponential': 'updates per generation grow by more than 10 percent per generation over the second half of the run',
                'polynomial': 'the rest; updates per generation grow as a power p of the generation, a volume dimension of about p + 1'},
    'runs': {'scan': 'up to 3,000 updates, 400 generations, 60,000 relations or 8 seconds per rule',
             'long': '27 rules up to 20,000 updates or 60 seconds, with a size series of intervals'},
    'intervals': 'root: from the root event to the six events with the largest pasts; mid: from the event 30 percent of the way up '
                 'with the largest interval; series: one top event, bottom tips at growing depth; at least 120 elements (60 in a series)',
    'long_runs_second_pass': [dict(r, group=g2.get(r['rule'])) for r in (load('long.jsonl') if os.path.exists('long.jsonl') else [])
                              for _ in [r.pop('per_gen', None)]],
    'desitter': json.load(open('desitter.json')) if os.path.exists('desitter.json') else None,
    'four_dimensional_criterion': 'Myrheim-Meyer, midpoint and longest-chain readings all between 3.5 and 4.5 on the same interval',
    'validation': json.load(open('validation.json')),
    'census': census, 'long_runs': deep}
json.dump(res, open(OUT + '/A42_Search_Results.json', 'w'), indent=1)
with open(OUT + '/A42_All_Rules.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['class_of_rules', 'rule', 'bits', 'behavior', 'updates', 'generations', 'exp_rate', 'power',
                'root_n', 'root_myrheim_meyer', 'root_midpoint', 'mid_n', 'mid_myrheim_meyer', 'mid_midpoint'])
    w.writerows(allrows)
print(len(allrows), 'rows;', len(deep), 'long runs')
