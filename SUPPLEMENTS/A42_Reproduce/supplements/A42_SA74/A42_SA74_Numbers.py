#!/usr/bin/env python3
"""Collect every number the SA74 additions to A42 quote, straight from the run files, into sa74_numbers.json."""
import json, sys, math, glob
sys.path.insert(0, '/home/claude/sa74/scripts')
import a42_classify as C

R = '/home/claude/sa74/runs/'
def jl(p): return [json.loads(l) for l in open(R + p)]
N = {}
pub = json.load(open(R + 'A42_Search_Results.json'))
bands = json.load(open(R + 'bands.json'))
N['bands'] = bands
m4 = 2.1443

# ---- near misses against the bands (published runs)
def near(runs):
    out = []
    for r in runs:
        for s in r.get('series', []):
            if 3.5 <= s['mm'] <= 4.5 and 3.5 <= s['mid'] <= 4.5:
                z = C.zscores(s, 4)
                out.append({'rule': r['rule'], 'n': s['n'], 'z_chain': z['chain_d'], 'z_mm': z['mm'], 'z_mid': z['mid'],
                            'ratio': s['chain'] / (m4 * s['n'] ** 0.25)})
    return out
n1 = near(pub['long_runs']); n2 = near(pub['long_runs_second_pass'])
N['near_first'] = {'count': len(n1), 'z_chain': [min(x['z_chain'] for x in n1), max(x['z_chain'] for x in n1)],
                   'ratio': [min(x['ratio'] for x in n1), max(x['ratio'] for x in n1)]}
N['near_second'] = {'count': len(n2), 'z_chain': [min(x['z_chain'] for x in n2), max(x['z_chain'] for x in n2)],
                    'ratio': [min(x['ratio'] for x in n2), max(x['ratio'] for x in n2)]}

# ---- two-dimensional pair against the 2D bands
pair = ['{{x,y},{x,z}}->{{x,y},{x,w},{z,x}}', '{{x,y},{z,y}}->{{x,y},{y,z},{w,y}}']
zs = []; ratios = []; firsts = []
for src in ('long_runs', 'long_runs_second_pass'):
    for r in pub[src]:
        if r['rule'] not in pair: continue
        for s in r['series']:
            z = C.zscores(s, 2); zs.append(z['chain_d'])
            mu, sd = C.band(2, s['n'], 'chain_d')
            L_sprk = 1.9544 * s['n'] ** (1 / mu)
            ratios.append(s['chain'] / L_sprk)
            if s['n'] > 20000: firsts += [s['mm'], s['mid']]
N['pair2d'] = {'z_chain': [min(zs), max(zs)], 'chain_longer': [min(ratios), max(ratios)],
               'first_two_largest': [min(firsts), max(firsts)]}

# ---- update orders
orders = ('oldest', 'newest', 'random1', 'random2')
runs = {}
for o in orders:
    for r in jl(f'orders_{o}.jsonl'):
        runs.setdefault(r['rule'], {})[o] = r
ci2 = {x['rule']: x for x in jl('ci_22_32.jsonl')}
for x in json.load(open(R + 'ci_long_sample.json')): ci2[x['rule']] = x
ci3 = {x['rule']: x for x in json.load(open(R + 'ci_long27_steps3.json'))}
def ci_status(rule):
    if ci2[rule]['causally_invariant_within_steps']: return '2'
    if ci3.get(rule, {}).get('pass'): return '3'
    return 'no'
def kind(r):
    if r.get('status') == 'stops or stays small': return 'stops'
    if r['gens'] >= 3000 and r['events'] <= r['gens'] * 2.01: return 'single chain'
    fit = r.get('fit', {})
    if fit.get('exp_rate', 0) > math.log(1.1): return 'exponential'
    return 'polynomial'
def largest(r):
    ser = r.get('series', [])
    if not ser: return None
    return max(ser, key=lambda s: s['n'])
ordtab = {}
reads4 = []; reads4_small = []; maxchain = (0, None)
for rule, d in runs.items():
    row = {'ci': ci_status(rule), 'group': d['oldest'].get('group')}
    for o in orders:
        if o not in d: continue
        r = d[o]
        L = largest(r)
        row[o] = {'kind': kind(r), 'events': r['events'], 'gens': r['gens'], 'links_mean': r.get('mean_links'),
                  'links_max': r.get('max_links'),
                  'largest': ({k: L[k] for k in ('n', 'mm', 'mid', 'chain_d')} if L else None),
                  'range': ({k: [min(s[k] for s in r['series']), max(s[k] for s in r['series'])]
                             for k in ('n', 'mm', 'mid', 'chain_d')} if r.get('series') else None)}
        for s in r.get('series', []) + [dict(r[l], src=l) for l in ('root', 'mid') if l in r]:
            if C.reads(s, 4):
                (reads4 if s['n'] >= 1000 else reads4_small).append((rule, o, s['n'], s['mm'], s['mid'], s['chain_d']))
            if s['n'] >= 120 and s['chain_d'] > maxchain[0]:
                maxchain = (s['chain_d'], (rule, o, s['n'], s['mm'], s['mid']))
        row[o]['candidate'] = C.candidate(r, 4)[0]
    ordtab[rule] = row
N['orders'] = ordtab
N['orders_reads4_big'] = reads4
N['orders_reads4_small'] = reads4_small
N['orders_maxchain'] = maxchain
# published runs: small intervals reading 4 within bands
pub_small = []
for src in ('long_runs', 'long_runs_second_pass'):
    for r in pub[src]:
        for s in r['series']:
            if C.reads(s, 4): pub_small.append((r['rule'], src, s['n'], s['mm'], s['mid'], s['chain_d']))
N['pub_reads4'] = pub_small

# ---- causal invariance, whole class
ci_all = jl('ci_22_32.jsonl')
pass2 = {x['rule'] for x in ci_all if x['causally_invariant_within_steps']}
try:
    s3 = jl('ci_steps3.jsonl')
    pass3 = {x['rule'] for x in s3 if x['causally_invariant_within_steps']}
    N['ci_steps3_done'] = len(s3)
except FileNotFoundError:
    pass3 = set(); N['ci_steps3_done'] = 0
import csv
beh = {row['rule']: row for row in csv.DictReader(open(R + 'A42_All_Rules.csv')) if row['class_of_rules'].startswith('two relations into three')}
from collections import Counter
N['ci'] = {'rules': len(ci_all), 'pass2': len(pass2), 'pass3_extra': len(pass3), 'pass_total': len(pass2 | pass3),
           'pass2_behavior': dict(Counter(beh[r]['behavior'] for r in pass2)),
           'pass3_behavior': dict(Counter(beh[r]['behavior'] for r in pass3)),
           'all_behavior': dict(Counter(b['behavior'] for b in beh.values()))}
two_d = [row['rule'] for row in beh.values() if row['behavior'] == 'polynomial' and row['root_myrheim_meyer'] and
         1.75 <= float(row['root_myrheim_meyer']) <= 2.25 and 1.75 <= float(row['root_midpoint']) <= 2.25]
N['ci']['two_d_rules'] = len(two_d)
N['ci']['two_d_pass'] = sorted(r for r in two_d if r in (pass2 | pass3))

# ---- links: sprinkled 2D and 4D
N['links_sprinkled'] = json.load(open(R + 'links_sprinkled.json'))

# ---- transitive percolation grid (protocol) and tuned trajectory, originary percolation, de Sitter
tp = jl('tp.jsonl')
grid = []
for r in tp:
    L = r.get('root')
    grid.append({'N': r['N'], 'p': r['p'], 'pN': r['pN'], 'seed': r['seed'],
                 'root': ({k: L[k] for k in ('n', 'mm', 'mid', 'chain_d')} if L else None),
                 'candidate': C.candidate(r, 4)[0]})
N['tp_grid'] = grid
N['tp_grid_candidates'] = sum(g['candidate'] for g in grid)
tuned = {}
for r in jl('tp_tuned.jsonl'):
    if r['N'] < 4000: continue
    t = tuned.setdefault(r['N'], {'p': r['p'], 'pN': r['pN'], 'seeds': []})
    t['seeds'].append({k: r['root'][k] for k in ('n', 'mm', 'mid', 'chain_d', 'chain')} | {'candidate': C.candidate(r, 4)[0],
                       'z': C.zscores(r['root'], 4)})
for Nn, t in tuned.items():
    for k in ('n', 'mm', 'mid', 'chain_d'):
        t[k] = sum(s[k] for s in t['seeds']) / len(t['seeds'])
N['tp_tuned'] = tuned
op = {}
for r in jl('op.jsonl'):
    op[f"{r['N']}_{r['p']}"] = {'N': r['N'], 'p': r['p'], 'root': {k: r['root'][k] for k in ('n', 'mm', 'mid', 'chain_d', 'chain')},
                                'series': [{k: s[k] for k in ('n', 'mm', 'mid', 'chain_d')} for s in r['series']],
                                'candidate': C.candidate(r, 4)[0]}
N['op'] = op
N['ds_table'] = json.load(open(R + 'ds_table.json'))
big = json.load(open(R + 'ds_big.json'))
dsb = {}
for x in big:
    g = dsb.setdefault(x['efolds'], [])
    g.append(x)
N['ds_big'] = {k: {'n': v[0]['n'], 'reps': len(v), **{key: [min(e[key] for e in v), max(e[key] for e in v), sum(e[key] for e in v) / len(v)] for key in ('mm', 'mid', 'chain_d', 'chain')}} for k, v in dsb.items()}
json.dump(N, open(R + 'sa74_numbers.json', 'w'), indent=1, default=str)
if __name__ == '__main__':
    print(json.dumps({k: N[k] for k in ('near_first', 'near_second', 'pair2d', 'orders_reads4_big', 'orders_maxchain', 'ci', 'tp_grid_candidates')}, indent=1, default=str))
    print('small 4D-reading intervals (orders):', len(N['orders_reads4_small']), 'published:', len(N['pub_reads4']))
    for x in N['orders_reads4_small'][:20]: print('   ', x)
    for x in N['pub_reads4'][:20]: print('   pub', x)
    for Nn, t in sorted(tuned.items()): print('TUNED', Nn, round(t['p'], 6), round(t['pN'], 1), round(t['n']), round(t['mm'], 2), round(t['mid'], 2), round(t['chain_d'], 2))
    print('DS big', json.dumps(N['ds_big']))
