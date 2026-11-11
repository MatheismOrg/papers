#!/usr/bin/env python3
"""Apply the SA74 protocol (A42_Protocol_SA74.txt, SHA-256 04eca208...) to run outputs: noise bands interpolated in
log N, each reading's distance from the d-dimensional sprinkled mean in standard deviations, and the candidate rule."""
import json, math, sys

BANDS = json.load(open('/home/claude/sa74/runs/bands.json'))
SIZES = (120, 250, 500, 1000, 2000, 4000)
TESTS = (('mm', 'mm'), ('mid', 'mid'), ('chain_d', 'chain_d'))

def band(d, n, key):
    """(mean, sd) of a test at size n, linear in log n between the protocol's sizes, nearest size outside."""
    if n <= SIZES[0]: r = BANDS[f'{d}_{SIZES[0]}']; return r[key + '_mean'], r[key + '_sd']
    if n >= SIZES[-1]: r = BANDS[f'{d}_{SIZES[-1]}']; return r[key + '_mean'], r[key + '_sd']
    for a, b in zip(SIZES, SIZES[1:]):
        if a <= n <= b:
            w = (math.log(n) - math.log(a)) / (math.log(b) - math.log(a))
            ra, rb = BANDS[f'{d}_{a}'], BANDS[f'{d}_{b}']
            return (ra[key + '_mean'] * (1 - w) + rb[key + '_mean'] * w, ra[key + '_sd'] * (1 - w) + rb[key + '_sd'] * w)

def zscores(iv, d):
    out = {}
    for name, key in TESTS:
        mu, sd = band(d, iv['n'], key)
        out[name] = (iv[key] - mu) / sd
    return out

def reads(iv, d, k=3.0):
    return all(abs(z) <= k for z in zscores(iv, d).values())

def window(iv, lo=3.5, hi=4.5):
    return all(lo <= iv[key] <= hi for _, key in TESTS)

def candidate(run, d=4):
    """Protocol item 5: the largest interval of at least 1,000 elements reads d-dimensional, and the two largest
    intervals of the size series with at least 1,000 elements both do."""
    ivs = []
    for label in ('root', 'mid'):
        if label in run:                       # medians over up to six intervals: use them as one interval reading
            r = run[label]; ivs.append({'n': r['n'], 'mm': r['mm'], 'mid': r['mid'], 'chain_d': r['chain_d']})
    ser = [s for s in run.get('series', []) if s['n'] >= 1000]
    big = [iv for iv in ivs + ser if iv['n'] >= 1000]
    if not big: return False, 'no interval of 1,000 elements'
    largest = max(big, key=lambda iv: iv['n'])
    if not reads(largest, d): return False, 'largest interval fails'
    ser = sorted(ser, key=lambda s: s['n'], reverse=True)[:2]
    if len(ser) < 2: return False, 'fewer than two series intervals of 1,000'
    if not all(reads(s, d) for s in ser): return False, 'series drifts'
    return True, 'candidate'

def best_interval(run, d=4):
    """The interval (any, any size) whose largest |z| across the three tests is smallest."""
    ivs = []
    for label in ('root', 'mid'):
        if label in run:
            r = run[label]; ivs.append(dict(n=r['n'], mm=r['mm'], mid=r['mid'], chain_d=r['chain_d'], src=label))
    for s in run.get('series', []): ivs.append(dict(s, src='series'))
    if not ivs: return None
    def worst(iv): return max(abs(z) for z in zscores(iv, d).values())
    b = min(ivs, key=worst)
    return dict(b, worst_z=round(worst(b), 2), z={k: round(v, 2) for k, v in zscores(b, d).items()})

if __name__ == '__main__':
    path = sys.argv[1]
    for line in open(path):
        run = json.loads(line)
        if run.get('status') == 'stops or stays small': continue
        print(run.get('rule') or (run.get('N'), run.get('p')), candidate(run), best_interval(run))
