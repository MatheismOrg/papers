#!/usr/bin/env python3
"""Writes PAPER_NUMBERS.json: the numbers Section 8.10 of A42 quotes for its larger grids, with their standard errors, from the result
files of the ORIGINAL chains (regge_checks.numbers). It was run once, in the tree the paper was built in (SA74 T32, November 10, 2026);
the zip holds its output, not those result files, which regenerate. The bootstrap error of the zero crossing (multihist.py, T32) comes
from a rerun of multihist.py on the original chains whose other numbers equal the paper's file exactly (checked here)."""
import json, shutil, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); import regge_checks as RC
SA = Path('/home/claude/sa74'); N = SA / 'regge/r2/night'; FS = SA / 'regge/r2/fss'
old = json.load(open(N / 'multihist_T30.json')); new = json.load(open(N / 'multihist_T32x.json'))
strip = lambda d: {k: ({L: {a: b for a, b in v.items() if a != 'crossing_boot_sd'} for L, v in x.items()} if k == 'per_L' else x) for k, x in d.items() if k != 'description'}
assert strip(old) == strip(new), 'the rerun with the crossing error must reproduce multihist_T30.json'
with tempfile.TemporaryDirectory() as t:
    t = Path(t); (t / 'n').mkdir(); (t / 'f').mkdir()
    for f in ('a1_results_T31.json', 'multihist_skip50_T30.json', 'au_results.json'): shutil.copy(N / f, t / 'n' / f)
    json.dump(new, open(t / 'n' / 'multihist_T30.json', 'w')); shutil.copy(FS / 'fss_results.json', t / 'f' / 'fss_results.json')
    P = RC.numbers(t / 'n', t / 'f')
for L, v in P['au_crossing'].items(): assert abs(v['v'] - v['v_check']) < 1e-9, (L, v)
P = dict(description=__doc__, **P)
json.dump(P, open(Path(__file__).resolve().parent / 'PAPER_NUMBERS.json', 'w'), indent=1)
r = lambda x, n=4: round(x, n)
for k in RC.KS:
    print('k', k, 'power %+.3f(%.3f) lowbias %+.3f(%.3f)' % (P['power'][k]['v'], P['power'][k]['se'], P['power_lowbias'][k]['v'], P['power_lowbias'][k]['se']),
          'chi', ' '.join('%s:%.2f(%.2f)' % (L, v['v'], v['se']) for L, v in P['fixed_k'][k].items()))
print('hamber power %.3f' % P['hamber_power'], 'peak power', P['peak_power'], P['peak_power_second_halves'])
print('peak', {L: (r(v['v'], 2), r(v['se'], 2)) for L, v in P['peak'].items()})
print('crossing', {L: (r(v['v'], 5), r(v['se'], 5)) for L, v in P['crossing'].items()})
print('start z', P['start_check_z'])
print('tau', {L: (round(v['median']), v['n'], r(v['se_ln'], 2)) for L, v in P['tau'].items()})
print('au crossing', {L: (r(v['v']), r(v['se'])) for L, v in P['au_crossing'].items()})
print('au chi 0.46', {L: (r(v['v'], 3), r(v['se'], 3)) for L, v in P['au_chi_046'].items()})
print('zmax', {m: round(RC.zmax(m), 2) for m in (1, 2, 3, 4, 5, 6, 8)})
