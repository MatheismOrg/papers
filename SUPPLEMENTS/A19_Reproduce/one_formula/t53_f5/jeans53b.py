"""SA73 T53: is Tian et al.'s acceleration radius r* (solved from their g_bar with the projected Sersic mass) also the radius of
their g_obs? Jeans rebuild (jeans53.py) at r*. Usage: python3 -I jeans53b.py <table1 txt>"""
import sys, re, json, numpy as np
from scipy.optimize import brentq
sys.argv, tab = sys.argv[:1] + ['/dev/null'], sys.argv[1]
src = open('jeans53.py').read().split('rows = []')[0]; exec(src)
rows = []
for l in open(tab):
    if not l.strip(): continue
    ident, rest = l.strip().split(maxsplit=1); v = [float(t) for t in re.findall(r'[-+]?\d+(?:\.\d+)?', rest)]
    rows.append(dict(id=ident.replace('∗', '').replace('*', ''), z=v[0], n=v[1], re=v[2], xl=v[3], m=v[4], b=v[5], lgM=v[6], sig=v[8], lgb=v[10], lgo=v[12]))
out = []
for q in rows:
    f = lambda x: np.log10(G * 10 ** q['lgM'] * f2d(x, q['n']) / (x * q['re']) ** 2 * ACC) - q['lgb']
    xs = brentq(f, 0.01, 50); r = xs * q['re']
    out.append((q['id'], q['xl'], xs, np.log10(gobs(r, q['n'], q['re'], q['sig'], q['m'], q['b'], q['re']) * ACC) - q['lgo']))
d = np.array([o[3] for o in out]); xs = np.array([o[2] for o in out]); xl = np.array([o[1] for o in out])
print('g_obs rebuilt at r* minus Tian: median %+.3f sd %.3f max|dev| %.3f' % (np.median(d), np.std(d), np.max(np.abs(d - np.median(d)))), ' corr with xl', round(float(np.corrcoef(xl, d)[0, 1]), 3))
json.dump([[o[0], o[1], round(o[2], 4), round(float(o[3]), 4)] for o in out], open('res53_jeans_rstar.json', 'w'))
