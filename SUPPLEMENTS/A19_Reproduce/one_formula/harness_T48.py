"""T48 (PROTOCOL_T48.txt): T45's formula with, for each CLASH cluster, M_p = all baryons (stars + hot gas) inside the
cluster's outermost benchmark radius, g_bar r^2 / G there, in place of the brightest galaxy's stars. MOGF0IESB."""
import sys
import numpy as np
import harness_T45 as P
H = P.H; Y = P.Y
MODEL = 'MOGF0IESB'; H.MODELS[MODEL] = H.MODELS['MOGF0']
_rows = H.cluster_rows
def cluster_rows():
    rows = [dict(r) for r in _rows()]; best = {}
    for r in rows:
        if r['kind'] == 'MaNGA_BCG': continue
        if r['name'] not in best or r['r'] > best[r['name']][0]: best[r['name']] = (r['r'], r['gb'] / H.ACC * r['r'] ** 2 / H.G)
    for r in rows:
        if r['kind'] != 'MaNGA_BCG': r['mp'] = max(r['ms'], best[r['name']][1])
    return rows
H.cluster_rows = cluster_rows
_fit, _gpred, _esd = H.fit, H.gpred, H.esd
def fit(model, *a): return _fit(P.MODEL if model == MODEL else model, *a)
def gpred(model, th, a): return _gpred(P.MODEL if model == MODEL else model, th, a)
def esd(model, th, st): return _esd(P.MODEL if model == MODEL else model, th, st)
H.fit, H.gpred, H.esd = fit, gpred, esd
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
