"""SA73 T45 exploration (all data): per-dwarf residuals of the frozen T44 formula."""
import numpy as np
import harness_T44 as M
H = M.H
B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
rows = H.dwarf_rows(); a = H.arr(rows); th = np.array([11.1182])
er = np.log10(a['go'] / H.gpred('MOGF0IESR', th, a))
keys = [k for k in rows[0].keys()]
print('row keys', keys)
o = np.argsort(er)
for i in o:
    r = rows[i]
    extra = ' '.join(f'{k}={r[k]}' for k in keys if k not in ('name', 'kind', 'r', 'gb', 'go', 'ms', 'mp', 're', 'z'))[:120]
    print(f"{r['name']:22s} res {er[i]:+.2f}  logM* {np.log10(r['ms']):5.2f}  Re {r['re']:6.3f}  r {r['r']:6.3f}  loggo {np.log10(r['go']):6.2f} {extra}")
print('rms', np.sqrt(np.mean(er**2)), 'n', len(er))
