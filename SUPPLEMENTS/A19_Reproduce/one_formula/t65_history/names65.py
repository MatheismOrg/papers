import numpy as np
import harness_T64 as M
H = M.H; B = M.B
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119]); TH = np.array([11.339])
d = H.arr(H.dwarf_rows()); m, c, rv = M.halo(TH, d)
print(sorted(d.keys())); print(len(d['name']))
o = np.argsort(-d['ms'])
for i in o: print('%-22s ms %9.3g re %6.3f r %6.3f r/re %5.2f logMh %5.2f c %5.1f fac %s' % (d['name'][i], d['ms'][i], d['re'][i], d['r'][i], d['r'][i] / d['re'][i], np.log10(m[i]), c[i], ('%.2f' % d['fac'][i]) if 'fac' in d else '-'))
