"""SA73 T56 look (all data): how the lensing chi-square of T48 and T56 depends on the cold gas assumed for the stacks' late
types (the SPARC HI relation, from HI-rich galaxies), scaled by f. Usage: python3 lensgas56.py (PYTHONPATH = dmfit)"""
import sys, json, numpy as np
import harness_T56 as M6
import harness_T48 as M8
P = M6.P; H = M6.H; B = M6.B
gs = H.disks(); (p1, p0), (q1, q0) = H.pop_coeff(gs)
print('SPARC relations: log Re = %.3f + %.3f logM*; log MHI = %.3f + %.3f logM*  -> MHI/M* at 1e9: %.2f, 1e10: %.2f, 1e11: %.2f' % (p0, p1, q0, q1, 10 ** (q0 + q1 * 9 - 9), 10 ** (q0 + q1 * 10 - 10), 10 ** (q0 + q1 * 11 - 11)))
for model, abc, mu, lmt in ((M8.MODEL, [11.7481, .8178, -.7119], 1., 11.339), (M6.MODEL, [11.602, 1.0003, -1.363], 4.41, 11.353)):
    B.FROZEN['abc'] = np.array(abc); M6.U['mu'] = mu; th = np.array([lmt])
    for f in (1., .5, .3, .1, 0.):
        row = []
        for drop in (False, True):
            L = H.Lens(gs, drop)
            for st in L.stacks: st.mp = st.ms + f * (st.mp - st.ms); st.bary = st.bary * 0 + st.mp[:, None] / np.pi * (st.bary / np.maximum((st.ms + (st.mp - st.ms) / max(f, 1e-9) if f > 0 else st.ms)[:, None] / np.pi, 1e-300)) if False else st.bary
            c2 = [round(float(np.sum(L.white(model, th, [j]) ** 2)), 1) for j in range(4)]; row.append((round(sum(c2), 1), c2))
        print(model, 'gas x', f, 'primary', row[0], 'strict', row[1], flush=True)
