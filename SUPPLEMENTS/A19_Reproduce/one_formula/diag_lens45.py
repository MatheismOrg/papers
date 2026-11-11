"""SA73 T45 exploration (all data, not held out): KiDS stack residuals of the frozen T44 formula, by bin, with the
early/late and dark/baryon split of the prediction. Usage: python3 diag_lens45.py"""
import json, numpy as np
import harness_T44 as M
H = M.H; S = M.S; K = H.K; T = M.T; Z = M.Z
B = __import__('harness_T15b'); B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
print('MOGF0 params', H.MODELS['MOGF0'][0])
th = np.array([11.1182]) if len(H.MODELS['MOGF0'][0]) == 1 else np.array([11.7481, .8178, -.7119, 11.1182])
gs = H.disks(); lens = H.Lens(gs, False)
kd = H.PKG / 'data/Lensing/KiDS'; C, _ = K.load_cov(kd / 'Fig-9_Massbins_covmatrix_compact.json')
def parts(st):
    m = S.mass(th, st.ms, st.mp, st.re, st.z, 0.); c = K.c_dm14(m, K.Z0) * M.ratio_all_relaxed(m)
    rv = (3 * m / (4 * np.pi * 200 * K.rho_crit(K.Z0))) ** (1 / 3)
    T.F, T.dF = Z.Fc, Z.dFc
    try: e = T.esd_dm(m, c, rv, st.mp, st.R, 0.)
    finally: T.F, T.dF = Z.F_ORIG, Z.DF_ORIG
    return m, c, np.sum(st.aw * e, axis=2), st.bary
tot = 0
for j, st in enumerate(lens.stacks):
    k = lens.masks[j]; y = lens.y[j]; sig = np.sqrt(np.diag(C)[k + j * 15])
    m, c, dm, by = parts(st); full = st.w @ (dm + by)
    early = st.ms == st.mp
    we = st.w * early; wl = st.w * ~early
    Rm = (st.w @ np.sqrt(st.R[:, :, 0] * st.R[:, :, -1])) * 1e3   # weighted bin-centre radius, kpc
    wh = lens.white('MOGF0IESR', th, [j]); tot += np.sum(wh ** 2)
    print(f'stack {j+1}: chi2 {np.sum(wh**2):.2f}  early weight {we.sum():.2f}  logM (w-mean) {np.log10(st.w @ m):.2f}  c {st.w @ c:.2f}')
    for i, b in enumerate(k):
        print(f'  R {Rm[b]:6.1f} kpc  data {y[i]:7.2f} +- {sig[i]:5.2f}  model {full[b]:7.2f}  (dm {st.w @ dm[:, b]:6.2f} bary {st.w @ by[:, b]:6.2f}; early-only {(we @ (dm + by)[:, b]) / max(we.sum(),1e-9):6.2f} late-only {(wl @ (dm + by)[:, b]) / wl.sum():6.2f})  pull {(y[i]-full[b])/sig[i]:+5.2f}')
print('total chi2', round(tot, 2))
