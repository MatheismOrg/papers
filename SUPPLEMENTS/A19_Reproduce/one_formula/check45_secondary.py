"""T45 secondary checks (PROTOCOL_T45, not pass criteria): KiDS early/late split and the Milky Way's switch, T45 against T44."""
import json, numpy as np
import harness_T45 as M
import kids_split as KS
H = M.H; B = M.B; S = M.S; Z = S.Z
B.FROZEN['abc'] = np.array([11.7481, .8178, -.7119])
gs = H.disks(); SS = KS.setup(gs)
lmt = H.fit(M.MODEL, .1, None, None, None, [0])['th'][0]
out = dict(Mt_T45=lmt, kids_T44=KS.check('MOGF0IESR', np.array([11.1182]), SS), kids_T45=KS.check(M.MODEL, np.array([lmt]), SS))
mstar = 5.43e10; a_ = np.linspace(1e-3, 60, 6000); discs = {'HI': (53.1, 7., 4., .085), 'H2': (2180., 1.5, 12., .045)}
mgas = sum(2 * np.pi * np.sum(a_ * p[0] * 1e6 * np.exp(-p[2] / a_ - a_ / p[1])) * (a_[1] - a_[0]) for p in discs.values())
args = (np.array([mstar]), np.array([mstar + mgas]), np.array([5.75]))
for nm, t in (('T44', 11.1182), ('T45', lmt)):
    out['MW_' + nm] = dict(q=float(Z.switch('MOG', np.array([t]), *args)[0]), M_halo=float(M.P.G.mass(np.array([t]), *args, np.array([0.]), 0.)[0]))
print(json.dumps(out, indent=0)); json.dump(out, open('res45_secondary.json', 'w'), indent=1)
