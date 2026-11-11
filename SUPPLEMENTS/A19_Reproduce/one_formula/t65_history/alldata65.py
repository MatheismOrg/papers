"""PROTOCOL_T65, Test 3: t and kappa of the core under the gate "any cold gas" (harness_T65.py), on all points of the 149
benchmark disks with T61's weights, A, B, C frozen at the all-disk values; and the same fit under the gas gate, as a check
that it returns T61's 0.1501 and 0.6229. Usage: PYTHONPATH=../dmfit python3 alldata65.py"""
import json, numpy as np
import harness_T65 as M5
M61 = M5.M61; H = M5.H; names = [g['name'] for g in H.disks()]; abc = np.array([11.7481, .8178, -.7119])
(t, k), loss, att, bh, npts = M5.fit_disks(names, abc); (t6, k6), loss6, att6, bh6, n6 = M61.fit_disks('MOGF0IESBC', names, abc)
out = dict(t=t, kappa=k, loss=loss, bound_hits=bh, points=npts, gas_gate=dict(t=t6, kappa=k6, loss=loss6), loss_over_gas_gate=loss / loss6)
print(json.dumps(out)); json.dump(out, open('res65_alldata_fit.json', 'w'), indent=1)
