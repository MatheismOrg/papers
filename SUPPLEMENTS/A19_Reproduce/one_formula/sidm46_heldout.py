"""SA73 T46: held-out scores for the SIDM runs (res46_sidm.json). In each of T45's 60 partitions the cross-section s0 is
chosen on the training disks alone (the halo masses need no fit: standard halos are published; the formula's use the
all-disk A, B, C, so its held-out disk score is slightly optimistic), and the held-out systems of every population are
scored as in T45. Also the published s0 = 2 (no fitting) and s0 = 0 (no self-interaction)."""
import json, numpy as np
res = json.load(open('res46_sidm.json'))
S0S = sorted({float(k.split('|')[1]) for k in res[0][2]})
mse = {}                                       # mse[model][s0][pop][host]
for pop, n, o in res:
    for key, v in o.items():
        model, s0 = key.split('|'); s0 = float(s0); er = np.array(v['er']); kinds = v['kind']
        d = mse.setdefault(model, {}).setdefault(s0, {})
        if pop == 'CLASH':
            for kk in ('CLASH_centre', 'CLASH_outer'):
                sel = np.array([x == kk for x in kinds]); d.setdefault(kk, {})[n] = float(np.mean(er[sel] ** 2))
        else: d.setdefault(pop, {})[n] = float(np.mean(er ** 2))
recs = json.load(open('res45_MOGF0IESM_primary.json'))['records']
POPS = ['disk', 'dwarfs', 'ETG', 'MaNGA', 'CLASH_centre', 'CLASH_outer']
alld = sorted(mse['std'][0.0]['disk'])
out = {}
for model in ('std', 'formula'):
    pooled = {p: {} for p in POPS}; chosen = []
    for r in recs:
        test = set(r['losses']['disk']); train = [h for h in alld if h not in test]
        s_best = min(S0S, key=lambda s: np.mean([mse[model][s]['disk'][h] for h in train])); chosen.append(s_best)
        for p in POPS:
            for h in r['losses'][p]: pooled[p].setdefault(h, []).append(mse[model][s_best][p][h])
    ho = {p: round(float(np.sqrt(np.mean([np.mean(v) for v in pooled[p].values()]))), 4) for p in POPS}
    alls = {s: {p: round(float(np.sqrt(np.mean(list(mse[model][s][p].values())))), 4) for p in POPS} for s in S0S}
    out[model] = dict(held_out_s0_chosen_on_disks=ho, s0_chosen=dict(zip(*np.unique(chosen, return_counts=True))) if False else
                      {str(s): chosen.count(s) for s in S0S}, all_data=alls)
    print(model, 'held out (s0 on training disks):', ho, 's0 chosen:', out[model]['s0_chosen'])
    for s in S0S: print('   all data s0 =', s, alls[s])
json.dump(out, open('res46_sidm_scores.json', 'w'), indent=1)
