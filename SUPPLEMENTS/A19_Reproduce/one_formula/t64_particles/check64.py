"""T64: every number quoted in the paper's new paragraph and in the round-2 report, recomputed from the results files.
Usage: python3 check64.py <dmfit folder> <paper .src> <report .src>"""
import sys, json, os, numpy as np
D, PAPER, REPORT = sys.argv[1], open(sys.argv[2]).read(), open(sys.argv[3]).read()
def J(p): return json.load(open(p))
ho = J('res64_heldout.json'); sh = J('res64_sidm_heldout.json'); ou = J('res64_outside.json'); de = J('res64_describe.json'); ck = J('res64_sidm_check.json'); h61 = J(os.path.join(D, 't61_core', 'res61_heldout.json'))
S = {k: J(os.path.join(D, f))['summary'] for k, f in dict(F='res64_MOGF0IESBF_primary.json', W='res64_MOGF0IESBW_primary.json', Ws='res64_MOGF0IESBW_strict.json', T48='res48_MOGF0IESB_primary.json', T61='res61_MOGF0IESBC_primary.json').items()}
sc = {k: list(v.values())[0] for k, v in S.items()}
pct = lambda a, b: 100 * (a / b - 1)
facts = [
 ('fuzzy bound-hit partitions', sc['F']['bound_hit_partitions'], '60'),
 ('warm lhm range -> 1e9 suns', [round(10 ** min(ho['fits']['lhm']) / 1e9, 2), round(10 ** max(ho['fits']['lhm']) / 1e9, 2)], '1.1 to 2.0'),
 ('warm keV range', [round(x, 2) for x in de['warm_relic']['keV_partitions']], '2.0 to 2.4'), ('warm keV all data', round(de['warm_relic']['keV_all_data'], 2), '2.2 keV'),
 ('warm all points', [ho['all']['T48'], ho['all']['W'], ho['all']['MOND'], ho['all']['T61']], '0.218 to 0.210 / 0.200 / 0.193'), ('warm all vs T48', ho['all']['W_minus_T48'], '-4 percent (-0.017 to -0.0003)'),
 ('warm dropped vs T48', ho['dropped']['W_minus_T48'], '-5 percent (-0.031 to -0.013)'), ('warm kept pct', ho['kept']['W_minus_T48']['pct'], '+0.5'),
 ('warm dwarfs', [sc['T48']['scores']['dwarfs'], sc['W']['scores']['dwarfs'], round(pct(sc['W']['scores']['dwarfs'], sc['T48']['scores']['dwarfs']), 1)], '0.370 to 0.769, +108'),
 ('warm lensing primary/strict', [sc['W']['lensing_chi2'], list(S['Ws'].values())[0]['lensing_chi2'], sc['T48']['lensing_chi2']], 'not quoted'),
 ('core vs T48 dropped', h61['dropped']['T61_minus_T48'], '-14 percent (-0.077 to -0.039)'), ('core vs T48 all', h61['all']['T61_minus_T48'], '-12 percent (-0.038 to -0.014)'), ('core kept pct', round(pct(h61['kept']['T61'], h61['kept']['T48']), 2), '-0.1'),
 ('SIDM chosen', sh['s0_counts'], '1.0 in 13 of 15'), ('SIDM dropped', sh['disk_dropped']['S_minus_S0'], '-14 percent (-0.082 to -0.033)'), ('SIDM all', sh['disk_all']['S_minus_S0'], '-9 percent (-0.034 to -0.005)'),
 ('SIDM kept', sh['disk_kept']['S_minus_S0'], '+5.7 percent (+0.003 to +0.015)'), ('SIDM dwarfs', sh['dwarfs']['S_minus_S0'], '+25 percent (+0.036 to +0.150)'), ('SIDM CLASH centre', sh['CLASH_centre']['S_minus_S0'], '+26 percent (-0.004 to +0.056)'),
 ('SIDM ETG, MaNGA pct', [sh['ETG']['S_minus_S0']['pct'], sh['MaNGA']['S_minus_S0']['pct']], '3 percent each'), ('SIDM vs MOND dropped', sh['disk_dropped']['S_minus_MOND'], '(-0.042 to -0.003)'),
 ('SIDM outside at s0 = 1', [sh['all_data']['1.0']['iorio_all'], sh['all_data']['1.0']['sparc_outside_all']], '0.322 and 0.401'), ('SIDM all-data choice', sh['s0_all_data'], '1'),
 ('outside Iorio', {k: ou['d_iorio_outside_all_points'][k] for k in ('T48', 'T61', 'W', 'MOND')}, '0.411 0.308 0.323 0.414'), ('outside SPARC', {k: ou['e_sparc_outside_benchmark_all_points'][k] for k in ('T48', 'T61', 'W', 'MOND')}, '0.435 0.392 0.407 0.415'),
 ('MOND vs T48 pct dropped, all, kept', [round(pct(ho[s]['MOND'], ho[s]['T48']), 1) for s in ('dropped', 'all', 'kept')], '-9, -8, +0.8'), ('MOND dwarfs pct', round(pct(0.633, sc['T48']['scores']['dwarfs']), 1), '+71'),
 ('fuzzy alone', {k.split()[1]: 10 ** v['best'] * 1e-22 for k, v in de['each_population_alone'].items() if 'BF' in k}, 'dwarfs 3e-20, ETG 2e-22, disks cold end'), ('fuzzy join radii best', {k: v['best_lm22'] for k, v in de['fuzzy_join_radius'].items()}, '6 for 2, 3, 4'),
 ('inner residuals', {k: {kk: (vv['mean'], vv['sem']) for kk, vv in v.items()} for k, v in de['inner_residuals_T48'].items()}, '-0.26+-0.10; -0.14+-0.07; -0.02+-0.05'),
 ('solver failures at 1.0', ck['failures'].get('1.0'), '9 ETG + 3 disks = 12 of 344'), ('r1 dwarfs median', ck['r1_kpc_at_s0_1']['dwarfs'][1], '0.4 kpc')]
for name, val, quoted in facts: print('%-36s %s   | quoted: %s' % (name, json.dumps(val), quoted))
need = {'paper': ['0.218 to 0.210', '0.370 to 0.769', '13 of 15', '14 percent (-0.082 to -0.033)', '9 percent (-0.034 to -0.005)', '5.7 percent (+0.003 to +0.015)', '25 percent (+0.036 to +0.150)', '26 percent (-0.004 to +0.056)', '0.322 and 0.401', '0.308 and 0.392', '1.1 to 2.0', '2.0 to 2.4 keV', '6.5 keV'],
        'report': ['0.218 to 0.210', '0.370 to 0.769', '13 of 15', '-14 percent (-0.082 to -0.033)', '-9 percent (-0.034 to -0.005)', '+5.7 percent (+0.003 to +0.015)', '+25 percent (+0.036 to +0.150)', '(-0.004 to +0.056)', '0.322 | 0.401', '0.308 | 0.392', '0.323 | 0.407', '0.411 | 0.435', '0.414 | 0.415', '2.0 to 2.4 keV', '40d5bbddf0cca4c61e85109b4f62b41858827ef7bc5dbb10cdc747df3487157a', '6cb4c1d89391fdddb2aa398f763762afcc987be640fa594046fdf6a561549320']}
for lab, txt in (('paper', PAPER), ('report', REPORT)):
    miss = [x for x in need[lab] if x not in txt]; print(lab, 'strings missing:', miss)
print('protocol hash file:', open(os.path.join(D, 'PROTOCOL_T64.sha256')).read().split()[0][:16], '| forecast:', open(os.path.join(D, 'FORECAST_T64.sha256')).read().split()[0][:16])
