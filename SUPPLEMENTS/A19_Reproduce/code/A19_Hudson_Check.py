"""Section 3.8: the frozen stellar rule (late types: the disks' size and HI scalings) against the CFHTLenS halo-model
masses of blue galaxies (Hudson et al. 2015, Table 2). The rule's comparable mass is its dark halo plus baryons,
(1 - f_b) M_ref + M_p, as Hudson et al. define M_h (baryons plus NFW M_200c).
Usage: python3 A19_Hudson_Check.py data blind_test"""
import sys, os, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import A19_KiDS_v9 as K
dd, bt = sys.argv[1:3]
rule = K.Rule(dd, bt); out = []
for r in csv.DictReader(open(os.path.join(dd, 'Lensing', 'hudson2015_table2_blue.csv'))):
    ms = float(r['Mstar_1e10']) * 1e10
    if ms < 1e10: continue
    m, mp = rule.mref(np.log10(ms), False, 'primary'); tot = (1 - K.FB) * m + mp; mh = float(r['Mh_1e11']) * 1e11
    out.append(dict(z=float(r['z_mean']), logMs=round(float(np.log10(ms)), 3), Mh_hudson=mh, Mh_err=float(r['Mh_err_1e11']) * 1e11,
                    rule_dark_plus_baryons=float('%.4g' % tot), dex_rule_minus_hudson=round(float(np.log10(tot / mh)), 3),
                    err_dex=round(float(float(r['Mh_err_1e11']) / float(r['Mh_1e11']) / np.log(10)), 3)))
print(json.dumps(out, indent=1))
