"""SA73 T59: parse the three downloaded sources into one JSON. Run with python3 -I parse59.py <bundle dir> <out.json>.
Sources (bundle59/SOURCES.txt): Iorio et al. (2017) rotation-curve tables, Dale et al. (2009) LVL Tables 1 and 2 (arXiv
source), Read et al. (2017) Table 2 (arXiv source). Nothing in the bundle is executed."""
import sys, os, re, json
B, OUT = sys.argv[1], sys.argv[2]
# file stem -> (display name, LVL name, SPARC name or None, Read et al. name)
MAP = {'cvidwa': ('CVnIdwA', 'UGCA292', None, 'CVnIdwA'), 'ddo101': ('DDO 101', 'UGC06900', None, 'DDO 101'),
       'ddo126': ('DDO 126', 'UGC07559', 'UGC07559', 'DDO 126'), 'ddo133': ('DDO 133', 'UGC07698', None, 'DDO 133'),
       'ddo154': ('DDO 154', 'UGC08024', 'DDO154', 'DDO 154'), 'ddo168': ('DDO 168', 'UGC08320', 'DDO168', 'DDO 168'),
       'ddo210': ('DDO 210', 'DDO210', None, 'DDO 210'), 'ddo216': ('DDO 216', 'UGC12613', None, 'DDO 216'),
       'ddo47': ('DDO 47', 'UGC03974', None, 'DDO 47'), 'ddo50': ('DDO 50', 'UGC04305', 'UGC04305', 'DDO 50'),
       'ddo52': ('DDO 52', 'UGC04426', None, 'DDO 52'), 'ddo53': ('DDO 53', 'UGC04459', None, 'DDO 53'),
       'ddo87': ('DDO 87', 'UGC05918', 'UGC05918', 'DDO 87'), 'ngc1569': ('NGC 1569', 'NGC1569', None, 'NGC 1569'),
       'ngc2366': ('NGC 2366', 'NGC2366', 'NGC2366', 'NGC 2366'), 'ugc8508': ('UGC 8508', 'UGC08508', None, 'UGC 8508'),
       'wlm': ('WLM', 'WLM', 'UGCA444', 'WLM')}
def lvl_name(s): return re.sub(r'\\[a-z]', '', s).strip()
f36, dlvl = {}, {}
for line in open(os.path.join(B, 'dale2009/table2.tex')):
    p = line.split('&')
    if len(p) < 11: continue
    v = p[4].strip(); m = re.match(r'^([0-9.]+)\\p([0-9.]+)E([+-][0-9]+)$', v)
    f36[lvl_name(p[0])] = [float(m.group(1)) * 10 ** int(m.group(3)), float(m.group(2)) * 10 ** int(m.group(3))] if m else None
for line in open(os.path.join(B, 'dale2009/table1.tex')):
    p = line.split('&')
    if len(p) < 11: continue
    try: dlvl[lvl_name(p[0])] = float(p[3])
    except ValueError: pass
rd = {}; flag = 'clean'
for line in open(os.path.join(B, 'read2017/readetal_table2_lines225_300.tex')):
    if 'Rogues' in line:
        flag = 'i-rogue' if '$i$' in line else 'disequilibrium' if 'Disequilibrium' in line else 'distance'; continue
    p = [x.strip() for x in line.split('&')]
    if len(p) < 14 or not p[0]: continue
    nm = re.sub(r'\\rowcolor\{gray!40\}|\$\^\{?\\bullet\}?\$', '', p[0]).strip()
    num = lambda s: float(re.match(r'^\$?([0-9.]+)', s).group(1))
    try: rd[nm] = dict(D=num(p[3]), Ms_1e7=num(p[4]), Mgas_1e7=num(p[5]), Rstar=num(p[6]), Rgas=num(p[7]), flag=flag)
    except AttributeError: pass
out = {}
for stem, (nm, lv, sp, rn) in sorted(MAP.items()):
    L = open(os.path.join(B, 'iorio2017', stem + '_onlinetab.txt')).read().split('\n')
    D = float(next(l for l in L if 'Distance' in l).split(':')[1].split()[0]); inc = float(next(l for l in L if 'Inclination' in l).split(':')[1].split()[0])
    t = [[float(x) for x in l.split()] for l in L if l.strip() and not l.startswith('#')]
    assert all(len(r) == 12 for r in t), stem
    out[nm] = dict(stem=stem, lvl=lv, sparc=sp, D=D, inc=inc, R=[r[1] for r in t], Vc=[r[6] for r in t], eVc=[r[7] for r in t], Vr=[r[2] for r in t],
                   Sd=[r[10] for r in t], eSd=[r[11] for r in t], f36=f36.get(lv), D_lvl=dlvl.get(lv), read=rd.get(rn))
json.dump(out, open(OUT, 'w'), indent=1)
for nm, g in out.items():
    print('%-9s D %5.2f inc %4.1f n %2d f36 %s Dlvl %s read %s' % (nm, g['D'], g['inc'], len(g['R']), g['f36'], g['D_lvl'], g['read']))
print(len(f36), 'LVL rows;', len(rd), 'Read rows')
