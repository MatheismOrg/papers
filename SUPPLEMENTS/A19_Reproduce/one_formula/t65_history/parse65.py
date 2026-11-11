"""SA73 T65: read two published tables from their arXiv LaTeX sources into one JSON. Nothing in the sources is executed.
 1. Read, Walker and Steger (2019, MNRAS 484, 1401; arXiv:1808.06634), Table 1: for 8 dwarf spheroidals and 8 dwarf
    irregulars, the star-formation truncation time t_trunc (Gyr ago) and the dark-matter density at 150 pc with its 68
    percent interval (1e8 Msun/kpc^3).
 2. Weisz et al. (2015, ApJ 804, 136; arXiv:1503.05195), Table 1: tau_90, log10 of the lookback time in years by which 90
    percent of the stars had formed, for 38 Local Group dwarfs.
Usage: python3 -I parse65.py <DM_heating.tex> <ms.tex> <out.json>"""
import sys, re, json
READ, WEISZ, OUT = sys.argv[1:4]
def num(s): return float(s.replace(',', ''))
rd = {}
for line in open(READ):
    p = [x.strip() for x in line.split('&')]
    if len(p) == 13 and p[1] in ('dSph', 'dIrr'):
        m = re.match(r'^\$([0-9.]+)_\{-([0-9.]+)\}\^\{\+([0-9.]+)\}\$$', p[10]); assert m, p[10]
        rd[p[0]] = dict(type=p[1], t_trunc=num(p[9]), rho150=num(m.group(1)), lo=num(m.group(2)), hi=num(m.group(3)), Mstar_1e6=num(re.match(r'^\$?([0-9.,]+)', p[3]).group(1)), M200_1e9=num(re.match(r'^\$?([0-9.,]+)', p[7]).group(1)))
assert len(rd) == 16 and sum(v['type'] == 'dSph' for v in rd.values()) == 8, len(rd)
wz = {}; on = False
for line in open(WEISZ):
    if '\\startdata' in line: on = True; continue
    if '\\enddata' in line:
        if wz: break
        continue
    p = [x.strip() for x in line.split('&')]
    if on and len(p) == 8:
        name = re.sub(r'\s+', ' ', re.sub(r'\{\\sc ([IVX]+)\}', r'\1', p[0]).replace('~', ' ')).strip()
        m = re.match(r'^([0-9.]+)\$\^\{\+([0-9.]+),([0-9.]+)\}_\{-([0-9.]+),([0-9.]+)\}\$', p[7]); assert m, p[7]
        wz[name] = dict(type=p[1], tau90=float(m.group(1)), err_up_total=float(m.group(3)), err_dn_total=float(m.group(5)), Mstar_1e6=float(p[3]), MHI_1e6=float(p[4]))
assert len(wz) == 38, len(wz)
json.dump(dict(read2019=rd, weisz2015=wz), open(OUT, 'w'), indent=1)
for k, v in rd.items(): print('%-9s %s t_trunc %5.2f rho150 %.2f -%.2f +%.2f' % (k, v['type'], v['t_trunc'], v['rho150'], v['lo'], v['hi']))
print(len(wz), 'Weisz rows:', ', '.join('%s %.2f' % (k, v['tau90']) for k, v in wz.items()))
