"""Tidal truncation of dwarf satellites at their Jacobi radius at pericentre.
Pericentres: Pace, Erkal and Li (2022) Table 3 (with the LMC; option without), for the 36 Milky Way and LMC satellites in
the sample they cover; other satellites use their present host distance (a lower bound on tides). Hosts as in the paper's
subhalo variant: NFW with M200 = 1.0e12 Msun (Milky Way, also for LMC satellites) or 1.5e12 (M31), Dutton-Maccio c at z = 0.
r_J^3 = r_p^3 m(<r_J) / ((3 - dlnM_h/dlnr) M_h(<r_p)); the satellite's dark mass beyond r_J is removed."""
import csv, json
import numpy as np
import harness_T15_hashed as H
K = H.K
PACE = {'Ant II': 'antlia_2', 'Aqu II': 'aquarius_2', 'Boo I': 'bootes_1', 'Boo II': 'bootes_2', 'Boo III': 'bootes_3',
        'CVn I': 'canes_venatici_1', 'CVn II': 'canes_venatici_2', 'Car': 'carina_1', 'Car II': 'carina_2', 'Car III': 'carina_3',
        'Col I': 'columba_1', 'CB': 'coma_berenices_1', 'Cra II': 'crater_2', 'Dra': 'draco_1', 'Dra II': 'draco_2',
        'Eri II': 'eridanus_2', 'For': 'fornax_1', 'Gru I': 'grus_1', 'Gru II': 'grus_2', 'Her': 'hercules_1', 'Hor I': 'horologium_1',
        'Hyd II': 'hydra_2', 'Hyi I': 'hydrus_1', 'Leo I': 'leo_1', 'Leo II': 'leo_2', 'Leo IV': 'leo_4', 'Leo V': 'leo_5',
        'Peg III': 'pegasus_3', 'Phx II': 'phoenix_2', 'Pis II': 'pisces_2', 'Ret II': 'reticulum_2', 'Ret III': 'reticulum_3',
        'Sgr II': 'sagittarius_2', 'Scl': 'sculptor_1', 'Seg 1': 'segue_1', 'Seg 2': 'segue_2', 'Sxt': 'sextans_1',
        'Tri II': 'triangulum_2', 'Tuc II': 'tucana_2', 'Tuc III': 'tucana_3', 'Tuc IV': 'tucana_4', 'Tuc V': 'tucana_5',
        'UMa I': 'ursa_major_1', 'UMa II': 'ursa_major_2', 'UMi': 'ursa_minor_1', 'Wil 1': 'willman_1'}
EXT = H.S / 'ext'
def pace(lmc=True):
    out = {}
    for r in csv.DictReader(open(EXT / 'pace2022_table3.csv')):
        out[PACE[r['name']]] = float(r['rperi_kpc' if lmc else 'rperi_noLMC_kpc'])
    return out
def dwarf_rows(lmc=True):
    """The 75 measured dwarfs (as harness_T16.dwarf_rows) with pericentre rp (kpc) and host mass mh (inf if no host)."""
    bt = H.PKG / 'blind_test'; rules = json.load(open(bt / 'frozen_disk_rules.json'))['frozen_on_all_disks']
    rr, src = H.V4.load_rows(str(H.PKG / 'data/LVDB/dwarf_all.csv'), str(bt))
    aa = H.V4.arrays(rr, rules, ml=2., major=False); pp = H.V3.predict(rr, rules, ml=2., major=False)
    P = pace(lmc); num = H.V3.num; out = []
    for i, x in enumerate(rr):
        if aa['ul'][i]: continue
        h = x['host']; rp, mh = np.inf, np.inf
        if h in H.V3.MW_HOSTS: mh = 1.0e12; rp = P.get(x['key'], num(x, 'distance_gc'))
        elif h == 'm_031': mh = 1.5e12; rp = num(x, 'distance_m31')
        out.append(dict(name=x['key'], kind='dwarfs', r=pp['r'][i], gb=pp['gb'][i], go=10 ** aa['lo'][i], ms=pp['Ms'][i], mp=pp['Ms'][i],
                        re=pp['Re'][i], z=0., fac=float(pp['fac'][i]), rp=float(rp), mh=float(mh), in_pace=x['key'] in P))
    assert len(out) == 75; return out
def r_jacobi(m, c, rv, rp, mh):
    """Jacobi radius (kpc) of an NFW satellite (total M200 m, concentration c, r200 rv in kpc; dark mass (1-fb) m) at rp."""
    ch = K.c_dm14(mh, 0.); rh = (H.G * mh / (100 * H.hz(0.) ** 2)) ** (1 / 3); xh = rp * ch / rh
    Mh = mh * H.ff(xh) / H.ff(ch); gam = xh ** 2 / ((1 + xh) ** 2 * H.ff(xh))
    msub = lambda r: (1 - H.FB) * m * H.ff(r * c / rv) / H.ff(c)
    lo, hi = np.log(1e-4) * np.ones_like(m), np.log(rp) * np.ones_like(m)
    for _ in range(80):
        mid = (lo + hi) / 2; r = np.exp(mid)
        f = r ** 3 * (3 - gam) * Mh - rp ** 3 * msub(r)
        lo = np.where(f < 0, mid, lo); hi = np.where(f < 0, hi, mid)
    return np.exp((lo + hi) / 2)
