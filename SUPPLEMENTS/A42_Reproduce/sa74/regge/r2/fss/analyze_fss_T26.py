"""T26 finite-size view of the curvature-squared ensemble (a = 0.005): average curvature R(k, L), its fluctuation
chi = 4 var(sum delta A) / <sum V> (Hamber's curvature fluctuation), and the zero crossing k0(L) by linear interpolation.
Pools every recorded chain with a = 0.005 at each (L, k): stage A/B/B2 audit chains, the 3^4 and 4^4 scans and the FSS stage.
Equilibration sweeps are dropped. Descriptive only: start independence was tested on 3^4 and 4^4 alone."""
import json, re, glob, sys
from pathlib import Path
import numpy as np
R2 = Path('/home/claude/sa74/regge/r2'); AUD = R2 / 'audit_a005'
sys.path.insert(0, str(AUD)); from analyze_r2 import batches, ess

def load(fn):
    if '/mac_pull/' in fn:
        a = np.loadtxt(fn, comments='#')
        return a[a[:, 1] == 1] if a.ndim == 2 and a.shape[1] >= 28 else None
    if '# done' not in open(fn).read()[-20000:]: return None   # skip chains still running
    a = np.loadtxt(fn)
    if a.ndim != 2 or a.shape[1] < 28: return None
    return a[a[:, 1] == 1]

runs = []
for fn in glob.glob(str(AUD / 'runs' / '*.dat')):
    name = Path(fn).stem
    if name.startswith('prepare'): continue
    m = re.match(r'(?:FSS_)?L(\d+)_k([0-9.]+)_(flat|rough)', name)
    if not m: continue
    runs.append((int(m.group(1)), float(m.group(2)), m.group(3), name, fn))
MAC = Path('/home/claude/sa74/NIGHT/mac_pull')   # Mac chains (plans M1b/N2): TAG.dat holds rows only, all parts appended
for fn in sorted(MAC.glob('N_L*_k*.dat')):
    m = re.match(r'N_L(\d+)_k([0-9.]+)$', fn.stem)
    if m: runs.append((int(m.group(1)), float(m.group(2)), 'flat', 'mac_' + fn.stem, str(fn)))
for fn in glob.glob(str(R2 / 'scan' / 'L3_k*.dat')) + glob.glob(str(R2 / 'scan4' / 'L4_k*.dat')):
    m = re.match(r'L(\d)_k([0-9.]+)_(flat|rough)', Path(fn).stem)
    runs.append((int(m.group(1)), float(m.group(2)), m.group(3), 'scan_' + Path(fn).stem, fn))

MINROWS = int(__import__('os').environ.get('MINROWS', '1000'))
table = {}
for L, k, st, name, fn in runs:
    a = load(fn)
    if a is None or len(a) < MINROWS: continue
    x = a[:, 4]; da = a[:, 3]; v = a[:, 2]; mq = a[:, 5]
    RH = float(mq.mean() * 2 * da.mean() / v.mean())   # Hamber's R(k) = <l^2> <2 sum delta A> / <sum V> (Hamber 2015, eq. 14)
    chi = 4 * da.var() / v.mean()
    blocks = np.array_split(np.arange(len(x)), 5)
    chib = [4 * da[b].var() / v[b].mean() for b in blocks]
    table.setdefault((L, round(k, 3)), []).append(dict(name=name, start=st, n=len(x), R=float(x.mean()),
        se=float(max(b['se'] for b in batches(x))), ess=float(ess(x)), chi=float(chi), chi_se=float(np.std(chib, ddof=1) / np.sqrt(5)),
        R_first=float(x[:len(x) // 2].mean()), R_last=float(x[len(x) // 2:].mean()), RH=RH,
        RH_se=float(max(b['se'] for b in batches(x)) * mq.mean())))

summary = []
for (L, k), rows in sorted(table.items()):
    Rm = np.mean([r['R'] for r in rows])
    se = (np.std([r['R'] for r in rows], ddof=1) / np.sqrt(len(rows))) if len(rows) > 2 else max(r['se'] for r in rows)
    RHm = np.mean([r['RH'] for r in rows])
    RHse = (np.std([r['RH'] for r in rows], ddof=1) / np.sqrt(len(rows))) if len(rows) > 2 else max(r['RH_se'] for r in rows)
    summary.append(dict(L=L, k=k, chains=len(rows), sweeps=int(sum(r['n'] for r in rows)), R=float(Rm), R_se=float(se), RH=float(RHm), RH_se=float(RHse),
                        chi=float(np.mean([r['chi'] for r in rows])), chi_se=float(np.sqrt(sum(r['chi_se'] ** 2 for r in rows)) / len(rows)),
                        ess_min=float(min(r['ess'] for r in rows)), starts=sorted(set(r['start'] for r in rows))))

def crossing(L):
    pts = sorted((s['k'], s['R'], s['R_se']) for s in summary if s['L'] == L)
    for (k1, r1, e1), (k2, r2, e2) in zip(pts, pts[1:]):
        if r1 < 0 <= r2:
            k0 = k1 + (k2 - k1) * (-r1) / (r2 - r1)
            # spread from moving each end by its standard error
            ks = []
            for d1 in (-e1, e1):
                for d2 in (-e2, e2):
                    a1, a2 = r1 + d1, r2 + d2
                    if a1 < 0 < a2: ks.append(k1 + (k2 - k1) * (-a1) / (a2 - a1))
            return dict(L=L, k0=float(k0), k_lo=float(min(ks)) if ks else None, k_hi=float(max(ks)) if ks else None, bracket=[k1, k2])
    return dict(L=L, k0=None, bracket=None)

cross = [crossing(L) for L in sorted(set(s['L'] for s in summary))]
out = dict(description=__doc__, summary=summary, crossings=cross, per_chain={f'{L}_{k}': v for (L, k), v in table.items()})
(R2 / 'fss' / 'fss_results.json').write_text(json.dumps(out, indent=1) + '\n')
print(f"{'L':>2} {'k':>6} {'ch':>3} {'R':>9} {'se':>7} {'chi':>8} {'chi_se':>7} {'essmin':>6} starts")
for s in summary:
    print(f"{s['L']:>2} {s['k']:>6.3f} {s['chains']:>3} {s['R']:>9.4f} {s['R_se']:>7.4f} {s['chi']:>8.3f} {s['chi_se']:>7.3f} {s['ess_min']:>6.0f} {','.join(s['starts'])}")
for c in cross: print('crossing', c)
