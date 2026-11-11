#!/usr/bin/env python3
"""Check the claims of A42, "Searching for a Short Discrete Theory of Everything: Rewriting Rules and Supplied-Lattice Tests"
(SA74 T32, November 10, 2026), from its code and input data. The package holds no outputs: every output is regenerated here
and compared with the number the paper states. Python 3.9 or later with NumPy (SciPy for some supplements); a C compiler.

    python3 verify.py                 files, engines, and the fast generators with the paper's headline numbers (about 15 minutes)
    python3 verify.py --regge         adds Section 8.10: analyzes the regenerated chains and compares every quoted number with the
                                      paper's, within the combined statistical errors (PAPER_NUMBERS.json, regge_checks.py)
                                      (first regenerate the chains: python3 regenerate.py cloud all; scan; run --all)
    python3 verify.py --only NAME     one step: files, engine, generators or regge

The steps work in _work/, a copy in which the scripts' original absolute paths (/home/claude/sa74 and /mnt/user-data/outputs/night,
where they ran) are replaced by this folder's path. The report goes to VERIFY_REPORT.txt."""
import argparse, hashlib, json, math, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WORK = ROOT / '_work'
MAP = {'/home/claude/sa74': 'sa74', '/mnt/user-data/outputs/night': 'ext/night'}
OUT = []; ROWS = []


def say(s=''):
    print(s, flush=True); OUT.append(s)


def record(step, name, ok, detail=''):
    st = {True: 'PASS', False: 'FAIL', None: 'SKIP'}[ok]
    ROWS.append((step, st, name, detail)); say(f'  [{st}] {name}' + (f': {detail}' if detail else ''))


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def run(cmd, cwd, env=None, timeout=3600):
    e = dict(os.environ); e.update(env or {})
    try:
        return subprocess.run(cmd, cwd=str(cwd), env=e, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as x:
        class R: returncode = -9; stdout = x.stdout or ''; stderr = f'timeout after {timeout} s'
        return R()


def tail(r, n=300):
    return ((r.stderr or '') + (r.stdout or ''))[-n:].replace('\n', ' | ')


def make_work():
    """Copy sa74/, ext/ and supplements/ to _work/, replacing the original absolute paths in every .py file; large chain files are linked."""
    if WORK.exists(): shutil.rmtree(WORK)
    for top in ('sa74', 'ext', 'supplements'):
        src = ROOT / top
        if not src.exists(): continue
        for p in src.rglob('*'):
            if not p.is_file(): continue
            q = WORK / p.relative_to(ROOT); q.parent.mkdir(parents=True, exist_ok=True)
            if p.suffix == '.py':
                t = p.read_text(encoding='utf-8')
                for old, new in MAP.items(): t = t.replace(old, (WORK / new).as_posix())
                q.write_text(t, encoding='utf-8')
            elif p.suffix in ('.dat', '.cfg', '.jsonl', '.log') and p.stat().st_size > 1_000_000:
                try: q.symlink_to(p)
                except OSError: shutil.copy2(p, q)
            else:
                shutil.copy2(p, q)
    (WORK / 'sa74' / 'out').mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- 1. files
def step_files():
    say('\n1. Files: SHA-256 of every file against MANIFEST.sha256')
    n = bad = 0
    for line in (ROOT / 'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith('#'): continue
        h, rel = line.split('\t')
        p = ROOT / rel
        if not p.exists(): bad += 1; say(f'    missing: {rel}'); continue
        n += 1
        if sha256(p) != h: bad += 1; say(f'    changed: {rel}')
    record('files', f'{n} files match their SHA-256', bad == 0, f'{bad} missing or changed' if bad else '')


# ---------------------------------------------------------------- 2. engines
def step_engine():
    say('\n2. Engines: build the Regge engines and check them')
    cc = os.environ.get('CC') or shutil.which('cc') or shutil.which('gcc') or shutil.which('clang')
    if not cc: record('engine', 'C compiler', None, 'no cc, gcc or clang found; set CC'); return
    d = WORK / 'engine'; d.mkdir(parents=True, exist_ok=True)
    P = WORK / 'sa74' / 'regge' / 'par'
    builds = [run([cc, '-O2', '-o', 'par', str(P / 'regge_r2_par.c'), '-lm', '-lpthread'], d), run([cc, '-O2', '-o', 'seq', str(P / 'regge_r2_big_x.c'), '-lm'], d),
              run([cc, '-O2', '-o', 'cloud', str(WORK / 'sa74' / 'regge' / 'r2' / 'regge_r2.c'), '-lm'], d), run([cc, '-O2', '-o', 'big', str(WORK / 'ext' / 'night' / 'regge_r2_big.c'), '-lm'], d)]
    ok = all(b.returncode == 0 for b in builds)
    record('engine', f'build the four engines (regge_r2_par.c, regge_r2_big_x.c, regge_r2.c, regge_r2_big.c) with {Path(cc).name} -O2', ok, '' if ok else ' '.join(tail(b) for b in builds))
    if not ok: return
    for exe in ('par', 'cloud'):
        r = run(['./' + exe, 'selftest', '4'], d, timeout=600); s = r.stdout
        def num(tag, key):
            for l in s.splitlines():
                if l.startswith(tag) and key in l:
                    try: return float(l.split(key)[1].split()[0])
                    except ValueError: return None
            return None
        checks = [r.returncode == 0, 'COUNTS_OK 1' in s, (num('FLAT', 'max|delta_h|=') or 1) < 1e-12, (num('REGULAR', 'max|diff|=') or 1) < 1e-12,
                  (num('RANDOM', 'max Schlaefli residual (rel)=') or 1) < 1e-4, (num('INCREMENTAL', 'max|dthsum|=') or 1) < 1e-9]
        record('engine', f'{exe} self-test on 4⁴: cell counts, flat deficits below 1e-12, the regular dihedral angle arccos(1/4), the Schlaefli identity, incremental against full recomputation',
               all(checks), '' if all(checks) else s[-400:])
    env = {'R2A': '0.005'}
    rows = lambda f: [l for l in open(d / f) if not l.startswith('#')]
    r1 = run(['./par', 'run', '3', '0.2475', '1', '3', '3', '4242', 'new.dat', '0.5', '-', 'new.cfg', '10', '0', '0'], d, env)
    r2 = run(['./seq', 'run', '3', '0.2475', '1', '3', '3', '4242', 'old.dat', '0.5', '-', 'old.cfg', '10', '0', '0'], d, env)
    record('engine', "without threads the multi-core engine gives the sequential engine's run (3⁴, a = 0.005, seed 4242: same rows, same final configuration)",
           r1.returncode == 0 and r2.returncode == 0 and (d / 'new.cfg').read_bytes() == (d / 'old.cfg').read_bytes() and rows('new.dat') == rows('old.dat'))
    r1 = run(['./par', 'run', '12', '0.248', '1', '1', '1', '777', 'p1.dat', '0.5', '-', 'p1.cfg', '10', '0', '0'], d, dict(env, PAR='1'), timeout=900)
    r2 = run(['./par', 'run', '12', '0.248', '1', '1', '1', '777', 'p2.dat', '0.5', '-', 'p2.cfg', '10', '0', '0'], d, dict(env, PAR='2'), timeout=900)
    record('engine', '1 and 2 threads give the same run (12⁴ block mode, k = 0.248: same final configuration)',
           r1.returncode == 0 and r2.returncode == 0 and (d / 'p1.cfg').read_bytes() == (d / 'p2.cfg').read_bytes())
    r1 = run(['./par', 'run', '8', '0.248', '1', '5', '5', '777', 'q1.dat', '0.5', '-', '-', '10', '0', '0'], d, dict(env, PAR='1'))
    r2 = run(['./par', 'run', '8', '0.248', '1', '5', '5', '777', 'q2.dat', '0.5', '-', '-', '10', '0', '0'], d, dict(env, PAR='2'))
    cut = lambda f: [[c for i, c in enumerate(l.split()) if i != 15] for l in rows(f)]
    record('engine', '1 and 2 threads give the same rows (8⁴ colored mode; all columns but the round-off diagnostic drift)', r1.returncode == 0 and r2.returncode == 0 and cut('q1.dat') == cut('q2.dat'))


# ---------------------------------------------------------------- 3. the fast generators and the paper's headline numbers
def getpath(obj, path):
    for k in path: obj = obj[k]
    return obj


def step_generators():
    say('\n3. Generators: rerun the fast scripts of the supplements from code and inputs, then compare their numbers with the paper')
    for job in json.load(open(ROOT / 'GENERATOR_JOBS.json')):
        folder, argv, tmo = job['folder'], job['argv'], job.get('timeout', 900)
        w = WORK / 'supplements' / folder
        for frm, to in job.get('prep_copy', {}).items(): shutil.copy2(w / frm, w / to)
        t0 = time.time(); r = run([sys.executable] + argv, w, job.get('env'), timeout=tmo)
        name = f'{folder}: {" ".join(argv[:4])}{" ..." if len(argv) > 4 else ""} ({time.time() - t0:.0f} s)'
        if r.returncode != 0: record('generators', name, False, tail(r)); continue
        if not job.get('claims'): record('generators', name + ' runs from code and inputs', True); continue
        for c in job['claims']:
            try:
                v = json.load(open(w / c['file']))
                v = len(v) if c.get('len') else getpath(v, c['path'])
                ok = abs(v - c['paper']) <= c.get('tol', 0)
            except Exception as x:
                v, ok = repr(x)[:80], False
            record('generators', f"{name}: {c['text']}", ok, f"regenerated {v}, paper {c['paper']}")


# ---------------------------------------------------------------- 4. Section 8.10 from regenerated chains
# The paper's numbers and their standard errors are in PAPER_NUMBERS.json; regge_checks.py computes the same numbers from the regenerated
# result files and gives the rule of comparison (z values against the combined standard errors; no hand-set ranges).
ANALYSES = [('regge/r2/night', 'multihist.py', dict(MINROWS='800', A1INC='1', LMAX='12', KPAD='0.002', SUFFIX='_T30')),
            ('regge/r2/night', 'multihist.py', dict(MINROWS='800', A1INC='1', LMAX='12', KPAD='0.002', SUFFIX='_T30', SKIPFRAC='0.5')),
            ('regge/r2/night', 'analyze_a1.py', dict(MHFILE='multihist_T30.json', OUTSUFFIX='_T37')),
            ('regge/r2/night', 'analyze_au.py', {}), ('regge/r2/fss', 'analyze_fss.py', dict(MINROWS='800'))]


def needed_chains():
    sys.path.insert(0, str(ROOT)); import regenerate as G
    jobs = G.plan_lines(); tags = [t for t in jobs if not t.startswith(('Q_L16_', 'Q_P8'))]
    cloud = []
    for p in G.CLOUD_ORDER:
        cloud += [j['tag'] for j in json.load(open(ROOT / 'sa74' / 'regge' / 'r2' / 'audit_a005' / p))['main'] if j['tag'] not in G.NOT_RUN]
    miss = [t for t in tags if not G.done(ROOT / 'sa74' / 'NIGHT' / 'mac_pull' / f'{t}.dat', jobs[t]['therm'] + jobs[t]['meas'])]
    miss += [t for t in cloud if not G.done(ROOT / 'sa74' / 'regge' / 'r2' / 'audit_a005' / 'runs' / f'{t}.dat')]
    for sub, L in (('scan', 3), ('scan4', 4)):   # the scans too: the 4^4 scan chain at k = 0.26 enters the reweighting (regenerate.py scan)
        d = ROOT / 'sa74' / 'regge' / 'r2' / sub
        for line in open(d / 'jobs.txt'):
            f = line.split()
            if f and not G.done(d / f'L{L}_k{f[0]}_{f[1]}.dat'): miss.append(f'{sub}/L{L}_k{f[0]}_{f[1]}')
    return miss


def step_regge():
    say('\n4. Section 8.10: analyze the regenerated chains and compare each quoted number with the paper')
    miss = needed_chains()
    if miss:
        record('regge', 'regenerated chains', None, f'{len(miss)} chains missing (for example {", ".join(miss[:4])}); run regenerate.py cloud all, scan and run --all first'); return
    for folder, script, env in ANALYSES:
        r = run([sys.executable, script], WORK / 'sa74' / folder, env, timeout=7200)
        record('regge', f'{" ".join(f"{k}={v}" for k, v in env.items())} {script}'.strip(), r.returncode == 0, '' if r.returncode == 0 else tail(r))
        if r.returncode: return
    sys.path.insert(0, str(ROOT)); import regge_checks as RC
    P = json.load(open(ROOT / 'PAPER_NUMBERS.json')); G = RC.numbers(WORK / 'sa74' / 'regge' / 'r2' / 'night', WORK / 'sa74' / 'regge' / 'r2' / 'fss')
    g = lambda L: f'{L}⁴'
    pm = lambda d, n=2: f"{d['v']:.{n}f}({d['se']:.{n}f})"

    def compare(title, new, paper, show, log=False):
        """One check: every number in paper against its regenerated value; passes when the largest |z| is within zmax(number of values)."""
        z = RC.zlist(new, paper, log); worst = max(abs(v) for v in z.values()); bound = RC.zmax(len(z))
        record('regge', title, worst <= bound, f'{show}; z = {", ".join(f"{v:+.1f}" for v in z.values())}; largest |z| {worst:.1f}, allowed {bound:.2f}')

    for k in RC.KS:
        kk = f'k = {float(k):.3f}'; pp, gp = P['power'][k], G['power'][k]
        compare(f'{kk}: the fluctuation on each grid ({", ".join(g(L) for L in P["fixed_k"][k])}; paper {min(v["v"] for v in P["fixed_k"][k].values()):.1f} to {max(v["v"] for v in P["fixed_k"][k].values()):.1f})',
                G['fixed_k'][k], P['fixed_k'][k], 'regenerated ' + ', '.join(f"{v['v']:.1f}" for v in G['fixed_k'][k].values()))
        compare(f'{kk}: fitted power of the fluctuation in L', {'p': gp}, {'p': pp}, f'regenerated {pm(gp)}, paper {pm(pp)}')
        record('regge', f'{kk}: no critical growth (a critical point needs a power of 0.94)', gp['v'] + 3 * gp['se'] < 0.94,
               f"regenerated power plus three standard errors {gp['v'] + 3 * gp['se']:+.2f}")
    k = '0.2440'; hp = G['hamber_power']; gap = lambda d: (hp - d[k]['v']) / d[k]['se']
    compare("k = 0.244: fitted power after the low-bias correction of the values from 16⁴ on", {'p': G['power_lowbias'][k]}, {'p': P['power_lowbias'][k]},
            f"regenerated {pm(G['power_lowbias'][k])}, paper {pm(P['power_lowbias'][k])}")
    record('regge', f"Hamber's fitted curve (a power of {hp:.2f} from L = 4 to 32) lies above the power found at his k_c", gap(G['power']) > 0 and gap(G['power_lowbias']) > 0,
           f"regenerated {gap(G['power']):.1f} standard errors above, {gap(G['power_lowbias']):.1f} after the low-bias correction (paper: {gap(P['power']):.1f} and {gap(P['power_lowbias']):.1f}; "
           'the size of the gap follows from the two powers checked above)')
    compare(f'largest fluctuation near the crossing on each grid ({", ".join(g(L) for L in P["peak"])})', G['peak'], P['peak'],
            'regenerated ' + ', '.join(f"{v['v']:.1f}" for v in G['peak'].values()) + '; paper ' + ', '.join(f"{v['v']:.1f}" for v in P['peak'].values()))
    compare('largest fluctuation: fitted power over 4⁴ to 12⁴', {'p': G['peak_power']}, {'p': P['peak_power']}, f"regenerated {pm(G['peak_power'])}, paper {pm(P['peak_power'])}")
    compare('largest fluctuation: fitted power from the second halves of the chains', {'p': G['peak_power_second_halves']}, {'p': P['peak_power_second_halves']},
            f"regenerated {pm(G['peak_power_second_halves'])}, paper {pm(P['peak_power_second_halves'])}")
    cv = [v['v'] for v in G['crossing'].values()]
    compare(f'zero crossing of the average curvature on each grid ({", ".join(g(L) for L in P["crossing"])}; paper: between 0.2516 and 0.2541)', G['crossing'], P['crossing'],
            f'regenerated {min(cv):.4f} to {max(cv):.4f}: ' + ', '.join(f"{v['v']:.4f}({v['se']:.4f})" for v in G['crossing'].values()))
    zs = list(G['start_check_z'].values())
    record('regge', "16⁴: chains started from another coupling's configuration agree with flat starts (paper: within about one standard error)",
           max(abs(x) for x in zs) <= RC.zmax(len(zs)), 'regenerated differences in standard errors: ' + ', '.join(f'{x:+.1f}' for x in zs) + f'; allowed {RC.zmax(len(zs)):.2f}')
    compare(f'autocorrelation time near the crossing, median over chains on each grid ({", ".join(g(L) for L in P["tau"])}; paper: 120 to 240 sweeps up to 20⁴; 385 on 24⁴ from three chains, 98 on 32⁴ from one)',
            G['tau'], P['tau'], 'regenerated ' + ', '.join(f"{v['median']:.0f}" for v in G['tau'].values()) + ' sweeps; paper ' + ', '.join(f"{v['median']:.0f}" for v in P['tau'].values()), log=True)
    compare(f'a = 0.02: zero crossing on {", ".join(g(L) for L in P["au_crossing"])} (paper: 0.479, 0.467, 0.468, 0.467)', G['au_crossing'], P['au_crossing'],
            'regenerated ' + ', '.join(f"{v['v']:.3f}({v['se']:.3f})" for v in G['au_crossing'].values()))
    far = {L: (v['v'] - 0.4244) / v['se'] for L, v in G['au_crossing'].items() if int(L) >= 6}
    record('regge', "a = 0.02: no drift toward Hamber's k_c = 0.4244 beyond 6⁴", all(x > 3 for x in far.values()),
           'the crossings lie above it by ' + ', '.join(f'{x:.0f}' for x in far.values()) + f' standard errors on {", ".join(g(L) for L in far)}')
    compare(f'a = 0.02: the fluctuation at k = 0.46 on {", ".join(g(L) for L in P["au_chi_046"])} (paper: 0.61, 0.55, 0.55, 0.62, 0.56)', G['au_chi_046'], P['au_chi_046'],
            'regenerated ' + ', '.join(f"{v['v']:.2f}" for v in G['au_chi_046'].values()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--regge', action='store_true'); ap.add_argument('--only', choices=['files', 'engine', 'generators', 'regge'])
    ap.add_argument('--keep', action='store_true', help='keep _work/')
    a = ap.parse_args()
    steps = [a.only] if a.only else ['files', 'engine', 'generators'] + (['regge'] if a.regge else [])
    say(f'A42 verification, {time.strftime("%Y-%m-%d %H:%M")}, Python {sys.version.split()[0]}, in {ROOT}')
    if 'files' in steps: step_files()
    if any(s in steps for s in ('engine', 'generators', 'regge')): make_work()
    if 'engine' in steps: step_engine()
    if 'generators' in steps: step_generators()
    if 'regge' in steps: step_regge()
    n = {s: sum(1 for r in ROWS if r[1] == s) for s in ('PASS', 'FAIL', 'SKIP')}
    say(f'\nSummary: {n["PASS"]} passed, {n["FAIL"]} failed, {n["SKIP"]} skipped.')
    (ROOT / 'VERIFY_REPORT.txt').write_text('\n'.join(OUT) + '\n', encoding='utf-8')
    if not a.keep and WORK.exists(): shutil.rmtree(WORK, ignore_errors=True)
    sys.exit(1 if n['FAIL'] else 0)


if __name__ == '__main__':
    main()
