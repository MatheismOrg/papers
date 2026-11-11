#!/usr/bin/env python3
"""Regenerate the lattice-gravity chains of A42 Sections 8.9 and 8.10 from their plans (SA74, November 2026).
Every chain is a Markov chain fixed by its plan line: grid L, coupling k, sweeps at the coupling (THERM) and measured (MEAS), seed,
curvature-squared weight a, and the chain whose final configuration it starts from (FROM; flat start otherwise). The chains are
written where the analyses read them: sa74/NIGHT/mac_pull/ (Mac plans) and sa74/regge/r2/audit_a005/runs/ (cloud plans).

    python3 regenerate.py list                       every chain of every plan, with its plan file
    python3 regenerate.py run TAG [TAG ...]          rerun these chains; a missing chain they start from runs first
    python3 regenerate.py run --plan queue3.txt      every chain of one Mac plan (jobs_native.txt, queue.txt, queue2.txt, queue3.txt)
    python3 regenerate.py run --all                  every Mac plan, in dependency order
    python3 regenerate.py cloud PLAN [WORKERS]       one cloud stage (plan.json, plan_B.json, ..., plan_hy1.json) with its own runner
    python3 regenerate.py cloud all [WORKERS]        every cloud stage, in the order they ran
    python3 regenerate.py scan                       the 3^4 pilot and the 3^4 and 4^4 scans (COMMANDS.txt of the T25 supplement)

Options: --threads N (threads of the multi-core engine, plan A5; default 8), --scale F (shorten every Mac chain by the factor F and
write into _smoke/: a quick test of the machinery, not of the physics), --force (rerun finished chains), --full (run the eight Mac
chains that stopped early in the original night to their planned lengths; by default they run as far as they ran, as the paper used them).

Runtimes: 4^4 to 8^4 chains take minutes, 12^4 a few hours, 16^4 about ten hours and 20^4 about 1.5 hours on eight cores per chain.
Exact digits: the cloud chains repeat byte for byte with gcc -O2 on x86-64 (the T25 supplement checked this). The Mac chains ran in
600-second segments, each reseeded (seed + 1000 x segment), with clang on Apple silicon; a rerun here follows another path of the same
Markov chain, so its averages agree within their errors and its digits differ."""
import os, sys, subprocess, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SA = ROOT / 'sa74'; MAC = SA / 'NIGHT' / 'mac_pull'; RUN = MAC / 'runner'; AUD = SA / 'regge' / 'r2' / 'audit_a005'
BIN = ROOT / '_bin'
ENGINES = {'big': (ROOT / 'ext' / 'night' / 'regge_r2_big.c', []), 'seq': (SA / 'regge' / 'par' / 'regge_r2_big_x.c', []),
           'par': (SA / 'regge' / 'par' / 'regge_r2_par.c', ['-lpthread'])}
# Mac plan files, the engine and environment their runner used (night_native.sh / night_queue.sh: regge_r2_big.c in random edge order;
# night_queue2.sh: regge_r2_big_x.c with SEQ=1 LINES=1; night_queue3.sh: regge_r2_par.c with PAR threads and LINES=1)
PLANS = {'jobs_native.txt': (ROOT / 'ext' / 'night' / 'jobs_native.txt', 'big', {}),
         'queue.txt': (RUN / 'queue.txt', 'big', {}),
         'queue2.txt': (RUN / 'queue2.txt', 'seq', {'SEQ': '1', 'LINES': '1'}),
         'queue3.txt': (RUN / 'queue3.txt', 'par', {'LINES': '1'})}
# starts that are not another Mac chain: Q_F8_k0.260 is the final configuration of the cloud chain FSS_L8_k0.26_flat (plan_fss8.json)
ALIAS = {'Q_F8_k0.260': AUD / 'runs' / 'FSS_L8_k0.26_flat.cfg'}
# Mac chains that stopped before their planned length in the original night (runner 1 stopped at 05:30 on Nov 7; the 16^4 Q chains
# were stopped at 01:31 and served only as starts for plan A1): the paper's analyses used them as they ran, so regenerate.py runs them
# to these lengths (sweeps at the coupling, measured sweeps); --full runs them to their planned lengths instead
AS_RUN = {'N_L12_k0.245': (1000, 1876), 'N_L12_k0.2525': (1000, 1880), 'Q_L16_k0.2525': (167, 0), 'Q_L16_k0.2475': (168, 0), 'Q_L12_k0.250': (1000, 1766),
          'Q_L8_k0.2525_r2': (1500, 1359), 'Q_L8_k0.2475_r3': (1432, 0), 'Q_L8_k0.2525_r3': (1435, 0)}
# plan A3's forks started from their parents' configurations after about 1,280 of the parents' 2,000 sweeps at k = 0.244
FORKS = {'A3_L16_k0.244_f1': ('A2_L16_k0.244_r1', 1280), 'A3_L16_k0.244_f2': ('A2_L16_k0.244_r2', 1280)}
CLOUD_ORDER = ['plan.json', 'plan_B.json', 'plan_B2.json', 'plan_fss.json', 'plan_fss2.json', 'plan_fss3.json', 'plan_fss4.json', 'plan_fss8.json',
               'plan_as.json', 'plan_as2.json', 'plan_as3.json', 'plan_as4.json', 'plan_as5.json', 'plan_au1.json', 'plan_au2.json', 'plan_hy1.json']
T19 = ROOT / 'supplements' / 'A42_Regge_R2_T25' / 'regge_audited_T19.c'
# FSS3's last 16 chains did not run in the original (the cloud runner stopped at its two-hour limit, NIGHT_LOG.md 06:22, Nov 7); the
# paper's analyses had none of them, so 'cloud' runs plan_fss3.json without them (as plan_fss3_asrun.json)
NOT_RUN = {f'FSS3_L4_k{k}_flat_r{r}' for k in ('0.2475', '0.2525') for r in (1, 2, 3)} | \
          {f'FSS3_L3_k{k}_flat_r{r}' for k, n in (('0.2450', 4), ('0.2475', 3), ('0.2525', 3)) for r in range(1, n + 1)}


def opt(name, default=None, flag=False):
    if name in sys.argv:
        i = sys.argv.index(name); sys.argv.pop(i)
        return True if flag else sys.argv.pop(i)
    return False if flag else default


THREADS = opt('--threads', '8'); SCALE = float(opt('--scale', '1')); FORCE = opt('--force', flag=True); FULL = opt('--full', flag=True)
OUTD = MAC if SCALE == 1 else ROOT / '_smoke'


def plan_lines():
    jobs = {}
    for name, (path, eng, env) in PLANS.items():
        for line in open(path):
            f = line.split()
            if not f or line.startswith('#'): continue
            jobs[f[0]] = dict(tag=f[0], L=int(f[1]), k=f[2], therm=int(f[3]), meas=int(f[4]), seed=int(f[5]), a=f[6],
                              frm=f[7] if len(f) > 7 else None, plan=name, eng=eng, env=env)
            if f[0] in AS_RUN and not FULL: jobs[f[0]]['therm'], jobs[f[0]]['meas'] = AS_RUN[f[0]]
    return jobs


def engine(eng):
    BIN.mkdir(exist_ok=True); exe = BIN / eng
    if not exe.exists():
        src, libs = ENGINES[eng]; cc = os.environ.get('CC') or shutil.which('cc') or shutil.which('gcc') or shutil.which('clang')
        subprocess.run([cc, '-O2', '-o', str(exe), str(src), '-lm'] + libs, check=True)
    return exe


def done(dat, rows=None):
    """A chain is complete when the engine wrote its '# done' line or, for the Mac chains (segments appended, rows only), when it has all its rows."""
    if not dat.exists(): return False
    t = dat.read_text(errors='ignore')
    return '# done' in t[-20000:] or (rows is not None and sum(1 for l in t.splitlines() if l.strip() and not l.startswith('#')) >= rows)


def run_chain(j, jobs, stack=()):
    OUTD.mkdir(parents=True, exist_ok=True)
    tag = j['tag']; dat = OUTD / f'{tag}.dat'; cfg = OUTD / f'{tag}.cfg'
    if done(dat, j['therm'] + j['meas']) and not FORCE:
        print(f'{tag}: finished already'); return
    load = '-'
    if tag in FORKS:   # the parent's start, run for the fork's head start at the parent's coupling, then the fork
        par, n = FORKS[tag]; p = jobs[par]; snap = OUTD / f'{tag}.start.cfg'
        pl = start_of(p, jobs, stack + (tag,))
        sc = lambda x: max(1, round(x * SCALE)) if x else 0
        env = dict(os.environ, R2A=p['a'], **p['env'], **({'PAR': THREADS} if p['eng'] == 'par' else {}))
        subprocess.run([str(engine(p['eng'])), 'run', str(p['L']), p['k'], '1', str(sc(n)), '0', str(p['seed']), str(OUTD / f'{tag}.start.dat'), '0.5',
                        pl, str(snap), '10', '0', '0'], check=True, env=env, stdout=subprocess.DEVNULL)
        load = str(snap)
    else:
        load = start_of(j, jobs, stack + (tag,))
    sc = lambda x: max(1, round(x * SCALE)) if x else 0
    env = dict(os.environ, R2A=j['a'], **j['env'], **({'PAR': THREADS} if j['eng'] == 'par' else {}))
    cmd = [str(engine(j['eng'])), 'run', str(j['L']), j['k'], '1', str(sc(j['therm'])), str(sc(j['meas'])), str(j['seed']), str(dat), '0.5', load, str(cfg), '10', '0', '0']
    print(f"{tag}: {j['L']}^4, k = {j['k']}, a = {j['a']}, {sc(j['therm'])} + {sc(j['meas'])} sweeps, start {'flat' if load == '-' else Path(load).name} ({j['plan']})", flush=True)
    with open(OUTD / f'{tag}.log', 'w') as lg:
        subprocess.run(cmd, check=True, env=env, stdout=lg, stderr=subprocess.STDOUT)


def start_of(j, jobs, stack):
    frm = j['frm']
    if not frm: return '-'
    if frm in ALIAS:
        if not ALIAS[frm].exists(): sys.exit(f"{j['tag']} starts from {ALIAS[frm].name}: run 'regenerate.py cloud plan_fss8.json' first")
        return str(ALIAS[frm])
    if frm not in jobs: sys.exit(f"{j['tag']}: unknown start {frm}")
    if frm in stack: sys.exit(f'start loop at {frm}')
    if not (OUTD / f'{frm}.cfg').exists() or not done(OUTD / f'{frm}.dat', jobs[frm]['therm'] + jobs[frm]['meas']): run_chain(jobs[frm], jobs, stack)
    return str(OUTD / f'{frm}.cfg')


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('list', 'run', 'cloud', 'scan'): print(__doc__); sys.exit(1)
    cmd = sys.argv[1]; jobs = plan_lines()
    if cmd == 'list':
        for j in jobs.values(): print(f"{j['plan']:16s} {j['tag']:26s} {j['L']:2d}^4 k={j['k']:7s} a={j['a']:6s} {j['therm']:5d}+{j['meas']:5d} seed {j['seed']} {'from ' + j['frm'] if j['frm'] else ''}")
        for p in CLOUD_ORDER: print(f'cloud            {p}')
        return
    if cmd == 'run':
        args = sys.argv[2:]
        if args[:1] == ['--all']: tags = list(jobs)
        elif args[:1] == ['--plan']: tags = [t for t, j in jobs.items() if j['plan'] == args[1]]
        else: tags = args
        for t in tags:
            if t not in jobs: sys.exit(f'unknown chain {t} (see: regenerate.py list)')
            run_chain(jobs[t], jobs)
        return
    if cmd == 'cloud':   # the stage's own runner: builds regge_r2.c from the T19 engine (make_r2.py), gcc -O2, runs the plan's preparations and chains
        which = sys.argv[2] if len(sys.argv) > 2 else 'all'; workers = sys.argv[3] if len(sys.argv) > 3 else '2'
        for p in (CLOUD_ORDER if which == 'all' else [which]):
            if p == 'plan_fss3.json':
                import json
                pl = json.load(open(AUD / p)); pl['main'] = [j for j in pl['main'] if j['tag'] not in NOT_RUN]
                (AUD / 'plan_fss3_asrun.json').write_text(json.dumps(pl, indent=1)); p = 'plan_fss3_asrun.json'
            subprocess.run([sys.executable, str(AUD / 'run_audit_r2.py'), str(T19), workers, p], check=True, cwd=str(AUD))
        return
    if cmd == 'scan':   # COMMANDS.txt: pilot preparation (a = 0), then the scans with R2A = 0.005
        d = SA / 'regge' / 'r2'; exe = BIN / 'cloud'
        if not exe.exists():
            BIN.mkdir(exist_ok=True)
            subprocess.run([sys.executable, str(d / 'make_r2.py'), str(T19), str(BIN / 'regge_r2.c')], check=True)
            subprocess.run(['gcc', '-O2', '-Wall', '-o', str(exe), str(BIN / 'regge_r2.c'), '-lm'], check=True)
        pil = d / 'pilot'; pil.mkdir(exist_ok=True)
        env0 = {k: v for k, v in os.environ.items() if k != 'R2A'}
        if not (pil / 'prep_L3_s0.cfg').exists():
            subprocess.run([str(exe), 'run', '3', '0.1', '1', '1000', '2000', '92030', str(pil / 'prep_L3_s0.dat'), '0.5', '-', str(pil / 'prep_L3_s0.cfg'), '10', '0', '0'], check=True, env=env0)
        for sub, L, nt, nm, seed, prep in (('scan', 3, 2000, 2000, 94300, pil / 'prep_L3_s0.cfg'), ('scan4', 4, 2000, 4000, 96400, AUD / 'runs' / 'prepare_L4_s0.cfg')):
            for line in open(d / sub / 'jobs.txt'):
                f = line.split()
                if not f: continue
                K, start = f[0], f[1]; out = d / sub / f'L{L}_k{K}_{start}'
                if done(Path(str(out) + '.dat')) and not FORCE: continue
                if start == 'rough' and not prep.exists(): sys.exit(f'{out.name} starts from {prep.name}: run the pilot or "regenerate.py cloud plan.json" first')
                print(out.name, flush=True)
                subprocess.run([str(exe), 'run', str(L), K, '1', str(nt), str(nm), str(seed), str(out) + '.dat', '0.5', str(prep) if start == 'rough' else '-', str(out) + '.cfg', '10', '0', '0'],
                               check=True, env=dict(env0, R2A='0.005'))
        return


if __name__ == '__main__':
    main()
