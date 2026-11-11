"""T25: rerun of the T19 Regge audit design with Hamber's curvature-squared term (a = 0.005).
Builds regge_r2.c from the unchanged T19 engine source by make_r2.py, then runs the T19 plan.
Python stdlib only. Usage: python3 run_audit_r2.py ENGINE_SOURCE.c [workers] [plan.json]"""
from pathlib import Path
import subprocess, json, time, hashlib, platform, sys, os
from concurrent.futures import ThreadPoolExecutor
D = Path(__file__).resolve().parent
SRC = Path(sys.argv[1]).resolve()
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
PLAN = sys.argv[3] if len(sys.argv) > 3 else 'plan.json'
R2C = D / 'regge_r2.c'; EXE = D / 'regge_r2'; RUNS = D / 'runs'; RUNS.mkdir(exist_ok=True)

def job(j):
    tag = j['tag']; load = RUNS / j['load'] if j.get('load') else '-'
    out = RUNS / (tag + '.dat')
    if out.exists() and '# done' in out.read_text()[-20000:]:
        return dict(job=j, skipped=True, returncode=0)
    cmd = [str(EXE), 'run', str(j['L']), str(j['k']), '1', str(j['thermal']), str(j['measure']), str(j['seed']),
           str(out), '0.5', str(load), str(RUNS / (tag + '.cfg')), '10', '0', '0']
    env = dict(os.environ); env.pop('R2A', None)
    if j.get('a', 0.0): env['R2A'] = repr(j['a'])
    start = time.time()
    with (RUNS / (tag + '.log')).open('w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    info = dict(job=j, command=cmd, R2A=env.get('R2A', 'unset'), returncode=r.returncode, wall_seconds=time.time() - start)
    (RUNS / (tag + '.execution.json')).write_text(json.dumps(info, indent=2))
    print(tag, r.returncode, round(info['wall_seconds'], 1), flush=True)
    return info

def batch(js):
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        return list(pool.map(job, js))

def main():
    plan = json.loads((D / PLAN).read_text())
    subprocess.run([sys.executable, str(D / 'make_r2.py'), str(SRC), str(R2C)], check=True)
    subprocess.run(['gcc', '-O2', '-Wall', '-o', str(EXE), str(R2C), '-lm'], check=True)
    for L in [3, 4]:
        for a in ['', '0.005']:
            env = dict(os.environ); env.pop('R2A', None)
            if a: env['R2A'] = a
            with (D / f'selftest_L{L}_a{a or 0}.txt').open('w') as f:
                subprocess.run([str(EXE), 'selftest', str(L)], stdout=f, check=True, env=env)
    meta = dict(python=platform.python_version(), platform=platform.platform(),
                compiler=subprocess.check_output(['gcc', '--version'], text=True).splitlines()[0],
                engine_source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),
                r2_source_sha256=hashlib.sha256(R2C.read_bytes()).hexdigest(),
                plan=PLAN, plan_sha256=hashlib.sha256((D / PLAN).read_bytes()).hexdigest(), workers=WORKERS)
    (D / f'environment_{Path(PLAN).stem}.json').write_text(json.dumps(meta, indent=2))
    prep = batch(plan['preparation'])
    if any(x['returncode'] for x in prep): raise RuntimeError('Preparation failed; no replacement starts selected.')
    # L=4 chains first: they are the long pole
    main_jobs = sorted(plan['main'], key=lambda j: -j['L'])
    results = batch(main_jobs)
    (D / f'executions_{Path(PLAN).stem}.json').write_text(json.dumps(prep + results, indent=2))
    print('ALL_DONE', PLAN, flush=True)

if __name__ == '__main__': main()
