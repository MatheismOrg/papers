"""Portable independent execution of unchanged SA74 Regge engine. Python stdlib only."""
from pathlib import Path
import subprocess,json,time,hashlib,platform,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
D=Path(__file__).resolve().parent
SRC=D.parent/'T18_ORIGINAL/SUPPLEMENTS/A42_SA74/A42_SA74_Regge.c'
EXE=D/'regge';RUNS=D/'runs';RUNS.mkdir(exist_ok=True)
def job(j):
    tag=j['tag'];load=RUNS/j['load'] if j.get('load') else '-'
    cmd=[str(EXE),'run',str(j['L']),str(j['k']),'1',str(j['thermal']),str(j['measure']),str(j['seed']),str(RUNS/(tag+'.dat')),'0.5',str(load),str(RUNS/(tag+'.cfg')),'10','0','0']
    start=time.time()
    with (RUNS/(tag+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    info=dict(job=j,command=cmd,returncode=r.returncode,wall_seconds=time.time()-start)
    (RUNS/(tag+'.execution.json')).write_text(json.dumps(info,indent=2))
    print(tag,r.returncode,round(info['wall_seconds'],1),flush=True)
    return info

def batch(js):
    with ThreadPoolExecutor(max_workers=8) as pool:
        return [f.result() for f in as_completed([pool.submit(job,j) for j in js])]

def main():
    plan=json.loads((D/'plan.json').read_text())
    subprocess.run(['gcc','-O2','-Wall','-o',str(EXE),str(SRC),'-lm'],check=True)
    for L in [3,4]:
        with (D/f'selftest_L{L}.txt').open('w') as f:subprocess.run([str(EXE),'selftest',str(L)],stdout=f,check=True)
    meta=dict(python=platform.python_version(),platform=platform.platform(),compiler=subprocess.check_output(['gcc','--version'],text=True).splitlines()[0],source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),plan_sha256=hashlib.sha256((D/'plan.json').read_bytes()).hexdigest(),workers=8)
    (D/'environment.json').write_text(json.dumps(meta,indent=2))
    prep=batch(plan['preparation'])
    # Never substitute a selected successful start for a failed preparation.
    if any(x['returncode'] for x in prep):raise RuntimeError('Preparation failed; inspect saved logs. No replacement starts selected.')
    results=batch(plan['main'])
    (D/'executions.json').write_text(json.dumps(prep+results,indent=2))
if __name__=='__main__':main()
