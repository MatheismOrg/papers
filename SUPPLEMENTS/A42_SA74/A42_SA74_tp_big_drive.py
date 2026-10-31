import subprocess, json, math, os
import os; B = os.environ.get('TP_BIG', './tp_big')
env = dict(os.environ, SKETCH='1', MIDALL='1')
def run(args): return json.loads(subprocess.run([B] + [str(a) for a in args], capture_output=True, text=True, check=True, env=env).stdout)
tuned = {json.loads(l)['tune']: json.loads(l)['p'] for l in open('trajectory.jsonl') if '"tune"' in l}
tuned[64000] = 0.00027221022950744544      # the Python run's tuned value
out = open('tp_big_trajectory.jsonl', 'w')
for N in (64000, 128000, 256000, 512000, 1024000, 2048000):
    if N not in tuned:
        lo, hi = math.log(8.0 / N), math.log(40.0 / N)
        for _ in range(11):
            mid = (lo + hi) / 2; r = run(['mm', N, math.exp(mid), 101, 1500])
            if r['mm'] > 4.0: lo = mid
            else: hi = mid
        tuned[N] = math.exp((lo + hi) / 2)
    for s in (1, 2, 3):
        r = run(['tp', N, tuned[N], s, 3000]); out.write(json.dumps(r) + '\n'); out.flush()
        print(N, s, r['n'], r['mm'], r['mid'], r['chain_d'], r['sec'], flush=True)
