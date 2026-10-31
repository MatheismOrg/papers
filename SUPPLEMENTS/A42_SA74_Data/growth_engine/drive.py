import subprocess, json, math, sys
B = '/home/claude/sa74/big/tp_big'
def run(args):
    return json.loads(subprocess.run([B] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout)
out = open('/home/claude/sa74/big/trajectory.jsonl', 'a')
for N in [int(x) for x in sys.argv[1:]]:
    lo, hi = math.log(8.0 / N), math.log(40.0 / N)        # pN from 8 (sparse, reads high) to 40 (dense, reads low)
    steps = []
    for _ in range(11):
        mid = (lo + hi) / 2; p = math.exp(mid)
        r = run(['mm', N, p, 101, 800]); steps.append((p * N, r['mm'], r['n']))
        if r['mm'] > 4.0: lo = mid
        else: hi = mid
    p = math.exp((lo + hi) / 2)
    out.write(json.dumps({'tune': N, 'p': p, 'pN': p * N, 'steps': steps}) + '\n'); out.flush()
    for s in (1, 2):
        r = run(['tp', N, p, s, 600]); out.write(json.dumps(r) + '\n'); out.flush()
        print(N, s, r['n'], r['mm'], r['mid'], r['chain_d'], r['sec'], flush=True)
