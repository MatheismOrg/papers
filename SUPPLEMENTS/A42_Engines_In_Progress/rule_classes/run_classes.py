#!/usr/bin/env python3
"""Run complete rule classes with the TOE@home client code (toehome.py, version 1), without a coordinator.
   python3 run_classes.py WORKER NWORKERS CLASS[,CLASS...]
Worker w takes the blocks whose position in the class's block list is w modulo NWORKERS, in order, and appends one line per
block to out/CLASS.wW.jsonl (block label, digest, seconds, records). A restart skips the blocks already in the file."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import toehome as T
w, n = int(sys.argv[1]), int(sys.argv[2]); classes = sys.argv[3].split(',')
os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
try: os.nice(19)
except Exception: pass
for cls in classes:
    path = os.path.join(HERE, 'out', '%s.w%d.jsonl' % (cls, w)); done = set()
    if os.path.exists(path):
        for line in open(path):
            try: done.add(json.loads(line)['uid'])
            except ValueError: pass
    units = [u for i, u in enumerate(T.unit_list(cls)) if i % n == w]
    t0 = time.time(); nr = 0
    with open(path, 'a') as fh:
        for u in units:
            if T.uid(u) in done: continue
            t1 = time.time(); res = T.run_unit(u); nr += len(res['rules'])
            fh.write(json.dumps({'uid': res['uid'], 'digest': T.digest(res), 'sec': round(time.time() - t1, 1), 'rules': res['rules']}, separators=(',', ':')) + '\n'); fh.flush()
    print('%s worker %d of %d: %d blocks, %d rules run now, %.0f s' % (cls, w, n, len(units), nr, time.time() - t0), flush=True)
