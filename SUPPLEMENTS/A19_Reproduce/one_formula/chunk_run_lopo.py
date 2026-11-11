"""T27 runner: chunk_run.py with one change. MaNGA and CLASH are left out of every training set, so their scores are
out-of-population predictions; everything else (partitions, disk rule, fits, losses, records, scoring) is as chunk_run.py.
The model is registered under the alias <model>X so results never mix with the trained runs.
Usage: python3 chunk_run_lopo.py <harness module> <mask> <cfg> <out.jsonl> <budget seconds>
       python3 chunk_run_lopo.py merge <harness module> <mask> <out.jsonl> <res.json>"""
import sys, json, csv, time, hashlib, importlib
from pathlib import Path
import numpy as np
def parts_done(path):
    done = set()
    if Path(path).exists():
        for l in open(path):
            r = json.loads(l); done.add((r['assignment'], r['fold']))
    return done
def main(modname, mask, cfg, out, budget):
    M = importlib.import_module(modname); Z = M
    F = __import__('harness_T21')
    for k in list(F.FLOOR): F.FLOOR[k + 'X'] = F.FLOOR[k]
    Hh = __import__('harness_T15_hashed')
    for k in list(Hh.MODELS): Hh.MODELS.setdefault(k + 'X', Hh.MODELS[k])
    H = __import__('harness_T15_hashed'); B = __import__('harness_T15b')
    t0 = time.monotonic()
    gs = H.disks(); folds = list(csv.DictReader(open(H.PKG / 'data/SPARC/A19_Disk_Folds.csv')))
    maps = json.load(open(H.AST / 'T56_fold_maps.json'))['external_fold_maps']; cl = H.cluster_rows()
    ext = {'dwarfs': H.arr(H.dwarf_rows()), 'ETG': H.arr(H.etg_rows()), 'MaNGA': H.arr([r for r in cl if r['kind'] == 'MaNGA_BCG']),
           'CLASH': H.arr([r for r in cl if r['kind'] != 'MaNGA_BCG'])}
    done = parts_done(out)
    for asg in sorted(maps):
        for fold in range(5):
            if (asg, fold) in done: continue
            if time.monotonic() - t0 > budget: return
            held = {r['galaxy'] for r in folds if r[asg] == str(fold)}
            tr = [g for g in gs if g['name'] not in held]; te = [g for g in gs if g['name'] in held]
            dtr, dte = H.arr(H.disk_rows(tr)), H.arr(H.disk_rows(te)); etr, ete = {}, {}
            for p, a in ext.items():
                k = np.array([maps[asg][p][n] == fold for n in a['name']]); etr[p] = H.sub(a, ~k); ete[p] = H.sub(a, k)
            for p in ('MaNGA', 'CLASH'): del etr[p]
            lens = H.Lens(tr, mask == 'strict')
            B.FROZEN['abc'] = np.array(H.fit('ORIG', None, dtr, etr, lens, [0, 1, 2])['th'])
            model, w = cfg.split('@'); recs = []
            for j in range(4):
                groups = [i for i in range(4) if i != j]
                f = H.fit(model, float(w), dtr, etr, lens, groups); th = np.array(f['th'])
                l = {'disk': H.losses(model, th, dte)}
                for p, a in ete.items(): l[p] = H.losses(model, th, a)
                for kind in ['CLASH_centre', 'CLASH_outer', 'ETG_inner', 'ETG_outer']:
                    a = ete['CLASH' if kind.startswith('CLASH') else 'ETG']; l[kind] = H.losses(model, th, H.sub(a, a['kind'] == kind))
                z = lens.white(model, th, [j])
                recs.append(dict(config=cfg, assignment=asg, fold=fold, lens=j, fit=f, losses=l, test_lensing_chi2=float(z @ z)))
            with open(out, 'a') as fh:
                fh.write(json.dumps(dict(assignment=asg, fold=fold, records=recs)) + '\n')
            print(mask, cfg, asg, fold, round(time.monotonic() - t0, 1), flush=True)
def merge(modname, mask, out, res):
    M = importlib.import_module(modname); H = __import__('harness_T15_hashed')
    recs = [r for l in open(out) for r in json.loads(l)['records']]
    assert len(recs) == 60, len(recs)
    cfgs = sorted(set(r['config'] for r in recs))
    json.dump(dict(mask=mask, configs=cfgs, summary=H.summarize(recs), records=recs, runner='chunk_run.py',
                   runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), open(res, 'w'))
    print(json.dumps(H.summarize(recs), indent=1))
if __name__ == '__main__':
    if sys.argv[1] == 'merge': merge(*sys.argv[2:6])
    else: main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
