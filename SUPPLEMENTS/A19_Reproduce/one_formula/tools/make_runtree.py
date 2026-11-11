"""Make a folder from which the one-formula scripts run, without changing any file of this repository.

Why this exists: the scripts, protocols and results in one_formula/ are byte-identical to the files that were hashed
before each run, so they still name the base folder of the session in which they ran. Two harness modules
(harness.py, harness_T15_hashed.py) and three scripts name it; everything else finds its inputs through them or
through the command line. This tool copies the code into a new folder and replaces that one base path in the copies.

Usage:  python3 tools/make_runtree.py <new folder>            (run from one_formula/, or give the tool's full path)

It builds, under <new folder>:
  dmfit/                   a copy of one_formula/ (code, protocols, saved results), base path replaced in the .py files
  pkg/A19_Data_Package/    a copy of A19_Reproduce/code, blind_test and data (the inputs of protocols 1 to 9)
  so73_t63/                T56_fold_maps.json (the fold assignments) and inputs/data/ETG (the early-type sample, copied
                           from A19_Original_Analysis/DM/A43_Dark_Matter_Fits/dm_broad_groups/data)
  ext/, ext_hecs/          the small external tables (see DATA_SOURCES.md)
Then, for example, the held-out benchmark of the final formula (about ten minutes on four cores):
  cd <new folder>/dmfit
  PYTHONPATH=. python3 chunk_split.py harness_T61 primary MOGF0IESBC@0.1 my_parts.jsonl 6000 0 1
  PYTHONPATH=. python3 chunk_split.py merge harness_T61 primary my_res.json my_parts.jsonl
  python3 ../compare_runs.py my_res.json res61_MOGF0IESBC_primary.json
Needs Python 3 with numpy and scipy; some scripts also need colossus (pip install colossus)."""
import sys, shutil, os, json
from pathlib import Path

OLD = '/tmp/claude-0/-home-claude/0a1f618d-0405-5db3-a0ab-84f1f9545369/scratchpad'
EXPECT = ['harness.py', 'harness_T15_hashed.py', 'blind_hecs_T40.py', 't59_outside/drift61.py', 't59_outside/validate59.py']

COMPARE = '''"""Compare the summary of a rerun with a saved result. Usage: python3 compare_runs.py <my_res.json> <saved_res.json>"""
import sys, json
a, b = (json.load(open(p)) for p in sys.argv[1:3])
def flat(d, pre=''):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict): out.update(flat(v, pre + k + '.'))
        elif isinstance(v, (int, float)) and not isinstance(v, bool): out[pre + k] = float(v)
    return out
fa, fb = flat(a.get('summary', {})), flat(b.get('summary', {}))
keys = sorted(set(fa) & set(fb)); worst = 0.0
for k in keys:
    d = abs(fa[k] - fb[k]); rel = d / max(abs(fb[k]), 1e-12); worst = max(worst, rel)
    print(f"{k:60s} {fa[k]:.6g}  {fb[k]:.6g}  {'same' if rel < 1e-6 else 'differs by %.2e' % rel}")
print(len(keys), 'numbers compared; largest relative difference', '%.2e' % worst)
'''

def main(target):
    here = Path(__file__).resolve().parent.parent              # one_formula/
    rep = here.parent                                          # A19_Reproduce/
    etg = rep.parent / 'A19_Original_Analysis/DM/A43_Dark_Matter_Fits/dm_broad_groups/data'
    for need in (rep / 'code', rep / 'blind_test', rep / 'data', etg / 'ETG_Lelli2017.mrt', etg / 'Rotmod_ETG.zip', here / 'astra_T63/T56_fold_maps.json'):
        if not need.exists(): sys.exit(f'missing in this copy of the repository: {need}')
    T = Path(target).resolve()
    if T.exists() and any(T.iterdir()): sys.exit(f'{T} exists and is not empty; give a new folder')
    if here == T or here in T.parents: sys.exit('give a folder outside the repository')
    (T / 'dmfit').mkdir(parents=True)
    patched = []
    for root, dirs, files in os.walk(here):
        rel = Path(root).relative_to(here)
        dirs[:] = sorted(d for d in dirs if d != '__pycache__' and not (rel == Path('.') and d in ('tools', 'ext', 'ext_hecs', 'astra_T63')))
        (T / 'dmfit' / rel).mkdir(parents=True, exist_ok=True)
        for f in sorted(files):
            src, dst = Path(root) / f, T / 'dmfit' / rel / f
            if f.endswith('.py'):
                s = src.read_text()
                if OLD in s: s = s.replace(OLD, str(T)); patched.append(str(rel / f) if str(rel) != '.' else f)
                dst.write_text(s)
            else: shutil.copy2(src, dst)
    for sub in ('code', 'blind_test', 'data'):
        shutil.copytree(rep / sub, T / 'pkg/A19_Data_Package' / sub, ignore=shutil.ignore_patterns('__pycache__'))
    (T / 'so73_t63/inputs/data/ETG').mkdir(parents=True)
    shutil.copy2(here / 'astra_T63/T56_fold_maps.json', T / 'so73_t63/T56_fold_maps.json')
    for f in ('ETG_Lelli2017.mrt', 'Rotmod_ETG.zip'): shutil.copy2(etg / f, T / 'so73_t63/inputs/data/ETG' / f)
    shutil.copytree(here / 'ext', T / 'ext'); shutil.copytree(here / 'ext_hecs', T / 'ext_hecs')
    (T / 'compare_runs.py').write_text(COMPARE)
    if sorted(patched) != sorted(EXPECT): print('note: base path replaced in', sorted(patched), 'expected', sorted(EXPECT))
    left = [str(p.relative_to(T)) for p in (T / 'dmfit').rglob('*.py') if OLD in p.read_text().replace(str(T), '')]
    if left: sys.exit(f'base path still present in {left}')
    json.dump(dict(base=str(T), patched=sorted(patched), source=str(here)), open(T / 'RUNTREE.json', 'w'), indent=1)
    print(f'run tree ready at {T}\n  base path replaced in {len(patched)} scripts: {", ".join(sorted(patched))}\n  next: cd {T}/dmfit  (see the usage text at the top of this tool, or README.md)')

if __name__ == '__main__':
    if len(sys.argv) != 2: sys.exit(__doc__)
    main(sys.argv[1])
