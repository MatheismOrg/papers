"""Build the SA74 night supplement zip: plans with their hashes, engines with theirs, the Mac queue and logs, every night chain
(cloud stages AS to AS5, AU1, AU2, HY1, FSS2 to FSS4 and the 8^4 FSS chains; Mac N2b, Q1-Q6, runner 2 and the VM slice runs), the analysis
scripts, their result files and figures, and the night log. Usage: night_supplement.py OUT.zip [main|chains] or night_supplement.py OUTDIR lean: main leaves out the cloud chain files, chains holds only them,\nlean writes a folder without any chain files (for the portfolio)."""
import sys, zipfile, glob, hashlib, json, os
from pathlib import Path
out = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
R2 = Path('/home/claude/sa74/regge/r2'); AUD = R2 / 'audit_a005'; NIGHT = Path('/home/claude/sa74/NIGHT'); MAC = NIGHT / 'mac_pull'
files = []   # (path on disk, path in zip)
def add(p, arc):
    p = Path(p)
    if p.is_file(): files.append((p, arc))
for f in ['plan_as.json', 'plan_as2.json', 'plan_as3.json', 'plan_as4.json', 'plan_hy1.json', 'plan_fss2.json', 'plan_fss3.json', 'plan_fss4.json', 'plan_as5.json', 'plan_fss8.json', 'plan_au1.json', 'plan_au2.json',
          'run_audit_r2.py', 'analyze_r2.py']:
    add(AUD / f, f'plans_cloud/{f}')
for pat in ['AS_*', 'AS2_*', 'AS3_*', 'AS4_*', 'AS5_*', 'AU1_*', 'AU2_*', 'HY1_*', 'FSS2_*', 'FSS3_*', 'FSS4_*', 'FSS_L8_*', 'prepare_L4_*', 'prepare_L6_*']:
    for f in sorted(glob.glob(str(AUD / 'runs' / pat))):
        if f.endswith(('.dat', '.log', '.execution.json')): add(f, f'chains_cloud/{Path(f).name}')
for f in sorted(MAC.glob('*.dat')) + sorted(MAC.glob('*.log')): add(f, f'chains_mac/{f.name}')
for f in sorted((MAC / 'slices').glob('*')): add(f, f'slices_vm/{f.name}')
for f in ['/mnt/user-data/outputs/night/' + x for x in ['regge_r2_big.c', 'regge_r2_big_s.c', 'regge_r2_big_x.c', 'night_queue.sh', 'start_queue.sh', 'switch_to_mac.sh',
                                                       'night_queue2.sh', 'start_queue2.sh']]:
    add(f, f'mac_runner/{Path(f).name}')
for f in ['queue.txt', 'queue2.txt', 'night_queue.log', 'night_queue2.log']:   # final versions, staged from the Mac at the end
    add(MAC / 'runner' / f, f'mac_runner/{f}')
add(R2 / 'regge_r2.c', 'engines/regge_r2.c'); add(R2 / 'make_r2.py', 'engines/make_r2.py')
for f in ['night/fig_night.py', 'night/hysteresis.py', 'night/shape_stats.py', 'night/analyze_as3.py', 'night/night_final.py', 'night/night_supplement.py',
          'fss/analyze_fss.py', 'fss/fit_delta_a.py', 'night/reweight.py', 'night/multihist.py', 'night/fig_reweight.py', 'night/defects.py', 'night/analyze_au.py', 'night/analyze_as5.py']:
    add(R2 / f, f'analysis/{Path(f).name}')
for f in ['night/night_summary.json', 'night/hysteresis.json', 'night/shape_stats.json', 'night/as3_results.json', 'fss/fss_results.json', 'fss/fit_delta_a.json',
          'night/chi_vs_L.png', 'night/R_vs_k.png', 'night/hysteresis_8.png', 'night/reweight_skip00.json', 'night/multihist.json',
          'night/chi_reweighted.png', 'night/defects.json', 'night/au_results.json', 'night/as5_results.json']:
    add(R2 / f, f'results/{Path(f).name}')
add(NIGHT / 'NIGHT_LOG.md', 'NIGHT_LOG.md'); add(NIGHT / 'NIGHT_PLAN.md', 'NIGHT_PLAN.md')
if mode == 'main': files = [(p, a) for p, a in files if not a.startswith('chains_cloud/')]
if mode == 'lean': files = [(p, a) for p, a in files if not a.startswith(('chains_cloud/', 'chains_mac/'))]
if mode == 'chains': files = [(p, a) for p, a in files if a.startswith('chains_cloud/')]
if mode.startswith('chains:'):   # chains:PREFIX,PREFIX,... selects cloud chain files by stage prefix (to keep each zip small)
    pre = tuple('chains_cloud/' + x for x in mode[7:].split(','))
    files = [(p, a) for p, a in files if a.startswith(pre)]
def _done(p):   # cloud runs still going (or cut off) are left out: a .dat counts only with its '# done' line, its .log/.execution.json with it
    d = Path(str(p).replace('.execution.json', '.dat').replace('.log', '.dat'))
    return d.exists() and '# done' in d.read_text(errors='ignore')[-20000:]
files = [(p, a) for p, a in files if not a.startswith('chains_cloud/') or a.startswith('chains_cloud/prepare') or _done(p)]
add(Path('/home/claude/sa74/regge/r2/night/SUPPLEMENT_README.md'), 'README_NIGHT.md') if not mode.startswith('chains') else None
man = {arc: hashlib.sha256(p.read_bytes()).hexdigest() for p, arc in files}
readme = ['# SA74 night of 6-7 November 2026: the curvature-squared Regge ensemble on larger grids', '',
          'Every plan was written and hashed before its runs (see NIGHT_LOG.md for the times); exploratory runs are labeled so in their plans.',
          'Cloud chains: run_audit_r2.py with engines/regge_r2.c (built by make_r2.py; R2A = a). Mac chains: mac_runner/regge_r2_big.c, the same engine',
          'with a time box (MAXSEC) and grids up to 16^4, run by night_queue.sh from queue.txt (plans Q1 to Q6, fixed in place) and by the native runner',
          '(plan N2b); runner 2 (night_queue2.sh, plan A1 in queue2.txt) runs regge_r2_big_x.c with SEQ=1 (edges in order) and LINES=1 (Wilson lines).',
          'Slice runs: regge_r2_big_s.c, identical rows plus per-slice curvature sums, run in the Mac VM in time-boxed chunks.',
          'Columns of every chain file are listed in its header; analysis/night_final.py rebuilds results/ from the chain files.', '',
          f'{len(files)} files; SHA-256 of each in MANIFEST.json.']
if mode == 'lean':   # a folder for the portfolio: everything except the per-sweep chain files, which go in the delivered zips
    import shutil
    if os.path.isdir(out): shutil.rmtree(out)
    for p, arc in files:
        os.makedirs(os.path.dirname(os.path.join(out, arc)), exist_ok=True); shutil.copy(p, os.path.join(out, arc))
    readme += ['', 'LEAN DISTRIBUTION NOTE: the per-sweep chain files (chains_cloud/, chains_mac/) are in the two night zips delivered with T28',
               '(SA74_T28_A42_REGGE_NIGHT_MAIN_*.zip with the Mac chains, SA74_T28_A42_REGGE_NIGHT_CHAINS_*.zip with the cloud chains).']
    open(os.path.join(out, 'MANIFEST.sha256.json'), 'w').write(json.dumps(man, indent=1) + '\n')
    open(os.path.join(out, 'README_FILES.md'), 'w').write('\n'.join(readme) + '\n')
    print(out, len(files), 'files', round(sum(p.stat().st_size for p, a in files) / 1e6, 2), 'MB'); sys.exit(0)
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p, arc in files: z.write(p, 'SA74_NIGHT/' + arc)
    z.writestr('SA74_NIGHT/MANIFEST.json', json.dumps(man, indent=1)); z.writestr('SA74_NIGHT/README.md', '\n'.join(readme) + '\n')
print(out, len(files), 'files', round(os.path.getsize(out) / 1e6, 1), 'MB')
