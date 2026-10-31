"""Run all independent audits from the supplied supplementary directory tree."""
from pathlib import Path
import subprocess,sys,json,hashlib
R=Path(__file__).resolve().parent
for name in ['verify_mond.py','verify_mond_bootstrap.py','verify_physics.py','verify_halo.py','verify_cluster.py','verify_cluster_bootstrap.py']:
 subprocess.run([sys.executable,str(R/name)],check=True)
reports={k:json.loads((R/(k+'_verification.json')).read_text()) for k in ['mond','halo','cluster','physics']}
assert all(v['status']=='passed' for v in reports.values())
main=[reports[k] for k in ['mond','halo','cluster']]
out={'status':'passed','prediction_tables':sum(v['prediction_tables'] for v in main),'prediction_rows':sum(v['prediction_rows'] for v in main),'fit_sets':sum(v['fit_sets'] for v in main),'optimizer_attempts':sum(v['optimizer_attempts'] for v in main),'failed_optimizer_attempts':sum(v['failed_optimizer_attempts'] for v in main),'descriptive_bootstrap_comparisons_checked':sum(v['descriptive_bootstrap_comparisons_checked'] for v in main),'individual_reports':reports,'audit_file_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(R.glob('*.py'))}}
(R/'verification.json').write_text(json.dumps(out,indent=2));print('ALL INDEPENDENT AUDITS PASS',out['prediction_tables'],'tables',out['prediction_rows'],'predictions')
