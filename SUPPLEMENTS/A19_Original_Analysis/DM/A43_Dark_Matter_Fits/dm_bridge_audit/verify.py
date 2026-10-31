"""Rerun independent audits and write one aggregate pass/fail report.

Run from any directory: python path/to/dm_bridge_audit/verify.py
Requires Python 3, NumPy and SciPy. The sibling analysis directories and their
frozen inputs must remain in place. No fitting code is imported or executed.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys,json,hashlib
HERE=Path(__file__).resolve().parent
NAMES=['shared_halo','cosmological_halo','gas_bridge','field_bridge']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name):
 script=HERE/('verify_'+name+'.py');p=subprocess.run([sys.executable,str(script)],cwd=HERE,capture_output=True,text=True)
 (HERE/(name+'_audit.log')).write_text(p.stdout+p.stderr)
 if p.returncode:
  return name,{'status':'FAIL','returncode':p.returncode,'log':name+'_audit.log','error_tail':(p.stdout+p.stderr)[-3000:]}
 report=json.loads((HERE/(name+'_verification.json')).read_text())
 if report.get('status','').lower() not in ['pass','passed'] or report['audit_script_sha256']!=sha(script):return name,{'status':'FAIL','reason':'Audit did not return a current passing report.'}
 report['status']='PASS';return name,report
with ThreadPoolExecutor(max_workers=4) as pool:results=dict(pool.map(run,NAMES))
status='PASS' if all(v['status']=='PASS' for v in results.values()) else 'FAIL'
report={'status':status,'audit_scope':'Independent arithmetic, density integrals, physical bounds, feature provenance and declared validation. Passing verifies reproduction and mathematical implementation, not physical truth, predictive universality, or independent-discovery status.','branches':results,'aggregate':{},'entrypoint_sha256':sha(Path(__file__))}
if status=='PASS':
 totals={
  'prediction_tables':sum(v.get('prediction_tables',v.get('tables',0)) for v in results.values()),
  'prediction_rows':sum(v['prediction_rows'] for v in results.values()),
  'fit_sets':sum(v['fit_sets'] for v in results.values()),
  'optimizer_attempts':sum(v['optimizer_attempts'] for v in results.values()),
  'descriptive_bootstrap_comparisons':sum(v.get('bootstrap_comparisons',v.get('descriptive_bootstrap_comparisons',0)) for v in results.values()),
  'full_disk_diagnostic_rows_separate':results['gas_bridge']['full_disk_diagnostic_rows'],
  'copied_comparator_tables_separate':results['cosmological_halo']['comparator_tables_verified'],
  'independent_cluster_density_integrals':results['gas_bridge']['cluster_density_integrals'],
  'independent_kernel_basis_and_midpoint_integrals':results['field_bridge']['independent_basis_integrals'],
  'max_recomputed_prediction_discrepancy_dex':max(v['max_prediction_discrepancy_dex'] for v in results.values())}
 assert totals['prediction_tables']==284 and totals['prediction_rows']==158920 and totals['fit_sets']==426 and totals['optimizer_attempts']==1362
 report['aggregate']=totals
(HERE/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({'status':status,'aggregate':report['aggregate'],'branches':{k:v['status'] for k,v in results.items()}},indent=2))
sys.exit(0 if status=='PASS' else 1)
