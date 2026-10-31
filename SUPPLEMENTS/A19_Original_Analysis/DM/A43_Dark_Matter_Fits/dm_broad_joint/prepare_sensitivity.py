from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parent
r=list(csv.DictReader((ROOT.parent/'dm_broad_dwarfs/results/dwarfs_benchmark.csv').open()))
allowed={v['system'] for v in r if max(float(v['sigma_err_plus_km_s']),float(v['sigma_err_minus_km_s']))/float(v['sigma_los_km_s'])<=.25 and v['possible_tidal_or_classification_caution']=='False'}
rows=list(csv.DictReader((ROOT/'combined_input.csv').open()))
selected=[v for v in rows if v['population']!='dSph' or v['system'] in allowed]
with (ROOT/'combined_high_precision.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(selected)
meta=json.loads((ROOT/'input_metadata.json').read_text())
meta['sensitivity']='Post-search sensitivity: retain only25dSph with maxquoteddispersionerror/sigma<=0.25 and no source caution flag; allotherpopulations unchanged. Not an independentconfirmatorytest.'
meta['total_rows']=len(selected);meta['total_physical_systems']=len({v['system'] for v in selected});meta['total_validation_groups']=len({v['group'] for v in selected})
meta['counts']['dSph']['rows']=len(allowed);meta['counts']['dSph']['systems']=len(allowed);meta['counts']['dSph']['validation_groups']=len(allowed)
(ROOT/'high_precision_metadata.json').write_text(json.dumps(meta,indent=2))
print('retaineddwarfs',len(allowed),'totalrows',len(selected))
