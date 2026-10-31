from pathlib import Path
import csv, json, hashlib, zipfile
import numpy as np

R=Path(__file__).resolve().parent;T=R.parent/'dm_unified_stress';C=R.parent/'dm_unified_limits/results/results.json'
D=R.parent/'dm_fixed_radius/fixed_radius_inputs.csv';J=json.loads((T/'results.json').read_text());F=json.loads(C.read_text())
assert hashlib.sha256(D.read_bytes()).hexdigest()==J['input_sha256']
assert hashlib.sha256(C.read_bytes()).hexdigest()==J['fit_result_sha256']
assert hashlib.sha256((T/'PROTOCOL.txt').read_bytes()).hexdigest()==J['protocol_sha256']
frozen=R.parent.parent/'output/dm_environment/A43_Code_and_Data.zip'
frozen_input_sha256='265d50f5a572d66f615ebb349abd4944ac92c89d400d3a65c94a22e0c96da04c'
assert hashlib.sha256(D.read_bytes()).hexdigest()==frozen_input_sha256
archive_checked=False
if frozen.exists():
    with zipfile.ZipFile(frozen) as z:
        assert z.read('A43_Dark_Matter_Fits/dm_fixed_radius/fixed_radius_inputs.csv')==D.read_bytes()
    archive_checked=True
raw=list(csv.DictReader(D.open()));lookup={(r['cluster'],float(r['radius_kpc'])):r for r in raw};assert len(lookup)==21
assert len({r['cluster'] for r in raw})==7 and {float(r['radius_kpc']) for r in raw}=={500.,750.,1000.}
pred=list(csv.DictReader((T/'predictions.csv').open()));assert len(pred)==1890
def direct(k,c,gb,r):
    x=np.log10(gb/1e-10);u=np.log10(r);q=x+u;m=x+2*u
    if k=='RAR':ae=10**c[0]
    elif k=='potential_quad':ae=10**(c[0]+c[1]*q+c[2]*q*q)
    elif k=='mass_bounded':
        z=np.tanh((m-c[3])/c[4]);ae=10**(c[0]+c[1]*z+c[2]*z*z)
    elif k=='power_sum':ae=10**c[0]+10**(c[1]-c[3]*q)+10**(c[2]+c[4]*q)
    elif k=='half_power':ae=10**c[0]+10**(c[1]-.5*q)+10**(c[2]+c[3]*q)
    else:
        t=q if k=='mixed_power' else x+c[5]*u if k=='index_power' else m
        ae=10**c[0]+10**(c[1]-c[3]*m)+10**(c[2]+c[4]*t)
    ae*=1e-10
    if k=='mass_additive':gp=gb+np.sqrt(gb*ae)
    elif k=='mass_screened':gp=gb+np.sqrt(gb*ae)*np.exp(-np.sqrt(gb/(1e-10*10**c[0])))
    else:gp=gb/(1-np.exp(-np.sqrt(gb/ae)))
    return gp
keys=set();maxerr=0.
for r in pred:
    case=r['case'];model=r['model'];cluster=r['cluster'];rad=float(r['radius_kpc']);orig=lookup[(cluster,rad)];target,baryon=case.split('__')
    key=(case,model,cluster,rad);assert key not in keys;keys.add(key)
    assert target in ['forward','nfw','einasto'] and baryon in ['gbar','gbar_no_stars','gbar_full_projected_stars']
    gb=float(r['gbar_m_s2']);go=float(r['gobs_m_s2']);gp=float(r['gpred_m_s2'])
    assert gb==float(orig[baryon]) and go==float(orig['gobs_'+target])
    expected=direct(model,F['models'][model]['coefficients'],gb,rad)
    err=abs(np.log10(gp/expected));maxerr=max(maxerr,err);assert err<1e-12
    assert abs(float(r['log10_residual'])-np.log10(gp/go))<1e-12
for case,v in J['cases'].items():
    assert len(v)==10
    for model,met in v.items():
        rr=[r for r in pred if r['case']==case and r['model']==model];assert len(rr)==21
        clusters=sorted({r['cluster'] for r in rr});assert len(clusters)==7
        e=np.log10([float(r['gpred_m_s2'])/float(r['gobs_m_s2']) for r in rr]);g=np.array([r['cluster'] for r in rr])
        rms=np.sqrt(np.mean([np.mean(e[g==c]**2) for c in clusters]));median=100*np.median(abs(10**e-1))
        assert abs(rms-met['RMS_dex'])<1e-12 and abs(median-met['median_acceleration_error_percent'])<1e-10
        for c in clusters:assert abs(np.sqrt(np.mean(e[g==c]**2))-met['per_cluster_RMS_dex'][c])<1e-12
assert len(J['cases'])==9
out={'status':'passed','prediction_rows_checked':len(pred),'cases_checked':len(J['cases']),'models_checked':len(F['models']),'maximum_direct_equation_error_dex':maxerr,'input_sha256':J['input_sha256'],'fit_result_sha256':J['fit_result_sha256'],'previous_archive_compared_byte_for_byte':archive_checked,'checks':['input matches frozen previous-delivery digest','21 fixed-radius inputs:7 hosts,500/750/1000 kpc','all nine model/systematic cases supplied','all ten models use frozen all-data five-population coefficients','all 1890 predictions and residuals independently recomputed','equal-host RMS,pooled median errors,per-host RMS independently recomputed','protocol/input/fit provenance hashes'],'interpretation':'No fitting to X-COP occurred in the frozen transfer. The already-inspected hydrostatic source estimates are correlated and model-dependent; nine cases are sensitivities, not nine independent datasets.'}
(R/'stress_verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
print('Primary X-COP transfer:')
for k,v in J['cases']['forward__gbar'].items():print(k,round(v['RMS_dex'],6),round(v['median_acceleration_error_percent'],3))
