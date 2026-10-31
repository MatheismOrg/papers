"""Run the actual massive local circuit on the geodesic packet family.
Uses the independently implemented spectral reference in geodesic_packet.py.
"""
from pathlib import Path
import importlib.util
import json
import numpy as np
from scipy.signal import resample
from geodesic_packet import initial_state, quantum, ray

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('massive',root/'massive/massive_curved.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def center(state):
    n=len(state);x=np.arange(n)*2*np.pi/n;x0=np.pi/2
    xx=x0+(x-x0+np.pi)%(2*np.pi)-np.pi
    rho=np.sum(abs(state)**2,axis=1)
    return float(np.sum(xx*rho))


def main():
    target_ray=ray([0.,1.])[-1]['x'];rows=[]
    for mass in [8.,16.,32.,64.]:
        qrecord,reference=quantum(2048,mass)
        reference=reference.reshape(2,2048).T
        for n in [512,1024,2048,4096,8192]:
            state,meta=module.local_circuit(n,.25,mass,.9,initial_state=initial_state(n,mass).T)
            reference_n=resample(reference,n,axis=0)*np.sqrt(2048/n)
            mean=center(state)
            row={'mass':mass,**meta,'state_error_vs_quantum_reference':float(np.linalg.norm(state-reference_n)),
                 'packet_center':mean,'quantum_reference_center':qrecord['trajectory'][-1]['mean_x'],
                 'point_ray_center':target_ray,'center_error_vs_quantum':mean-qrecord['trajectory'][-1]['mean_x'],
                 'center_error_vs_ray':mean-target_ray}
            rows.append(row);print(json.dumps(row),flush=True)
    bad,bmeta=module.local_circuit(2048,.25,64.,.9,wrong_mass=True,initial_state=initial_state(2048,64.).T)
    badrow={'mass':64.,**bmeta,'packet_center':center(bad),'initial_center':np.pi/2,
            'point_ray_center_for_correct_geometry':target_ray}
    ratios=[]
    for mass in [8.,16.,32.,64.]:
        r=[x for x in rows if x['mass']==mass]
        ratios.append({'mass':mass,'state_error_ratios':[r[i]['state_error_vs_quantum_reference']/r[i+1]['state_error_vs_quantum_reference'] for i in range(len(r)-1)]})
    assert all(3.5<ratio<4.5 for r in ratios for ratio in r['state_error_ratios'])
    out={'scope':'Actual calibrated alternating local unitary circuit, including mass coins, evaluated on the same narrowing positive-band packets as the independent GR ray/quantum reference. Finite mass quantum spread is separated from discretization error.',
         'rows':rows,'refinement':ratios,'wrong_mass_coupling':badrow,'all_assertions_passed':True}
    Path(__file__).with_name('circuit_geodesic_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
