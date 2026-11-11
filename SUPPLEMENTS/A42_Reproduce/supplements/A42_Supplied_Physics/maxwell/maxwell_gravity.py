"""Same local rule applied to transverse Maxwell fields on supplied geometry.
Metric: -N(x)^2 dt^2+dx^2+dy^2+dz^2; fields depend only on t,x.
Eulerian components E_y,B_z. Classical waves, not quantized photons.
"""
from pathlib import Path
import importlib.util
import json
import numpy as np

spec=importlib.util.spec_from_file_location('massive',Path(__file__).resolve().parents[1]/'massive/massive_curved.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def packet(y,Y,offset):
    j=np.arange(1,6);c=np.exp(-.5*((j-3)/1.2)**2)
    return np.sqrt(2/(Y*np.sum(c*c)))*np.sum(c[:,None]*np.cos((2*np.pi*j/Y)[:,None]*(y[None,:]-offset)),axis=0)


def analytic(n,t):
    dx,x,N,Nmid,y,Y=module.geometry(n,.25)
    chi=np.column_stack([.8*packet(y-t,Y,.9),.6*packet(y+t,Y,3.5)])/np.sqrt(N[:,None])
    return chi,dx,x,N


def fields(chi,N):
    return (chi[:,0]+chi[:,1])/(2*np.sqrt(N)),(chi[:,0]-chi[:,1])/(2*np.sqrt(N))


def main():
    rows=[]
    for n in [64,128,256,512,1024,2048]:
        chi0,dx,x,N=analytic(n,0.)
        initial=np.sqrt(dx)*chi0.astype(complex);normalizer=np.linalg.norm(initial);initial/=normalizer
        out,meta=module.local_circuit(n,.25,0.,.9,initial_state=initial)
        exact,_,_,_=analytic(n,1.);exact/=normalizer
        chi=out/np.sqrt(dx);e,b=fields(chi,N);ee,bb=fields(exact,N)
        norm=lambda z:float(np.sqrt(dx*np.sum(abs(z)**2)))
        energy=float(.5*dx*np.sum(N*(abs(e)**2+abs(b)**2)))
        field_relative=float(np.sqrt(dx*np.sum(abs(e-ee)**2+abs(b-bb)**2))/np.sqrt(dx*np.sum(abs(ee)**2+abs(bb)**2)))
        row={'sites':n,**meta,'characteristic_state_error':norm(chi-exact),
             'relative_maxwell_field_error':field_relative,'killing_energy':energy,
             'killing_energy_drift':abs(energy-.25),'max_imaginary_field':float(max(abs(e.imag).max(),abs(b.imag).max()))}
        rows.append(row);print(json.dumps(row),flush=True)
    ratios=[rows[i]['relative_maxwell_field_error']/rows[i+1]['relative_maxwell_field_error'] for i in range(len(rows)-1)]
    assert all(3.7<r<4.3 for r in ratios[1:])
    assert max(x['killing_energy_drift'] for x in rows)<1e-11
    assert max(x['max_imaginary_field'] for x in rows)<1e-12
    out={'scope':'Two counterpropagating classical transverse plane Maxwell waves on supplied 3+1 static metric with one-dimensional dependence. No quantized-photon or general 3D Maxwell claim.',
         'metric':'ds²=-N(x)²dt²+dx²+dy²+dz²; N=1+.25cosx',
         'eulerian_field_equations':'partial_t E_y=-partial_x(N B_z); partial_t B_z=-partial_x(N E_y)',
         'same_operator_mapping':'chi_plus/minus=sqrt(N)(E_y +/- B_z), partial_t chi_plus/minus=minus/plus[N partial_x+Nprime/2]chi_plus/minus',
         'energy':'Killing energy per unit transverse area = integral N(E_y²+B_z²)/2 dx = integral(|chi_plus|²+|chi_minus|²)/4 dx',
         'energy_flux':'N² E_y B_z',
         'gauss_constraints':'div E=div B=0 identically for E=(0,E_y(t,x),0), B=(0,0,B_z(t,x))',
         'rows':rows,'field_error_ratios':ratios,'all_assertions_passed':True}
    Path(__file__).with_name('maxwell_gravity_results.json').write_text(json.dumps(out,indent=2)+'\n')
    paragraph=('The same calibrated massless gates also pass a transverse Maxwell-wave check on the supplied 3+1 plane metric. The combinations sqrt(N)(E_y +/- B_z) obey the earlier transport operator. Two counterpropagating real waves show about fourfold field-error reduction per grid doubling; at 2,048 sites the relative field error is %.3g. The conserved Killing energy is preserved to %.1g. This covers this one-dimensional plane-wave sector, not general 3D electromagnetism or quantized photons.'%(rows[-1]['relative_maxwell_field_error'],max(x['killing_energy_drift'] for x in rows)))
    Path(__file__).with_name('MAXWELL_SUMMARY.json').write_text(json.dumps({'report_paragraph':paragraph},indent=2)+'\n')


if __name__=='__main__':main()
