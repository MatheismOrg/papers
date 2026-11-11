"""Audit the independently solved continuum oscillators and largest fermion run."""
from pathlib import Path
import json
import numpy as np
import bosonization_test as b
import lattice_scaling as s
root=Path(__file__).parent;r=json.loads((root/'bosonization_results.json').read_text());assert r['complete']
reference=r['continuum'][-1];tight=b.continuum_modes(80,tol=2e-13)
control=b.continuum_modes(8,coupling=0.)
sudden=b.continuum_modes(4,tau=1e-4,tol=2e-13)
k=np.arange(1,5)*np.pi/8;mu=1/np.sqrt(np.pi)
wi=np.sqrt(k*k+(mu*s.profile(-.001,1e-4))**2)
wf=np.sqrt(k*k+(mu*s.profile(.001,1e-4))**2)
sudden_exact=(wf-wi)**2/(4*wi*wf)
massless_free=s.make_basis(10,boundary=-.25)
free_result,_=s.evolve(massless_free,length=8,mass=0,coupling=0,tau=.15,improved=True)
variants=[]
for flux in [0.,-.25]:
    chain=s.make_basis(12,boundary=flux)
    for improved in [False,True]:
        z,_=s.evolve(chain,length=8,mass=0,coupling=1,tau=.15,improved=improved)
        variants.append(dict(boundary=flux,improved=improved,survival=z['ground_survival'],excess_energy=z['final_excess_energy']))
chain=s.make_basis(20,boundary=-.25)
z,_=s.evolve(chain,length=8,mass=0,coupling=1,tau=.15,improved=True,tol=2e-11)
largest=r['rows'][-1]
out=dict(continuum_tight_survival_difference=abs(reference['ground_survival']-tight['ground_survival']),
    continuum_tight_energy_difference=abs(reference['final_excess_energy']-tight['final_excess_energy']),
    mode80_to160_survival_change=abs(reference['ground_survival']-r['continuum'][-2]['ground_survival']),
    mode80_to160_energy_change=abs(reference['final_excess_energy']-r['continuum'][-2]['final_excess_energy']),
    free_boson_control=control,free_fermion_control=free_result,
    sudden_occupation_max_error=float(np.max(abs(np.array(sudden['occupations'])-sudden_exact))),
    largest_fermion_survival_tolerance_difference=abs(largest['ground_survival']-z['ground_survival']),
    largest_fermion_energy_tolerance_difference=abs(largest['final_excess_energy']-z['final_excess_energy']),
    largest_tight_result=z,boundary_mass_controls=variants)
assert out['continuum_tight_survival_difference']<1e-9
assert out['continuum_tight_energy_difference']<1e-8
assert abs(control['ground_survival']-1)<1e-12
assert abs(free_result['ground_survival']-1)<1e-10
assert abs(free_result['final_excess_energy'])<1e-9
assert out['sudden_occupation_max_error']<1e-7
assert out['largest_fermion_survival_tolerance_difference']<1e-7
assert out['largest_fermion_energy_tolerance_difference']<1e-7
out['all_assertions_passed']=True
(root/'bosonization_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
