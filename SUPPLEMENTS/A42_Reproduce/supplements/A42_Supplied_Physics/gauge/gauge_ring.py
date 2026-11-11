"""Supplied U(1) links for the pass-7 local bond circuit.

Tests exact static/time-dependent gauge covariance, spectral flux shifts,
and a one-loop two-arm phase comparison against an analytic massless target.
Python + NumPy only. Units hbar=c=q=1, circumference 2*pi.
"""
from pathlib import Path
import json
import numpy as np


def bond(state, theta, links, parity, fraction=1.0):
    n = len(state)
    left = np.arange(parity, n, 2)
    right = (left + 1) % n
    c = np.cos(fraction * theta[left])
    s = np.sin(fraction * theta[left])
    u = links[left]
    if state.ndim == 2:
        c, s, u = c[:, None], s[:, None], u[:, None]
    a, b = state[left].copy(), state[right].copy()
    state[left] = c * a - s * u.conj() * b
    state[right] = s * u * a + c * b


def step(state, theta, links, parity=0):
    result = state.copy()
    bond(result, theta, links, parity, .5)
    bond(result, theta, links, 1-parity)
    bond(result, theta, links, parity, .5)
    return result


def gauged_links(links, phases):
    return np.roll(phases, -1) * links * phases.conj()


def covariance_test(n=128, steps=80):
    dx = 2*np.pi/n
    dt = .8*dx
    x = dx*np.arange(n)
    theta = np.arcsin(.8*(.9+.2*np.cos(x+dx/2))/2)
    links = np.exp(1j*(1.1/n + .2*(np.sin(x+dx)-np.sin(x))))
    rng = np.random.default_rng(74281)
    initial = rng.normal(size=n)+1j*rng.normal(size=n)
    initial /= np.linalg.norm(initial)
    gauge = np.exp(1j*(.7*np.sin(x)+.2*np.cos(3*x)))
    transformed_links = gauged_links(links, gauge)
    base = initial.copy()
    changed = gauge*initial
    wrong = changed.copy()
    max_error = 0.
    for k in range(steps):
        base = step(base, theta, links, k%2)
        changed = step(changed, theta, transformed_links, k%2)
        wrong = step(wrong, theta, links, k%2)
        max_error = max(max_error, float(np.linalg.norm(changed-gauge*base)))
    static = {
        'sites':n, 'steps':steps, 'max_step_state_covariance_error':max_error,
        'final_density_l1_difference':float(np.sum(abs(abs(changed)**2-abs(base)**2))),
        'wilson_loop_difference':float(abs(np.prod(links)-np.prod(transformed_links))),
        'wrong_untransformed_links_state_error':float(np.linalg.norm(wrong-gauge*base)),
        'norm_error':float(abs(np.linalg.norm(changed)-1)),
    }
    def phases(t):
        return np.exp(1j*(.7*np.sin(x)*np.cos(.9*t)+.2*np.cos(2*x)*np.sin(.6*t)))
    base = initial.copy()
    changed = phases(0)*initial
    wrong = changed.copy()
    max_error = 0.
    max_plaquette_error = 0.
    for k in range(steps):
        before, after = phases(k*dt), phases((k+1)*dt)
        links_now = gauged_links(links, before)
        links_after = gauged_links(links, after)
        temporal_link = after*before.conj()
        # Spatial forward, temporal forward at its right end, spatial
        # backward on the next slice, temporal backward at its left end.
        plaquette = links_now*np.roll(temporal_link,-1)*links_after.conj()*temporal_link.conj()
        max_plaquette_error = max(max_plaquette_error,float(np.max(abs(plaquette-1))))
        base = step(base, theta, links, k%2)
        changed = temporal_link*step(changed, theta, links_now, k%2)
        wrong = step(wrong, theta, links_now, k%2)
        max_error = max(max_error, float(np.linalg.norm(changed-after*base)))
    dynamic = {
        'sites':n,'steps':steps,'dt':dt,
        'max_step_state_covariance_error':max_error,
        'max_spacetime_plaquette_difference_from_static_base':max_plaquette_error,
        'final_density_l1_difference':float(np.sum(abs(abs(changed)**2-abs(base)**2))),
        'wrong_omitted_temporal_link_state_error':float(np.linalg.norm(wrong-phases(steps*dt)*base)),
        'wrong_omitted_temporal_link_density_l1_difference':float(np.sum(abs(abs(wrong)**2-abs(base)**2))),
        'norm_error':float(abs(np.linalg.norm(changed)-1)),
    }
    assert static['max_step_state_covariance_error'] < 2e-12
    assert dynamic['max_step_state_covariance_error'] < 2e-12
    assert dynamic['max_spacetime_plaquette_difference_from_static_base'] < 2e-12
    assert static['wrong_untransformed_links_state_error'] > .1
    assert dynamic['wrong_omitted_temporal_link_state_error'] > .1
    return {'static': static, 'time_dependent': dynamic}


def bloch_spectrum(n, flux, mode, nu=1.):
    dx = 2*np.pi/n
    dt = nu*dx
    theta = np.full(n, np.arcsin(nu/2))
    links = np.full(n, np.exp(1j*flux/n))
    cell_phase = np.exp(1j*mode*dx*np.arange(0,n,2))
    basis = np.zeros((n,2),dtype=complex)
    basis[::2,0] = cell_phase/np.sqrt(n/2)
    basis[1::2,1] = cell_phase/np.sqrt(n/2)
    evolved = step(basis, theta, links)
    matrix = basis.conj().T @ evolved
    values, vectors = np.linalg.eig(matrix)
    smooth = np.array([1, np.exp(1j*mode*dx)])/np.sqrt(2)
    chosen = np.argmax(abs(smooth.conj()@vectors)**2)
    omega = -np.angle(values[chosen])/dt
    p = mode-flux/(2*np.pi)
    discrete_exact = 2*np.arcsin(np.sin(theta[0])*np.sin(p*dx))/dt
    return {'sites':n,'flux_phase':float(flux),'mode':mode,'frequency':float(omega),
            'continuum_frequency':float(p),'finite_step_exact_frequency':float(discrete_exact),
            'continuum_error':float(abs(omega-p)),
            'finite_step_formula_error':float(abs(omega-discrete_exact)),
            'bloch_subspace_leakage':float(np.linalg.norm(evolved-basis@matrix))}


def flux_spectrum_tests():
    rows = [bloch_spectrum(n, .7*np.pi, 3) for n in [32,64,128,256,512,1024]]
    scan = [bloch_spectrum(256,flux,mode) for flux in [0.,.7*np.pi,1.3*np.pi]
            for mode in [-2,0,1,3]]
    # Whole-spectrum periodicity under one flux quantum, independent of labels.
    n = 64
    theta = np.full(n,np.arcsin(.5))
    phases = []
    for flux in [.7*np.pi,2.7*np.pi]:
        matrix = step(np.eye(n,dtype=complex),theta,np.full(n,np.exp(1j*flux/n)))
        phases.append(np.sort(np.angle(np.linalg.eigvals(matrix))))
    periodicity = float(np.max(abs(phases[0]-phases[1])))
    max_formula = max(r['finite_step_formula_error'] for r in rows+scan)
    no_links = bloch_spectrum(1024,0.,3)
    assert max_formula < 1e-11
    assert periodicity < 1e-12
    assert rows[-2]['continuum_error']/rows[-1]['continuum_error'] > 3.9
    return {'refinement':rows,'flux_and_mode_scan':scan,
            'whole_spectrum_2pi_flux_periodicity_error':periodicity,
            'max_finite_step_formula_error':max_formula,
            'dropped_links_frequency_at_1024_sites':no_links['frequency'],
            'dropped_links_error_against_flux_reference':abs(no_links['frequency']-2.65),
            'dropped_links_limiting_error':.35}


def one_loop_interference(n, flux=.7*np.pi):
    dx = 2*np.pi/n
    # T=2pi, dt=dx, exactly n steps with alternating parity.
    x = dx*np.arange(n)
    modes = np.arange(1,6)
    coeff = np.exp(-.5*((modes-3)/1.)**2)
    initial = np.sum(coeff[:,None]*np.exp(1j*modes[:,None]*(x[None,:]-.7)),axis=0)
    initial /= np.linalg.norm(initial)
    theta = np.full(n,np.arcsin(.5))
    links = np.full(n,np.exp(1j*flux/n))
    zero = np.ones(n,dtype=complex)
    reference, threaded = initial.copy(), initial.copy()
    for k in range(n):
        reference = step(reference,theta,zero,k%2)
        threaded = step(threaded,theta,links,k%2)
    overlap = np.vdot(reference,threaded)
    measured = .5*(1+overlap.real)
    exact = .5*(1+np.cos(flux))
    # A coherent path ancilla controls reference/threaded evolution; its plus
    # measurement gives this probability. Density on one chiral arm alone
    # cannot reveal the spatially uniform flux-dependent phase.
    row = {'sites':n,'flux_phase':float(flux),'time':2*np.pi,'dt_over_dx':1.,
           'measured_plus_probability':float(measured),'exact_plus_probability':float(exact),
           'probability_error':float(abs(measured-exact)),
           'relative_arm_state_error':float(np.linalg.norm(threaded-np.exp(1j*flux)*reference)),
           'threaded_raw_continuum_state_error':float(np.linalg.norm(threaded-np.exp(1j*flux)*initial)),
           'arm_overlap_real':float(overlap.real),'arm_overlap_imag':float(overlap.imag),
           'max_norm_error':float(max(abs(np.linalg.norm(reference)-1),abs(np.linalg.norm(threaded)-1))),
           'dropped_links_plus_probability':1.,'dropped_links_probability_error':float(1-exact)}
    return row


def main():
    covariance = covariance_test()
    spectrum = flux_spectrum_tests()
    interference = [one_loop_interference(n) for n in [64,128,256,512,1024]]
    assert interference[-2]['probability_error']/interference[-1]['probability_error'] > 3.8
    assert max(r['max_norm_error'] for r in interference) < 2e-12
    out = {'scope':'Supplied classical U(1) link field; one-particle local unitary circuit. No dynamical Maxwell field, many-body QED, charge derivation, or gravity derivation.',
           'conventions':'hbar=c=q=1, ring length 2pi, periodic scalar spatial boundary; psi_prime=exp(i chi) psi, forward u_j=exp(i integral A dx), H=-i(partial_x-i A) on smooth right-moving branch.',
           'covariance':covariance,'spectrum':spectrum,'one_loop_interference':interference,
           'all_assertions_passed':True}
    Path(__file__).with_name('gauge_ring_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
