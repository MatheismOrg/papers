#!/usr/bin/env python3
"""Independent analytical/numerical audit of the RST endpoint baseline.

This verifies published charge/flux conventions and derived geometry identities.
It is not a simulation of quantum gravity or evidence for a new endpoint law.
Run: python tmp/toe_pass18/audit_rst.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.optimize import brentq


KAPPA = 12.0
LAMBDA = 0.7
X0 = 1.3
A = KAPPA * LAMBDA / 4.0
MUS = [0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]


def fraction(q):
    return q / (1.0 + q)


def old_charge(q, mu):
    return A * (mu - math.log1p(q) - fraction(q))


def corrected_charge(q, mu):
    return A * (mu - math.log1p(q) + fraction(q))


def old_flux(q):
    # Stable algebraic form of A*lambda*[1-(1+q)^(-2)].
    r = fraction(q)
    return A * LAMBDA * r * (2.0 - r)


def corrected_flux(q):
    return A * LAMBDA * fraction(q) ** 2


def correction(q):
    return 2.0 * A * fraction(q)


def retarded_derivative(function, q):
    # Fourth-order centered finite difference in affine time u,
    # independently using q(u+du)=q(u)*exp(lambda*du).
    h = 1.0e-3 / LAMBDA
    fm2, fm1, fp1, fp2 = [function(q * math.exp(LAMBDA * h * n))
                          for n in (-2, -1, 1, 2)]
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h)


def integrate_flux(flux, upper_q, old):
    # Integrate in t=log(1+q), a bounded interval even for very large q_s.
    # Integrand is the supplied flux times the full du/dt Jacobian.
    def integrand(t):
        q = math.expm1(t)
        if q == 0.0:
            return 2.0 * A if old else 0.0
        return flux(q) * (1.0 + q) / (LAMBDA * q)
    return quad(integrand, 0.0, math.log1p(upper_q),
                epsabs=1e-11, epsrel=1e-12, limit=100)


def fold_gap(s):
    # s=z-1; stable close to the critical field value.
    if s < 1e-3:
        return sum(((-1.0) ** n) * s**n / n for n in range(2, 16))
    return s - math.log1p(s)


def invert_gap(gap):
    if gap < 0:
        raise ValueError("No real dilaton inverse exists for W<1.")
    if gap == 0:
        return 1.0
    upper = max(1.0, 2.0 * gap + 2.0 * math.sqrt(2.0 * gap))
    while fold_gap(upper) < gap:
        upper *= 2.0
    s = brentq(lambda v: fold_gap(v) - gap, 0.0, upper,
               xtol=1e-15, rtol=1e-14)
    return 1.0 + s


def normalized_error(a, b, scale=1.0):
    return abs(a - b) / scale


def audit():
    records = []
    derivative_errors = []
    quadrature_errors = []
    endpoint_errors = []
    geometry_errors = []
    clip_threshold = brentq(
        lambda mu: mu - math.log1p(mu) - mu / (1.0 + mu), 0.01, 5.0)

    for mu in MUS:
        mass = A * mu
        qs = math.expm1(mu)
        d = -math.expm1(-mu)
        xs = qs / mu
        ys = mu / d
        # Retarded rays through the full pre-endpoint history.
        qvalues = [math.expm1(mu * t) for t in (0.001, 0.01, 0.1, 0.5, 0.9, 1.0)]
        max_old_deriv = max_corrected_deriv = max_delta_deriv = 0.0
        for q in qvalues:
            old_d = retarded_derivative(lambda z: old_charge(z, mu), q)
            new_d = retarded_derivative(lambda z: corrected_charge(z, mu), q)
            delta_d = retarded_derivative(correction, q)
            max_old_deriv = max(max_old_deriv, normalized_error(old_d, -old_flux(q), A*LAMBDA))
            max_corrected_deriv = max(max_corrected_deriv, normalized_error(new_d, -corrected_flux(q), A*LAMBDA))
            max_delta_deriv = max(max_delta_deriv, normalized_error(delta_d, old_flux(q)-corrected_flux(q), A*LAMBDA))
            # Verify finite correction, without deriving one charge from the other in code.
            endpoint_errors.append(normalized_error(
                corrected_charge(q, mu)-old_charge(q, mu), correction(q), A))

        integral_old, quad_old_err = integrate_flux(old_flux, qs, True)
        integral_new, quad_new_err = integrate_flux(corrected_flux, qs, False)
        quadrature_errors += [abs(integral_old - A*(mu+d))/A,
                              abs(integral_new - A*(mu-d))/A]
        endpoint_errors += [abs(old_charge(0, mu)-mass)/A,
                            abs(corrected_charge(0, mu)-mass)/A,
                            abs(old_charge(qs, mu)+A*d)/A,
                            abs(corrected_charge(qs, mu)-A*d)/A]

        # Endpoint consistency in dimensionless geometry.
        w_endpoint = xs*ys - math.log(xs*ys) - mu*(xs-1.0)
        geometry_errors += [abs((ys-mu)*xs-1), abs(w_endpoint-1)]

        qclip = brentq(lambda q: math.log1p(q)+fraction(q)-mu,
                      0, qs, xtol=1e-14, rtol=1e-14)
        clip_gap = mu-math.log1p(qclip)
        physical_clip_horizon = qclip/mu >= 1.0
        geometry_errors.append(abs(clip_gap-fraction(qclip)))
        clip_corrected = corrected_charge(qclip, mu)
        endpoint_errors.append(abs(clip_corrected-2*A*fraction(qclip))/A)

        derivative_errors += [max_old_deriv, max_corrected_deriv, max_delta_deriv]
        records.append({
            "dimensionless_mass_mu": mu,
            "incoming_mass_m": mass,
            "endpoint_q": qs,
            "endpoint_X": xs,
            "endpoint_Y": ys,
            "old_endpoint_charge": old_charge(qs, mu),
            "corrected_endpoint_charge": corrected_charge(qs, mu),
            "conventional_matter_thunderpop_energy": -A*d,
            "integrated_old_flux": integral_old,
            "integrated_corrected_flux": integral_new,
            "quadrature_reported_errors": [quad_old_err, quad_new_err],
            "old_flux_balance_scaled_error": abs(old_charge(qs,mu)+integral_old-mass)/A,
            "corrected_flux_balance_scaled_error": abs(corrected_charge(qs,mu)+integral_new-mass)/A,
            "old_charge_derivative_scaled_error": max_old_deriv,
            "corrected_charge_derivative_scaled_error": max_corrected_deriv,
            "correction_derivative_scaled_error": max_delta_deriv,
            "naive_old_charge_zero_rule": {
                "q_stop": qclip,
                "formal_horizon_X": qclip/mu,
                "horizon_after_incoming_shock": physical_clip_horizon,
                "formal_W_minus_1": clip_gap,
                "actual_postshock_horizon_gap": clip_gap if physical_clip_horizon else None,
                "corrected_charge_at_stop": clip_corrected,
                "meaning": "Negative control tied to old charge, not universal physical mass."
            }
        })

    # Conditional divergence as a horizon approaches z=1 from above.
    curvature = []
    for gap in (1e-2,1e-4,1e-6,1e-8,1e-10):
        z = invert_gap(gap)
        r = 4*LAMBDA**2*z/(z-1)
        curvature.append({"W_minus_1":gap, "z_minus_1":z-1,
                          "R_over_lambda_squared":r/LAMBDA**2,
                          "R_sqrt_gap_over_lambda_squared":r*math.sqrt(gap)/LAMBDA**2,
                          "gap_reconstruction_error":abs(fold_gap(z-1)-gap)})

    checks = {
        "both_charge_flux_derivative_laws": max(derivative_errors) < 1e-8,
        "both_charge_flux_integral_laws": max(quadrature_errors) < 1e-10,
        "initial_and_endpoint_charge_and_correction": max(endpoint_errors) < 1e-10,
        "endpoint_geometry_and_clip_identity": max(geometry_errors) < 1e-7,
        "horizon_curvature_increases_as_fold_is_approached": all(
            curvature[i+1]["R_over_lambda_squared"] > curvature[i]["R_over_lambda_squared"]
            for i in range(len(curvature)-1)),
        "horizon_curvature_asymptotic_coefficient": abs(
            curvature[-1]["R_sqrt_gap_over_lambda_squared"] - 2*math.sqrt(2)) < 5e-5,
    }
    result = {
        "scope": "Reproduction and audit of reduced semiclassical RST baseline; no new endpoint dynamics.",
        "parameters": {"kappa":KAPPA,"lambda":LAMBDA,"incoming_x_plus_0":X0,"energy_unit_A":A},
        "source_mapping": {
            "original_RST": "https://arxiv.org/pdf/hep-th/9206070 equations (2.1)-(2.9), (3.2)-(3.3)",
            "conventional_flux": "https://arxiv.org/pdf/1401.0590 page 6 before equation (10)",
            "corrected_charge": "https://arxiv.org/pdf/2510.18298v2 equations (33), (43), (A4)-(A6)",
            "notation_warning": "Our mu=4m/(kappa*lambda); 2026 paper's mu=kappa/2. They are different."
        },
        "max_errors": {"derivative_scaled":max(derivative_errors),
                       "quadrature_scaled":max(quadrature_errors),
                       "charge_scaled":max(endpoint_errors),
                       "geometry_absolute":max(geometry_errors)},
        "old_zero_rule_postshock_horizon_mass_threshold":clip_threshold,
        "checks":checks,"all_pass":all(checks.values()),
        "records":records,"conditional_curvature_sequence":curvature,
        "interpretation_limits": [
            "Different charge definitions require different fluxes; their mixing is inconsistent.",
            "Positive corrected mass does not remove the RST critical fold or prove a regular endpoint.",
            "The curvature obstruction assumes unchanged RST field map, bulk equation and Kruskal gauge.",
            "The 2026 corrected flux is a full charge-balance flux, not the conventional matter flux alone.",
            "The last-ray corrected charge is not itself a prediction of a stable remnant.",
            "The apparent horizon formula after the shock requires X>=1."
        ]
    }
    return result


if __name__ == "__main__":
    result = audit()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"output":str(output),"all_pass":result["all_pass"],
                      "checks":result["checks"],"max_errors":result["max_errors"]},indent=2))
    if not result["all_pass"]:
        raise SystemExit(1)
