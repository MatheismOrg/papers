#!/usr/bin/env python3
"""
Minimal Physics Specifications
Symbolic summaries for QM, GR, and the minimal Standard Model.
These strings require an external mathematical interpretation, not supplied here.
Each function returns descriptive data, not a numerical solver or complete physical model.
"""

def quantum_mechanics():
    return {
        "state_space": "complex Hilbert space H",
        "state": "normalized vector |psi> in H, or density matrix rho",
        "evolution": "i hbar d|psi>/dt = H |psi>",
        "observable": "self-adjoint operator A on H",
        "measurement": "eigenvalues a_i, probability |<a_i|psi>|^2",
        "density_evolution": "drho/dt = -i/hbar [H, rho]",
        "free_inputs": ["hbar (Planck constant)"],
    }

def general_relativity():
    return {
        "manifold": "4D pseudo-Riemannian (M, g), signature (-,+,+,+)",
        "field": "metric tensor g_{mu nu}(x)",
        "action": "S = (1/16piG) integral sqrt(-g) (R - 2 Lambda) d4x + S_matter",
        "equation": "G_{mu nu} + Lambda g_{mu nu} = 8piG T_{mu nu}",
        "geodesic": "d2x^mu/dtau2 + Gamma^mu_{alpha beta} dx^alpha/dtau dx^beta/dtau = 0",
        "free_inputs": ["G (Newton constant)", "Lambda (cosmological constant)"],
    }

def standard_model():
    gauge = "SU(3)_c x SU(2)_L x U(1)_Y"
    fermions = []
    for gen in range(1, 4):
        fermions.append({
            "generation": gen,
            "quarks": {
                "Q_L": {"SU3": "3", "SU2": "2", "Y": "1/6"},
                "u_R": {"SU3": "3", "SU2": "1", "Y": "2/3"},
                "d_R": {"SU3": "3", "SU2": "1", "Y": "-1/3"},
            },
            "leptons": {
                "L_L": {"SU3": "1", "SU2": "2", "Y": "-1/2"},
                "e_R": {"SU3": "1", "SU2": "1", "Y": "-1"},
            },
        })
    higgs = {"SU3": "1", "SU2": "2", "Y": "1/2", "potential": "V = -mu^2 |phi|^2 + lambda |phi|^4"}
    lagrangian_terms = [
        "-1/4 F^a_{mu nu} F^{a mu nu}  (gauge kinetic, each factor)",
        "i psi_bar gamma^mu D_mu psi  (fermion kinetic + gauge coupling)",
        "|D_mu phi|^2 - V(phi)  (Higgs kinetic + potential)",
        "-y_ij psi_bar_L phi psi_R + h.c.  (Yukawa couplings)",
    ]
    free_inputs = [
        "g1, g2, g3 (gauge couplings)",
        "mu^2, lambda (Higgs potential parameters)",
        "Y_u (complex up Yukawa matrix; physical counts require all shared field redefinitions)",
        "Y_d (3x3 down Yukawa matrix)",
        "Y_e (3x3 charged lepton Yukawa matrix)",
        "theta_QCD (CP-violating phase)",
    ]
    return {
        "gauge_group": gauge,
        "fermion_content": fermions,
        "higgs": higgs,
        "lagrangian": lagrangian_terms,
        "free_inputs": free_inputs,
        "model_variant": "minimal Standard Model with massless neutrinos, including theta_QCD",
        "total_free_parameters": 19,
    }

if __name__ == "__main__":
    import json, sys
    for name, func in [("QM", quantum_mechanics), ("GR", general_relativity), ("SM", standard_model)]:
        spec = func()
        txt = json.dumps(spec, indent=None, separators=(',',':'))
        print(f"\n{'='*60}")
        print(f"{name}: {len(txt)} chars, {len(txt.encode('utf-8'))} bytes, {len(txt.encode('utf-8'))*8} bits")
        print(f"{'='*60}")
        print(json.dumps(spec, indent=2))
