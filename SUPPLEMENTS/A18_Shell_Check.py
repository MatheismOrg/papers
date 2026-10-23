#!/usr/bin/env python3
"""A18 supplement (SO70 T166): checks Appendix A of A18, equations (1) to (4), with the standard library only.

(1) The interior Yukawa potential of a thin shell, by averaging the kernel -G exp(-m s)/s over the shell numerically
    and comparing with the closed form -G Ms exp(-m a) sinh(m r)/(m a r).
(2) Its radial derivative r Y'(r), by central finite differences of (1), and its sign (negative inside the shell).
    The point-source factor F = 1 + f(m0 r)/3 - 4 f(m2 r)/3, with f(x) = (1 + x) exp(-x), never exceeds 4/3 (checked on a grid).
(3) and (4) The central mass plus shell circular speed at r = 1, a = 2, m0 = 5, m2 = 1, Mc = 1, Ms = 100, G = 1, against the
    Newtonian speed of the same masses (inside the shell, the central mass alone): the ratio, about 3.307, exceeds 4/3.

usage: python3 A18_Shell_Check.py [OUT.json]   (default: A18_Shell_Check_Results.json beside this file)"""
import json, math, os, sys

G = 1.0


def f(x):
    return (1 + x) * math.exp(-x)


def F_point(r, m0, m2):
    return 1 + f(m0 * r) / 3 - 4 * f(m2 * r) / 3


def Y_closed(r, a, m, Ms):
    return -G * Ms * math.exp(-m * a) * math.sinh(m * r) / (m * a * r)


def rY_closed(r, a, m, Ms):
    return -G * Ms * math.exp(-m * a) * (m * r * math.cosh(m * r) - math.sinh(m * r)) / (m * a * r)


def Y_shell_average(r, a, m, Ms, n=200000):
    # Average of -G Ms exp(-m s)/s over the shell's directions: u = cos(angle) is uniform on [-1, 1],
    # s^2 = r^2 + a^2 - 2 r a u. Midpoint rule in u.
    total = 0.0
    for i in range(n):
        u = -1 + (i + 0.5) * 2 / n
        s = math.sqrt(r * r + a * a - 2 * r * a * u)
        total += math.exp(-m * s) / s
    return -G * Ms * total / n


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'A18_Shell_Check_Results.json')
    res = {'paper': 'A18, Appendix A', 'units': 'G = 1'}

    # (1) closed form against the numerical shell average
    eq1 = []
    for (r, a, m, Ms) in [(1.0, 2.0, 5.0, 100.0), (1.0, 2.0, 1.0, 100.0), (0.5, 2.0, 1.0, 1.0), (1.5, 2.0, 3.0, 1.0)]:
        num, cl = Y_shell_average(r, a, m, Ms), Y_closed(r, a, m, Ms)
        rel = abs(num - cl) / abs(cl)
        eq1.append({'r': r, 'a': a, 'm': m, 'Ms': Ms, 'numerical': num, 'closed_form': cl, 'relative_error': rel})
        assert rel < 1e-8, (r, a, m, rel)
    res['eq1_shell_average'] = eq1

    # (2) derivative by finite differences, and its sign inside the shell
    eq2 = []
    for (r, a, m, Ms) in [(1.0, 2.0, 5.0, 100.0), (1.0, 2.0, 1.0, 100.0), (0.3, 2.0, 2.0, 1.0)]:
        h = 1e-5
        fd = r * (Y_closed(r + h, a, m, Ms) - Y_closed(r - h, a, m, Ms)) / (2 * h)
        cl = rY_closed(r, a, m, Ms)
        rel = abs(fd - cl) / abs(cl)
        eq2.append({'r': r, 'a': a, 'm': m, 'Ms': Ms, 'finite_difference': fd, 'closed_form': cl, 'relative_error': rel})
        assert rel < 1e-6, (r, m, rel)
    neg = all(rY_closed(0.01 + 1.98 * k / 400, 2.0, m, 1.0) < 0 for k in range(401) for m in (0.1, 1.0, 5.0))
    assert neg
    res['eq2_derivative'] = eq2
    res['eq2_negative_inside_shell'] = neg

    # Point-source ceiling: f is nonincreasing from f(0) = 1, so F <= 4/3
    xs = [0.002 * k for k in range(5001)]
    f_nonincreasing = all(f(xs[k + 1]) <= f(xs[k]) for k in range(len(xs) - 1))
    Fmax = max(1 + f(x) / 3 - 4 * f(y) / 3 for x in xs[::10] for y in xs[::10])
    assert f_nonincreasing and Fmax <= 4 / 3 + 1e-12
    res['point_source'] = {'f_nonincreasing_on_0_to_10': f_nonincreasing, 'max_F_on_grid': Fmax, 'ceiling': 4 / 3}

    # (3) and (4): central mass plus outer shell
    r, a, m0, m2, Mc, Ms = 1.0, 2.0, 5.0, 1.0, 1.0, 100.0
    central = G * Mc / r * F_point(r, m0, m2)
    scalar = rY_closed(r, a, m0, Ms) / 3
    spin2 = -4 * rY_closed(r, a, m2, Ms) / 3
    vN2 = G * Mc / r
    ratio = (central + scalar + spin2) / vN2
    assert abs(ratio - 3.307) < 5e-4 and abs(spin2 - 3.319) < 5e-4 and abs(scalar + 0.045) < 5e-4 and ratio > 4 / 3
    res['eq4_instance'] = {'r': r, 'a': a, 'm0': m0, 'm2': m2, 'Mc': Mc, 'Ms': Ms,
                           'central_mass_term': central, 'shell_scalar_term': scalar, 'shell_spin2_term': spin2,
                           'v_c2_over_v_N2': ratio, 'exceeds_point_source_ceiling': ratio > 4 / 3,
                           'paper_values': {'ratio': 3.307, 'spin2_term': 3.319, 'scalar_term': -0.045}}
    res['all_checks_passed'] = True
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1)
    print(f'ratio {ratio:.6f}; spin-2 term {spin2:+.6f}; scalar term {scalar:+.6f}; central {central:.6f}; all checks passed')


if __name__ == '__main__':
    main()
