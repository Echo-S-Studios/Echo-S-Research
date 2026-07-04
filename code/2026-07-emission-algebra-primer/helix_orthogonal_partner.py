"""Producer: the three self-coupling obstructions and the orthogonal partner K.

Source paper: papers/2026-07-emission-algebra-primer/emission_algebra_primer.tex
Produces (Sec. 6 -- the three questions):
  * Prop 6.1-6.3  the three obstructions to self-coupling the golden mode
        (1) generative grade vanishes: Phi_R(R)=R^2-R-I=0, [R,R]=0
        (2) charge sublattice {0,2} closed, never reaches {1,3}
        (3) the Salem slot is empty (trace-down regimes disjoint)
  * Thm 6.6       the orthogonal partner K(x) = x^4 + 5x^2 - 5:
        real roots +/-K with K = 5^{1/4}/phi (terrain, inside unit circle),
        imaginary roots +/- i*beta (rotation), magnitude M(K) = beta^2 = phi^2 sqrt5,
        charge chi(K) = {0,1,2,3} -- completes Z/4Z
  * Prop 6.7      the parity criterion (even quartic w/ negative x^2-root)

Outputs:
  data/.../obstructions.json
  data/.../partner_K.json

Unlike tests/test_questions.py (asserts), this script solves the quartic
symbolically, splits the roots, and emits K's magnitude/charge record.
"""
from __future__ import annotations

import sympy as sp

from eap_core import (R, I2, phi, psi, sqrt5, comm, is_zero, mat_eq,
                      charge_one, charge_set, outside_unit)
from eap_io import write_json


def obstructions_payload():
    PhiR_at_R = R * R - sp.trace(R) * R + R.det() * I2
    # trace-down excess (theta - 1)^2/theta for real theta>1
    th = sp.symbols("theta", positive=True)
    excess = sp.simplify((th + 1 / th) - 2 - (th - 1) ** 2 / th)
    return {
        "obstruction_1_generative_grade": {
            "Phi_R(R)_eq_R^2-R-I": mat_to_rows_zero(PhiR_at_R),
            "is_zero": bool(mat_eq(PhiR_at_R, sp.zeros(2))),
            "[R,R]_is_zero": bool(mat_eq(comm(R, R), sp.zeros(2))),
            "reading": "self-application returns the keystone and annihilates; "
                       "no new direction is generated",
        },
        "obstruction_2_charge_closed": {
            "union_{0,2}": sorted({0, 2} | {0, 2}),
            "sum_mod4_{0,2}": sorted({(a + b) % 4 for a in (0, 2) for b in (0, 2)}),
            "double_{0,2}": sorted({(2 * a) % 4 for a in (0, 2)}),
            "reaches_{1,3}": bool(({0, 2} & {1, 3})),
            "reading": "the golden charge is trapped in the index-2 subgroup {0,2}",
        },
        "obstruction_3_salem_slot_empty": {
            "trace_down_excess_identity": "theta+1/theta-2 = (theta-1)^2/theta",
            "excess_is_zero_diff": bool(is_zero(excess)),
            "on_unit_circle": "t = 2 cos(a) in [-2, 2]",
            "real_growth": "t > 2 for real theta > 1 (disjoint regimes)",
            "reading": "every golden combination stays real; no complex "
                       "unit-circle (Salem) eigenvalue can arise",
        },
    }


def mat_to_rows_zero(M):
    """String form of a matrix that should be the zero matrix (for the record)."""
    M = sp.Matrix(M)
    return [[sp.sstr(sp.simplify(M[i, j])) for j in range(M.cols)] for i in range(M.rows)]


def partner_K_payload():
    x, y = sp.symbols("x y")
    # y = x^2 : y^2 + 5y - 5
    yroots = sp.solve(y ** 2 + 5 * y - 5, y)
    yp = sp.nsimplify((-5 + 3 * sqrt5) / 2)     # > 0
    ym = sp.nsimplify((-5 - 3 * sqrt5) / 2)     # < 0
    K = sp.sqrt(yp)
    beta = sp.sqrt(-ym)
    roots = [K, -K, sp.I * beta, -sp.I * beta]

    # magnitude: only the two imaginary roots lie outside the unit circle
    out = outside_unit(roots)
    M_exact = sp.nsimplify(sp.simplify(sp.prod([sp.Abs(v) for v in out])))

    # parity criterion: several even quartics with a negative y-root
    parity_examples = []
    for (b, c) in [(5, -5), (1, -3), (2, -7), (4, -1)]:
        yr = sp.solve(sp.symbols("y") ** 2 + b * sp.symbols("y") + c, sp.symbols("y"))
        negs = [r for r in yr if sp.im(r) == 0 and r < 0]
        charges = sorted({charge_one(sp.I * sp.sqrt(-yv)) for yv in negs})
        parity_examples.append({"quartic": f"x^4+{b}x^2+({c})",
                                "imaginary_root_charges": charges})

    return {
        "K": "x^4 + 5 x^2 - 5",
        "substitution": "y = x^2  ->  y^2 + 5y - 5",
        "y_roots_straddle_zero": {
            "y_plus": sp.sstr(yp), "y_plus_decimal": round(float(yp), 4),
            "y_minus": sp.sstr(ym), "y_minus_decimal": round(float(ym), 4),
        },
        "real_roots_terrain": {
            "value_pm_K": sp.sstr(sp.simplify(K)),
            "equals_5^(1/4)/phi": bool(is_zero(sp.simplify(K - 5 ** sp.Rational(1, 4) / phi))),
            "K_decimal": round(float(K), 4),
            "inside_unit_circle": bool(float(K) < 1),
        },
        "imaginary_roots_rotation": {
            "value_pm_i_beta": sp.sstr(sp.simplify(beta)),
            "beta_decimal": round(float(beta), 4),
        },
        "magnitude": {
            "num_outside_unit_circle": len(out),
            "M_K": sp.sstr(M_exact),                     # (5+3sqrt5)/2
            "equals_phi^2_sqrt5": bool(is_zero(sp.simplify(M_exact - phi ** 2 * sqrt5))),
            "equals_beta^2": bool(is_zero(sp.simplify(M_exact - beta ** 2))),
            "M_K_decimal": round(float(M_exact), 4),
        },
        "charge_completion": {
            "chi_K": charge_set(roots),                  # {0,1,2,3}
            "real_roots_supply": charge_set([K, -K]),    # {0,2}
            "imaginary_roots_supply": charge_set([sp.I * beta, -sp.I * beta]),  # {1,3}
            "completes_Z4": bool(charge_set(roots) == [0, 1, 2, 3]),
        },
        "parity_criterion": {
            "statement": "an even real quartic (polynomial in x^2) with a negative "
                         "x^2-root places its complex roots on the imaginary axis "
                         "(charges {1,3}); it is PARITY, not the field, that singles out K",
            "even_quartic_examples": parity_examples,
            "non_even_example": "x^4 + x + 1 places a complex root off-axis (charge-lattice miss)",
        },
    }


def main():
    p1 = write_json("obstructions.json", obstructions_payload(), __file__)
    p2 = write_json("partner_K.json", partner_K_payload(), __file__)
    print(f"wrote {p1}")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
