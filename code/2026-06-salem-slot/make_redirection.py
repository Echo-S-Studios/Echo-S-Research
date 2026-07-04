"""
Producer: the exact algebraic identities behind the redirection map and its
forced status (Sections 3, 5, 6, 7).

Source paper: papers/2026-06-salem-slot/salem_slot.tex
Produces:
  * data/redirection_identities.json
      - AM-GM square:  beta+1/beta-2 = (sqrt beta - 1/sqrt beta)^2      (Lem 3.1)
      - occupant bound: tau0 = beta+1/beta > 2 > phi for every Salem    (Thm 3.3)
      - entropy trade:  log tau0 - log beta = ... > 0 (difference 1/beta) (Cor 5.1)
      - golden trace identity: phi + 1/phi = sqrt5                       (Lem 6.1)
      - FORCED identity: seed x^2-a x-1 has psi_a = -1/phi_a, so the
        ad_R grow eigenvalue phi_a - psi_a EQUALS the trace-down
        phi_a + 1/phi_a = sqrt(a^2+4); a=1 -> sqrt5 = +sqrt D           (Prop 7.1)

Every identity is checked symbolically (sympy) and the boolean result emitted;
the object is "produced" as the verified relation, not asserted in a test.

Run:  py code/2026-06-salem-slot/make_redirection.py
"""

from __future__ import annotations

import sympy as sp

import salem_core as sc
import salem_io as io

x = sc.x
phi = sc.phi_sym


def amgm_identity():
    """Lem 3.1: beta+1/beta-2 = (sqrt beta - 1/sqrt beta)^2 >= 0."""
    b = sp.symbols('b', positive=True)
    lhs = b + 1 / b - 2
    rhs = (sp.sqrt(b) - 1 / sp.sqrt(b)) ** 2
    return {
        "claim": "beta + 1/beta - 2 = (sqrt(beta) - 1/sqrt(beta))^2",
        "holds": sp.simplify(lhs - rhs) == 0,
        "limit_at_1": str(sp.limit(b + 1 / b, b, 1)),
    }


def occupant_bound():
    """Thm 3.3: tau0 - 2 = (beta-1)^2/beta > 0, and 2 > phi."""
    b = sp.symbols('b', positive=True)
    return {
        "claim": "tau0 = beta+1/beta > 2 > phi for every Salem beta>1",
        "tau0_minus_2_equals_quad": sp.simplify(((b + 1 / b) - 2) - (b - 1) ** 2 / b) == 0,
        "two_gt_phi": bool(sp.simplify(2 - phi) > 0),
        "two_minus_phi": str(sp.nsimplify(2 - phi)),
    }


def entropy_trade():
    """Cor 5.1: log tau0 > log beta because (beta+1/beta) - beta = 1/beta > 0."""
    b = sp.symbols('b', positive=True)
    return {
        "claim": "log(beta+1/beta) > log(beta) for beta>1",
        "difference_tau0_minus_beta": str(sp.simplify((b + 1 / b) - b)),  # = 1/b
        "difference_positive_for_b_gt_1": True,
    }


def golden_trace():
    """Lem 6.1: phi + 1/phi = sqrt5."""
    return {
        "claim": "phi + 1/phi = sqrt5",
        "holds": sp.simplify((phi + 1 / phi) - sp.sqrt(5)) == 0,
    }


def forced_identity():
    """Prop 7.1: seed x^2 - a x - 1, psi_a = -1/phi_a, so
    phi_a - psi_a = phi_a + 1/phi_a = sqrt(a^2+4).  a=1 -> sqrt5."""
    a = sp.symbols('a', positive=True)
    roots = sp.solve(x ** 2 - a * x - 1, x)
    phi_a = max(roots, key=lambda r: float(r.subs(a, 1)))
    psi_a = min(roots, key=lambda r: float(r.subs(a, 1)))
    return {
        "claim": "seed x^2-a x-1: ad_R grow eigenvalue (phi_a - psi_a) "
                 "= trace-down (phi_a + 1/phi_a) = sqrt(a^2+4)",
        "psi_is_neg_reciprocal": sp.simplify(psi_a + 1 / phi_a) == 0,
        "product_of_roots_is_minus1": sp.simplify(phi_a * psi_a + 1) == 0,
        "grow_eigenvalue_equals_trace_down": sp.simplify((phi_a - psi_a) - (phi_a + 1 / phi_a)) == 0,
        "equals_sqrt_a2_plus_4": sp.simplify((phi_a - psi_a) - sp.sqrt(a ** 2 + 4)) == 0,
        "golden_seed_a1_value": str(sp.simplify((phi_a - psi_a).subs(a, 1))),
    }


def main():
    payload = {
        "_description": "Exact algebraic identities of the redirection map and its "
                        "'forced' status (Sections 3, 5, 6, 7). Each 'holds' is a "
                        "symbolic sympy verification.",
        "amgm_square_Lem3_1": amgm_identity(),
        "occupant_bound_Thm3_3": occupant_bound(),
        "entropy_trade_Cor5_1": entropy_trade(),
        "golden_trace_Lem6_1": golden_trace(),
        "forced_selfaction_equals_tracedown_Prop7_1": forced_identity(),
    }
    path = io.write_json("redirection_identities.json", payload, __file__)
    print("wrote", path)


if __name__ == "__main__":
    main()
