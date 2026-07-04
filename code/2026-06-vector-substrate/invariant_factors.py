"""
PRODUCER: invariant factors, Smith normal form, and the exact seed bridge.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex
  - Sec. 2.3: Cor. 2.9, Thm. 2.10, Ex. 2.11 (charpoly insufficient),
              Ex. 2.13 (Jordan needs full list)
  - Sec. 4.3: Thm. 4.2 / Ex. 4.4 (coordinates -> minimal polynomial bridge:
              m_alpha = largest invariant factor of xI - rho(alpha))

Invariant factors are produced independently of any library SNF, via
determinantal divisors (monic gcds of i x i minors of xI - A).

Emits:
  data/invariant_factors.json -- the bridge for 2sqrt6, the two non-similarity
                                  separations, and the completeness property.
"""
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "invariant_factors.py"


def bridge_2sqrt6():
    """Thm. 4.2 / Ex. 4.4: rho(2sqrt6)=rho(theta)^2-5I, its invariant factors are
    (x^2-24, x^2-24), largest = m_{2sqrt6} = x^2-24 (monic over Z => admitted)."""
    C = vc.companion_from_poly(x**4 - 10 * x**2 + 1)   # theta = sqrt2 + sqrt3
    R = C**2 - 5 * sp.eye(4)                            # = rho(2sqrt6), coords (-5,0,1,0)
    ifs = vc.invariant_factors(R)
    minp = sp.expand(vc.matrix_minpoly(R))
    lib_minp = sp.minimal_polynomial(2 * sp.sqrt(6), x)
    return {
        "coord": [-5, 0, 1, 0],
        "matrix": vc.mat_to_list(R),
        "matches_rho_of_coord": R == vc.rho((-5, 0, 1, 0), C),
        "invariant_factors": [vc.sval(f) for f in ifs],
        "minpoly_from_matrix": vc.sval(minp),
        "minpoly_from_sympy": vc.sval(sp.expand(lib_minp)),
        "largest_invariant_factor": vc.sval(ifs[-1]),
        "product_if_equals_charpoly": sp.expand(sp.prod(ifs)) == sp.expand(R.charpoly(x).as_expr()),
        "square_is_24": bool(sp.simplify((2 * sp.sqrt(6))**2 - 24) == 0),
        "monic_integer_admissible": True,
        "n_over_d": 2,
    }


def charpoly_insufficient():
    """Ex. 2.11: rho(phi)(+)rho(phi) vs C((x^2-x-1)^2): same charpoly (x^2-x-1)^2
    and trace 2, but invariant factors (x^2-x-1, x^2-x-1) vs single (x^2-x-1)^2
    -> NOT similar.  Independent witness: annihilation by x^2-x-1."""
    Cphi = vc.companion_from_poly(x**2 - x - 1)
    blk = sp.diag(Cphi, Cphi)
    comp = vc.companion_from_poly((x**2 - x - 1)**2)
    return {
        "block_sum": {
            "matrix": vc.mat_to_list(blk),
            "charpoly_factored": sp.sstr(sp.factor(blk.charpoly(x).as_expr())),
            "trace": vc.sval(blk.trace()),
            "invariant_factors": [vc.sval(f) for f in vc.invariant_factors(blk)],
            "annihilated_by_x2mxm1": blk * blk - blk - sp.eye(4) == sp.zeros(4, 4),
        },
        "companion_of_square": {
            "matrix": vc.mat_to_list(comp),
            "charpoly_factored": sp.sstr(sp.factor(comp.charpoly(x).as_expr())),
            "trace": vc.sval(comp.trace()),
            "invariant_factors": [vc.sval(f) for f in vc.invariant_factors(comp)],
            "annihilated_by_x2mxm1": comp * comp - comp - sp.eye(4) == sp.zeros(4, 4),
        },
        "similar": False,
    }


def jordan_full_list():
    """Ex. 2.13: J2(+)J2 vs J2(+)J1(+)J1 share charpoly x^4 and minpoly x^2 but
    have invariant factors (x^2,x^2) vs (x,x,x^2).  Independent witness: rank."""
    J2 = sp.Matrix([[0, 1], [0, 0]])
    J1 = sp.Matrix([[0]])
    A1 = sp.diag(J2, J2)
    A2 = sp.diag(J2, J1, J1)
    return {
        "A1_J2plusJ2": {
            "charpoly": vc.sval(sp.expand(A1.charpoly(x).as_expr())),
            "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(A1))),
            "invariant_factors": [vc.sval(f) for f in vc.invariant_factors(A1)],
            "rank": int(A1.rank()),
        },
        "A2_J2plusJ1plusJ1": {
            "charpoly": vc.sval(sp.expand(A2.charpoly(x).as_expr())),
            "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(A2))),
            "invariant_factors": [vc.sval(f) for f in vc.invariant_factors(A2)],
            "rank": int(A2.rank()),
        },
        "similar": False,
    }


def completeness_derogatory():
    """Thm. 2.10: product of invariant factors = charpoly; largest = minpoly.
    Shown on the derogatory rho(3) = 3I in Q(sqrt5)."""
    C = vc.companion_from_poly(x**2 - x - 1)
    R = vc.rho((3, 0), C)
    ifs = vc.invariant_factors(R)
    return {
        "matrix": vc.mat_to_list(R),
        "invariant_factors": [vc.sval(f) for f in ifs],
        "product": vc.sval(sp.expand(sp.prod(ifs))),
        "charpoly": vc.sval(sp.expand(R.charpoly(x).as_expr())),
        "largest_if": vc.sval(ifs[-1]),
        "minpoly": vc.sval(sp.expand(vc.matrix_minpoly(R))),
    }


def main():
    payload = {
        "bridge_2sqrt6": bridge_2sqrt6(),
        "charpoly_insufficient": charpoly_insufficient(),
        "jordan_full_list": jordan_full_list(),
        "completeness_derogatory": completeness_derogatory(),
    }
    p = vc.write_json("invariant_factors.json", payload, SCRIPT)
    print(f"wrote {p}")
    print("  bridge 2sqrt6 invariant factors =",
          payload["bridge_2sqrt6"]["invariant_factors"],
          "-> minpoly", payload["bridge_2sqrt6"]["largest_invariant_factor"])
    print("  Jordan IFs:", payload["jordan_full_list"]["A1_J2plusJ2"]["invariant_factors"],
          "vs", payload["jordan_full_list"]["A2_J2plusJ1plusJ1"]["invariant_factors"])


if __name__ == "__main__":
    main()
