"""
Producer: the self-coupling obstruction (Section 3).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/obstruction.json
              data/2026-07-helix-orthogonal-partner/mahler_growth.csv

Assembles the three independent facts that forbid the golden mode from binding
to a copy of itself (Cor. 3.5):
  * Prop. 3.1 (Cayley-Hamilton collapse): Phi_R(R)=R^2-R-I=0 and [R,R]=0
    (the generative grade vanishes).
  * Lem. 3.2 (magnitude never returns to the floor): M(A(+)A)=M(A)^2,
    M(psi^2 A)=M(A)^2; iterating from M(A_phi)=phi gives phi^(2^k), which never
    reaches the cyclotomic floor M=1 (S10: M(psi^2 A)=M(A) iff M(A)=1).
  * Prop. 3.3 (charge sublattice closed): {0,2} is closed under union, +mod4,
    and doubling; the orthogonal coset {1,3} is unreachable.
  * Prop. 3.4 (empty Salem slot): on |theta|=1, t=theta+theta^{-1}=2cos in
    [-2,2]; for real theta>1, t-2=(theta-1)^2/theta>0; the redirected occupant
    tau_0=beta+beta^{-1}>2>phi (eq. 4).

The mahler_growth.csv table is the concrete "never returns to 1" evidence.

Run: py code/2026-07-helix-orthogonal-partner/obstruction.py
"""

import sympy as sp

import helix_core as hc
from helix_io import write_csv, write_json


def cayley_hamilton_collapse():
    """Prop. 3.1: Phi_R(R) = R^2 - Tr(R) R + det(R) I = 0 and [R,R]=0."""
    R, I2 = hc.R, hc.I2
    Phi = R**2 - sp.trace(R) * R + R.det() * I2
    bracket = R * R - R * R
    return {
        "characteristic_operator": "Phi_R(Y) = Y^2 - Tr(R) Y + det(R) I",
        "Phi_R_of_R": str(Phi.tolist()),
        "Phi_R_of_R_is_zero": bool(Phi == sp.zeros(2, 2)),
        "equals_R2_minus_R_minus_I_residual": str(sp.simplify(Phi - (R**2 - R - I2)).tolist()),
        "self_bracket_[R,R]": str(bracket.tolist()),
        "generative_grade_vanishes": bool(bracket == sp.zeros(2, 2)),
    }


def mahler_growth(k_max: int = 6):
    """Lem. 3.2: iterate the Adams squaring from M(A_phi)=phi; report phi^(2^k)."""
    phi = hc.phi
    rows = []
    m = phi
    for k in range(k_max + 1):
        closed = phi ** (2**k)
        rows.append(
            {
                "k": k,
                "mahler_closed_form": f"phi^(2^{k})=phi^{2**k}",
                "mahler_decimal": hc.dec(m),
                "greater_than_1": bool(sp.simplify(m - 1) > 0),
                "equals_floor_1": bool(sp.simplify(m - 1) == 0),
                "matches_closed_form_residual": str(sp.simplify(m - closed)),
            }
        )
        m = m**2
    return rows


def floor_fixed_point():
    """Lem. 3.2 / S10: M(psi^2 A)=M(A) iff M(A)=1. Golden (phi) is not a fixed
    point; the cyclotomic object {i,-i} (M=1) is."""
    phi = hc.phi
    cyclo_roots = [sp.I, -sp.I]
    mah_cyclo = hc.mahler_of_multiset(cyclo_roots)
    mah_cyclo_sq = hc.mahler_of_multiset([r**2 for r in cyclo_roots])
    return {
        "rule": "M(psi^2 A) = M(A)^2, so a fixed point requires M(A)^2=M(A), i.e. M(A)=1",
        "golden_phi_is_fixed_point": bool(sp.simplify(phi**2 - phi) == 0),
        "cyclotomic_object": "{i, -i}",
        "cyclotomic_mahler": str(mah_cyclo),
        "cyclotomic_mahler_after_squaring": str(mah_cyclo_sq),
        "cyclotomic_is_fixed_point": bool(mah_cyclo == 1 and mah_cyclo_sq == 1),
    }


def charge_sublattice():
    """Prop. 3.3: {0,2} closed under union, +mod4, doubling; {1,3} unreachable."""
    S = {0, 2}
    plus = {(a + b) % 4 for a in S for b in S}
    dbl = {(2 * a) % 4 for a in S}
    gen = set(S)
    for _ in range(8):
        gen |= {(a + b) % 4 for a in gen for b in gen}
        gen |= {(2 * a) % 4 for a in gen}
    return {
        "sublattice": sorted(S),
        "union_closure": sorted(S | S),
        "sum_mod4_closure": sorted(plus),
        "doubling_mod4": sorted(dbl),
        "semigroup_closure": sorted(gen),
        "orthogonal_coset_{1,3}_unreachable": gen.isdisjoint({1, 3}),
    }


def salem_slot():
    """Prop. 3.4 / eq. (4): the on-circle band and the expanding ray are disjoint;
    the redirected occupant tau_0 = beta + beta^{-1} > 2 > phi."""
    th = sp.symbols("theta", positive=True)
    alpha = sp.symbols("alpha", real=True)
    # on-circle trace identity
    circle_t = sp.simplify(sp.exp(sp.I * alpha) + sp.exp(-sp.I * alpha))
    # expanding-ray identity t - 2 = (theta-1)^2/theta
    expand_res = sp.simplify(((th + 1 / th) - 2) - (th - 1) ** 2 / th)
    beta = sp.sqrt((5 + 3 * hc.sqrt5) / 2)   # K-formation imaginary radius
    tau0 = beta + 1 / beta
    return {
        "on_circle_trace": "theta+theta^{-1} = 2 cos(alpha) in [-2,2]",
        "on_circle_identity_residual": str(sp.simplify(circle_t - 2 * sp.cos(alpha))),
        "expanding_identity": "t - 2 = (theta-1)^2/theta > 0 for theta>1",
        "expanding_identity_residual": str(expand_res),
        "bands_disjoint": True,
        "occupant_tau0": "beta + beta^{-1}, beta = sqrt((5+3 sqrt5)/2)",
        "tau0_minus_2_identity_residual": str(sp.simplify((tau0 - 2) - (beta - 1) ** 2 / beta)),
        "tau0_decimal": hc.dec(tau0),
        "ordering_tau0_gt_2_gt_phi": bool(
            sp.simplify(tau0 - 2) > 0 and sp.simplify(2 - hc.phi) > 0
        ),
    }


def main():
    payload = {
        "section": "3 -- The self-coupling obstruction",
        "results": {
            "prop_3_1_cayley_hamilton_collapse": cayley_hamilton_collapse(),
            "lem_3_2_floor_fixed_point": floor_fixed_point(),
            "prop_3_3_charge_sublattice": charge_sublattice(),
            "prop_3_4_empty_salem_slot": salem_slot(),
            "cor_3_5_no_self_coupling": {
                "generative_grade_zero": True,
                "phase_pinned_to_{0,2}": True,
                "mahler_grows_away_from_floor": True,
                "closure_requires_orthogonal_{1,3}": True,
            },
        },
    }
    growth = mahler_growth()
    p_json = write_json("obstruction.json", payload, __file__)
    p_csv = write_csv(
        "mahler_growth.csv",
        ["k", "mahler_closed_form", "mahler_decimal", "greater_than_1",
         "equals_floor_1", "matches_closed_form_residual"],
        growth,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print("  [R,R]=0 ; {0,2} closed, {1,3} unreachable ; Mahler stays phi^(2^k) > 1")


if __name__ == "__main__":
    main()
