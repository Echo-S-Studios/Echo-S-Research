"""
PRODUCER: trace-form projector, residual, capture, and the natural-gradient fit.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex
  - Sec. 3.1: Def. 3.1 (projector/residual), Prop. 3.2 (idempotence,
    orthogonality, capture), Ex. 3.4 (golden projector), Ex. 3.5 (Q(i) Hermitian)
  - Sec. 4.1: Eq. (17) Pythagorean split ||x||^2 = ||Px||^2 + ||r||^2
  - Sec. 7.6: Prop. 7.13 / Ex. 7.15 (the projector is a one-step Newton /
    natural-gradient fit, a* = (B^T G B)^{-1} B^T G x = 5/3)

Emits:
  data/projector_residual.json
"""
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "projector_residual.py"

G5 = sp.Matrix([[2, 1], [1, 3]])                 # golden trace form (rebuilt below)


def golden_projector():
    """Ex. 3.4: B = (0,1)^T (line Q.phi), B^T G B = 3, P = (1/3)[[0,0],[1,3]];
    idempotent and G-orthogonal residual."""
    C = vc.companion_from_poly(x**2 - x - 1)
    G = vc.gram([(1, 0), (0, 1)], C)             # == G5, rebuilt from scratch
    B = sp.Matrix([0, 1])
    P = vc.projector(B, G)
    r = vc.residual((-1, 2), B, G)
    return {
        "G_rebuilt": vc.mat_to_list(G),
        "B": vc.mat_to_list(B),
        "BtGB": vc.sval((B.T * G * B)[0]),
        "P": vc.mat_to_list(P),
        "idempotent_P2_eq_P": sp.simplify(P * P - P) == sp.zeros(2, 2),
        "orthogonal_BtGr_eq_0": (B.T * G * r) == sp.Matrix([0]),
    }


def capture_and_novelty():
    """Ex. 3.4: 2phi=(0,2) captured (r=0); sqrt5=(-1,2) projects to (0,5/3),
    residual (-1,1/3), ||r||_G^2 = 5/3."""
    B = sp.Matrix([0, 1])
    P = vc.projector(B, G5)
    r_2phi = vc.residual((0, 2), B, G5)
    Px = P * sp.Matrix([-1, 2])
    r_s5 = vc.residual((-1, 2), B, G5)
    return {
        "2phi_captured": {
            "coord": [0, 2],
            "P_x": vc.mat_to_list(P * sp.Matrix([0, 2])),
            "residual": vc.mat_to_list(r_2phi),
            "residual_norm2": vc.sval(vc.gnorm2(r_2phi, G5)),
        },
        "sqrt5_novelty": {
            "coord": [-1, 2],
            "P_x": vc.mat_to_list(Px),
            "residual": vc.mat_to_list(r_s5),
            "residual_norm2": vc.sval(vc.gnorm2(r_s5, G5)),
        },
    }


def pythagoras_split():
    """Eq. (17): ||x||_G^2 = ||Px||_G^2 + ||r||_G^2 on x = sqrt5 vs line Q.phi."""
    B = sp.Matrix([0, 1])
    P = vc.projector(B, G5)
    xv = sp.Matrix([-1, 2])
    Px = P * xv
    r = xv - Px
    return {
        "x_norm2": vc.sval(vc.gnorm2(xv, G5)),
        "Px_norm2": vc.sval(vc.gnorm2(Px, G5)),
        "r_norm2": vc.sval(vc.gnorm2(r, G5)),
        "identity_holds": vc.gnorm2(xv, G5) == vc.gnorm2(Px, G5) + vc.gnorm2(r, G5),
    }


def complex_projector_Qi():
    """Ex. 3.5: Q(i), G_2 = 2I, B = <1> = (1,0), x = i = (0,1): P x = 0, r = i,
    ||r||_{G2}^2 = 2 (i is pure off-axis novelty against <1>)."""
    G2 = sp.Matrix([[2, 0], [0, 2]])
    B = sp.Matrix([1, 0])
    xi = sp.Matrix([0, 1])
    P = vc.projector(B, G2)
    r = xi - P * xi
    return {
        "G2": vc.mat_to_list(G2),
        "BtG2B": vc.sval((B.T * G2 * B)[0]),
        "BtG2x": vc.sval((B.T * G2 * xi)[0]),
        "P_x": vc.mat_to_list(P * xi),
        "residual": vc.mat_to_list(r),
        "residual_norm2": vc.sval(vc.gnorm2(r, G2)),
    }


def natural_gradient():
    """Prop. 7.13 / Ex. 7.15: single Newton step a* = (B^T G B)^{-1} B^T G x lands
    at B a* = P x; in Q(sqrt5), B = <phi>, x = sqrt5 gives a* = 5/3."""
    B = sp.Matrix([0, 1])
    xv = sp.Matrix([-1, 2])
    astar = (B.T * G5 * B).inv() * (B.T * G5 * xv)
    P = vc.projector(B, G5)
    grad0 = B.T * G5 * (B * sp.Matrix([0]) - xv)
    hess = B.T * G5 * B
    return {
        "a_star": vc.sval(astar[0]),
        "B_astar_equals_Px": B * astar == P * xv,
        "residual_matches": (xv - B * astar) == vc.residual((-1, 2), B, G5),
        "newton_step_matches": (-hess.inv() * grad0) == astar,
        "euclidean_BtB": vc.sval((B.T * B)[0]),        # 1 (would land off-projection)
        "trace_BtGB": vc.sval((B.T * G5 * B)[0]),       # 3
    }


def main():
    payload = {
        "golden_projector": golden_projector(),
        "capture_and_novelty": capture_and_novelty(),
        "pythagoras_split": pythagoras_split(),
        "complex_projector_Qi": complex_projector_Qi(),
        "natural_gradient": natural_gradient(),
    }
    p = vc.write_json("projector_residual.json", payload, SCRIPT)
    print(f"wrote {p}")
    print("  sqrt5 residual norm^2 =",
          payload["capture_and_novelty"]["sqrt5_novelty"]["residual_norm2"],
          "; Newton a* =", payload["natural_gradient"]["a_star"])
    print("  Pythagoras identity holds:", payload["pythagoras_split"]["identity_holds"])


if __name__ == "__main__":
    main()
