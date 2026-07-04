r"""
Producer: Section 2 "Learning as exact basis growth" of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits the concrete Section-2 artifacts:
  * eq (Gpower): the power-basis trace-form Gram of K = Q(sqrt2+sqrt3)
      -> data/trace_gram_powerbasis.csv
  * Prop 2.1: the exact rational projector P = B(B^T G B)^{-1} B^T G for B={1,theta}
      -> data/projector_B_one_theta.json
  * Ex 2.4 (a worked episode): w = 2 sqrt6 = theta^2 - 5, residual, ||r||_G^2 = 96,
      seed minpoly x^2-24
      -> data/worked_episode.json
  * Ex 2.6 (a hand-checkable witness link): the genesis SHA-256 witness digest
      31f1f1e05ac9a35a and the canonical body
      -> data/witness_digest.json

Run:  py code/2026-06-residual-return-learning/learning_and_witness.py

The trace Gram is REBUILT from the regular representation of theta (power sums),
and the seed minimal polynomial is computed by the paper's largest-invariant-factor
bridge -- the producers compute, they do not assert.
"""
import sympy as sp
from sympy import Matrix, eye

import rrl_core as core

SCRIPT = "learning_and_witness.py"
M_THETA = [1, 0, -10, 0, 1]          # theta = sqrt2+sqrt3, m_theta = x^4-10x^2+1
x = core.X


def produce_gram():
    """eq (Gpower): G_{ij} = Tr(theta^{i-1}theta^{j-1}), rebuilt from power traces."""
    G = core.trace_form_powerbasis(M_THETA)
    minors = [int(G[:k, :k].det()) for k in range(1, 5)]   # Sylvester: all > 0
    rows = [[i + 1] + [int(G[i, j]) for j in range(4)] for i in range(4)]
    core.write_csv(SCRIPT, "trace_gram_powerbasis.csv",
                   ["row", "c1", "c2", "c3", "c4"], rows)
    return G, minors


def produce_projector(G):
    """Prop 2.1: exact projector onto col(B) for B = {1, theta}."""
    B = Matrix([[1, 0], [0, 1], [0, 0], [0, 0]])          # columns 1, theta
    P = core.projector(B, G)
    payload = {
        "field": "Q(sqrt2+sqrt3) = Q[x]/(x^4-10x^2+1)",
        "forced_basis_B": "columns {1, theta}",
        "projector_P": core.mat_to_strlist(P),
        "idempotent_P2_eq_P": bool(P * P == P),
        "fixes_columns_PB_eq_B": bool(P * B == B),
        "largest_numerator": int(max(abs(sp.fraction(P[i, j])[0])
                                     for i in range(4) for j in range(4))),
        "largest_denominator": int(max(abs(sp.fraction(P[i, j])[1])
                                       for i in range(4) for j in range(4))),
    }
    core.write_json(SCRIPT, "projector_B_one_theta.json", payload)
    return P


def produce_episode(G, P):
    """Ex 2.4: w = 2 sqrt6 = theta^2 - 5, coords (-5,0,1,0); x = 1 + w; r = w;
    ||r||_G^2 = Tr(24) = 96; seed minpoly x^2-24 by the bridge."""
    xvec = Matrix([-4, 0, 1, 0])                 # x = 1 + w = theta^2 - 4
    r = xvec - P * xvec
    score = int((r.T * G * r)[0])
    seed_coords = [-5, 0, 1, 0]                  # w = theta^2 - 5
    minpoly, is_int = core.coords_to_minpoly(seed_coords, M_THETA)
    # G-orthogonality of w to col(B) = span{1,theta}
    w = Matrix(seed_coords)
    orth = [int((Matrix(col).T * G * w)[0]) for col in
            ([1, 0, 0, 0], [0, 1, 0, 0])]
    payload = {
        "off_axis_element": "w = 2*sqrt6 = theta^2 - 5",
        "observation_x": "x = 1 + w = theta^2 - 4",
        "x_coords": [-4, 0, 1, 0],
        "residual_r_coords": [int(v) for v in r],
        "residual_equals_w": bool(r == w),
        "residual_norm_G": score,                # == 96
        "residual_norm_via_trace_Tr24": int((24 * eye(4)).trace()),
        "G_orthogonal_to_col_B": orth,           # [0, 0]
        "seed_coords": seed_coords,
        "seed_minpoly": core.poly_str(minpoly),  # x^2 - 24
        "seed_minpoly_coeffs": core.poly_coeffs(minpoly),
        "seed_is_algebraic_integer": bool(is_int),
    }
    core.write_json(SCRIPT, "worked_episode.json", payload)
    return payload


def produce_witness():
    """Ex 2.6: the genesis growth body and its SHA-256 witness digest."""
    body = {
        "coords": ["-5", "0", "1", "0"],
        "event": "basis_growth",
        "index": 0,
        "min_poly": [1, 0, -24],
        "num_seeds": 3,
        "prev_hash": "genesis",
        "snap": "exact",
        "streak": 4,
    }
    canonical, digest = core.witness_digest("genesis", body)
    payload = {
        "chain_rule": "h_k = SHA-256(h_{k-1} || canonical_json(body))[0:16]",
        "prev_hash": "genesis",
        "body": body,
        "canonical_json": canonical,
        "digest": digest,                        # == 31f1f1e05ac9a35a
    }
    core.write_json(SCRIPT, "witness_digest.json", payload)
    return payload


def main():
    G, minors = produce_gram()
    P = produce_projector(G)
    ep = produce_episode(G, P)
    wt = produce_witness()
    print(f"[{SCRIPT}] Gram Sylvester minors (all>0): {minors}")
    print(f"[{SCRIPT}] residual norm ||r||_G^2 = {ep['residual_norm_G']} "
          f"(seed {ep['seed_minpoly']})")
    print(f"[{SCRIPT}] witness digest = {wt['digest']}")
    print(f"[{SCRIPT}] wrote 4 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
