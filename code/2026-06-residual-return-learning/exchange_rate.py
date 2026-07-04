r"""
Producer: Section 4.3 "The exchange rate, derived: lambda = 2c" and the
self-action / frame-shift material of

    "Residual Return: Exact Learning Dynamics and Language over the Vector
     Substrate" (papers/2026-06-residual-return-learning/residual_return_learning.tex)

Emits:
  * Thm 4.5 (lambda=2c) + Rem 4.6/4.10 (the two floors from mu_S): the Smyth
      plastic constant mu_S (root of x^3-x-1) and the floors
      2 log mu_S = 0.5624, 4x = 2.2496, 8x = 4.4992
      -> data/smyth_floors.json
  * Ex 4.3b (certified GROW): rigorous mpmath verified-interval enclosures of
      log 7 and log 24, and the certified verdicts 56 >= 2 log7, 96 >= 2 log24
      -> data/certified_grow.json
  * Prop 4.7 (self-action spectrum of the C-ladder): the 4x4 adjoint [R_C,.]
      spectrum {-sqrt(1+4C), 0, +sqrt(1+4C)} for C in {1/4,1/2,1,3}
      -> data/selfaction_spectrum.csv
  * Prop 4.7 aside (corrected 2026-07-04): R_1 = [[0,1],[1,-1]] is conjugate to
      -R (shared charpoly x^2+x-1), NOT to the keystone R (charpoly x^2-x-1);
      both share the self-action gap sqrt5 (discriminant 5)
      -> data/R1_similarity.json
  * Def 4.8 / Table 4 (frame-shift canonicalization): gap sqrt(1+4C) and
      c = sqrt(1+4C)/(2C) on the gates C in {1/4,1/2,1}
      -> data/canonicalization_table.csv

Run:  py code/2026-06-residual-return-learning/exchange_rate.py

Everything algebraic is exact sympy; only the transcendental log/mu_S use mpmath
(mu_S) and mpmath.iv verified interval arithmetic (the certified enclosures).
"""
import mpmath as mp
import sympy as sp
from sympy import Matrix, Rational as Q, sqrt, symbols

import rrl_core as core

SCRIPT = "exchange_rate.py"
x = core.X


def adjoint_matrix(Rmat):
    """4x4 matrix of the Lie self-action X -> [R,X] = RX - XR on M_2 in the basis
    E11, E12, E21, E22."""
    basis = [Matrix([[1, 0], [0, 0]]), Matrix([[0, 1], [0, 0]]),
             Matrix([[0, 0], [1, 0]]), Matrix([[0, 0], [0, 1]])]
    cols = []
    for B in basis:
        img = Rmat * B - B * Rmat
        cols.append([img[0, 0], img[0, 1], img[1, 0], img[1, 1]])
    return Matrix(cols).T


def produce_smyth_floors():
    """Thm 4.5 + Rem 4.6/4.10: lambda=2c, and the floors 2c log mu_S at c in {1,n}."""
    mp.mp.dps = 50
    muS = mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))
    base = 2 * mp.log(muS)
    payload = {
        "identity": "lambda = 2c  (from D_KL = (1/2c)||r||_G^2 >= log Mah)",
        "lambda_at_c_1_gaussian": 2,
        "lambda_at_c_n_maxentropy": "2n",
        "smyth_constant_mu_S": mp.nstr(muS, 25),
        "mu_S_is_plastic_root_of_x3_minus_x_minus_1":
            bool(abs(muS**3 - muS - 1) < mp.mpf("1e-40")),
        "floor_2logmuS_c1": mp.nstr(base, 12),                 # 0.5623991486
        "floor_n2logmuS_c_n4": mp.nstr(4 * base, 12),          # 2.2496...
        "floor_n2logmuS_c_n8": mp.nstr(8 * base, 12),          # 4.4992...
        "note": "the two shipped floors are the c=1 and c=n readings of one identity",
    }
    core.write_json(SCRIPT, "smyth_floors.json", payload)
    return payload


def produce_certified_grow():
    """Ex 4.3b: certified GROW via mpmath.iv verified interval arithmetic; the
    exact rational gain must dominate a rigorous UPPER bound on the cost."""
    mp.iv.dps = 40
    records = []
    for name, Mah, gain in [("sqrt7", 7, 56), ("2*sqrt6", 24, 96)]:
        enc = mp.iv.log(mp.iv.mpf(Mah))            # rigorous enclosure of log(Mah)
        cost_hi = float((2 * enc).b)               # lambda=2 -> upper bound on cost
        records.append({
            "seed": name, "mahler": Mah,
            "log_enclosure_lo": mp.nstr(enc.a, 10),
            "log_enclosure_hi": mp.nstr(enc.b, 10),
            "lambda": 2,
            "cost_upper_bound_2logMah": mp.nstr((2 * enc).b, 10),
            "exact_gain": gain,
            "certified_grow": bool(gain >= cost_hi),
        })
    # STOP: lattice noise 1/10 < Smyth constant floor 0.5624 (exact rational compare)
    stop = {
        "lattice_noise_magnitude": "1/10",
        "smyth_floor_rational_bracket": "5624/10000",
        "stopped": bool(Q(1, 10) < Q(5624, 10000)),
    }
    payload = {"certified_grow": records, "certified_stop": stop}
    core.write_json(SCRIPT, "certified_grow.json", payload)
    return payload


def produce_selfaction_spectrum():
    """Prop 4.7: adjoint [R_C,.] spectrum {0,0,+-sqrt(1+4C)}, kernel span{I,R_C}."""
    rows = []
    for C in (Q(1, 4), Q(1, 2), sp.Integer(1), sp.Integer(3)):
        RC = Matrix([[0, C], [1, -1]])
        A = adjoint_matrix(RC)
        eigs = A.eigenvals()
        gap = sp.sqrt(1 + 4 * C)
        zero_mult = int(eigs.get(sp.Integer(0), 0))
        ker_dim = len(A.nullspace())
        rows.append([str(C), core.poly_str(sp.expand(RC.charpoly(x).as_expr())).replace(" ", ""),
                     str(sp.nsimplify(gap)), zero_mult, ker_dim])
    core.write_csv(SCRIPT, "selfaction_spectrum.csv",
                   ["C", "charpoly_R_C", "self_action_gap", "zero_eig_mult",
                    "centralizer_kernel_dim"], rows)
    return rows


def produce_R1_similarity():
    """Prop 4.7 aside (corrected): R_1 ~ -R (not R); shared self-action gap sqrt5."""
    R1 = Matrix([[0, 1], [1, -1]])
    Rk = Matrix([[0, 1], [1, 1]])        # keystone
    negRk = -Rk
    cp = lambda M: core.poly_str(sp.expand(M.charpoly(x).as_expr()))
    # self-action gaps
    g1 = {sp.simplify(sp.Abs(e)) for e in adjoint_matrix(R1).eigenvals() if e != 0}
    gk = {sp.simplify(sp.Abs(e)) for e in adjoint_matrix(Rk).eigenvals() if e != 0}
    disc = lambda M: int(M.trace()**2 - 4 * M.det())
    payload = {
        "R_1": [[0, 1], [1, -1]], "charpoly_R_1": cp(R1),            # x^2+x-1
        "keystone_R": [[0, 1], [1, 1]], "charpoly_R": cp(Rk),        # x^2-x-1
        "neg_keystone": [[0, -1], [-1, -1]], "charpoly_negR": cp(negRk),  # x^2+x-1
        "R_1_conjugate_to_neg_R": bool(cp(R1) == cp(negRk)),         # True
        "R_1_conjugate_to_R": bool(cp(R1) == cp(Rk)),               # False
        "reason": "characteristic polynomial is a similarity invariant; "
                  "R_1 and -R share x^2+x-1 with distinct eigenvalues => similar over Q",
        "discriminant_R_1": disc(R1), "discriminant_R": disc(Rk),   # 5, 5
        "shared_self_action_gap": str(list(g1)[0]),                 # sqrt5
        "gaps_coincide": bool(g1 == gk == {sqrt(5)}),
    }
    core.write_json(SCRIPT, "R1_similarity.json", payload)
    return payload


def produce_canonicalization():
    """Def 4.8 / Table 4: gap sqrt(1+4C), c = sqrt(1+4C)/(2C), lambda=sqrt(1+4C)/C."""
    rows = []
    for C in (Q(1, 4), Q(1, 2), sp.Integer(1)):
        gap = sp.sqrt(1 + 4 * C)
        c = gap / (2 * C)
        lam = gap / C
        rows.append([str(C), str(sp.nsimplify(gap)), str(sp.nsimplify(c)),
                     str(sp.nsimplify(lam)), mp.nstr(mp.sqrt(float(1 + 4 * C)), 10)])
    core.write_csv(SCRIPT, "canonicalization_table.csv",
                   ["C_gate", "gap_sqrt(1+4C)", "c=gap/(2C)", "lambda=gap/C",
                    "gap_numeric"], rows)
    return rows


def main():
    sf = produce_smyth_floors()
    cg = produce_certified_grow()
    ss = produce_selfaction_spectrum()
    r1 = produce_R1_similarity()
    ct = produce_canonicalization()
    print(f"[{SCRIPT}] floors: {sf['floor_2logmuS_c1']}, "
          f"{sf['floor_n2logmuS_c_n4']}, {sf['floor_n2logmuS_c_n8']}")
    print(f"[{SCRIPT}] certified GROW: "
          f"{[(r['seed'], r['certified_grow']) for r in cg['certified_grow']]}")
    print(f"[{SCRIPT}] self-action gaps: {[(r[0], r[2]) for r in ss]}")
    print(f"[{SCRIPT}] R_1 ~ -R: {r1['R_1_conjugate_to_neg_R']}, "
          f"R_1 ~ R: {r1['R_1_conjugate_to_R']}, gap {r1['shared_self_action_gap']}")
    print(f"[{SCRIPT}] canon c-values: {[(r[0], r[2]) for r in ct]}")
    print(f"[{SCRIPT}] wrote 5 data files -> {core.data_dir()}")


if __name__ == "__main__":
    main()
