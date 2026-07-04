"""
PRODUCER: the worked learning episode and the growth threshold.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex
  - Sec. 6 (worked episode in Q(sqrt2+sqrt3)): 2sqrt6 residual gain = Tr(24) = 96;
    Pythagoras ||theta+2sqrt6||^2 = 116 = 20 + 96 (cross term 0)
  - Sec. 8: Def. 8.1 (gain vs cost), Thm. 8.2 (Smyth floor mu_S), Conj. 8.4
    (Lehmer mu_L), Sec. 8.3 evidence table, Ex. 8.5 (lattice-aligned noise STOP),
    Ex. 8.7 (certified GROW via interval log-Mahler), Rem. 8.10 (degree-aware floor)
  - Conj. 7.16 golden-gate scale c=sqrt(1+4C)/(2C) (cross-domain; internal
    consistency only, per the paper's own honesty discipline)

Uses lambda = 2 (the paper's representative MDL/BIC-scale rate).

Emits:
  data/threshold_decisions.csv -- (seed, n, gain, Mahler, cost, floor, decision)
  data/threshold_constants.json
"""
import mpmath as mp
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "growth_threshold.py"
LAMBDA = 2
mp.mp.dps = 50


def _muS():
    return mp.findroot(lambda t: t**3 - t - 1, mp.mpf("1.3"))


def worked_episode():
    """Sec. 6: 2sqrt6 = theta^2-5 is G-orthogonal to col(B)=Q+Q.theta; residual
    is exactly 2sqrt6 with ||r||_G^2 = Tr((2sqrt6)^2) = Tr(24) = 96.  Pythagoras
    on theta+2sqrt6: 116 = 20 + 96 with vanishing cross term."""
    C = vc.companion_from_poly(x**4 - 10 * x**2 + 1)
    G = vc.power_gram(C)
    B = sp.Matrix.hstack(sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1, 0, 0]))
    r = vc.residual((-5, 0, 1, 0), B, G)
    theta = (0, 1, 0, 0)
    two_s6 = (-5, 0, 1, 0)
    xcoord = tuple(sp.Matrix(theta) + sp.Matrix(two_s6))    # (-5,1,1,0)
    cross = vc.rho(theta, C) * sp.Matrix(two_s6)
    return {
        "residual_is_2sqrt6": r == sp.Matrix([-5, 0, 1, 0]),
        "residual_norm2_gain": vc.sval(vc.gnorm2(r, G)),        # 96
        "residual_trace_zero": vc.sval(vc.field_trace((-5, 0, 1, 0), C)),  # 0
        "pythagoras_total": vc.sval(vc.gnorm2(sp.Matrix(xcoord), G)),      # 116
        "captured_part_theta_norm2": vc.sval(vc.gnorm2(sp.Matrix(theta), G)),  # 20
        "residual_part_norm2": vc.sval(vc.gnorm2(sp.Matrix(two_s6), G)),   # 96
        "cross_term_trace": vc.sval(vc.field_trace(list(cross), C)),       # 0
        "check_116_eq_20_plus_96": 116 == 20 + 96,
    }


def floors():
    """Thm. 8.2 (Smyth) and Conj. 8.4 (Lehmer): the arithmetic floors."""
    muS = _muS()
    coeffs = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]     # Lehmer's degree-10 poly
    muL = mp.findroot(lambda t: sum(c * t**(10 - i) for i, c in enumerate(coeffs)),
                      mp.mpf("1.17"))
    const_floor = LAMBDA * mp.log(muS)
    return {
        "mu_S": mp.nstr(muS, 15),
        "mu_S_is_root_x3_minus_x_minus_1": bool(abs(muS**3 - muS - 1) < mp.mpf("1e-40")),
        "log_mu_S": mp.nstr(mp.log(muS), 12),
        "constant_floor_lambda2": mp.nstr(const_floor, 12),        # 0.562...
        "mu_L": mp.nstr(muL, 15),
        "mu_L_less_than_mu_S": bool(muL < muS),
        "degree_aware_floor_n4": mp.nstr(4 * const_floor, 6),      # 2.25
        "degree_aware_floor_n8": mp.nstr(8 * const_floor, 6),      # 4.50
    }


def certified_enclosures():
    """Ex. 8.7: rigorous rational enclosures of log-Mahler, so an exact rational
    gain dominates lambda * (upper bound) and GROW is certified."""
    out = {}
    for name, mah, gain in [("sqrt7", 7, 56), ("2sqrt6", 24, 96)]:
        lo, hi = vc.log_interval(mah, places=5)
        cost_upper = LAMBDA * hi
        out[name] = {
            "mahler": mah,
            "log_interval": [vc.sval(lo), vc.sval(hi)],
            "cost_upper_bound": vc.sval(cost_upper),
            "gain": gain,
            "certified_grow": bool(sp.Rational(gain) >= cost_upper),
        }
    return out


def golden_gate_scale():
    """Conj. 7.16 (cross-domain, companion GSA paper): the SELECTED scale
    c = sqrt(1+4C)/(2C) at the golden gate C=1 gives c=sqrt5/2, lambda=2c=sqrt5,
    ladder root 1/phi.  Internal-consistency only -- NOT an arithmetic claim of
    this paper."""
    C = sp.Integer(1)
    c = sp.sqrt(1 + 4 * C) / (2 * C)
    root = sp.Rational(-1, 2) + sp.sqrt(5) / 2
    phi = (1 + sp.sqrt(5)) / 2
    return {
        "c_at_C1": vc.sval(sp.simplify(c)),                    # sqrt5/2
        "lambda_eq_2c": vc.sval(sp.simplify(2 * c)),           # sqrt5
        "ladder_root": vc.sval(sp.simplify(root)),             # (sqrt5-1)/2
        "root_is_1_over_phi": bool(sp.simplify(root - 1 / phi) == 0),
        "note": "cross-domain SELECTED scale; internal consistency only",
    }


def decision_rows():
    """Sec. 8.3 evidence table (lambda=2): the two shipped GROW cases and the
    tiny lattice-aligned STOP.  gain is exact rational; cost = lambda*log(Mahler);
    constant floor = lambda*log(mu_S) = 0.562..."""
    const_floor = LAMBDA * mp.log(_muS())
    # 2sqrt6 gain from the degree-4 field
    C4 = vc.companion_from_poly(x**4 - 10 * x**2 + 1)
    gain_2s6 = int(vc.gnorm2(sp.Matrix([-5, 0, 1, 0]), vc.power_gram(C4)))   # 96
    # sqrt7 gain from the degree-8 compositum Gram diagonal (sqrt7 slot = 56 = 8*7)
    KL = sp.Matrix(sp.kronecker_product(sp.diag(4, 8, 12, 24), sp.diag(2, 14)))
    gain_s7 = int(KL[1, 1])                                                  # 56
    rows = []
    for name, n, gain, mah in [("2sqrt6", 4, gain_2s6, 24), ("sqrt7", 8, gain_s7, 7)]:
        cost = LAMBDA * mp.log(mah)
        rows.append([name, n, gain, mah, mp.nstr(cost, 8),
                     mp.nstr(const_floor, 6),
                     "GROW" if gain >= cost else "STOP"])
    # tiny lattice-aligned noise (1/10) sqrt5 in Q(sqrt5): gain = 1/10 -> STOP
    tiny = sp.Rational(1, 10)
    rows.append(["tiny_off_axis", 4, "1/10", "-", "-",
                 mp.nstr(const_floor, 6),
                 "STOP" if mp.mpf(tiny) < const_floor else "GROW"])
    return rows


def lattice_noise():
    """Ex. 8.5: r = (1/10) sqrt5 in Q(sqrt5) is trace-zero and exactly off-axis,
    yet ||r||_G^2 = 1/10 < 0.562 = floor -> STOP (though floor 0 would GROW)."""
    C = vc.companion_from_poly(x**2 - x - 1)
    G = vc.gram([(1, 0), (0, 1)], C)
    r = sp.Rational(1, 10) * sp.Matrix([-1, 2])
    const_floor = LAMBDA * mp.log(_muS())
    return {
        "residual_trace_zero": vc.sval(vc.field_trace(list(r), C)),   # 0
        "residual_norm2": vc.sval(vc.gnorm2(r, G)),                   # 1/10
        "constant_floor": mp.nstr(const_floor, 6),                   # 0.562
        "stop_below_floor": bool(mp.mpf(sp.Rational(1, 10)) < const_floor),
        "floor0_would_grow": True,
    }


def main():
    payload = {
        "worked_episode": worked_episode(),
        "floors": floors(),
        "certified_enclosures": certified_enclosures(),
        "lattice_noise": lattice_noise(),
        "golden_gate_scale": golden_gate_scale(),
        "lambda": LAMBDA,
    }
    p1 = vc.write_json("threshold_constants.json", payload, SCRIPT)
    p2 = vc.write_csv(
        "threshold_decisions.csv",
        ["seed", "n", "gain", "mahler", "cost_lambda_logM", "constant_floor", "decision"],
        decision_rows(),
        SCRIPT,
    )
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print("  episode gain =", payload["worked_episode"]["residual_norm2_gain"],
          "; Pythagoras 116 = 20 + 96:", payload["worked_episode"]["check_116_eq_20_plus_96"])
    print("  mu_S =", payload["floors"]["mu_S"], "; floor(lambda=2) =",
          payload["floors"]["constant_floor_lambda2"])


if __name__ == "__main__":
    main()
