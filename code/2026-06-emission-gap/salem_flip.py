"""Producer: Salem = flip-straddle, the emission delta, and its boundary cases.

Source paper: papers/2026-06-emission-gap/emission_gap_paper.tex
Produces:
  data/2026-06-emission-gap/salem_flip_tests.csv        (Lemma 8.2, Thm. 8.3, Prop. 11.1)
  data/2026-06-emission-gap/tracedown_polynomials.json  (Lemma 8.2, Prop. 8.4, App. A)
Paper results: a reciprocal integer polynomial is Salem iff its trace-down
T (R(x)=x^m T(x+1/x)) is totally real with exactly one root in (2, inf) and the
rest in (-2, 2) (Lemma 8.2); the emission delta fires on Salem carriers (Lehmer,
beta_4) and not on cyclotomic Phi_10 or a non-reciprocal quartic (Thm. 8.3); the
fourth roots of unity have trace-downs {2, 0, -2}, never interior (Prop. 8.4);
the traceless reciprocal quartic x^4+bx^2+1 has trace-down t^2+(b-2) with
symmetric roots, so it never straddles the flip (Prop. 11.1).
"""
import mpmath as mp
import sympy as sp

import emgap_core as C

SCRIPT = "salem_flip.py"
x, t = C.x, C.t
PHI = (1 + mp.sqrt(5)) / 2

# named test polynomials (highest-first), with the paper location they illustrate
CASES = {
    "Lehmer (deg 10)":        ([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1], "Lemma 8.2 / Thm 8.3"),
    "beta_4 (deg 4)":         ([1, -1, -1, -1, 1], "Lemma 8.2 / Cor 10.4"),
    "Phi_10 cyclotomic":      ([1, -1, 1, -1, 1], "Thm 8.3 negative (cyclotomic)"),
    "x^4 - 2 (non-recip)":    ([1, 0, 0, 0, -2], "Thm 8.3 negative (non-reciprocal)"),
    "traceless x^4+x^2+1":    ([1, 0, 1, 0, 1], "Prop 11.1 (symmetric trace-down)"),
    "traceless x^4-x^2+1":    ([1, 0, -1, 0, 1], "Prop 11.1 (symmetric trace-down)"),
}


def salem_flip_rows():
    rows = []
    for name, (coeffs, loc) in CASES.items():
        Q = sp.Poly(coeffs, x)
        deg = Q.degree()
        palindromic = C.is_palindromic(coeffs)
        if deg % 2 == 0 and palindromic:
            T = C.trace_down(Q)
            straddle = C.flip_straddle(T)
            tstr = str(T.as_expr())
        else:
            straddle, tstr = False, "n/a (not even-degree palindromic)"
        rows.append({
            "polynomial": name,
            "coeffs": " ".join(str(c) for c in coeffs),
            "degree": deg,
            "palindromic": palindromic,
            "trace_down_T": tstr,
            "flip_straddle": straddle,
            "is_salem": C.has_salem_factor(Q),
            "mahler_measure": C.s(C.mahler(coeffs), 12),
        })
    return rows


def tracedown_details():
    # Lehmer trace-down and its five real roots
    T_lehmer = C.trace_down(sp.Poly(C.LEHMER, x))
    lehmer_roots = sorted(float(r) for r in C.mp_roots(
        [int(c) for c in T_lehmer.all_coeffs()]) if abs(r.imag) < 1e-9)
    # beta_4 trace-down should equal t^2 - t - 3
    T_b4 = C.trace_down(sp.Poly([1, -1, -1, -1, 1], x))
    # fourth roots of unity trace-downs rho(z) = z + 1/z
    fourth = {}
    for label, z in {"1": 1, "i": 1j, "-1": -1, "-i": -1j}.items():
        rho = z + 1 / z
        fourth[label] = round(rho.real, 12)
    # traceless reciprocal quartic family x^4 + b x^2 + 1 -> t^2 + (b-2)
    family = []
    for b in range(-8, 9):
        T = C.trace_down(sp.Poly([1, 0, b, 0, 1], x))
        family.append({
            "b": b, "trace_down_T": str(T.as_expr()),
            "equals_t2_plus_b_minus_2": bool(sp.expand(T.as_expr() - (t**2 + (b - 2))) == 0),
            "flip_straddle": C.flip_straddle(T),
        })
    return {
        "lehmer_trace_down": {
            "T": str(T_lehmer.as_expr()),
            "roots": [round(r, 4) for r in lehmer_roots],
            "n_roots_gt_2": sum(1 for r in lehmer_roots if r > 2),
            "n_roots_in_minus2_2": sum(1 for r in lehmer_roots if -2 < r < 2),
            "flip_straddle": C.flip_straddle(T_lehmer),
        },
        "beta4_trace_down": {
            "T": str(T_b4.as_expr()),
            "equals_t2_minus_t_minus_3": bool(sp.expand(T_b4.as_expr() - (t**2 - t - 3)) == 0),
            "flip_straddle": C.flip_straddle(T_b4),
        },
        "fourth_roots_of_unity_trace_downs": fourth,
        "fourth_root_trace_down_set": sorted(set(fourth.values())),
        "traceless_reciprocal_quartic_family": family,
        "any_traceless_quartic_straddles": any(f["flip_straddle"] for f in family),
    }


def main():
    rows = salem_flip_rows()
    C.write_csv("salem_flip_tests.csv",
                ["polynomial", "coeffs", "degree", "palindromic", "trace_down_T",
                 "flip_straddle", "is_salem", "mahler_measure"], rows, SCRIPT)
    details = tracedown_details()
    C.write_json("tracedown_polynomials.json", details, SCRIPT)

    print("salem_flip_tests.csv:")
    for r in rows:
        print(f"  {r['polynomial']:24s} salem={r['is_salem']!s:5s} straddle={r['flip_straddle']}")
    print(f"Lehmer T = {details['lehmer_trace_down']['T']}, "
          f"roots {details['lehmer_trace_down']['roots']}")
    print(f"fourth-root trace-down set = {details['fourth_root_trace_down_set']}")


if __name__ == "__main__":
    main()
