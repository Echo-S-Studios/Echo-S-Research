"""
Producer: the canonical flip D = 1 + 4C (Section 4).

Source paper: papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex
Produces    : data/2026-07-helix-orthogonal-partner/flip.json
              data/2026-07-helix-orthogonal-partner/flip_family.csv

For the gate family g_C(x) = x^2 + x - C (Def. 4.1, discriminant D = 1+4C) this
recomputes Theorem 4.2:
  * roots (-1 +/- sqrt D)/2, with the double root -1/2 at the flip C=-1/4;
  * the trace-form Gram G = [[2,-1],[-1,1+2C]] in the power basis {1,theta},
    with det G = 1+4C = D (eq. 5), built from the field trace over Q(theta);
  * the signature flip: positive-definite (Riemannian, (2,0)) for D>0,
    indefinite (Lorentzian, (1,1)) for D<0, degenerate at D=0;
  * the D<0 rotation channel spec(ad_M) = {0,0,+i sqrt|D|,-i sqrt|D|};
  * the golden face C=1 gives D=5=(sqrt5)^2 (connects to spec(ad_R)).

flip_family.csv sweeps representative C on both sides of the fold; each row is a
self-contained "which face" record. The symbolic det G = D is in flip.json.

Run: py code/2026-07-helix-orthogonal-partner/flip.py
"""

import sympy as sp

import helix_core as hc
from helix_io import write_csv, write_json

C = sp.symbols("C")


def field_trace(expr_of_theta, theta_sym):
    """Field trace over Q(theta)/Q for theta a root of x^2+x-C: sum over the two
    conjugate roots (-1 +/- sqrt(1+4C))/2."""
    D = 1 + 4 * C
    t1 = (-1 + sp.sqrt(D)) / 2
    t2 = (-1 - sp.sqrt(D)) / 2
    return sp.simplify(expr_of_theta.subs(theta_sym, t1) + expr_of_theta.subs(theta_sym, t2))


def symbolic_flip():
    """Def. 4.1 + Thm. 4.2: discriminant, root formula, symbolic Gram, det G = D."""
    x = hc.x
    th = sp.symbols("theta")
    g = x**2 + x - C
    a, b, c = sp.Poly(g, x).all_coeffs()
    D = sp.simplify(b**2 - 4 * a * c)
    tr1 = field_trace(sp.Integer(1) + 0 * th, th)
    trt = field_trace(th, th)
    trt2 = field_trace(th**2, th)
    G = sp.Matrix([[tr1, trt], [trt, trt2]])
    return {
        "gate": "g_C(x) = x^2 + x - C",
        "discriminant_D": str(D),
        "D_equals_1_plus_4C": bool(sp.simplify(D - (1 + 4 * C)) == 0),
        "roots": "(-1 +/- sqrt(D))/2",
        "double_root_at_flip": str(
            sp.solve((x**2 + x - C).subs(C, sp.Rational(-1, 4)), x)
        ),
        "flip_location_C": "-1/4",
        "trace_form_Gram": [[str(tr1), str(trt)], [str(trt), str(trt2)]],
        "gram_det": str(sp.simplify(G.det())),
        "gram_det_equals_D": bool(sp.simplify(G.det() - (1 + 4 * C)) == 0),
        "golden_face_C1_D": int((1 + 4 * sp.Integer(1))),
        "golden_face_sqrtD_is_sqrt5": bool(
            sp.simplify(sp.sqrt(1 + 4 * sp.Integer(1)) - hc.sqrt5) == 0
        ),
    }


def gram_signature(cval):
    """Signature of the trace form G(C)=[[2,-1],[-1,1+2C]] from eigenvalue signs."""
    G = sp.Matrix([[2, -1], [-1, 1 + 2 * cval]])
    evs = [sp.nsimplify(e) for e in G.eigenvals().keys()]
    pos = sum(1 for e in evs if e > 0)
    neg = sum(1 for e in evs if e < 0)
    zero = sum(1 for e in evs if e == 0)
    return G, (pos, neg, zero)


def rotation_channel(cval):
    """Thm. 4.2: for D<0 the ad-channel of the companion M_c is {0,0,+i sqrt|D|,-i sqrt|D|}."""
    D = 1 + 4 * cval
    Mc = sp.Matrix([[0, cval], [1, -1]])          # companion of x^2 + x - C
    eig = hc.eig_multiset(hc.ad_matrix(Mc))
    root = sp.sqrt(-D)
    return {
        "C": str(cval),
        "D": int(D),
        "companion": str(Mc.tolist()),
        "ad_spectrum": "{0, 0, +i sqrt|D|, -i sqrt|D|}",
        "mult_zero": int(eig.get(sp.Integer(0), 0)),
        "mult_plus_i_root": int(eig.get(sp.I * root, 0)),
        "mult_minus_i_root": int(eig.get(-sp.I * root, 0)),
        "sqrt_abs_D_decimal": hc.dec(root),
    }


def family_rows():
    """Sweep C on both sides of the fold C=-1/4; one 'which face' record per C."""
    x = hc.x
    values = [sp.Integer(1), sp.Rational(1, 4), sp.Integer(0),
              sp.Rational(-1, 4), sp.Rational(-1, 2), sp.Integer(-1), sp.Integer(-2)]
    rows = []
    for cval in values:
        D = 1 + 4 * cval
        Dsign = "positive" if D > 0 else ("zero" if D == 0 else "negative")
        if D > 0:
            root_type, field, face = "two real", "real quadratic Q(sqrt D)", "terrain (growth)"
        elif D == 0:
            root_type, field, face = "double root -1/2", "Q (degenerate)", "the fold"
        else:
            root_type, field, face = "complex conjugate pair", "imaginary quadratic", "rotation"
        G, (pos, neg, zero) = gram_signature(cval)
        if zero:
            sig = "degenerate (det 0)"
        elif neg == 0:
            sig = "(2,0) Riemannian"
        else:
            sig = "(1,1) Lorentzian"
        rows.append(
            {
                "C": str(cval),
                "D=1+4C": str(D),
                "D_sign": Dsign,
                "root_type": root_type,
                "field": field,
                "channel_face": face,
                "gram_leading_minor": int(G[0, 0]),
                "gram_det": str(sp.simplify(G.det())),
                "trace_form_signature": sig,
            }
        )
    return rows


def main():
    payload = {
        "section": "4 -- The flip: D = 1 + 4C",
        "results": {
            "def_4_1_thm_4_2_flip": symbolic_flip(),
            "thm_4_2_rotation_channel_C_minus_1": rotation_channel(sp.Integer(-1)),
            "signature_examples": {
                "D_positive_C1": {"signature": "(2,0)", "detail": str(gram_signature(sp.Integer(1))[1])},
                "D_zero_Cminus_quarter": {"signature": "degenerate", "detail": str(gram_signature(sp.Rational(-1, 4))[1])},
                "D_negative_Cminus1": {"signature": "(1,1)", "detail": str(gram_signature(sp.Integer(-1))[1])},
            },
        },
    }
    rows = family_rows()
    p_json = write_json("flip.json", payload, __file__)
    p_csv = write_csv(
        "flip_family.csv",
        ["C", "D=1+4C", "D_sign", "root_type", "field", "channel_face",
         "gram_leading_minor", "gram_det", "trace_form_signature"],
        rows,
        __file__,
    )
    print(f"wrote {p_json}")
    print(f"wrote {p_csv}")
    print("  det G = 1+4C = D ; signature flips (2,0)->degenerate->(1,1) at C=-1/4")


if __name__ == "__main__":
    main()
