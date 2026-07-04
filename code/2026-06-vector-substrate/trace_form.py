"""
PRODUCER: trace-form Grams, discriminants, the field catalog, and trace duality.

Source paper: papers/2026-06-vector-substrate/vector_substrate.tex
  - Thm. 2.14 (G = M^T M, det G = d_K), Ex. 2.15 (golden), Ex. 2.16 (biquadratic
    diag(4,8,12,24), index 2), Ex. 2.18 (Q(i): trace form diag(2,-2), Hermitian
    G_2 = 2I), Ex. 2.19 (Q(cbrt2): det -108, signature (2,1))
  - Table 2 catalog: det G = d_K, covol = sqrt|d_K|, N(different)=N(m'(theta))=|d_K|
  - Ex. 3.9 trace duality: dual basis, different (sqrt5)

Emits:
  data/trace_forms.json  -- Gram matrices, dets, signatures, G=M^T M numeric
                            cross-check, and the golden trace-duality data.
  data/field_catalog.csv -- Table 2: (r1,r2), det G, covol, N(different).
"""
import mpmath as mp
import sympy as sp

import vsub_core as vc

x = vc.x
SCRIPT = "trace_form.py"

# coords of sqrt2, sqrt3, sqrt6 in the power basis {1,theta,theta^2,theta^3},
# theta = sqrt2 + sqrt3 (derived once, reused):
S2 = (0, sp.Rational(-9, 2), 0, sp.Rational(1, 2))
S3 = (0, sp.Rational(11, 2), 0, sp.Rational(-1, 2))
S6 = (sp.Rational(-5, 2), 0, sp.Rational(1, 2), 0)
ONE4 = (1, 0, 0, 0)


def _signature(G):
    """(#positive eigenvalues, #negative eigenvalues) of an exact symmetric G,
    counted with multiplicity."""
    eigs = G.eigenvals()
    pos = sum(m for e, m in eigs.items() if e > 0)
    neg = sum(m for e, m in eigs.items() if e < 0)
    return int(pos), int(neg)


def gram_matrices():
    out = {}

    # Ex. 2.15 golden field, basis {1, phi}
    Cphi = vc.companion_from_poly(x**2 - x - 1)
    Gg = vc.gram([(1, 0), (0, 1)], Cphi)
    out["golden_Q_sqrt5"] = {
        "basis": "{1, phi}",
        "gram": vc.mat_to_list(Gg),
        "det": vc.sval(Gg.det()),
        "d_K": 5,
        "traces": [vc.sval(vc.field_trace((1, 0), Cphi)),
                   vc.sval(vc.field_trace((0, 1), Cphi)),
                   vc.sval(vc.field_trace((1, 1), Cphi))],   # Tr(1),Tr(phi),Tr(phi^2)
    }

    # Ex. 2.16 biquadratic Q(sqrt2,sqrt3), product basis {1,sqrt2,sqrt3,sqrt6}
    C4 = vc.companion_from_poly(x**4 - 10 * x**2 + 1)
    Gb = vc.gram([ONE4, S2, S3, S6], C4)
    dK_biquad = 8 * 12 * 24                               # 2304 = 2^8 3^2
    out["biquadratic_Q_sqrt2_sqrt3"] = {
        "basis": "{1, sqrt2, sqrt3, sqrt6}",
        "gram": vc.mat_to_list(Gb),
        "diagonal": [vc.sval(Gb[i, i]) for i in range(4)],
        "det": vc.sval(Gb.det()),
        "d_K": dK_biquad,
        "index": vc.sval(sp.sqrt(sp.Rational(int(Gb.det()), dK_biquad))),
    }

    # Ex. 2.18 Q(i): power-basis trace form (indefinite) + Hermitian companion
    Ci = vc.companion_from_poly(x**2 + 1)
    Gi = vc.power_gram(Ci)
    Mi = sp.Matrix([[1, sp.I], [1, -sp.I]])
    G2i = sp.simplify(Mi.conjugate().T * Mi)
    pos, neg = _signature(Gi)
    out["gaussian_Q_i"] = {
        "basis": "{1, i}",
        "trace_form": vc.mat_to_list(Gi),
        "trace_form_signature": [pos, neg],
        "hermitian_G2": vc.mat_to_list(G2i),
        "det_G2": vc.sval(G2i.det()),
        "abs_d_K": 4,
        "covol": vc.sval(sp.sqrt(G2i.det())),
        "geometric_rescaled_covol": vc.sval(sp.Rational(1, 2) * sp.sqrt(G2i.det())),
    }

    # Ex. 2.19 Q(cbrt2): power-basis trace form, det -108, signature (2,1)
    Cc = vc.companion_from_poly(x**3 - 2)
    Gc = vc.power_gram(Cc)
    posc, negc = _signature(Gc)
    out["cubic_Q_cbrt2"] = {
        "basis": "{1, cbrt2, cbrt4}",
        "trace_form": vc.mat_to_list(Gc),
        "det": vc.sval(Gc.det()),
        "d_K": -108,
        "signature": [posc, negc],
        "eigenvalues": sorted(vc.sval(e) for e in Gc.eigenvals().keys()),
    }
    return out


def g_equals_mtm_check():
    """Thm. 2.14: cross-check G = M^T M in the power basis against numeric
    Minkowski embeddings for six fields; emit the max abs discrepancy."""
    cases = [x**2 - x - 1, x**2 - 2, x**2 - 7, x**3 - 2, x**2 + 1, x**4 - 10 * x**2 + 1]
    rows = {}
    mp.mp.dps = 40
    for poly in cases:
        C = vc.companion_from_poly(poly)
        Gexact = vc.power_gram(C)
        M = vc.embedding_matrix_power(poly, dps=40)
        MtM = M.T * M                                    # unconjugated -> trace form
        n = C.shape[0]
        err = mp.mpf(0)
        for i in range(n):
            for j in range(n):
                err = max(err, abs(mp.mpf(MtM[i, j].real) - mp.mpf(int(Gexact[i, j]))))
                err = max(err, abs(mp.mpf(MtM[i, j].imag)))
        rows[sp.sstr(poly)] = {
            "det_G_exact": vc.sval(Gexact.det()),
            "max_abs_error_MtM_vs_G": mp.nstr(err, 5),
        }
    return rows


def field_catalog_rows():
    """Table 2: (r1,r2), det G = d_K (signed), covol = sqrt|d_K|, and
    N(different) = N(m'(theta)) = |d_K|."""
    catalog = [
        # name, min_poly, (r1,r2), signed d_K, covol_str
        ("Q(sqrt5)",   x**2 - x - 1, (2, 0), 5,    "sqrt(5)"),
        ("Q(sqrt2)",   x**2 - 2,     (2, 0), 8,    "2*sqrt(2)"),
        ("Q(sqrt3)",   x**2 - 3,     (2, 0), 12,   "2*sqrt(3)"),
        ("Q(sqrt7)",   x**2 - 7,     (2, 0), 28,   "2*sqrt(7)"),
        ("Q(i)",       x**2 + 1,     (0, 1), -4,   "2"),
        ("Q(cbrt2)",   x**3 - 2,     (1, 1), -108, "6*sqrt(3)"),
    ]
    rows = []
    for name, poly, (r1, r2), dK, covol in catalog:
        C = vc.companion_from_poly(poly)
        G = vc.power_gram(C)
        detG = int(G.det())
        # different generator = m'(theta); N(different) = |N(m'(theta))| = |d_K|
        mprime = sp.Poly(poly, x).diff(x)
        n = C.shape[0]
        deriv_coords = [mprime.coeff_monomial(x**i) for i in range(n)]
        Nd = int(vc.field_norm(deriv_coords, C))
        rows.append([name, f"({r1},{r2})", detG, dK, covol, abs(Nd),
                     "OK" if abs(detG) == abs(dK) == abs(Nd) else "MISMATCH"])
    return rows


def trace_duality_golden():
    """Ex. 3.9: G^{-1} = (1/5)[[3,-1],[-1,2]]; dual basis 1^v=(3-phi)/5,
    phi^v=(2phi-1)/5; Tr(1*1^v)=1, Tr(phi*1^v)=0; different (2phi-1)=(sqrt5)."""
    C = vc.companion_from_poly(x**2 - x - 1)
    G = vc.gram([(1, 0), (0, 1)], C)
    Ginv = G.inv()
    one_dual = [Ginv[0, 0], Ginv[0, 1]]
    phi_dual = [Ginv[1, 0], Ginv[1, 1]]
    prod = vc.rho((0, 1), C) * sp.Matrix(one_dual)       # coords of phi * 1^v
    return {
        "G_inverse": vc.mat_to_list(Ginv),
        "dual_of_1": [vc.sval(v) for v in one_dual],     # (3/5, -1/5)
        "dual_of_phi": [vc.sval(v) for v in phi_dual],   # (-1/5, 2/5)
        "Tr_1_times_1dual": vc.sval(vc.field_trace(one_dual, C)),   # 1
        "Tr_phi_times_1dual": vc.sval(vc.field_trace(list(prod), C)),  # 0
        "different_generator": "2*phi - 1 = sqrt(5)",
        "norm_different": vc.sval(sp.Abs(vc.field_norm((-1, 2), C))),   # 5
    }


def main():
    payload = {
        "gram_matrices": gram_matrices(),
        "G_equals_MtM_numeric_check": g_equals_mtm_check(),
        "trace_duality_golden": trace_duality_golden(),
    }
    p1 = vc.write_json("trace_forms.json", payload, SCRIPT)
    rows = field_catalog_rows()
    p2 = vc.write_csv(
        "field_catalog.csv",
        ["field", "(r1,r2)", "det_G", "d_K", "covol", "N_different", "check"],
        rows,
        SCRIPT,
    )
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print("  golden Gram det =", payload["gram_matrices"]["golden_Q_sqrt5"]["det"],
          "biquad index =", payload["gram_matrices"]["biquadratic_Q_sqrt2_sqrt3"]["index"])
    print("  catalog checks:", [r[-1] for r in rows])


if __name__ == "__main__":
    main()
