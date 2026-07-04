"""Shared producer library for *The Pisot Cross-Shell Residue* whitepaper.

Source paper: papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex
    (AceTheDactyl, Echo S Studios Research Developments, July 2026).

This module is the clean, reusable engine that the producer scripts in this
folder import.  It rebuilds the paper's objects FROM SCRATCH -- companion
matrices, resultants, high-precision roots, exact cyclotomic factor matching --
so every emitted datum is computed, not restated.  It is intentionally
independent of the sibling verifier module in tests/2026-07-pisot-residue/
(the tests assert; these producers emit).

Engines: sympy 1.14 (exact resultants / factorization / irreducibility),
mpmath 1.3 (high-precision roots and unimodularity), numpy 2.4 (poly-from-roots).
Run under the Windows ``py`` launcher (Python 3.12).
"""
import json
import os

import sympy as sp
from sympy import symbols, Poly, resultant, cyclotomic_poly, totient
import mpmath as mp

x, y = symbols('x y')

# ---- provenance constants (embedded into every emitted artifact) -------------
SOURCE_PAPER = "papers/2026-07-pisot-residue/pisot_residue_whitepaper.tex"
DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "2026-07-pisot-residue")


# =============================================================================
# Ratio object  Rat_p = prim  Res_y( p(y), p(x*y) )        (Proposition 2.1)
# =============================================================================
def rat_object(coeffs):
    """Primitive part of Res_y(p(y), p(x*y)) as a sympy Poly in x.

    ``coeffs`` is a monic integer polynomial, high->low degree, e.g.
    [1, a, b, c, d, e] for a monic quintic.  Returns Rat_p, whose degree is
    exactly n^2 and whose root multiset is the n^2 ordered ratios a_j/a_i.
    """
    p_y = Poly(coeffs, y)
    n = p_y.degree()
    expr = sum(coeffs[k] * (x * y) ** (n - k) for k in range(len(coeffs)))
    p_xy = Poly(sp.expand(expr), y)
    R = resultant(p_y.as_expr(), p_xy.as_expr(), y)
    return Poly(sp.expand(R), x).primitive()[1]


def cyclotomic_scan(poly_in_x):
    """Exact cyclotomic 'contact scan': {m: multiplicity of Phi_m} in a poly.

    Factor over Z, then match each irreducible factor against Phi_m for the
    finitely many m with phi(m)=deg(factor)  (Lemma 2.3 makes this a finite,
    all-integer decision procedure, since Phi_m | P forces m <= 2 (deg P)^2).
    """
    P = Poly(poly_in_x, x)
    out = {}
    for fac, mult in P.factor_list()[1]:
        d = fac.degree()
        for m in range(1, 2 * d * d + 11):           # phi(m)=d  =>  m <= 2 d^2
            if int(totient(m)) == d and Poly(cyclotomic_poly(m, x), x) == fac:
                out[m] = out.get(m, 0) + mult
                break
    return out


def multiplicity_of_factor(poly_in_x, lin):
    """Exact multiplicity of a linear factor `lin` (e.g. x-1) in a polynomial."""
    q = Poly(poly_in_x, x)
    L = Poly(lin, x)
    m = 0
    while q.rem(L).is_zero:
        q = q.quo(L)
        m += 1
    return m


# =============================================================================
# Companion matrix (last-column form): eigenvalues = roots of the polynomial.
# =============================================================================
def companion(coeffs):
    n = len(coeffs) - 1
    C = sp.zeros(n, n)
    for i in range(1, n):
        C[i, i - 1] = 1
    for i in range(n):
        C[i, n - 1] = -coeffs[n - i]          # -a0 .. -a_{n-1}
    return C


# =============================================================================
# High-precision root machinery (mpmath).
# =============================================================================
def hp_roots(coeffs, dps=60):
    with mp.workdps(dps):
        return mp.polyroots([mp.mpf(int(c)) for c in coeffs],
                            maxsteps=600, extraprec=4 * dps)


def is_pisot(coeffs, dps=50):
    """Certify a Pisot number by root geometry.

    Returns (is_pisot, n_real_strict_inside, n_nonreal_pairs), the counts taken
    over the roots other than the dominant theta.  Requires exactly one root
    outside the closed unit disk (that root real > 1) and NO root on the unit
    circle -- a genuine modulus gap, sound because the certified inputs are
    irreducible and non-reciprocal.
    """
    with mp.workdps(dps):
        tol = mp.mpf(10) ** (-dps // 4)
        rts = hp_roots(coeffs, dps)
        if any(abs(abs(r) - 1) < tol for r in rts):
            return (False, None, None)
        outside = [r for r in rts if abs(r) > 1]
        if len(outside) != 1:
            return (False, None, None)
        th = outside[0]
        if abs(th.imag) > tol or th.real <= 1:
            return (False, None, None)
        inside = [r for r in rts if abs(r) < 1]
        nonreal = [r for r in inside if abs(r.imag) > tol]
        real_in = [r for r in inside if abs(r.imag) <= tol]
        return (True, len(real_in), len(nonreal) // 2)


def dominant_root(coeffs, dps=50):
    with mp.workdps(dps):
        rts = hp_roots(coeffs, dps)
        return max((r.real for r in rts
                    if abs(r.imag) < mp.mpf(10) ** (-dps // 4)), default=None)


def is_root_of_unity(z, Mmax=1680, dps=80):
    """True iff z is a root of unity of order <= Mmax: |z|=1 and z^m=1 for some
    1<=m<=Mmax.  A unimodular z at an irrational angle returns False (its powers
    never return near 1).  This is the C_2 negative-certificate test."""
    with mp.workdps(dps):
        eps = mp.mpf(10) ** (-dps // 3)
        if abs(abs(z) - 1) > eps:
            return False
        p = mp.mpf(1) * z
        for _ in range(1, Mmax + 1):
            p = p * z
            if abs(p - 1) < eps:
                return True
        return False


def phi(dps=50):
    with mp.workdps(dps):
        return (1 + mp.sqrt(5)) / 2


def is_pm_reciprocal(coeffs):
    """Coefficient vector equals its reversal (palindromic) or minus-reversal
    (anti-palindromic): the +-reciprocal carriers of unit-circle roots."""
    rev = coeffs[::-1]
    return coeffs == rev or coeffs == [-c for c in rev]


# =============================================================================
# Totient / candidate-count helpers (Lemma 2.3).
# =============================================================================
def totient_candidates(K):
    """Sorted list of all m>=1 with phi(m) <= K (finite: phi(m)>=sqrt(m/2))."""
    return [m for m in range(1, 4 * K * K + 10) if int(totient(m)) <= K]


# =============================================================================
# Formatting / serialization helpers.
# =============================================================================
def poly_str(coeffs):
    """Human-readable monic polynomial from a high->low coefficient vector,
    e.g. [1,-2,-2,-2,-2,-2] -> 'x^5 - 2x^4 - 2x^3 - 2x^2 - 2x - 2'."""
    n = len(coeffs) - 1
    terms = []
    for i, c in enumerate(coeffs):
        c = int(c)
        p = n - i
        if c == 0:
            continue
        mag = abs(c)
        if p == 0:
            body = str(mag)
        elif p == 1:
            body = ("" if mag == 1 else str(mag)) + "x"
        else:
            body = ("" if mag == 1 else str(mag)) + f"x^{p}"
        sign = "-" if c < 0 else "+"
        terms.append((sign, body))
    if not terms:
        return "0"
    s0, b0 = terms[0]
    out = ("-" if s0 == "-" else "") + b0
    for sign, body in terms[1:]:
        out += f" {sign} {body}"
    return out


def scan_str(scan):
    """Render a cyclotomic-scan dict {m: mult} as '{Phi_1^10}' etc."""
    if not scan:
        return "{}"
    parts = [f"Phi_{m}^{mult}" if mult != 1 else f"Phi_{m}"
             for m, mult in sorted(scan.items())]
    return "{" + ", ".join(parts) + "}"


def scan_json(scan):
    """JSON-friendly form of a scan: list of {m, multiplicity} plus a string."""
    return {"factors": [{"m": int(m), "multiplicity": int(mult)}
                        for m, mult in sorted(scan.items())],
            "pretty": scan_str(scan)}


def write_json(filename, payload, generated_by):
    """Write a JSON artifact with embedded provenance to data/<paper>/."""
    os.makedirs(DATA_DIR, exist_ok=True)
    obj = {
        "_source_paper": SOURCE_PAPER,
        "_generated_by": f"code/2026-07-pisot-residue/{generated_by}",
    }
    obj.update(payload)
    path = os.path.join(DATA_DIR, filename)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)
        fh.write("\n")
    return path


def csv_header(generated_by):
    """Leading provenance comment line for a CSV artifact."""
    return (f"# source: {SOURCE_PAPER}; "
            f"generated by: code/2026-07-pisot-residue/{generated_by}")


def write_csv(filename, header_cols, rows, generated_by):
    """Write a CSV artifact with a leading provenance comment line."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, filename)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(csv_header(generated_by) + "\n")
        fh.write(",".join(header_cols) + "\n")
        for row in rows:
            fh.write(",".join(str(v) for v in row) + "\n")
    return path
