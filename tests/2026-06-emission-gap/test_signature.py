"""Field signatures, the trace form, and the signature lattice
(Prop. 9.1; Prop. 12.3 circulant; Thm. 10.3; App. A).
"""
import mpmath as mp
import sympy as sp

import emgap_util as U

mp.mp.dps = 40

x = U.x


def test_salem_field_signature_and_trace_form():
    """Prop. 9.1: a Salem field of degree 2m has signature (2, m-1) and trace
    form signature (m+1, m-1). Checked on beta_4 (deg 4, m=2) and Lehmer
    (deg 10, m=5)."""
    # beta_4 : x^4 - x^3 - x^2 - x + 1, m = 2
    assert U.signature_from_minpoly([1, -1, -1, -1, 1]) == (2, 1)
    assert U.trace_form_signature([1, -1, -1, -1, 1]) == (3, 1)   # (m+1, m-1)
    # Lehmer : degree 10, m = 5
    L = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    assert U.signature_from_minpoly(L) == (2, 4)
    assert U.trace_form_signature(L) == (6, 4)                     # (m+1, m-1)


def test_K_seed_signature_and_trace_form():
    """Prop. 9.1 / App. A: the catalog K-formation x^4+5x^2-5 has field
    signature (2,1) and trace form (3,1) -- the one Lorentzian generator; its
    complex place is off the circle (modulus 2.4195 >> 1)."""
    assert U.signature_from_minpoly([1, 0, 5, 0, -5]) == (2, 1)
    assert U.trace_form_signature([1, 0, 5, 0, -5]) == (3, 1)
    rts = U.mp_roots([1, 0, 5, 0, -5])
    cplx = [r for r in rts if abs(r.real) < 1e-30]
    assert abs(abs(cplx[0]) - mp.mpf("2.4195")) < 1e-4


def test_totally_real_seeds_are_definite():
    """Prop. 9.1 / App. A: sqrt2, sqrt5, phi seeds have signature (2,0) and a
    positive-definite trace form (2,0)."""
    for coeffs in ([1, 0, -2], [1, 0, -5], [1, -1, -1]):
        assert U.signature_from_minpoly(coeffs) == (2, 0)
        pos, neg = U.trace_form_signature(coeffs)
        assert (pos, neg) == (2, 0)


def test_signature_lattice_computed_values():
    """Thm. 10.3 / App. A: representative subfield signatures inside
    K = Q(sqrt2, sqrt3, 5^{1/4}):
      Q(5^{1/4}) = (2,1), Q(sqrt2,5^{1/4}) = (4,2),
      Q(sqrt2,sqrt3,sqrt5) = (8,0), K = (8,4)."""
    assert U.field_signature(sp.root(5, 4))[0] == (2, 1)
    assert U.field_signature(sp.sqrt(2) + sp.root(5, 4))[0] == (4, 2)
    assert U.field_signature(sp.sqrt(2) + sp.sqrt(3) + sp.sqrt(5))[0] == (8, 0)
    assert U.field_signature(sp.sqrt(2) + sp.sqrt(3) + sp.root(5, 4)) == ((8, 4), 16)


def test_signature_lattice_shape_only_salem_at_degree_four():
    """Thm. 10.3: every non-totally-real subfield of K has signature (2k, k);
    the Salem shape (2, m-1) forces k=1, m=2, i.e. degree 4 (only Q(5^{1/4}))."""
    non_real = [
        U.field_signature(sp.root(5, 4))[0],                 # (2,1)
        U.field_signature(sp.sqrt(2) + sp.root(5, 4))[0],    # (4,2)
        U.field_signature(sp.sqrt(2) + sp.sqrt(3) + sp.root(5, 4))[0],  # (8,4)
    ]
    for (r1, r2) in non_real:
        assert r1 == 2 * r2                                  # shape (2k, k)
    # Salem shape (2, m-1) intersect (2k,k): r1=2 => k=1 => r2=1 => m=2 => deg 4
    salem_shaped = [(r1, r2) for (r1, r2) in non_real if r1 == 2 and r2 == r1 // 2]
    assert salem_shaped == [(2, 1)]


def test_cyclotomic_fields_are_totally_complex_not_salem():
    """Prop. 12.3: circulant eigenvalues lie in Q(zeta_n), which is totally
    complex (CM) -- signature (0, phi(n)/2) -- while a Salem field is (2, m-1);
    the two signatures are disjoint, so no circulant emits a Salem number."""
    for n, phi_n in ((5, 4), (7, 6), (8, 4), (12, 4)):
        cp = sp.Poly(sp.cyclotomic_poly(n, x), x)
        sig = U.signature_from_minpoly([int(c) for c in cp.all_coeffs()])
        assert sig == (0, phi_n // 2)          # totally complex
        assert sig[0] == 0                      # no real embedding => not Salem
