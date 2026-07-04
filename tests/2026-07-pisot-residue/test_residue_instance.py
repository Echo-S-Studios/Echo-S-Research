"""Section 5: the residue instruments, executed on the first live cross-shell
instance  p = x^5 - 2x^4 - 2x^3 - 2x^2 - 2x - 2.

Re-derives:
  * Rat_p^o = Rat_p/(x-1)^5 is degree 20, squarefree, irreducible (Sec. 6 claim);
  * shell detector = 4 unimodular roots of Rat_p^o (distinct shells, Prop. 5.4);
  * Prop. 5.2 modulus multiset structure of Rat_p^o;
  * Prop. 5.3 the composed-square C_2 has degree d^2 = 400, its roots are the
    ordered products of Rat_p^o-roots, and the negative certificate holds:
    the ONLY products that are roots of unity are the 20 that equal 1
    (=> scan {Phi_1^20}, zero mirrored cross-shell classes, Theorem 6.1).

The C_2 certificate is checked at high precision on the 400 ordered products,
with an exact root-of-unity test (order <= 1680, the totient bound for deg 400).
"""
import os
import sys
from collections import Counter

import sympy as sp
from sympy import symbols, Poly, div, expand, gcd, diff

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pisot_residue_lib import rat_object, hp_roots, is_root_of_unity

x, y = symbols('x y')
P = [1, -2, -2, -2, -2, -2]                       # x^5-2x^4-2x^3-2x^2-2x-2


def _rat_circ():
    R = rat_object(P)
    q, r = div(R.as_expr(), (x - 1) ** 5, x)
    assert sp.simplify(r) == 0
    return Poly(expand(q), x)


_RC = None
_ROOTS = None


def rat_circ():
    global _RC
    if _RC is None:
        _RC = _rat_circ()
    return _RC


def rc_roots():
    global _ROOTS
    if _ROOTS is None:
        Rc = rat_circ()
        _ROOTS = hp_roots([int(c) for c in Rc.all_coeffs()], 80)
    return _ROOTS


def test_rat_circ_degree_squarefree_irreducible():
    """Sec. 6: Rat_p^o is degree 20, squarefree, and irreducible (so S*=Rat_p^o,
    deg C_2 = 400)."""
    Rc = rat_circ()
    assert Rc.degree() == 20
    g = gcd(Rc.as_expr(), diff(Rc.as_expr(), x))
    assert Poly(g, x).degree() == 0                # squarefree
    assert Rc.is_irreducible


def test_shell_detector_reads_four():
    """Prop. 5.4: for the two-pair pattern the number of unimodular roots of
    Rat_p^o is 4 iff |lambda_1| != |lambda_2| (distinct shells)."""
    import mpmath as mp
    rts = rc_roots()
    with mp.workdps(80):
        eps = mp.mpf(10) ** -25
        unimod = [r for r in rts if abs(abs(r) - 1) < eps]
    assert len(unimod) == 4


def test_prop52_modulus_multiset_structure():
    """Prop. 5.2: modulus multiset of Rat_p^o is
      { theta/|l1|, theta/|l2|, |l1|/theta, |l2|/theta   (x2 each),
        |l1|/|l2|, |l2|/|l1|                              (x4 each),
        1                                                 (x4) },
    i.e. 7 distinct values with multiplicity multiset [2,2,2,2,4,4,4] and
    modulus 1 of multiplicity 4; closed under reciprocal."""
    import mpmath as mp
    rts = rc_roots()
    with mp.workdps(80):
        mods = [abs(r) for r in rts]
        # bucket by value at high precision
        buckets = []
        for m in mods:
            for b in buckets:
                if abs(b[0] - m) < mp.mpf(10) ** -20:
                    b[1] += 1
                    break
            else:
                buckets.append([m, 1])
        mults = sorted(b[1] for b in buckets)
        assert len(buckets) == 7
        assert mults == [2, 2, 2, 2, 4, 4, 4]
        one = [b for b in buckets if abs(b[0] - 1) < mp.mpf(10) ** -20]
        assert one and one[0][1] == 4
        # reciprocal closure of the distinct moduli
        vals = [b[0] for b in buckets]
        for v in vals:
            assert any(abs(v * w - 1) < mp.mpf(10) ** -18 for w in vals)


def test_prop53_C2_degree_and_root_products():
    """Prop. 5.3(b): C_2 = Res_y(S*(y), y^d S*(x/y)) is monic of degree d^2 with
    root multiset ALL ordered products of S*-roots.  Verified on the small
    self-reciprocal S* = Z* = x^4-3x^2+1: recompute C_2 exactly, check it is
    monic of degree 16, and confirm its coefficient vector equals that of
    prod_{i,j}(x - r_i r_j) built independently from the S*-roots (this pins the
    full multiset, including the repeated roots C_2 genuinely has)."""
    import numpy as np
    import mpmath as mp
    S = [1, 0, -3, 0, 1]                                              # Z* = x^4-3x^2+1
    assert S == S[::-1]                                               # self-reciprocal
    d = 4
    Sy = Poly(S, y)
    Sxy = sum(S[k] * (x / y) ** (d - k) for k in range(len(S)))       # S(x/y)
    C2 = Poly(expand(sp.resultant(Sy.as_expr(), expand(y ** d * Sxy), y)), x)
    assert C2.degree() == d * d == 16
    lead = C2.all_coeffs()[0]
    assert abs(int(lead)) == 1                                        # monic (up to sign)
    c2_norm = [complex(c) / complex(lead) for c in C2.all_coeffs()]
    with mp.workdps(50):
        sroots = [complex(r) for r in hp_roots(S, 50)]
    prods = [a * b for a in sroots for b in sroots]                   # 16 ordered products
    assert len(prods) == 16
    prod_poly = np.poly(np.array(prods))                             # monic prod(x - r_i r_j)
    assert len(prod_poly) == 17
    for a, b in zip(c2_norm, prod_poly):
        assert abs(a - b) < 1e-6


def test_C2_negative_certificate_zero_mirrored_classes():
    """Prop. 5.3(c) + Theorem 6.1: the ONLY ordered products r_i r_j that are
    roots of unity are the 20 equal to 1 (Phi_1^20).  Every other product that
    lands on the unit circle sits at an IRRATIONAL angle -> zero mirrored
    cross-shell classes."""
    import mpmath as mp
    rts = rc_roots()
    with mp.workdps(80):
        eps1 = mp.mpf(10) ** -25
        epsc = mp.mpf(10) ** -14
        equal_one = 0
        oncircle_not_one = 0
        rogue = 0                       # on-circle, != 1, yet a root of unity
        for a in rts:
            for b in rts:
                z = a * b
                if abs(z - 1) < eps1:
                    equal_one += 1
                elif abs(abs(z) - 1) < epsc:
                    oncircle_not_one += 1
                    if is_root_of_unity(z, Mmax=1680, dps=80):
                        rogue += 1
    assert equal_one == 20               # Phi_1^20
    assert oncircle_not_one > 0          # there genuinely ARE unimodular ratios
    assert rogue == 0                    # ... but none is torsion
    assert rat_circ().degree() ** 2 == 400   # deg C_2
