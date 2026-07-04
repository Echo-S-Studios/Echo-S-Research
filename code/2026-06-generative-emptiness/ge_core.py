"""
ge_core -- producer-side computational core for

    papers/2026-06-generative-emptiness/generative_emptiness.tex
    "The Generative Content of a Conserved Emptiness:
     Kinematic Voids as Superselection Generators --
     the Salem Slot and the Five Objects Its Charge Produces"

This module is the shared engine for the object1..object5 / minimality PRODUCER
scripts in this folder.  It re-implements, independently of the test suite under
tests/2026-06-generative-emptiness, the three operators of the emission algebra
(Def. 2.1), the Mahler measure, and the Z/4Z angle charge (Def. 2.1 / Thm 2.2).

The producers CALL these primitives and EMIT machine-readable data files into
data/2026-06-generative-emptiness/.  Nothing here asserts; the verification role
belongs to the independent tests.

Operators (Def. 2.1 "emission algebra"):
  tensor  (x)   : poly whose roots are pairwise products  lambda_i * mu_j
  sq      ( )^2 : poly whose roots are the squares         lambda_i^2
  dsum    (+)   : direct sum = polynomial product          union of root multisets

Angle charge (Def. 2.1): multiset of  arg(root) / (pi/2)  mod 4, flagged 'OFF'
the Z/4Z lattice when an argument is not an integer multiple of pi/2.

Canonical seeds (Def. 2.1 / Thm 2.2):
  phi seed  x^2 - x - 1        (minimal polynomial of the golden ratio)
  K   seed  x^4 + 5 x^2 - 5    (the quartic "Lorentzian" generator)
"""
from __future__ import annotations

import sympy as sp
from sympy import symbols, Poly, resultant, expand
import mpmath as mp

mp.mp.dps = 50
x, y = symbols("x y")

# ----------------------------------------------------------------------------
# polynomial hygiene
# ----------------------------------------------------------------------------


def as_poly(e):
    """Coerce an expression/Poly into a Poly in x."""
    return e if isinstance(e, Poly) else Poly(e, x)


def normalize(p):
    """Return the primitive integer polynomial with positive leading coeff."""
    p = Poly(p, x).primitive()[1]
    return -p if p.LC() < 0 else p


def coeff_list(P):
    """Descending integer coefficient list, e.g. x^2-x-1 -> [1, -1, -1]."""
    return [int(c) for c in as_poly(P).all_coeffs()]


def pretty(P):
    """Compact algebraic string, e.g. 'x**2 - x - 1'."""
    return str(as_poly(P).as_expr())


# ----------------------------------------------------------------------------
# the three operators of the emission algebra (Def. 2.1)
# ----------------------------------------------------------------------------


def tensor(P, Q):
    """Tensor (x): roots {lambda_i * mu_j}.

    (P (x) Q)(x) = Res_y( Q(y),  y^{deg P} P(x/y) ).
    """
    P, Q = as_poly(P), as_poly(Q)
    dP = P.degree()
    Pxy = expand(y**dP * P.as_expr().subs(x, x / y))
    R = expand(resultant(Q.as_expr().subs(x, y), Pxy, y))
    return normalize(R)


def sq(P):
    """Squaring ( )^2: roots {lambda_i^2}.

    (P^2)(x) = Res_y( P(y),  x - y^2 ).
    """
    P = as_poly(P)
    R = expand(resultant(P.as_expr().subs(x, y), x - y**2, y))
    return normalize(R)


def dsum(P, Q):
    """Direct sum (+): union of root multisets = polynomial product."""
    return normalize(as_poly(P).as_expr() * as_poly(Q).as_expr())


# ----------------------------------------------------------------------------
# factorization, roots, Mahler measure
# ----------------------------------------------------------------------------


def irred_factors(P):
    """List of (primitive irreducible factor over Q, multiplicity)."""
    _, facs = sp.factor_list(as_poly(P).as_expr(), x)
    return [(normalize(f), m) for f, m in facs]


def _factor_roots(fp):
    coeffs = [mp.mpf(int(c)) for c in fp.all_coeffs()]
    return mp.polyroots(coeffs, maxsteps=1000, extraprec=500)


def roots_mp(P):
    """All complex roots (with multiplicity), computed factor-by-factor."""
    out = []
    for f, m in irred_factors(P):
        rr = _factor_roots(f)
        for _ in range(m):
            out.extend(rr)
    return out


def mahler(P):
    """Mahler measure |lead| * prod_{|root|>1} |root|, exact via factors.

    Computed factor-by-factor so repeated roots (e.g. (x^2+5)^4) are exact.
    """
    m = mp.mpf(1)
    for f, mult in irred_factors(P):
        fm = mp.mpf(abs(int(f.LC())))
        for r in _factor_roots(f):
            if abs(r) > 1:
                fm *= abs(r)
        m *= fm ** mult
    return m


# ----------------------------------------------------------------------------
# the Z/4Z angle charge (Def. 2.1)
# ----------------------------------------------------------------------------

_HALF_PI = mp.pi / 2
_LATTICE_TOL = mp.mpf(10) ** (-20)


def charge_of_root(r):
    """Angle charge arg(r)/(pi/2) mod 4, or 'OFF' if off the Z/4Z lattice."""
    a = mp.arg(r) % (2 * mp.pi)
    q = a / _HALF_PI
    nq = mp.nint(q)
    if abs(q - nq) > _LATTICE_TOL:
        return "OFF"
    return int(nq) % 4


def charges(P):
    """Sorted charge multiset of P (integers 0..3, with any 'OFF' last)."""
    out = [charge_of_root(r) for r in roots_mp(P)]
    return sorted(out, key=lambda v: (v == "OFF", v))


def on_circle_roots(P, tol=mp.mpf(10) ** (-22)):
    """Roots with modulus 1 (within tol)."""
    return [r for r in roots_mp(P) if abs(abs(r) - 1) < tol]


def is_root_of_unity(r, maxord=24, tol=mp.mpf(10) ** (-17)):
    """True if r is a modulus-1 root of unity of order <= maxord."""
    if abs(abs(r) - 1) > tol:
        return False
    return any(abs(r**k - 1) < tol for k in range(1, maxord + 1))


# ----------------------------------------------------------------------------
# canonical seeds and constants
# ----------------------------------------------------------------------------


def phi_seed():
    """Minimal polynomial of the golden ratio: x^2 - x - 1."""
    return as_poly(x**2 - x - 1)


def K_seed():
    """The quartic 'Lorentzian' generator K = x^4 + 5 x^2 - 5."""
    return as_poly(x**4 + 5 * x**2 - 5)


def PHI():
    """The golden ratio phi = (1 + sqrt 5)/2, at mp.dps precision."""
    return (1 + mp.sqrt(5)) / 2


# ----------------------------------------------------------------------------
# small provenance helpers used by every producer
# ----------------------------------------------------------------------------

SOURCE_PAPER = "papers/2026-06-generative-emptiness/generative_emptiness.tex"


def provenance(script_name):
    """Return the two JSON provenance fields for a producer script."""
    return {
        "_source_paper": SOURCE_PAPER,
        "_generated_by": f"code/2026-06-generative-emptiness/{script_name}",
    }


def csv_header(script_name):
    """Return the leading provenance comment line for a CSV file."""
    return (
        f"# source: {SOURCE_PAPER}; "
        f"generated by: code/2026-06-generative-emptiness/{script_name}"
    )


def mpf_str(v, n=30):
    """Fixed-precision string for an mpmath value (stable across platforms)."""
    return mp.nstr(mp.mpf(v), n)


if __name__ == "__main__":
    # tiny self-check that the engine loads and the seeds behave; emits nothing
    print("phi seed  :", pretty(phi_seed()), "  charges", charges(phi_seed()))
    print("K   seed  :", pretty(K_seed()), "  charges", charges(K_seed()))
    print("M(phi)    :", mpf_str(mahler(phi_seed()), 12))
