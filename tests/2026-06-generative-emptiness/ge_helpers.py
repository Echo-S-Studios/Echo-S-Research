"""
Shared re-derivation helpers for the independent verification of
    papers/2026-06-generative-emptiness/generative_emptiness.tex
    "The Generative Content of a Conserved Emptiness".

Nothing here restates a paper value; these are neutral primitives (operator
constructions on Z[x], Mahler measure, the angle charge) that the test files
use to REBUILD every claimed object from the paper's stated premises.

Operators (Def. 2.1 "emission algebra"):
  * tensor  (x)  : polynomial whose roots are the pairwise products  lambda*mu
  * squaring ( )^2: polynomial whose roots are the squares            lambda^2
  * dsum    (+)  : direct sum = polynomial product (union of root sets)

Angle charge (Def. 2.1): the multiset of  arg(root)/(pi/2) mod 4.
"""
import sympy as sp
from sympy import symbols, Poly, resultant, expand
import mpmath as mp

mp.mp.dps = 50
x, y = symbols('x y')


def as_poly(e):
    return e if isinstance(e, Poly) else Poly(e, x)


def _norm(p):
    """primitive integer poly, positive leading coeff"""
    p = Poly(p, x).primitive()[1]
    return -p if p.LC() < 0 else p


def tensor(P, Q):
    """Roots {lambda_i * mu_j}.  (P (x) Q)(x) = Res_y( Q(y), y^{deg P} P(x/y) )."""
    P, Q = as_poly(P), as_poly(Q)
    dP = P.degree()
    Pxy = expand(y**dP * P.as_expr().subs(x, x / y))
    R = expand(resultant(Q.as_expr().subs(x, y), Pxy, y))
    return _norm(R)


def sq(P):
    """Roots {lambda_i^2}.  (P^2)(x) = Res_y( P(y), x - y^2 )."""
    P = as_poly(P)
    R = expand(resultant(P.as_expr().subs(x, y), x - y**2, y))
    return _norm(R)


def dsum(P, Q):
    """Direct sum = product of polynomials (union of root multisets)."""
    return _norm(as_poly(P).as_expr() * as_poly(Q).as_expr())


def irred_factors(P):
    """List of (integer-primitive irreducible factor, multiplicity) over Q."""
    _, facs = sp.factor_list(as_poly(P).as_expr(), x)
    return [(_norm(f), m) for f, m in facs]


def _factor_roots(fp):
    coeffs = [mp.mpf(int(c)) for c in fp.all_coeffs()]
    return mp.polyroots(coeffs, maxsteps=1000, extraprec=500)


def roots_mp(P):
    """All complex roots (with multiplicity), computed factor-by-factor for accuracy."""
    out = []
    for f, m in irred_factors(P):
        rr = _factor_roots(f)
        for _ in range(m):
            out.extend(rr)
    return out


def mahler(P):
    """Mahler measure = |lead| * prod_{|root|>1} |root|, exact via irreducible factors."""
    m = mp.mpf(1)
    for f, mult in irred_factors(P):
        fm = mp.mpf(abs(int(f.LC())))
        for r in _factor_roots(f):
            if abs(r) > 1:
                fm *= abs(r)
        m *= fm ** mult
    return m


def charges(P):
    """Multiset (sorted list) of angle charges arg/(pi/2) mod 4; 'OFF' if a root
    argument is not an integer multiple of pi/2 (i.e. off the Z/4Z lattice)."""
    out = []
    half = mp.pi / 2
    for r in roots_mp(P):
        a = mp.arg(r) % (2 * mp.pi)
        q = a / half
        nq = mp.nint(q)
        if abs(q - nq) > mp.mpf(10) ** (-20):
            out.append('OFF')
        else:
            out.append(int(nq) % 4)
    return sorted(out, key=lambda v: (v == 'OFF', v))


def on_circle_roots(P, tol=mp.mpf(10) ** (-25)):
    return [r for r in roots_mp(P) if abs(abs(r) - 1) < tol * 1000]


def is_root_of_unity(r, maxord=24, tol=mp.mpf(10) ** (-20)):
    if abs(abs(r) - 1) > tol * 1000:
        return False
    return any(abs(r ** k - 1) < tol * 1000 for k in range(1, maxord + 1))


# canonical seeds (Def. 2.1 / Thm 2.2)
def phi_seed():
    return as_poly(x**2 - x - 1)          # minimal poly of the golden ratio


def K_seed():
    return as_poly(x**4 + 5 * x**2 - 5)   # the quartic "Lorentzian" generator K


def PHI():
    return (1 + mp.sqrt(5)) / 2
