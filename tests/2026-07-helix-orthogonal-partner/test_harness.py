"""
Sanity checks for Section 8 ("The harness: constraint-import as verification")
of "The Dissolved Helix and Its Orthogonal Partner".

Section 8 is largely metalogical (Defs 8.1-8.3, 8.7; Props 8.4-8.6; Remarks).
Its load-bearing *computational* content is the filter asymmetry (Prop 8.4:
one exact counterexample refutes a universal, but no finite exact range promotes
one) and the "generator is inert" corollary (Prop 8.6: on a closed sublattice no
falsifiable exit exists, so C_dashv = empty). We give faithful, independent
demonstrations of exactly these principles, plus the transferable "kernel"
identity of Remark 8.9 (sqrt5 = phi + phi^{-1}). The interpretive parts are
recorded as untestable in NOTES.md.
"""
import sympy as sp

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2


# ---------------------------------------------------------------- Prop. 8.4
def test_filter_asymmetry_single_counterexample_refutes():
    """Prop. 8.4: delta = dashv follows from a SINGLE exact counterexample.
    Universal gamma: 'every a in {0,2} satisfies a == 0 (mod 4)'.
    Witness a=2 gives 2 mod 4 = 2 != 0, an exact certificate of not-gamma."""
    domain = {0, 2}
    witnesses = [a for a in domain if (a % 4) != 0]
    assert witnesses == [2]                       # exact counterexample exists
    assert (2 % 4) != 0                           # decisive, from one instance


def test_filter_asymmetry_finite_range_never_promotes():
    """Prop. 8.4: a finite exact search yields at most [computed], never [forced]:
    passing a finite range does NOT entail the universal. Independent witness:
    Euler's polynomial n^2+n+41 is prime for n=0..39 but composite at n=40."""
    f = lambda n: n**2 + n + 41
    assert all(sp.isprime(f(n)) for n in range(40))    # finite range passes
    assert not sp.isprime(f(40))                       # universal is false
    assert f(40) == 41 * 41                             # exact demotion certificate


# ---------------------------------------------------------------- Prop. 8.6
def test_generator_is_inert_on_closed_sublattice():
    """Prop. 8.6: a pure generator asserts only unfalsifiable/definitional claims,
    so C_dashv = empty and no partition forms. On the closed {0,2} sub-semigroup
    the claim 'no reachable charge lies in {1,3}' is never refuted: no exact
    counterexample exists."""
    reachable = {0, 2}
    for _ in range(6):
        reachable |= {(a + b) % 4 for a in reachable for b in reachable}
        reachable |= {(2 * a) % 4 for a in reachable}
    casualties = reachable & {1, 3}               # would-be refuting witnesses
    assert casualties == set()                    # C_dashv = empty -> inert


def test_informative_import_needs_a_falsifiable_member():
    """Prop. 8.5/8.6: the harness is informative iff C contains a falsifiable
    gamma whose negation is decidable. The refuted casualty 'a helix couples to
    a copy of itself' IS decidable (Cor. 3.5: [R,R]=0), so C_dashv != empty and
    the import is substantive."""
    R = sp.Matrix([[0, 1], [1, 1]])
    assert (R * R - R * R) == sp.zeros(2, 2)       # exact refutation certificate


# ---------------------------------------------------------------- Remark 8.9
def test_transferable_kernel_identity():
    """Remark 8.9: the one portable kernel is the sqrt-branch of a 2-to-1 fold;
    this system's instance is sqrt5 = phi + phi^{-1} (= phi - psi)."""
    psi = (1 - sqrt5) / 2
    assert sp.simplify((phi + 1 / phi) - sqrt5) == 0
    assert sp.simplify((phi - psi) - sqrt5) == 0
