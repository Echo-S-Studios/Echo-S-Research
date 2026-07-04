"""
Object I -- the Z/4Z grading (Thm 2.2 `thm:charge`, Def 2.1 `def:charge`).

Paper claims, verified here by INDEPENDENTLY computing root arguments:
  * every emittable object carries charge in Z/4Z (arg in (pi/2)Z);
  * the operators realise the group law:  (x) = add, ()^2 = double, (+) = union;
  * Z/4Z is closed under +, x2, union;
  * charge multisets:  phi:{0,2},  phi(x)phi:{0,0,2,2},  phi^2:{0,0},  K:{0,1,2,3}
    (the quartic generator K realises the FULL group; its imaginary pair at +-pi/2).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ge_helpers as H
from ge_helpers import phi_seed, K_seed, tensor, sq, dsum, charges


# ----- the four stated charge multisets (Thm 2.2 parenthetical) -----

def test_charge_multiset_phi():
    """Thm 2.2: 'phi:{0,2}'.  phi = x^2-x-1 has a positive root (arg 0, q=0)
    and a negative root (arg pi, q=2)."""
    assert charges(phi_seed()) == [0, 2]


def test_charge_multiset_phi_tensor_phi():
    """Thm 2.2: 'phi(x)phi:{0,0,2,2}'."""
    assert charges(tensor(phi_seed(), phi_seed())) == [0, 0, 2, 2]


def test_charge_multiset_phi_squared():
    """Thm 2.2: 'phi^2:{0,0}' -- squaring sends both real roots to positives."""
    assert charges(sq(phi_seed())) == [0, 0]


def test_charge_multiset_K_full_group():
    """Thm 2.2: 'K:{0,1,2,3}'.  K=x^4+5x^2-5 realises the FULL Z/4Z:
    a real pair at args 0,pi (charges 0,2) and an imaginary pair +-i*beta at
    args +-pi/2 (charges 1,3)."""
    assert charges(K_seed()) == [0, 1, 2, 3]


# ----- the operators act as add / double / union ON the charge -----

def test_tensor_acts_as_addition_on_charge():
    """Thm 2.2: '(x): add (arg lambda*mu = arg lambda + arg mu)'.
    The charge multiset of phi(x)K must equal the pairwise (a+b) mod 4."""
    a = charges(phi_seed())
    b = charges(K_seed())
    predicted = sorted((i + j) % 4 for i in a for j in b)
    assert charges(tensor(phi_seed(), K_seed())) == predicted


def test_squaring_acts_as_doubling_on_charge():
    """Thm 2.2: '()^2: double (arg lambda^2 = 2 arg lambda)'."""
    a = charges(K_seed())
    predicted = sorted((2 * i) % 4 for i in a)
    assert charges(sq(K_seed())) == predicted


def test_dsum_acts_as_union_on_charge():
    """Thm 2.2: '(+): union'.  charge(phi (+) K) = charge(phi) U charge(K)."""
    predicted = sorted(charges(phi_seed()) + charges(K_seed()))
    assert charges(dsum(phi_seed(), K_seed())) == predicted


# ----- Z/4Z is closed under the three operations (group-theoretic fact) -----

def test_Z4_closed_under_add_double_union():
    """Thm 2.2: 'internal to Z/4Z since it is closed under +, x2, and union'."""
    G = {0, 1, 2, 3}
    assert {(a + b) % 4 for a in G for b in G} == G     # addition closed
    assert {(2 * a) % 4 for a in G} <= G                # doubling maps in
    assert (G | G) == G                                 # union closed


# ----- no finite operator word escapes the lattice (charge is conserved) -----

def test_no_word_leaves_Z4_lattice():
    """Thm 2.2 discussion: 'no finite word in the operators leaves Z/4Z.'
    Generate a two-generation orbit from {phi,K}; assert NO object ever
    acquires an off-lattice ('OFF') root argument."""
    seeds = [phi_seed(), K_seed()]
    gen1 = []
    for P in seeds:
        gen1.append(sq(P))
        for Q in seeds:
            gen1.append(tensor(P, Q))
            gen1.append(dsum(P, Q))
    orbit = seeds + gen1
    for P in gen1[:6]:            # one more generation (kept small for speed)
        orbit.append(sq(P))
        orbit.append(tensor(P, phi_seed()))
    for P in orbit:
        assert 'OFF' not in charges(P)


# ----- K is a genuine irreducible quartic (needed to realise charges 1,3) -----

def test_K_is_irreducible_quartic():
    """Def/Thm: K is 'the quartic generator'.  A reducible K could not supply
    the +-i*beta place carrying charges {1,3}."""
    import sympy as sp
    from ge_helpers import x
    facs = sp.factor_list(K_seed().as_expr(), x)[1]
    assert len(facs) == 1 and facs[0][1] == 1
    assert K_seed().degree() == 4
