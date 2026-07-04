"""
Object V -- the conserved current (Prop 6.1 `prop:current` block-diagonal
evolution, Prop 6.2 `prop:clean` clean radial growth).

Claims verified independently:
  * block-diagonal / conservation: no operator word carries amplitude into the
    off-charge (Salem) sector -- every orbit object stays on the Z/4Z lattice;
  * clean radial growth: under iteration every Mahler measure lands in
    {1} U [phi,infty) -- never the forbidden band (1,phi) -- and every on-circle
    root stays a root of unity;
  * the specific listed orbit measures reproduce:
        phi, phi^2, phi^4, 46.98, 76.63, 122.99, 8049.92.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import ge_helpers as H
from ge_helpers import (phi_seed, K_seed, tensor, sq, dsum, mahler, charges,
                        on_circle_roots, is_root_of_unity, PHI)

PHI_V = PHI()


def _orbit():
    """Two-generation orbit generated from the seeds {phi, K}."""
    seeds = [phi_seed(), K_seed()]
    objs = list(seeds)
    gen1 = []
    for P in seeds:
        gen1.append(sq(P))
        for Q in seeds:
            gen1.append(tensor(P, Q))
            gen1.append(dsum(P, Q))
    objs += gen1
    for P in gen1:
        objs.append(sq(P))
        objs.append(tensor(P, phi_seed()))
    return objs


def test_block_diagonal_no_offcharge_transition():
    """Prop 6.1: 'no admissible dynamics transitions into the off-charge sector.'
    Every object in the orbit keeps all root arguments on the Z/4Z lattice."""
    for P in _orbit():
        assert 'OFF' not in charges(P)


def test_clean_radial_growth_measures_avoid_forbidden_band():
    """Prop 6.2: 'every Mahler measure stays in {1} U [phi,infty) -- never the
    forbidden band (1,phi).'"""
    for P in _orbit():
        m = mahler(P)
        assert abs(m - 1) < mp.mpf(10)**(-20) or m >= PHI_V - mp.mpf(10)**(-20), \
            f"measure {m} fell in the forbidden band"


def test_clean_radial_growth_on_circle_roots_are_roots_of_unity():
    """Prop 6.2: '...and every on-circle root stays a root of unity.'"""
    for P in _orbit():
        for r in on_circle_roots(P):
            assert is_root_of_unity(r), f"on-circle root {r} is not a root of unity"


# ---- the seven listed orbit measures, each rebuilt from an operator word ----

def test_listed_measure_phi():
    """Prop 6.2 list, value 1: phi = M(phi seed)."""
    assert abs(mahler(phi_seed()) - PHI_V) < mp.mpf(10)**(-30)
    assert round(float(mahler(phi_seed())), 3) == 1.618


def test_listed_measure_phi_squared():
    """Prop 6.2 list, value 2: phi^2 = M(square(phi))."""
    assert abs(mahler(sq(phi_seed())) - PHI_V**2) < mp.mpf(10)**(-30)


def test_listed_measure_phi_fourth():
    """Prop 6.2 list, value 3: phi^4 = M(square^2(phi))."""
    assert abs(mahler(sq(sq(phi_seed()))) - PHI_V**4) < mp.mpf(10)**(-30)


def test_listed_measure_46_98():
    """Prop 6.2 list, value 4: 46.98 = phi^8 = M(phi^4 (x) phi^4)."""
    phi4 = sq(sq(phi_seed()))
    m = mahler(tensor(phi4, phi4))
    assert abs(m - PHI_V**8) < mp.mpf(10)**(-25)
    assert round(float(m), 2) == 46.98


def test_listed_measure_76_63():
    """Prop 6.2 list, value 5: 76.63 = M(phi (x) K)."""
    m = mahler(tensor(phi_seed(), K_seed()))
    assert round(float(m), 2) == 76.63


def test_listed_measure_122_99():
    """Prop 6.2 list, value 6: 122.99 = phi^10, realised as
    M( square(phi) (+) [phi^4 (x) phi^4] ) = phi^2 * phi^8."""
    phi4 = sq(sq(phi_seed()))
    obj = dsum(sq(phi_seed()), tensor(phi4, phi4))
    m = mahler(obj)
    assert abs(m - PHI_V**10) < mp.mpf(10)**(-20)
    assert round(float(m), 2) == 122.99


def test_listed_measure_8049_92():
    """Prop 6.2 list, value 7: 8049.92 = M(phi (x) phi (x) K)."""
    m = mahler(tensor(tensor(phi_seed(), phi_seed()), K_seed()))
    assert round(float(m), 2) == 8049.92
