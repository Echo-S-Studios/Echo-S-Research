"""Pure arithmetic / counting consistency of the tallies stated across the paper.

These are the mechanical box-size, Burnside, ledger, and next-step arithmetic
claims -- each derived independently and checked against the paper's number.
"""


def test_quintic_box_size():
    """Sec. 6: [-2,2]^5 has 5^5 = 3125 candidates."""
    assert 5 ** 5 == 3125


def test_quintic_reject_tally_sums_to_box():
    """Sec. 6 table: 625+50+638+1318+411 rejects + 83 Pisot = 3125."""
    assert 625 + 50 + 638 + 1318 + 411 + 83 == 3125


def test_quintic_pattern_split():
    """Sec. 6: real5=0, mixed=16, two-pair=67, total 83."""
    assert 0 + 16 + 67 == 83


def test_quintic_scan_sizes():
    """Sec. 6: deg Rat = 25 = 5^2, bound 1250 = 2*5^4, C_2 degree 400 = 20^2."""
    assert 5 ** 2 == 25
    assert 2 * 5 ** 4 == 1250
    assert 20 ** 2 == 400


def test_salem_burnside_arithmetic():
    """Sec. 7.1: 729 = 3^6 vectors, 27 = 3^3 twist-fixed, orbits (729+27)/2 = 378."""
    assert 3 ** 6 == 729
    assert 3 ** 3 == 27
    assert (729 + 27) // 2 == 378
    assert (729 + 27) % 2 == 0


def test_salem_cascade_sum():
    """Sec. 7.1: 39 + 257 + 45 + 37 = 378."""
    assert 39 + 257 + 45 + 37 == 378


def test_salem_scan_sizes():
    """Sec. 7.1: deg Rat = 144 = 12^2, bound 41472 = 2*12^4."""
    assert 12 ** 2 == 144
    assert 2 * 12 ** 4 == 41472


def test_quartic_box_and_split():
    """Sec. 7.2: [-3,3]^4 = 7^4 = 2401 candidates; 102 complex + 1 real = 103."""
    assert 7 ** 4 == 2401
    assert 102 + 1 == 103


def test_ledgerM_signature_and_scan_counts():
    """Sec. 7.3: 14 canonical + 37 census + 103 quartic + 83 quintic = 237
    signatures; two decision paths each => 474 = 2*237 scan executions."""
    assert 14 + 37 + 103 + 83 == 237
    assert 2 * 237 == 474


def test_lehmer_scan_sizes():
    """Sec. 2: Lehmer deg Rat = 100 = 10^2, bound 20000 = 2*10^4."""
    assert 10 ** 2 == 100
    assert 2 * 10 ** 4 == 20000


def test_N3_census_extension_burnside():
    """N3 / P5: c in {-2,..,2}^6 has 5^6 = 15625 vectors, 5^3 = 125 twist-fixed,
    Burnside (15625+125)/2 = 7875 twist-classes (forced arithmetic)."""
    assert 5 ** 6 == 15625
    assert 5 ** 3 == 125
    assert (15625 + 125) // 2 == 7875
    assert (15625 + 125) % 2 == 0
