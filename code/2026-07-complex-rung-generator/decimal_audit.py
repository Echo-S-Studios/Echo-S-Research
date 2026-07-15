#!/usr/bin/env python3
# decimal_audit.py -- handoff SS1 decimal audit of register literals.
# Fresh mpmath at 50 dps. Matcher spec (handoff SS1): a printed literal MATCHES iff
#   (a) |literal - value| <= 0.51 ulp at the printed precision, OR
#   (b) the literal is the exact floor-truncation of the higher-precision value.
# Both directions: paper->recompute (literal matches fresh value) and
# recompute->paper (round-half-up and truncation renderings reproduce the literal
# or are reported for inspection).

import sys
from mpmath import mp, mpf, sqrt, log, pi, ellipe, acos, degrees, floor

mp.dps = 50

phi = (1 + sqrt(5)) / 2
tau = 1 / phi
gap = phi ** -4
K = sqrt(1 - gap)
kap = pi / (2 * log(phi))

LITERALS = [
    # (name, printed literal string, fresh 50-dps value)
    ("L_res", "3.7312193125",
     sqrt(1 + kap ** 2) * ellipe(kap ** 2 / (1 + kap ** 2))),
    ("d_1 (deg)", "43.4026803146",
     degrees(acos(5 ** mpf("0.25") * tau ** mpf("1.5")))),
    ("pi*K/2", "1.4516928502",
     pi * K / 2),
    ("K*ln(phi)/(4*pi)", "0.0353900591",
     K * log(phi) / (4 * pi)),
]

def audit(name, lit_str, val):
    dec = len(lit_str.split(".")[1])
    ulp = mpf(10) ** -dec
    lit = mpf(lit_str)
    diff = abs(lit - val)
    within = diff <= mpf("0.51") * ulp
    trunc = floor(val * mpf(10) ** dec) / mpf(10) ** dec
    is_trunc = (lit == trunc)
    rnd = floor(val * mpf(10) ** dec + mpf("0.5")) / mpf(10) ** dec
    is_round = (lit == rnd)
    verdict = "PASS" if (within or is_trunc) else "FAIL"
    print(f"{verdict} {name}")
    print(f"    literal        = {lit_str}  ({dec} dp)")
    print(f"    fresh 50 dps   = {mp.nstr(val, 30)}")
    print(f"    |diff| in ulp  = {mp.nstr(diff / ulp, 6)}  (0.51-ulp: {'yes' if within else 'no'})")
    print(f"    trunc @ {dec}dp   = {mp.nstr(trunc, dec + 4, strip_zeros=False)}  "
          f"(floor-truncation match: {'yes' if is_trunc else 'no'})")
    print(f"    round @ {dec}dp   = {mp.nstr(rnd, dec + 4, strip_zeros=False)}  "
          f"(round-half match: {'yes' if is_round else 'no'})")
    return verdict == "PASS"

ok = True
for name, lit_str, val in LITERALS:
    ok = audit(name, lit_str, val) and ok
    print()

# cross-form corroboration for d_1 (both exact forms, as in NX-H2)
d1a = degrees(acos(5 ** mpf("0.25") * tau ** mpf("1.5")))
d1b = degrees(acos(sqrt(tau) * K))
agree = abs(d1a - d1b) < mpf(10) ** -45
print(("PASS" if agree else "FAIL"),
      "d_1 cross-form: arccos(5^(1/4) tau^(3/2)) vs arccos(sqrt(tau) K) agree < 1e-45")
ok = ok and agree

print()
print("DECIMAL AUDIT:", "ALL PASS" if ok else "FAILURES PRESENT")
sys.exit(0 if ok else 1)
