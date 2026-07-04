# Code — The Dissolved Helix and Its Orthogonal Partner

Producer scripts for **`papers/2026-07-helix-orthogonal-partner/helix_orthogonal_partner.tex`**
(*The Dissolved Helix and Its Orthogonal Partner — Terrain, rotation, and the completion of Z/4Z*).

Each script recomputes one section's results **from the paper's stated premises**
(the Fibonacci companion `R = [[0,1],[1,1]]`; `phi, psi` the roots of `x^2-x-1`;
the quartic `Kf = x^4+5x^2-5`; the registry quartics `cons`, `res`) using exact
`sympy` arithmetic over Q, Q(√5), Q(5^{1/4}), with `mpmath` for decimal display
and independent high-precision Mahler cross-checks. **No float crosses a decision
boundary.** Every run **emits** machine-readable data (it does not assert): the
scripts are producers, the files in `tests/` are independent verifiers, and the
two never share code. Emitted `*_residual` fields are the produced evidence — an
exact quantity minus the paper's closed form, which is `0` when they agree.

Run any script from the repository root; outputs land in
`data/2026-07-helix-orthogonal-partner/`.

## Helpers (not run directly)

| module | role |
|--------|------|
| `helix_core.py` | shared kernel: golden data (`R, phi, psi, sqrt5`), the two characters `M` / `chi` (Def. 2.4), the `ad_M` builder, exact Mahler / reciprocity / parity helpers |
| `helix_io.py` | provenance-stamped CSV/JSON writers (leading `# source:` comment; JSON `_source_paper` / `_generated_by`) |

## Producers

| script | run command | paper result produced | data artifact(s) |
|--------|-------------|-----------------------|------------------|
| `keystone.py` | `py code/2026-07-helix-orthogonal-partner/keystone.py` | §2 Prop. 2.2 (keystone: charpoly, Tr=1, det=−1, spec{φ,ψ}, R²=R+I, eq. 1 √5=φ+φ⁻¹), Prop. 2.3 (spec(ad_R)={0,±√5}), Def. 2.4 (golden object M=φ, χ={0,2}) | `keystone.json`, `ad_spectrum.csv` |
| `obstruction.py` | `py code/2026-07-helix-orthogonal-partner/obstruction.py` | §3 Prop. 3.1 ([R,R]=0), Lem. 3.2 (M grows as φ^(2^k), never 1), Prop. 3.3 ({0,2} closed, {1,3} unreachable), Prop. 3.4 (empty Salem slot, τ₀), Cor. 3.5 | `obstruction.json`, `mahler_growth.csv` |
| `flip.py` | `py code/2026-07-helix-orthogonal-partner/flip.py` | §4 Def. 4.1 / Thm. 4.2 (D=1+4C; Gram [[2,−1],[−1,1+2C]], det G=D; signature flip (2,0)→(1,1); rotation channel {0,±i√\|D\|}; golden face C=1⇒D=5) | `flip.json`, `flip_family.csv` |
| `kformation.py` | `py code/2026-07-helix-orthogonal-partner/kformation.py` | §5 Prop. 5.1 / eqs. 6–8 (Kf; y±=(−5±3√5)/2; K=√y₊=5^{1/4}/φ; β=√\|y₋\|; Mah(Kf)=β²=φ²√5) | `kformation.json` |
| `partner.py` | `py code/2026-07-helix-orthogonal-partner/partner.py` | §6 Prop. 6.1 (χ(Kf)={0,1,2,3}, imag. roots Re=0), Prop. 6.2 (parity criterion; Q(√5) factorisations; Re parts φ², −φ; Q(5^{1/4}); sig (2,1); Galois D₄ order 8; non-reciprocal; Mahler φ²√5, 11+5√5, (43+19√5)/2) | `partner.json`, `registry_quartics.csv` |
| `synthesis.py` | `py code/2026-07-helix-orthogonal-partner/synthesis.py` | §7 Table 1 (the [posited] helix dictionary over the [forced] substrate), Thm. 7.2 ({0,2}⊔{1,3}=Z/4Z) | `synthesis.json`, `helix_dictionary.csv` |
| `harness.py` | `py code/2026-07-helix-orthogonal-partner/harness.py` | §8 Prop. 8.4 (filter asymmetry), Prop. 8.6 (generator inert), Rem. 8.9 (kernel √5=φ+φ⁻¹), Table 2 (survivors/casualties/residue ledger) | `harness.json`, `harness_ledger.csv` |

`partner.py` is the slowest (three `sympy.galois_group` calls); the others are
near-instant. Regenerate everything with:

```
for s in keystone obstruction flip kformation partner synthesis harness; do
  py code/2026-07-helix-orthogonal-partner/$s.py
done
```

## How these differ from `tests/2026-07-helix-orthogonal-partner/`

The tests **assert** each paper value (`assert ... == 0`) and stop; these
producers **compute and serialise** the same objects as reusable data — spectra,
closed forms, decimals, invariant tables, the two paper tables, and the exact
residuals that certify the match. The two trees share no code (producers import
only `helix_core` / `helix_io`, never `tests/`), so they are genuinely
independent reconstructions of the same substrate.
