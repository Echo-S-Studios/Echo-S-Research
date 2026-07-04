# Tests — mechanical verification of the papers' derivations

This tree independently re-derives and checks the mathematical claims of every
paper under [`../papers/`](../papers). Each paper has a sibling folder here with
the same `YYYY-MM-shortname` name, containing its test modules (`test_*.py`),
any helper module, and a `NOTES.md` recording exactly what was verified, what was
flagged, and what could not be mechanically checked.

```
tests/<same-shortname-as-paper>/
    test_*.py     # independent re-derivation tests (pytest)
    NOTES.md      # per-paper claim ledger (verified / flagged / untestable)
    _*.py         # optional per-paper helper module (uniquely named)
```

## How to run

From the repository root (Windows; use the `py` launcher):

```bash
py -m pip install --user sympy mpmath numpy pytest   # one-time: install deps
py -m pytest                     # whole suite (tests/ is the configured testpath)
py -m pytest tests/2026-06-salem-slot -v             # a single paper, verbose
py -m pytest -rx                 # also print the reason for each xfail (the flags)
```

`pytest.ini` (repo root) sets `testpaths = tests` and `--import-mode=importlib`
so identically-named test files coexist across folders.

**Current status: `1041 passed, 6 xfailed, 0 failed, 0 errors`** (~3.3 min).

## Methodology

- **Independent re-derivation, not restatement.** Every test rebuilds a result
  from the paper's *stated premises* — companion matrices from their definitions,
  Mahler measures from 40–60-digit root products, invariant factors from
  determinantal divisors, Fisher matrices from log-likelihood Hessians, Galois
  groups from scratch — and only then compares to the value the paper reports.
  Numbers are never copied from the paper and asserted against themselves.
- **Tools.** `sympy` for exact symbolic algebra, `mpmath` (`dps ≥ 40`) for
  high-precision transcendental/root computations, `numpy` for linear algebra.
  Equality is exact (`simplify(lhs-rhs)==0`) where possible, else within a
  tolerance matched to the paper's stated digits.
- **Cross-domain / deliberately-derived constants** (e.g. a golden-gate
  `c = √5/2`, externally-seeded `φ`/`√5` values, `[computed]`/`[posited]`
  quantities) are **not** asserted as arithmetic; they are tested only for
  internal consistency with the paper's own stated procedure, per the paper's
  own epistemic tags.
- **Proofs / arguments** (universal statements, cited theorems, open problems)
  are exercised by their computable core — limiting/boundary cases, small-N
  instances, algebraic identities the proof relies on — and the non-computable
  remainder is documented as *untestable* in the paper's `NOTES.md`.
- **Flagged discrepancies are isolated, not fixed.** Where an independent
  re-derivation genuinely disagrees with the paper, the test is marked
  `@pytest.mark.xfail(reason=...)` (so the suite stays green and the discrepancy
  is visible) and written up in `NOTES.md` for human review. **No paper was
  modified. No tolerance was loosened to hide a discrepancy.** The 6 xfails below
  are these flags.

## Summary table

Counts are **claim-level** (one mathematical claim may map to several pytest
assertions). "Flagged" = a genuine discrepancy isolated as `xfail` for review.

| Paper (folder) | Claims tested | Passed | Flagged¹ | Untestable² | pytest (pass / xfail) |
|----------------|:-------------:|:------:|:--------:|:-----------:|:---------------------:|
| `2026-06-charge-measure-coupling` | 58 | 58 | 0 | 10 | 58 / 0 |
| `2026-06-emission-gap` | 42 | 42 | 0 | 5 | 47 / 0 |
| `2026-06-generative-emptiness` | 37 | 37 | 0 | 6 | 44 / 0 |
| `2026-06-lambda-2c` | 57 | 57 | 0 | 8 | 72 / 0 |
| `2026-06-lehmers-box` | 44 | 43 | **1** | 6 | 45 / 1 |
| `2026-06-operator-algebra` | 45 | 45 | 0 | 4 | 57 / 0 |
| `2026-06-residual-return-learning` | 72 | 71 | **1** | 10 | 71 / 1 |
| `2026-06-salem-slot` | 39 | 38 | **1** | 5 | 46 / 1 |
| `2026-06-vector-substrate` | 62 | 62 | 0 | 7 | 63 / 0 |
| `2026-06-z5-no-salem-dichotomy` | 50 | 49 | **1** | 5 | 54 / 1 |
| `2026-07-emission-algebra-primer` | 50 | 50 | 0 | 6 | 104 / 0 |
| `2026-07-helix-orthogonal-partner` | 45 | 45 | 0 | 4 | 61 / 0 |
| `2026-07-pisot-residue` | 42 | 41 | **1** | 9 | 51 / 1 |
| `2026-07-relational-charge` | 85 | 84 | **1** | 5 | 268 / 1 |
| **Total** | **728** | **722** | **6** | **90** | **1041 / 6** |

¹ Isolated as `xfail`; see below and the paper's `NOTES.md`. ² Cited external
theorems, open problems, companion-paper results, implementation guarantees, and
interpretive/`[posited]` framing — documented per paper.

## Flagged for human review (the 6 xfails)

All are **wording / bucketing-convention level** — every one has its load-bearing
mathematical content independently verified as *correct*; only a label, sign, or
intermediate tally differs. None affects a paper's conclusions.

1. **`relational-charge` & `pisot-residue` — same degree‑12 Salem census split.**
   Both papers report the intermediate cascade as `257` trace‑Sturm rejects /
   `45` reducible; two independent enumerations reproduce `256 / 46` (relational‑
   charge pins it to exactly one polynomial `x¹²−x¹¹−x¹⁰−x⁹−x⁷−x⁶−x⁵−x³−x²−x+1`,
   which has a `Φ₁₀` factor yet a coincidental `(1,0,5)` trace pattern). **Every
   robust invariant matches** (378 orbits, 37 census Salems, all inert `{Φ₁^12}`,
   combined `302`). Two independent agents flagging the same split is the
   strongest signal here — likely a single shared counting-convention choice.
2. **`z5-no-salem-dichotomy` — Sec. 7 reciprocal count.** Table says `5`
   reciprocal charge‑`Z/5Z` objects; independent enumeration of the paper's exact
   window finds `4`. The distinct measures `{1, φ², 2+√3}` **and** the
   non‑reciprocal count `13` reproduce exactly (which is what makes the off‑by‑one
   notable but benign).
3. **`salem-slot` — Prop 6.4 quadratic sign.** The `O((φ−β)²)` term of
   `√5 − τ₀(β)` is displayed with `+`; independent Taylor expansion gives `−`
   (`= −φ⁻³`). The curvature *magnitude* (`√5−2 = φ⁻³`) and the linear slope
   `φ⁻¹` are correct; the geometric-rate results are unaffected. Typo‑level.
4. **`lehmers-box` — Thm 5.3 note.** The note says "all `27` subgroups" of
   `G = C₂×C₂×D₄`; `G` actually has `158` subgroups. The `27` correctly counts the
   subgroups fixing a subfield of `K` (i.e. the `27` subfields), whose signature
   census is reproduced exactly — only the word "subgroups" is imprecise.
5. **`residual-return-learning` — Prop 4.7 aside.** Calls `R₁ = [[0,1],[1,−1]]`
   "conjugate to" the keystone `R = [[0,1],[1,1]]`, but they have different
   characteristic polynomials (`x²+x−1` vs `x²−x−1`) so are not similar
   (`R₁ ∼ −R`). The load-bearing shared self-action gap `√5` (both discriminant 5)
   is correct.

*(6 flags across 5 issues — items 1 counts twice because the same census split is
flagged independently in both `pisot-residue` and `relational-charge`.)*

## What "untestable" means here

The 90 untestable claims are **not** failures — they are claims a finite
computation cannot settle: cited external theorems (Kronecker, Smyth, Dobrowolski,
Cencov, Yuzvinskii–Bowen, Siegel), acknowledged open problems (Lehmer's problem
and conjecture, general Pisot inertness), results proved in companion papers,
implementation/protocol guarantees, and explicitly `[posited]`/interpretive
framing. Where such a claim has a computable shadow (a specific instance, a
boundary case, an algebraic identity), that shadow **was** verified; each paper's
`NOTES.md` records which.
