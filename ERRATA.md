# Errata register

Adjudications of errata reported against the corpus. Each entry records what
was alleged, what an audit of this repository found, and where the erratum
actually applies. Corrections whose targets exist in this lineage are fixed
in the papers directly (see `tests/README.md` for the six adjudicated-and-
corrected flags, and the git history for the fix commits); this file tracks
reports whose targets do **not** exist here.

An external review batch (2026-07) raised four items, D1 and E1–E3. D1 — a
stale "twenty-three cross-engine signatures" count in the relational-charge
Methodology replication list — was confirmed in-repo and corrected (the bullet
now reads "all twenty-six cross-engine items of ledger entry M (twenty-four
signatures, two factorization checks)", matching the paper's own Second-engine
paragraph and Appendix A entry M). E1–E3 were adjudicated as follows.

## E1 — degree-2 census 6→7; prediction P4 12→13; prediction P3 2661→2662

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

E1 alleged three wrong numerals in the pisot/relational lineage: (a) a
degree-2 census count printed 6 that should be 7; (b) a prediction "P4" count
12 that should be 13; (c) a prediction "P3" count 2661 that should be 2662.
Audit findings:

- Neither `2661` nor `2662` appears in any `papers/**/*.tex`, in the committed
  deposit builds (`pisot_residue_whitepaper` v13, `relational_charge_paper`
  rev. l), or in `tests/`.
- Both P3 and P4 in this lineage (`pisot_residue_whitepaper.tex` §"Predictions",
  `relational_charge_paper.tex` open-problems list) are count-free statements.
- The only degree-2 census count anywhere in the archive
  (`lehmers_box.tex`, subfield-signature census, the (2,2,0) column) already
  reads 7.

E1 therefore targets an external Zenodo deposit revision that was never
committed to this repository. Action: none in-repo; if the affected deposit
version is identified, file the erratum against that DOI version record.

## E2 — degree-7 Pisot polynomial display vs. its c-vector

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

No degree-7 polynomial display and no "c-vector" notation exists in any
committed `.tex`, PDF, test, or data file. Every `x^7` occurrence is an
interior term of a higher-degree polynomial (Lehmer's degree-10 polynomial;
the degree-12 reducible twist-class). Degree 7 appears only as prose future
work (pisot-residue item N4). If this erratum is real it targets an external
Zenodo deposit revision only.

## E3 — "both-straddle" count 163 vs. 160

**Not applicable to this lineage** (adjudicated 2026-07-14 at HEAD `2a6b4be`).

The token "both-straddle" appears nowhere in the repository; every "straddle"
in this lineage is a boolean predicate (the Salem flip-straddle; K-formation
y-roots straddling zero), never a counted category. 163 occurs only as an RGB
color component in paper preambles; no 160/163 tally pair exists in any
`.tex`, committed PDF, test, or data file. If this erratum is real it targets
an external Zenodo deposit revision only.

## D5 — pisot-residue §6: `Rat°` irreducible on 66/67, not all 67 (in-lineage)

**Applies to this lineage; correction pending an owner decision** (finding
verified 2026-07-27, from the Exact Arithmetic booklet's fourth audit; see
`papers/2026-07-exact-arithmetic/`).

`pisot_residue_whitepaper.tex` §"Decision (stage 2)" states that on the 67
two-pair instances the cross-shell residue `Rat°_p` was "squarefree **and in
fact irreducible of degree 20 on all 67**", and Observation `obs:generic` tags
the Galois-genericity of that as `[PLAUSIBLE]`. The squarefree-on-all-67 claim
(and Theorem `thm:box`, which asserts only squarefree `67/67`) is **correct**.
The **irreducibility** claim is off by one: `Rat°_p` is irreducible on **66 of
67**, the sole exception being

    p = x^5 - x^3 - 2x^2 - 2x - 1,   Gal(p) = D_5,

whose `Rat°_p` splits into two self-reciprocal degree-10 factors. Certified by
exact re-multiplication (engine-independent), and independently reproduced here:
`sympy` gives `Gal(p) = D_5` (order 10; discriminant 2209 = 47², so `Gal ⊆ A_5`;
`p` irreducible over ℚ), and the shipped `code/2026-07-exact-arithmetic/
lesson6_census.py --full` re-derives it (its D3 check now encodes the corrected
count rather than the source note's "all 67").

**Theorem `thm:box` is untouched.** Both degree-10 factors carry 2 unimodular
roots, so `S* = Rat°_p` remains the full degree-20 product, `deg C_2 = 400`, and
the complete cyclotomic scan still returns `{Φ_1^20}`: the `D_5` instance is
relationally inert, exactly as the theorem concludes.

**The upgrade worth folding in** (supersedes `obs:generic`'s `[PLAUSIBLE]` with a
`[FORCED]` criterion): for squarefree `Rat°_p`, irreducibility over ℚ is
transitivity on the `n(n−1)` ordered distinct root pairs, i.e. **iff `Gal(p)` is
2-transitive**. The Galois census of the 67 is **S₅ × 65, M₂₀ × 1, D₅ × 1**;
`M₂₀ = AGL(1,5)` (order 20) is sharply 2-transitive → irreducible, while `D₅`
(order 10) has two length-10 orbits → the observed split — reproducing 66/67 with
no residue. Consequence for the degree 6–7 frontier (pisot-residue item N4):
irreducible `Rat°` needs `n(n−1) ∣ |G|`, i.e. `30 ∣ |G|` at degree 6 (only
PSL(2,5), PGL(2,5), A₆, S₆) and `42 ∣ |G|` at degree 7.

**Action:** this is the *erratum-note* form. Whether to instead **fold the
correction into the paper** — rewrite the §"Decision (stage 2)" clause to
"squarefree 67/67, irreducible 66/67 (exception `x^5−x^3−2x^2−2x−1`, `Gal = D_5`)"
and replace `obs:generic` with the 2-transitivity criterion, then recompile the
pisot PDF and capture it in the next tagged release — is the maintainer's call
(the note's own N6 anticipates it). Until then the finding lives in the booklet's
`lesson6_census.py` D3 check.
