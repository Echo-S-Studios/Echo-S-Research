# HANDOFF — Exact Arithmetic booklet → Echo-S repos

**To:** a local Claude session with write access to `math-research-pipelines` (private) and
`Echo-S-Studios/Echo-S-Research` (public).
**From:** a web session that produced the thirteen-lesson booklet and four revision passes.
**Date of handoff:** 2026-07-27.

You are inheriting a finished artifact and an unfinished queue. Most of the queue consists of
items that were **hard-blocked for me and are trivial for you**: the web session could not
`grep` the archive repo (project knowledge indexes only the pipelines repo) and could not fetch
`github.com/.../tree/...` (robots). You have local write access. Three of the four remaining
open questions are one shell command each.

---

## 0. Set your paths first

```bash
export PIPELINES=~/src/math-research-pipelines
export ARCHIVE=~/src/Echo-S-Research
export HANDOFF=~/handoff            # where this bundle was unpacked
```

Everything below uses those. Adjust to your actual checkout locations before running anything.

---

## 1. Files handed to you

| Path | What it is | Use it for |
|---|---|---|
| `$HANDOFF/exact-arithmetic-thirteen-lessons.html` | the booklet, 417 KB, self-contained, 15 pages | the artifact to file (Phase 2) |
| `$HANDOFF/checks/*.py` | **14 verification scripts + `run_all.py`** | re-run them; they become the contribution's `code/` tree |
| `$HANDOFF/checks/run_all.py` | reproduction harness: env check, digests, quick/full modes | the one command to run first |
| `$HANDOFF/checks/REPRODUCE.md` | environment pins, expected runtimes, what `--quick` does not claim | ships with `code/` |
| `$HANDOFF/checks/SHA256SUMS` | digests of all 15 | must match the booklet's audit manifest |
| `$HANDOFF/reviews/Remaining_Revisions.txt` | third-pass external audit | context for why the prose reads as it does |
| `$HANDOFF/reviews/Final_Revisions.txt` | fourth-pass external audit | ditto; also names the open frontier |
| `$HANDOFF/lessons/lesson-*.md` | the ten original lesson transcripts | source text if you need to re-derive prose |

**First thing to do, before reading further:**

```bash
cd $HANDOFF/checks
python3 run_all.py --list      # environment, runtimes, what each script decides
python3 run_all.py --full      # everything, ~3.5 min, exits non-zero on any failure
```

Expected: `checksums: all match`, then `ALL GREEN`. Counts are 16, 32, 33, 25, 23, 20, 6, 39, 20,
32, 26, 13, 21, and 22 executable certificates + 1 cited theorem. **Do not trust this document's
summaries of what the scripts prove — trust the scripts.** That is the house rule, and it applies
to me.

Environment is pinned in `REPRODUCE.md`: CPython 3.12.3, SymPy 1.14.0, mpmath 1.3.0, with
3.12–3.13 / 1.13–1.14 known-compatible. The SymPy floor matters — `Poly.intervals(all=True)`
changed shape and silently broke the census's disk certificate; the shipped code now tolerates
both. `python3 run_all.py` without `--full` is a ~50 s smoke run that deliberately does **not**
assert the published census tallies.

---

**Three repairs were made to this bundle after the final audit, and they matter to you:**

1. `lesson6_checks.py` carried a **dead census section** that crashed on load under current
   SymPy. Excised; the file is checks A–C, the 20 its badge names.
2. `revision_checks.py` **hard-coded `/home/claude`**. Now resolves siblings via `__file__`.
3. `lesson6_census.py` check **D3 asserted the source note's claim** (`Rat°` irreducible on all
   67) and therefore *failed* when run, because the truth is 66. **That failure was the D₅
   finding.** Reporting the artifact as 6/6 had hidden a finding inside a passing badge. D3 now
   asserts the corrected counts and names the exception, so the finding lives in the check.

That third one is the template for how the D₅ erratum should read when you file it: the check
encodes the corrected state, not the state being corrected.

## 2. House discipline — non-negotiable

These are repo law, not preferences. Several are already enforced by
`.github/pull_request_template.md` and CI.

1. **No uncertified approximation crosses a decision boundary.** Not "no numbers" — the
   corrected form. Symbolic equality where possible; exhaustive enumeration where the set is
   finite; **directed interval** exclusion where a transcendental comparison enters. Floats may
   display, never decide.
2. **Every document patch is anchored and all-or-nothing.** Anchor on markup or distinctive
   prose, never line numbers. Validate *every* anchor before writing *any* of them. This is the
   `docs/HTML_EDITS.md` idiom; it caught four real errors during these revisions.
3. **Function → mechanism → implementation must match.** A check's label may not claim more
   than its predicate decides. The fourth audit found several that did; they are fixed in
   `certificates.py` and you should hold new checks to the same bar.
4. **A summary may be shorter than its source but never stronger.** Titles, figures, headings,
   deks, tags, handoffs, and ledger rows are all copies of a claim. When one changes, the others
   go stale. This is the single most common defect found across four audits.
5. **Edit `.tex`, never PDFs.** Regenerate `SHA256SUMS` and manifests after any artifact change.
6. **Determinism, for the archive's drift gate:** seeded RNG, rounded floats, sorted writes.
7. **Tags are pinned to their checking file.** A tag with no check ID behind it is decoration.

---

## 3. Phase 0 — settle the four open questions (do this first; ~15 minutes)

Each of these has been carried as "pending" across multiple sessions purely because the web
session lacked repo access. Resolve them, then update the ledgers.

### 0.1 — D1: the cross-engine signature count

```bash
grep -rn "twenty-three" $ARCHIVE/papers/ $PIPELINES/papers/
grep -rn "twenty-six\|twenty-four" $ARCHIVE/papers/*relational*
```

**Expected:** empty for the first. The served PDF reads *twenty-six items* / *twenty-four contact
signatures*, and the numeral 237 belongs to the Pisot note §7.3 run 3 (237/237 signatures, 474
dual-path executions) — a different ledger entirely.
**Done when:** either D1 is struck from the blocker list as never-shipped, or you have the line
number and can patch it.

### 0.2 — E2 and E3: numerals absent from every indexed source

```bash
grep -rn "163\|160" $ARCHIVE/papers/*.tex | grep -i "straddle"
grep -rn -A3 -B3 "degree-7\|degree 7" $ARCHIVE/papers/*pisot* | grep -i "c-vector\|cvector"
```

- **E2:** a displayed degree-7 Pisot polynomial contradicts its own c-vector. The **c-vector is
  machine-confirmed correct**, so the displayed polynomial is the thing to fix.
- **E3:** the both-straddle count should read **163, not 160** — the difference is 3
  twist-fixed single-member classes.

**Done when:** located and patched in `.tex`, or confirmed absent (in which case they were
already integrated pre-tag and the blocker list is stale).

### 0.3 — E1: confirm the deposit agrees with the served compendium

E1 is **closed on the pipelines side** — the compendium's subfield census table reads
`(2,2,0) → 7`, which is correct. The seven quadratic subfields of
ℚ(√2, √3, 5^(1/4)) are ℚ(√d) for d ∈ {2,3,5,6,10,15,30}, and the one an earlier count omitted
is ℚ(√5), presented as `x² − 3x + 1` (roots φ^±2, discriminant 5). Independently reproduced by
group theory in `lesson8_checks.py` (B1–B3) and `lesson9_checks.py` (C1–C3).

**Remaining:** diff the Zenodo v1.0.0 deposit's TeX against the current archive TeX to confirm
the sealed deposit also reads 7. If it reads 6, that is a versioned erratum, not a live defect.

### 0.4 — G11: does any module guard total reality?

```bash
grep -rn "count_roots\|totally.real\|is_totally_real" $PIPELINES/*/L00M/ $PIPELINES/*/src/
```

**Expected:** no hit. Across the modules I read (`integral_basis.py`, `residual_learner.py`,
`field_extension.py`, `compositum.py`, `capacity.py`) the guards are G8 (float refusal) and G10
(monic-integer) only. If that holds, see Phase 3.1.

---

## 4. Phase 1 — file the errata

### 1.1 — The D₅ erratum + the orbit criterion (highest value)

**Target:** the Pisot cross-shell residue note, §6 stage-2 sentence and Observation 6.2.

**The defect.** The note states Rat°_p is *"in fact irreducible of degree 20 on all 67"*. Truth:
**squarefree on all 67, irreducible on 66 of 67.** The exception is

```
p = x^5 - x^3 - 2x^2 - 2x - 1,     Gal(p) = D_5
```

whose Rat° splits into two irreducible self-reciprocal degree-10 factors. Certified by exact
re-multiplication, so the factorization *is* the proof — engine-independent.

**The theorem survives.** Both factors carry 2 unimodular roots, so S\* is still the full degree-20
product, deg C₂ = 400, and the complete scan returns {Φ₁²⁰}. The instance is inert; Theorem 6.1's
conclusion is untouched.

**Why it escaped.** The observables — `deg C2 == 400` and `contacts == {Φ₁^20}` — are byte-identical
whether Rat° is irreducible or splits with both factors admissible. The missing line is one
`assert Rat0.is_irreducible` per instance.

**The upgrade (this is the part worth publishing).** Observation 6.2's `[PLAUSIBLE]` P7 should be
replaced by an iff-criterion:

> **Lemma.** If Rat°_p is squarefree, then Rat°_p is irreducible over ℚ **iff** Gal(p) acts
> 2-transitively on the roots of p.

Irreducibility is transitivity on the n(n−1) ordered distinct root pairs, which *is* 2-transitivity.
Galois census of the 67: **S₅ × 65, M₂₀ × 1, D₅ × 1**. M₂₀ = AGL(1,5) has order 20 and is sharply
2-transitive → irreducible. D₅ has order 10 < 20 → two orbits of 10 → exactly the observed split.
Reproduces 66/67 with no residue.

**Consequence for N4 (degree 6–7), worth stating in the same note.** Irreducibility requires
n(n−1) | |G|, so **30 | |G| at degree 6** and **42 | |G| at degree 7**. At degree 6 only the four
2-transitive groups — PSL(2,5), PGL(2,5), A₆, S₆ (orders 60, 120, 360, 720) — can give an
irreducible Rat°; every other transitive group *forces* a split S\* and a cheaper C₂. The next
census can pre-sort by Galois group instead of stumbling on the exceptions.

**Evidence:** `checks/lesson6_checks.py`, `checks/lesson6_census.py`, and
`checks/lesson9_checks.py` C5–C8.

**Filing.** Ace's call between an erratum note and a v1.1 fold-in (the note's own N6 anticipates
this). Draft manifest, to be conformed to the real schema — read `$ARCHIVE/CONTRIBUTING.md` and
`$ARCHIVE/scripts/new_contribution.py` first, and prefer scaffolding over hand-writing:

```yaml
title: "Erratum and sharpening: Rat° irreducibility in the quintic census"
targets: [papers/<pisot-folder>/pisot_residue.tex]     # §6 stage-2 sentence; Obs 6.2
kind: erratum+sharpening
claims:
  - id: E-D5-1
    text: "Rat° irreducible on 66/67; the D5 instance splits 10+10."
    tag: FORCED          # exhibited factorization, certified by re-multiplication
  - id: E-D5-2
    text: "Rat°_p irreducible <=> Gal(p) 2-transitive (given Rat° squarefree)."
    tag: FORCED          # orbit-stabilizer; supersedes Obs 6.2's PLAUSIBLE P7
  - id: E-D5-3
    text: "Inertness certificate for the D5 instance is clean: deg C2 = 400, {Phi_1^20}."
    tag: FORCED
checks:
  - factor_list(Rat0) yields two degree-10 irreducible self-reciprocal factors
  - exact re-multiplication of the factors equals Rat0        # engine-independent
  - S* = Rat0 (per-factor unimodular count = 2 each, trace-fold Sturm on (-2,2))
  - complete cyclotomic scan of C2 (deg 400, 790 candidates, max m = 1680) == {Phi_1^20}
  - orbit of an ordered pair: |S5| 20, |M20| 20, |D5| 10  ->  66/67 reproduced
```

### 1.2 — F1: a display slip in the ATR paper

Three occurrences of **0.046** should read **0.045**. The exact tube quantity at k = 10, n = 3 is
1/(2 + k(n−1)) = **1/22 = 0.045454…**. Mechanism identified: `0.167 × 6/22 = 0.045545 → 0.046`,
i.e. the value was scaled from the already-rounded first figure rather than computed from 1/22.

```bash
grep -rn "0\.046" $ARCHIVE/papers/*transport* $ARCHIVE/papers/*accumulated*
```

Expected in the §2.1 text, Table 1, and §14 (C4). Display-only: no decision depends on it, and the
paper's own monotonicity claims are exact and unaffected (∂ₖ[1/(2+k(n−1))] = −(n−1)/(2+k(n−1))² < 0
for all n ≥ 2). Evidence: `checks/lesson9_checks.py` B3, B4, B7.

### 1.3 — Schinzel: a citation the corpus is missing, and a tag it earns

Zero mentions across the corpus. Verified against the literature this session:

> **A. Schinzel**, *On the product of the conjugates outside the unit circle of an algebraic
> number*, **Acta Arith. 24 (1973), 385–399**; addendum ibid. **26 (1974/75), 329–331**.
> For totally real α ≠ 0, ±1 the Weil height satisfies h(α) ≥ h(φ) = ½ log φ. Since
> h = (1/d)·log M, that is exactly **M ≥ φ^(d/2)** — tight at x² − x − 1, and growing with degree.
> **No unit hypothesis is required.**

**What it buys.** The universal (1, φ) gap, currently tagged *forced for quadratics, plausible in
general*, is **forced unconditionally on the totally real sector** — and there the bound is not a
flat φ but φ^(d/2).

**What it does not buy, and this is the part to state carefully.** Wall 2 permits arguments in
(π/2)ℤ — real **or purely imaginary** — so the emission image is *not* totally real, and the
floor-attaining witness lies outside Schinzel's hypothesis: q₂ = x⁴ + x² − 1 carries its whole
measure on ±i√φ. **The classical theorem and the corpus's contribution partition the problem
rather than duplicate it.** Say that explicitly wherever the tag is promoted, or the next reader
will conclude the whole wall is a 53-year-old theorem re-derived.

**A route left open.** ψ² sends real-or-imaginary to real, so it maps the whole image into
Schinzel's sector. Verified as a positive control on q₂ (M(q₂)² = φ² recovers M(q₂) ≥ φ exactly);
`[OPEN]` as a general argument pending degree bookkeeping.

**Next names to read against,** since this puts the corpus in their territory: **Smyth (1980)** and
**Flammang** on the totally real spectrum.

Evidence: `checks/synthesis_checks.py` S1–S8, `checks/certificates.py` C21.

---

## 5. Phase 2 — file the booklet as a contribution

This is the item that closes the custody graph, and the genealogy argues it is the natural next
row of the same invariant: *stage 0 was a single self-contained HTML with pinned provenance; so is
this.* By the corpus's own definition — **a paper is the prose face of a green test suite** — the
booklet already qualifies. It just hasn't been filed.

```bash
cd $ARCHIVE && python3 scripts/new_contribution.py   # scaffold; do not hand-write the manifest
# proposed folder: papers/2026-07-exact-arithmetic/
```

**Mapping onto the contract:**

| Contract slot | Content |
|---|---|
| paper artifact | `exact-arithmetic-thirteen-lessons.html` (`kind: doc`, self-contained, 429 KB) |
| `code/` | the 15 files from `$HANDOFF/checks/`, `REPRODUCE.md` included |
| `tests/` | `run_all.py` is the entry point and already exits non-zero on any failure |
| `data/` | none — every script is self-generating |
| `figures/` | none — the booklet's figures are inline CSS/SVG |
| `checks:` | `python3 code/run_all.py --full` — one entry, deterministic, ~3.5 min |

**Two things the manifest must say honestly, or don't file it:**

1. **The Lesson 11–13 authoring scripts are unrecoverable** (32 + 26 + 30 = 88 checks, container
   gone). `certificates.py` is a **selected headline recheck** — 22 checks, stronger per check and
   weaker in coverage. It does **not** supersede them. The booklet's ledger already says this;
   the manifest must too.
2. **Check counts:** 354 numbered checks across the thirteen lessons, plus 6 (census) + 13
   (revision pass) + 21 (synthesis) + 22 (certificates) as separate artifacts. Do not merge them
   into one headline number.

**Before filing, add a card.** `$PIPELINES/papers/catalog.json` drives the pipelines site;
`kind: doc` routes it to *Additional papers*. The archive gets a card per folder automatically.

---

## 6. Phase 3 — hardening and doc-sync

### 3.1 — G11: the learner's capture predicate

`Prop exactcapture`'s proof states the hypothesis outright — *"positive-definiteness of G on a
totally real field gives ‖r‖²_G = 0 ⟺ r = 0"* — and the paper prescribes the Hermitian form
G₂ = M\*M when complex places appear. The constructor does not check it.

**Exact witness**, on the paper's own example field:

```
K = Q(cbrt 2), forced basis {1}.  B^T G B = 3 != 0, so no NO_PROJECTION sentinel fires.
Observe theta = cbrt 2:  r = (0,1,0) != 0  but  ||r||^2_G = Tr(theta^2) = 0 exactly.
=> the shipped predicate `rn == 0` reports CAPTURED for the generator of its own field.
Second mode: over Q(i) with basis {1}, the residual of 1+i has norm -2 — a NEGATIVE gain,
which the capacity gate compares against a floor.
```

**Fix:** one exact Sturm count at construction — `count_roots(-oo, oo) == degree`.

**Scope it correctly in whatever you write.** Total reality is a **conservative sufficient**
guard, *not* a characterization of every rational failure mode. r₂ > 0 gives *real*
indefiniteness, which yields negative gains; a *rational* isotropic residual is a separate
arithmetic question (2a² − 4b² is indefinite over ℝ with no nonzero rational zero). Both of this
corpus's witnesses happen to supply one, so the conclusion stands — but the argument
"indefinite, therefore isotropic" does not.

Evidence: `checks/lesson10_checks.py` A8–A10, `checks/revision_checks.py` R2–R2d.

### 3.2 — ℚ-independence is discharged; the README hasn't noticed

The field-surprisal repo README still carries ℚ-independence of (log 2, log 3, log 5, log φ) as a
standing open item. It is **proved** — norms in ℚ(√5): N(2), N(3), N(5) = 4, 9, 25 against
N(φ) = −1, so 2^a 3^b 5^c φ^d = 1 forces (−1)^d = 1 and 2^2a 3^2b 5^2c = 1, hence a = b = c = 0 by
unique factorization, then d = 0. Move it from `[declared]`/standing-input to `[forced]` in the
README and in any outlook section that still conditions on it.

### 3.3 — Filter-side guard *(reported, NOT verified by me — verify before acting)*

An external audit reports that a **block-3 `a_in_span` guard omits the four-block affine relation
v₁₁ = v₁₀ + v₀₁ − v₀₀**. If real, it would *silently shrink a census* rather than leave a visible
gap. It shares its trigger — multiplicity ≥ 3 — with the fat-level frontier that Lesson 12's Q3
names on the classification side. I could not check it; you can.

### 3.4 — Doc-sync

- `$PIPELINES/papers/README.md` lists a smaller companion stack than the archive's current **18**
  papers. The archive has outrun the pipelines repo's self-description.
- `A3_DESIGN.md` still lists **F4 as OPEN** though the λ = 2c revision closed it.
- Release hygiene from an older review — **verify before trusting, some may be stale**: filenames
  carrying spaces or `.l` / `.v13` suffixes; `papers/COVERAGE.md` containing local Windows paths;
  a stale Lehmer's Box PDF (drop-zone card 14 pp vs archive 15 pp); MIT vs CC-BY-4.0 licence
  mismatch; no cross-links between the two repos; four dual-homed `.tex` sources with no declared
  canonical direction.

---

## 7. Phase 4 — the standing proposal

The corpus has **drift CI for outputs** and **nothing for prose**. Four consecutive audits found
the same failure shape: *a correction landed in one representation and failed to propagate into
another* — the title said ten lessons, a heading said "two open problems closed" while its own
paragraph said otherwise, a manifest described 16-character hashes while printing 64.

The proposal is a **claim registry**: canonical statement, scope, evidence tag, certificate ID,
and *approved summary wording* — with a CI gate that fails when a claim's short forms drift from
its long form. `contribution.yml` already does this for artifacts. This is the same contract for
sentences.

The house rule that falls out, and which is now in the booklet's method paragraph:
**a summary may be shorter than its source but never stronger.**

---

## 8. Do not

- **Do not** re-run the four revision passes against the booklet. It is at its final state; the
  patch scripts were one-shot and are not in this bundle deliberately.
- **Do not** promote a tag without a check ID behind it, or a published theorem covering exactly
  the stated sector.
- **Do not** let a finite search promote a universal claim. It may only demote one. (A *published
  unconditional theorem* may promote one — that asymmetry is the Schinzel case, and it is the only
  exception.)
- **Do not** describe `certificates.py` as replacing the lost 88 checks anywhere.
- **Do not** edit PDFs. Sources are canonical.

---

## 9. Provenance: what I verified versus what I am relaying

| Status | Items |
|---|---|
| **Verified this session, scripts included** | D₅ erratum and the orbit criterion; F1's exact value and mechanism; the Schinzel translation and its scope; G11's witness; E1's census by two disjoint methods; ℚ-independence; every mathematical claim in the booklet |
| **Verified against the literature** | The Schinzel citation, including the addendum |
| **Evidence gathered, conclusion pending your `grep`** | D1 (served PDF is clean; TeX unchecked) |
| **Unlocated, numerals absent from every source I could index** | E2, E3 |
| **Relayed from an external audit, not verified by me** | the filter-side `a_in_span` guard |
| **From older sessions, may be stale** | the release-hygiene list in §3.4 |

One process note worth carrying forward, because it recurred four times across these revisions and
is the reason the discipline exists: **when the engine and the hand disagreed, the hand was wrong
every time.** A hand-transcribed 2020 against a computed 1940; a `·` typed where the file held
`&#183;`; an anchor assuming a line break that wasn't there; a dropped `</span>` caught at 440/441
by a tag-balance count. Each was found by a check that had no reason to exist except that the
protocol required it.
