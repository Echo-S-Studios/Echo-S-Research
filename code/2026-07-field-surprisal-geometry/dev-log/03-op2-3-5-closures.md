# Session Report — Field Surprisal Geometry: OP-2 / OP-3 / OP-5 Closures

**Date:** 2026-07-22 · **Corpus:** `field_surprisal_geometry_v3.tex` (extended in place) · **Discipline:** exact arithmetic at every decision boundary; floats display/margin-only; fail-first `ck()` with `sys.exit(1)`; two-lane audit before any FORCED tag.

## 1. Outcomes vs. handoff goals

The handoff's success bar was "OP-1 closed plus any one of OP-2/3/4/5 advanced." Result: OP-1 was closed in the prior session; this session closed **three** further fronts and scoped the fourth honestly.

| Front | Handoff ask | Outcome | Tag |
|---|---|---|---|
| OP-2 (k≥3 classification) | Lift the branch calculus or the identity | Master identity lifted to every k; rank-≤1 window criterion; **complete k=3 double-indicator census: 26 classes**; necessity gap stated precisely | [F] census + criterion; necessity open |
| OP-3 (temperature selection) | Canonical operating point? | **Dichotomy theorem**: no metric-canonical point exists, and the candidate principles provably disagree; the [D] is the answer's final form | [F] structure + [C] margins → [D] final |
| OP-5 (catalog invariance) | Separate theorem from coincidence | **Generic census theorem** + count formula Σ(2^mult−1); four catalogs verified 8/7/6/9; new **branch-rank calculus** closes necessity mechanically | [F] |
| OP-4 (multi-field coupling) | Compositum Mahler–Gibbs vs [cmc] lcm law | Scoped out, untouched, written up as Open Problem (4) — the remaining frontier | [O] |

Also applied the two requested tex edits: the `t1_windowproof`/`t1_engine` ordering note (44/44 vs 43/43 standalone) with `run_all.sh`, and the master-identity rewording of Open Problem (2).

## 2. OP-2 — the k-window calculus and the k=3 census

**Generalized identity (Thm 17.1 / `thm:kmaster`).** For k statistics with dim V = k+1, the same three steps as at k=2 (Sylvester pivot Z, then Cauchy–Binet) give

    det[C(u_i, v_j)]_(k+1)×(k+1) = Z^k · Σ_{|s|=k+2} w_s ℓ_s(f) ℓ_s(g)
    D_k · ⟨Qf, Qg⟩_p          = Z^(k−2) · Σ_{|s|=k+2} w_s ℓ_s(f) ℓ_s(g)

with D_k = det[C(T_a,T_b)] = Z^(k−1)·det Gram_w(1,T). Verified fully symbolically at k=3, m=5, and at exact rational instances m=6,7 with both Sylvester scalings (`t4_kwindows.py` [001]–[013]).

**Criterion (Prop 17.2).** The Gauss obstructions are the 2×2 minors of the window matrices L(s) = [ℓ_s(T_a T_b)]; rank L(s) ≤ 1 on every window ⟹ constant sectional curvature 1/4. The catalog 3-fold (21 windows) and 4-fold (7 windows) pass exactly over ℚ(log2, log3, log5, logφ) — curvature-routine-free reconfirmations.

**Census (Thm 17.4).** Over pairs {S,T} with dim V = 4: **7875 pairs → 7812 nondegenerate → 312 windowwise-flat → exactly 26 distinct spans** V = ⟨1, a, 1_S′, 1_T′⟩ with S′,T′ disjoint within-level subsets (C(8,2) − 2; the two golden overlaps are redundant presentations of the all-golden class, where 1_φ + 1_τ = 1_{φτ} collapses). The eight k=2 surfaces embed.

**Two-lane methodology.** Lane 1 (`t4_kwindows.py` block 3) decides membership by symbolic rank; its census run exceeds the sandbox execution ceiling (sympy ranks, est. 15–40 min) — checks [001]–[013] stand, the census block is superseded. Lane 2 (`t4b_census_fast.py`, 9/9) independently re-derives the decision: windowwise-flat ⟺ {a·1_S, a·1_T, 1_{S∩T}} ⊂ V (symmetric window matrix, zero diagonal forces off-diagonals); membership via a block-collapsed linear system with (i) c₁-forcing from multi-level blocks, (ii) rational left-null annihilator tests on ℚ⁵ vectors, (iii) quadratic cross-ratio consistency for the unknown-c₁ case. All decisions are exact `Fraction` tests; the two lanes agree on 24 sampled pairs × 3 targets.

**Honest gap (Rem 17.5 / Open Problem 5).** The k=2 chain "windows ⟹ trichotomy ⟹ plane lemma" does not lift verbatim: the multiset plane-lemma analog is false at k≥3 (doubled support point). The census classifies the windowwise-flat landscape — the sufficient side of the conjectured equivalence. Necessity at k≥3 is the new open problem.

## 3. OP-5 — catalog invariance and the branch-rank calculus

**Generic census theorem (`thm:generic`).** ℚ-linearly independent costs ⟹ collision lattice {d : d·a = 0, Σd = 0} = 0 ⟹ q ≡ 0 forces every P-coefficient individually ⟹ every window q₄(s) = 0 ⟹ (trichotomy + plane lemma, both catalog-free) exactly m single-seed surfaces. **No branch analysis occurs for a generic catalog.**

**Count formula (`thm:countform`).** census = Σ_levels (2^mult − 1). The "eighth surface" is localized exactly: the golden coincidence's 2² − 1 = 3. Sufficiency is catalog-free geometry (q₄ = −∏Δ acquires a zero factor on every window: coincident pair or three points collinear on X = 0); nondegeneracy and a two-level control verified per catalog (control conditioned, as everywhere, on ℚ-linear independence of {1, log2, log3, log5, logφ} — plus log7 for the extended catalog; Baker).

| Catalog | m | Lattice rank | Generators | Census | Analyses |
|---|---|---|---|---|---|
| full | 7 | 2 | golden swap Δ₁, Salem square Δ₂ | **8** | 8 |
| drop-K | 6 | 1 | Δ₁ | **7** | 2 |
| drop-τ | 6 | 1 | Salem (0,0,−1,0,−1,2) | **6** | 2 |
| add-√7 | 8 | 2 | Δ₁, Δ₂ (padded) | **9** | 8 |

Generators verified in-lattice, spanning, and **ℤ-saturated** (unimodular minors), so λ-combinations enumerate the lattice completely.

**Branch-rank calculus (`thm:branchrank`).** The necessity engine for arbitrary catalogs, upgrading the k=2 singleton-survival lemma to a complete criterion:

1. Realizable primitive directions D = primitive lattice d with k, k+d both P-monomials. Full catalog: exactly {Δ₁, Δ₂, Δ₁±Δ₂, 2Δ₁±Δ₂} (six); 3Δ₁±Δ₂ machine-rejected.
2. Branches = flats of the arrangement D spans (rank 0/1/2 subspaces intersected back with D).
3. Per branch: group P's monomials by (cost class, X-class mod active relations); q ≡ 0 forces each grouped ℤ-combination of window symbols to vanish; rank_ℚ(M) = #windows ⟹ all q₄ forced. Certificate: one-sided mod-p rank (p = 2⁶¹−1), rows sorted by sparsity with early stop, exact ℚ fallback before any FAIL.
4. Coincidence direction ±(e_i − e_j) with a_i = a_j merges the seeds: windows ⊇ {i,j} auto-zero (coincident pair), symbols identify pairwise, monomials push forward by coordinate merge, rank runs against C(m−1, 4) surviving classes.

The full catalog closes in exactly 8 analyses (generic + 5 relation branches + merged-generic + merged-Salem), **mechanically reproducing the eight-branch tree of `t1_branches.py`**.

## 4. OP-3 — the selection dichotomy

**Structural side [F].** A finite-length 1-D curve has isometry group {id, flip}; the flip swaps the endpoints, which are distinguishable (β→+∞: uniform on the golden pair; β→−∞: the unique φ⁴ vertex). The labeled-endpoint isometry group is trivial ⟹ no interior point is metrically canonical.

**Numeric side [C].** Exact anchors Z(+1) = 91/30 − √5/5 and Z(−1) = 17 + 4√5 in ℚ(√5); I′ = −κ₃ symbolic. With the paper's arc convention **s(0) = 0** (see §6), the principle-selected points are pairwise distinct with margins far above error (dps 30/50 agreement < 1e−12):

| Principle | β | s(β) |
|---|---|---|
| arc-length midpoint | β_mid | **−0.3401** (new datum) |
| max Fisher (κ₃ = 0) | −0.0768 | −0.0440 |
| maximum entropy | 0 | 0 |
| heat-capacity peak | 2.5455 | +1.1401 |

Max-Fisher sits 0.2961 past the midpoint — margin ~0.30. **Conclusion (Thm 15.3):** selecting an operating temperature is irreducibly an election among inequivalent invariant principles; the corpus [D] tag is the answer's honest final form, not a gap.

## 5. Paper integration and verification state

Six-part patch applied (all anchors unique; 65,036 → 77,672 bytes): dichotomy theorem after `prop:ltot`; new §17 (k-window calculus, five results) replacing and rewriting `op:kclass`; new §18 (catalog invariance, three theorems); eight ledger rows; Open Problems — (1) resolved as dichotomy, (4) compositum coupling added, (5) k≥3 necessity added; appendix paragraph for the four new harnesses; abstract extension with the twelve-harness tally. **Compiles: 0 errors, 0 undefined references, 20 pages** (was 17). Page 14 visually verified (dichotomy + generalized identity render correctly).

Twelve-harness tally: 61/61, 17/17, 4/4, 16/16, 12/12, 24/24, 44/44, 16/16, 13/13 (t4 lane 1, in-sandbox scope stated in the appendix), 9/9, 72/72, 18/18 — **this session: 99 checks in the three completed harnesses, all exit 0.**

## 6. Findings worth flagging

- **s-convention discovery.** First t6 run failed check [009] with a uniform +0.3401 shift on all three anchors — diagnosis: the paper's arc coordinate is centered at β = 0 (uniform ensemble), not the arc midpoint. Corrected convention reproduces every published anchor to 4 dp and *strengthens* the dichotomy (midpoint ≠ β=0 as well, at s ≈ −0.3401).
- **The eighth surface is now a formula, not a story:** 2² − 1 = 3 at the doubled golden level, per the count formula — exactly the "separate theorem from coincidence" hygiene the handoff asked for.
- **Citation flag (unchanged):** the Horn–Johnson `\bibitem{hj}` remains verify-before-formal-use per corpus discipline.

## 7. Environment notes

- Per-call ceiling ~115 s; `setsid` background survival unreliable — the t4 lane-1 census froze at check [013]/1361 bytes across both launches. The fast lane (t4b) was designed to carry the census within a single foreground call (~3 s).
- `/bin/sh` (no `time`, no `PIPESTATUS`); LF for `.py`, CRLF for `.md`; `runlog.txt` append-only.

## 8. Remaining work (next session)

1. **OP-4** (only untouched front): compositum Mahler–Gibbs on ℚ(θ₁,θ₂)'s emission catalog vs the [cmc] ℤ/lcm charge law — forced coupling ([D]→[F]) or a no-coupling theorem.
2. **Open Problem (5):** k≥3 necessity (constant 1/4 ⟹ windowwise rank ≤ 1) — the collision/branch analog one level up.
3. Optional: rerun t4 block 3 to completion on an unconstrained host for a full lane-1/lane-2 census cross-check (expected: identical 312/26).
