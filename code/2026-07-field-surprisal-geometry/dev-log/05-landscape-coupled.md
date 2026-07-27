# SESSION REPORT — The Landscape Completed; The Coupled Family Curves (v3_4 → v3_5)

Scope: close the last two enumerated open items of the field-surprisal-geometry corpus:
(A) Open Problem 17.6(ii) — the k ≥ 3 landscape of windowwise-flat tuples; (B) the coupled
compositum family's full geometry and the Δ-vs-charge question. Both closed. The paper's
enumerated open-problem ledger is now fully resolved.

## (A) The landscape: partitioned-affine classification — conjecture REFUTED

Chain (all [FORCED], hand + machine):

1. Reformulation. For a spanning (k+2)-window, ℓ_s is the unique affine dependency λ of the
   points P_i = (T_1..T_k)(i), supported on the window's unique circuit, and
   L(s) = Σ λ_i P_i P_iᵀ; rank is affine-invariant.
2. Moment trichotomy (thm:momtri). Circuit size 2 ⟹ L = 0; size 3 ⟹ rank exactly 1
   (moment −∏(t_i−t_j) ≠ 0); size ≥ 4 ⟹ rank ≥ 2, by normalized minors: mixed minor
   x_j x_l x_r (d ≥ 3), det = xy(1−x−y) (d = 2) — products of circuit weights.
   Hence constant-¼ ⟺ every circuit has size ≤ 3.
3. Structure (thm:pastruct). Flat ⟹ rich lines pairwise disjoint + skew; no transversal
   triples; V = V(π) = {f affine-in-a per block} via the deficient-transversal circuit
   argument. Conversely dim-matching V(π) is flat (functional-evaluation proof).
4. Canonical form. Blocks with ≤ 2 distinct cost values ≡ within-level clusters; line blocks
   need ≥ 3 distinct cost values; dim V(π) = c + 2b₂.
5. Classification + counts (thm:paclass): constant-¼ = canonical V(π), c + 2b₂ = k+1.
   Catalog counts k = 2..6: **8, 56, 95, 31, 1**.
   - k=2 recovers the eight surfaces in one line.
   - k=3: 26 double-indicator + 30 split-affine ⟨1, a, 1_B, a·1_B⟩; the five (3,4)-splits
     with 3-side {φ, 1/φ, x} coincide with ⟨1, a, 1_golden, 1_x⟩ (no double count).
   - k=5 contains the golden merge ℝ[log M] (the totally geodesic family) among 31.
6. Corollary (cor:landres): the within-level indicator-join conjecture of op:kclass is
   **false** — split-affine classes contain no indicators beyond 1_B, 1_{B^c}. Thm 17.4
   stands exactly as the double-indicator stratum.
7. Scope (rem:fatlevels): verbatim for any catalog with level multiplicities ≤ 2 (drop-K,
   drop-1/φ, add-√7); multiplicity ≥ 3 admits vertical line blocks (outlook).

Session correction trail: the "exotic T5" incidence configuration explored mid-derivation is
DEGENERATE (its constraints force coplanarity) — consistent with, and superseded by, the skew-
line lemma; and the earlier within-session count 61 was corrected to 56 (five golden overlaps).
Nothing wrong was ever shipped: v3_4 stated the landscape as open with a conjecture; v3_5
resolves and refutes.

## (B) The coupled family curves; Δ is not charge-determined

- dim V = 4; a genuine 5-circuit at cells (φ,√5),(√3,φ),(φ⁴,√3),(√2,φ),(√2,φ⁴): all five
  dependency weights nonzero, interval-certified at the true logs ⟹ not partitioned-affine;
  the window matrix has rank ≥ 2 ⟹ not windowwise flat. [FORCED]
- Gauss obstructions at the uniform point nonzero (rigorous 60-digit enclosures) ⟹ NOT
  constant ¼. Sectional values at uniform:
    sec(cost₁,cost₂) = ¼ − 0.248531… ≈ 0.0015   (ghost of the product geometry)
    sec(cost_i, c)   = ¼ − 0.083802… ≈ 0.1662
  [FORCED for ≠ ¼; values displayed with certified enclosures]
- Cov_uniform(cost₁, cost₂) = 0 exactly, yet curvature already deviates: the coupling is
  invisible to the cross-covariance at uniform and visible to curvature. [FORCED]
- Δ is NOT determined by the charge data: (√2,√3) and (√2,√5) share charge pair (ℤ/2,ℤ/2) but
  Δ = (5/49, −5/49, 0, 20/49) ≠ (5/49, 2/49, −1/7, 20/49) over (log2, log3, log5, logφ),
  exactly. [FORCED]

## Machine lanes (both exit 0)

| harness            | checks | content                                                        |
|--------------------|--------|----------------------------------------------------------------|
| t9_landscape.py    | 21/21  | trichotomy symbolic d=2,3,4 + exact instances; exhaustive       |
|                    |        | window verification of all 8+56+95+31 families over ℚ(L);      |
|                    |        | 5 overlap equalities (symbolic rank 4); 1540 distinctness       |
|                    |        | certificates at the true logs (interval rank-5 minors);         |
|                    |        | 20 exact random flatness-violating controls                     |
| t10_coupled.py     | 9/9    | dim 4; 5-circuit + rank≥2 window (interval-certified);          |
|                    |        | 3/3 Gauss obstructions ≠ 0 at uniform; sectional enclosures;    |
|                    |        | exact Cov = 0; Δ-charge witness                                 |

Discipline: flatness checks are polynomial identities over ℚ(L2,L3,L5,Lφ) (unconditional);
all nonvanishing claims at the true catalog certified by mpmath.iv interval arithmetic
(dps 60); rational specializations used only to locate witnesses, never as certificates.
Session totals: 172 (v3_4 suite re-run) + 28 (t7) + 13 (t8) + 21 (t9) + 9 (t10) = 243, exit 0.

## Paper delta (v3_4 → v3_5, 12 anchored edits, 26 pp, 0 errors, 0 undefined refs)

- NEW §21 "The landscape completed: the partitioned-affine classification"
  (thm:momtri, thm:pastruct, thm:paclass, cor:landres, rem:planebypass, rem:fatlevels).
- NEW §22 "The coupled family curves"
  (thm:coupcurv, prop:covinvis, prop:deltacharge, rem:closed).
- op:kclass → resolved in both parts (conjecture refuted); rem:kobstruction → both halves
  closed, plane lemma bypassed; §sec:knec transition updated; rem:coupwarp tail → residuals
  settled; Open problems (2) resolved with counts, (4) residuals closed; Outlook paragraph
  added; abstract extended ("Sixteen harnesses … 21/21, 9/9"); +4 ledger rows; appendix
  paragraph for t9/t10.

## Outlook (honest, outward-facing)

Fat-level catalogs (vertical line blocks); composita of two distinct catalogs; the coupled
family's curvature landscape away from the uniform point; the standing ℚ-linear-independence
conditioning of the collision analysis. The temperature selection remains [DECLARED]-final by
the selection dichotomy — an answer, not a gap.
