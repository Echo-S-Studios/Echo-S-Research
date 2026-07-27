# MANIFEST — the verification suite

Twenty harnesses, **422 checks**, all exact-arithmetic and fail-first (each exits
non-zero on the first failed assertion). Exact decisions run over
`ℚ(log2, log3, log5, logφ)` (plus `log7` for the extended catalog); floats are
display-only; nonvanishing at the true catalog is certified by `mpmath.iv`
interval arithmetic. Run everything with `./run_all.sh`.

## Foundational (the five forced layers)

| harness | checks | proves | paper |
|---|---|---|---|
| `field_surprisal_v2.py` | 12 | Z=17+4√5; surprisal affinity; Fisher = Var(logM); dual-flat Bregman–KL; C=β²I; Bhattacharyya distances; embedding obstruction k₁k₂=−¼; Fisher monotonicity | §§2–9, 19 |
| `suspension_theorem.py` | 9 | the indicator-suspension theorem: any single-outcome-indicator surface has **constant curvature ¼** (round-sphere metric, base-measure invariant); generic controls give 4/21, 857/16928 | §10 |
| `field_surprisal_tier2.py` | 15 | multi-statistic landscape (only degree is constant-¼ among forced invariants); product geometry (block-diagonal Fisher, zero cross-curvature); charge grading (parity = 𝟙_{φ⁴} a second constant-¼ surface) | §§11, 12, 14 |
| `field_surprisal_classification.py` | 9 | **no (logM,X) surface is totally geodesic** (value-count d∈{6,7}); the ruled 2nd fundamental form (II_ac=II_cc=0, II_aa≠0); the **8** constant-¼ indicator surfaces among all 127; ~5000-sample non-indicator exclusion (later superseded by proof) | §§10, 12 |

## The k=2 classification (proven, then made conceptual)

| harness | checks | proves | paper |
|---|---|---|---|
| `t1_core.py` | 61 | cost-value structure; the **ruled system** {aX∈V, X²∈V} ⇔ 8 families; invariant-module dims N(X)=2 on families / 1 on controls; Möbius-elimination steps | §13 (thm:ruled) |
| `t1_reduction.py` | 17 | the reduction `Z²Dq = P(w)`; `q = Tr(ΓN)`; four-point rank-one law `q=κq₄`; the **rank-2 collision lattice** (golden swap + Salem square) | §13 |
| `t1_engine.py` | 4 | the full **462-monomial** coefficient dictionary of P; support-≤3 vanishing; the (3,1,1,1)=q₄ window identity at all 140 heavy placements. **Writes `P_coeffs.pkl`.** | §13 (prop) |
| `t1_branches.py` | 16 | branch survival across all 8 branches; merged six-point recursion; the 8 families kill all 35 windows identically (retained as an independent lane) | §13 |
| `t1_windowproof.py` | 44 | the **master window identity** `P = Z²·Σ_{\|s\|=4} q₄(s) w_s` by **Sylvester + Cauchy–Binet** — the conceptual proof; coefficient law + 462-census; `q₄=−∏Δ`; compound-ratio form. **Reads `P_coeffs.pkl`** (43/43 standalone if absent). | §13 (thm:master) |
| `t3c_partC_exact.py` | 16 | exact certificates (no rank heuristics) for the geodesic threshold; the level-set annihilator ∏(a−vᵣ𝟙)=0; a²∉V for the catalog 3-fold | §16 |

## Temperature, higher-k, census, selection

| harness | checks | proves | paper |
|---|---|---|---|
| `t2_temperature.py` | 12 | forced Gibbs form (measurable multiplicative Cauchy); the lattice-Cauchy pathology; exact anchors Z(∓1); I′=−κ₃; finite information length | §15 |
| `t3_suspension.py` | 24 | the **iterated join** — round ¼-spheres in every dimension (k=3,4 exact); curvature routine validated on sphere/hyperbolic/flat/simplex; catalog 3-fold; geodesic threshold k=5 | §16 |
| `t4_kwindows.py` | 13† | the **generalized window identity** (Sylvester + Cauchy–Binet over (k+2)-subsets); windowed Gauss obstructions; rank-≤1 ⇒ constant-¼; k=3 census, lane 1 | §17 |
| `t4b_census_fast.py` | 9 | the k=3 double-indicator census, lane 2 (exact ℚ⁵-annihilator): 7812 dim-4 pairs → 312 flat → **26 classes**; two-lane agreement | §17 (thm:kcensus) |
| `t5_catalog_census.py` | 72 | the **generic census** theorem; the count `Σ_ℓ(2^{m_ℓ}−1)`; four catalogs (full/drop-K/drop-τ/add-√7 → **8/7/6/9**); the branch-rank necessity closed per branch | §18 |
| `t6_selection.py` | 18 | the **temperature dichotomy**: no metric-canonical operating point; the four candidate principles select pairwise-distinct arc positions | §15 (thm:dichotomy) |

## The last two fronts

| harness | checks | proves | paper |
|---|---|---|---|
| `t7_knecessity.py` | 28 | **k-necessity** (2≤k≤5): constant-¼ ⇔ every window matrix rank ≤1, via the **squarefree-collision** mechanism (the Salem square can't act at window level) + twin symmetry; all four catalogs | §sec:knec (thm:knec) |
| `t8_compositum.py` | 13 | **OP-4 first horn**: the compositum cost is non-separable — interaction quanta **log2** and **−6logφ**; rational sector separable; Δ rank 5; lcm charge law on all 28 pairs | §sec:compositum |
| `t9_landscape.py` | 21 | **OP-17.6(ii)**: the partitioned-affine classification; moment trichotomy; exhaustive window check of all **8/56/95/31** families; the **conjecture is refuted** by 30 split-affine classes | §21 |
| `t10_coupled.py` | 9 | the coupled compositum family **curves** (nonzero Gauss obstruction at uniform; genuine 5-circuit; rank-≥2 window); Cov_uniform=0 yet curvature deviates; **Δ not charge-determined** | §22 |

† `t4_kwindows.py`: the 13 counted checks (lane 1) run quickly; its block-3 *symbolic*
census over 7812 pairs exceeds a small sandbox's time ceiling. That census is the
same result `t4b_census_fast.py` proves via the fast exact lane, with sampled
two-lane agreement — so `t4b` is the practical census lane.

## Independent-audit status

Every load-bearing claim above was independently re-derived and confirmed during
review (not merely re-run): the master identity `P=Z²Σq₄(s)w_s` (symbolic + random
rationals), the N(X) module dims, the four-point trichotomy, the squarefree
collision lattice ({0,±golden swap}), the twin symmetry, the compositum contrasts
(log2, −6logφ), the census counts (8/56/95/31 — including an independent partition
enumeration reproducing **56** at k=3), the dichotomy arc positions, and the
split-affine refutation (a constant-¼ family containing no indicators beyond
1_B/1_{B^c}). See `dev-log/` for the per-session record.
