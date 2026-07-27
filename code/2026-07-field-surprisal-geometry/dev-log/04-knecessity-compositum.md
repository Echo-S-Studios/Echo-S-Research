# Session report — Harnesses 2 verification + closure of the last two open fronts

Edition: field_surprisal_geometry **v3_4** (23 pp, compiles 0 errors, 0 undefined refs).
New harnesses: `t7_knecessity.py` (28/28, exit 0), `t8_compositum.py` (13/13, exit 0).

## 1. Harnesses 2 — independent re-execution (all green)

| harness | result | note |
|---|---|---|
| t6_selection | 18/18 | exact Q(sqrt5) anchors, dichotomy margins |
| t4b_census_fast | 9/9 | census lane 2, two-lane agreement |
| t3c_partC_exact | 16/16 | explicit certificates |
| t1_windowproof | 44/44 | full count (t1_engine artifact present) |
| t5_catalog_census | 72/72 | four catalogs, branch-rank necessity |
| t4_kwindows | 13/13 lane-1 | block-3 symbolic census exceeds sandbox ceiling, exactly as shipped runlog records; covered by t4b |

Rerun total: 172 checks, every executed run exit 0. Matches shipped `runlog.txt`.

## 2. Closure A — k-necessity (OP 17.6(i) / final open (5)) [FORCED]

**Theorem (new, sec:knec).** On the catalog, for 2 <= k <= 5 and statistics
(log M, T_2..T_k) with dim V = k+1: constant sectional curvature 1/4 <=> every
window matrix L(s) has rank <= 1.

Mechanism (the new idea): **window collisions are squarefree.** 1_s - 1_s' has
entries in {-1,0,1}; the collision lattice's general element j*d1 + k*d2 carries
entry 2k at K, so cube points are only +-d1 (golden swap). The Salem square is a
multiplicity-two relation and cannot act at window level — the same mechanism
that made the heavy cofactor exactly 1. A golden-swap collision requires every
statistic to tie on {phi,tau}; then the pair are twins, ell_swap(s) = ell_s
exactly (identical rows in the same sorted slot), so the grouped coefficient is
2*M(s), forcing each; twin-containing windows have L(s) = 0 outright. The Gauss
coefficients are exactly the 2x2 minors of L(s), and all minors vanish iff
rank <= 1.

Consequences now in the paper:
- The 26-class k=3 census is the **complete** constant-1/4 double-indicator
  landscape; the catalog 3- and 4-folds need no curvature routine.
- The k=2 eight-branch tree collapses to the golden dichotomy (the branches
  lived on the expanded monomials of P, where Z^2 reintroduces multiplicities);
  branch machinery retained as an independent lane and general-catalog tool.
- General catalogs: the theorem holds for transposition-tame catalogs (all four
  studied qualify; drop-tau has no window collisions at all); genuine {-1,0,1}
  trades fall to the window-level branch-rank test.
- Honest remainder: only OP 17.6(ii), the windowwise-flat landscape beyond
  double indicators — now pure linear algebra, decoupled from curvature.

Machine lane t7 (28/28): collision census all 4 catalogs at every window size;
cube fact; twin symmetry k=2 (10+10, by determinants) and k=3 (entrywise rows +
det spot-checks); minor/rank equivalence; the k=3 windowed identity
D_k G = Z^{k-2} sum_s w_s [minor] at an exact rational m=7 instance (9 quads);
3-fold (21 windows) and 4-fold (7 windows) rank <= 1 exact + failing control.

## 3. Closure B — OP (4) compositum coupling, resolved on the first horn [FORCED]

Construction: tensor polynomial p (x) q, roots all conjugate products —
gauge-free. Every conjugate log-magnitude is a Q-vector over
(log2, log3, log5, log phi); |K_out|^2 = sqrt5*phi^2, |K_in|^2 = sqrt5/phi^2.

**Theorem (sec:compositum).** The compositum cost is not additively separable:
- contrast({sqrt2,phi} x {sqrt3,phi}) = **log 2 exactly**
- golden sector contrast({phi,phi4}^2) = **-6 log phi exactly**

So the compositum's Mahler-Gibbs law factors through the marginals for no
beta != 0: the coupling the product construction left [declared] is, in the
compositum reading, forced. Structure: rational sector {sqrt2,sqrt3,sqrt5}
exactly separable; golden tie persists (c(phi,.) = c(tau,.)); double-centered
interaction tensor Delta has **exact rank 5** over Q(L) — five independent
coupling modes. Cross-Fisher Cov_beta(cost1,cost2) != 0 at beta in {1,-1,sqrt5}
[computed]: 4.74e-3, 2.19e-4, -3.58e-4. Charge tensor law n(p(x)q) =
lcm(n_p,n_q) verified exactly on all 28 unordered pairs (every tensor conjugate
is real or purely imaginary).

Machine lane t8 (13/13): exact 7x7 cost matrix (10 on-circle products excluded
symbolically; all signs certified by rigorous intervals, min margin > 1/100),
hand anchors, both witnesses, separability, tie, Delta rank, Cov, lcm law.

## 4. Paper delta (v3_3 -> v3_4)

Inserted sec:knec (thm:goldencollide, thm:knec, cor:knecclosed, rem:squarefree,
thm:tame) and sec:compositum (thm:coupled, prop:coupstruct, prop:chargelcm,
rem:coupwarp) before the trace-form section. Rewrote rem:kobstruction (failure
blocks only the landscape) and op:kclass (part (i) closed). Open problems: (2)
reduced to the landscape; (4) and (5) marked resolved with mechanisms. Ledger:
+4 rows. Abstract: closure paragraph + "Fourteen harnesses (..., 28/28, 13/13)".
Appendix: t7/t8 paragraph + note that t4b is the census's practical lane.

## 5. Discipline

Exact arithmetic over Q(L2,L3,L5,L7,Lphi) at every decision boundary; floats
display-only (Cov values tagged [computed]); on-circle exclusions symbolic;
sign decisions by rigorous interval arithmetic; fail-first ck() harnesses;
two-lane status: hand proofs in the paper, machine lanes t7/t8.

Session totals: 172 rerun + 41 new = 213 checks executed, all exit 0.

## 6. Remaining open (honest)

- OP 17.6(ii): windowwise-flat landscape beyond double indicators (linear
  algebra; multiset plane lemma false).
- Compositum: full geometry of the coupled family; whether Delta is determined
  by the charge data.
