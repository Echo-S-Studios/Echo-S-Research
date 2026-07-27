# Field Surprisal Geometry -- Session Report: OP-1 closure (window identity)
Date: 2026-07-22 | Echo-Squirrel Research | exact-arithmetic discipline throughout

## Outcome vs handoff success bar
- **OP-1 (flagship): CLOSED [F].** The window identity is now the master theorem
  P = Z^2 SUM_{|s|=4} q4(s) w_s, proved by Sylvester's determinant identity
  (pivot Z on the mixed 4x4 Gram) plus Cauchy-Binet over 4-subsets -- exactly
  the route named in the handoff. Exceeds the success bar: not just the heavy
  cofactor-1 case but the full coefficient law, with new corollaries.
- Baseline reproduction: complete (see harness table); t3 Part C closed by a
  supplementary exact harness after hitting the sandbox execution ceiling.

## The theorem and its corollaries (paper: Thm 13.4, Cor 13.5, Rem 13.6)
1. (S1) R(f,g) = det3 of the C-covariance matrix, rows (f,a,X) x cols (g,a,X).
2. (S2) C(u,v) = Z<u,v> - <u,1><1,v> = bordered 2x2 of the mixed Gram B;
   Sylvester with pivot Z:  det3[C] = Z^2 det4[B].
3. (S3) Cauchy-Binet:  det4[B] = SUM_{|s|=4} w_s l_s(f) l_s(g).
=> P = R(F,G) - R(H,H) = Z^2 SUM_s w_s [l_s(F)l_s(G) - l_s(H)^2] = Z^2 SUM_s q4(s) w_s.

Corollaries: (a) full coefficient law, census 140+210+105+7 = 462 forced;
(b) heavy cofactor exactly 1 (w_h^2 unique in Z^2) -- the 140 machine divisions
explained; (c) q == 0 <=> all 35 q4(s) = 0 with NO collision analysis -- the
branch calculus is retired from the critical path (kept as verification lane);
(d) q4(s) = - prod of the four triple minors => the trichotomy is immediate;
(e) D = Z SUM_t w_t Delta_t^2, so q = compound-sum ratio; explicit face kappa.

## Two-lane adversarial audit (FORCED tag justification)
- Lane 1 (prior session, files synced into container): authored
  t1_windowproof.py (44 checks: symbolic chain m=4,5; q4=-prodDelta symbolic;
  exact instances m=6,7; catalog-ring dictionary equality, all 462 coefficients).
- Lane 2 (this session): independent hand derivation of the identical chain
  before reading lane 1's file; hand verification of (d) on two exact
  configurations (q4=-1 vs prod=+1; q4=-40 vs prod=40); full code audit
  (sound; block-4 dictionary equality is the load-bearing m=7 identity);
  independent execution: **44/44 PASSED** (t1_wp_out.txt).

## Harness status (in-container this session)
t1_core.py           61/61  EXIT 0   (reproduced pre-cutoff)
t1_reduction.py      17/17  EXIT 0
t1_engine.py          4/4   EXIT 0   (P artifact regenerated, 462 monomials)
t1_branches.py       16/16  EXIT 0
t2_temperature.py    12/12  EXIT 0
t3_suspension.py     [001]-[020] reproduced twice; [021]-[024] blocked by the
                     per-call execution ceiling of this sandbox
t3c_partC_exact.py   16/16  EXIT 0   (NEW: certifies Part C [016]-[024] by
                     explicit certificates: level-set annihilator, closure
                     recursion with elementary-symmetric coefficients,
                     Vandermonde lower bound; conditioning stated: Q-linear
                     independence of log2, log3, log5, logphi)
t1_windowproof.py    44/44  PASSED   (OP-1; audit-lane execution)

## Paper (field_surprisal_geometry_v3.tex -> 17 pp, compiled 2x pdflatex, 0 errors)
- NEW Theorem [master window identity] + Corollary [coefficient law] + Remark
  [architecture; Open Problem 3 resolved] inserted after Prop [window].
- Abstract: fifth-front sentence; harness count 6 -> 8.
- Ledger row upgraded to the master form. Open Problems (3) marked resolved.
- Appendix: both new harnesses documented. New bibitem: Horn & Johnson,
  Matrix Analysis, 2nd ed., CUP 2013 (**verify edition/page before formal
  submission** -- citation discipline; the two identities are classical).
- Backup of the pre-edit tex: field_surprisal_geometry_v3.tex.bak.

## Environment notes (operationally relevant for the next session)
- Per-tool-call execution ceiling ~115 s; background processes may survive
  call boundaries but unreliably; file writes sync into the container with
  visible lag (files from the prior session arrived DURING this one).
- Consequence: long harnesses must be chunked or given explicit certificates
  (as done for t3 Part C).

## OP-2 status (not this session's claim)
t4_kwindows.py (k>=3 generalization) and a partial t4_out.txt (through check
[013], mid census) synced in from the prior lane. NOT audited or completed
here. Next session: audit t4 with the same two-lane discipline before trusting
any of its checks; the k=3 double-indicator census was mid-run at cutoff.

## Epistemic notes
- All decision boundaries exact over QQ / QQ[L]; floats display-only.
- The generic-m symbolic proof is at m=4,5; the m=7 catalog case is closed
  identically (all 462 coefficients as ring elements); the prose proof is
  m-uniform via the two classical identities.
- Tag movements: window identity [C] -> [F]. No other tag changed. The branch
  calculus results remain [F]; only their architectural role changed.
