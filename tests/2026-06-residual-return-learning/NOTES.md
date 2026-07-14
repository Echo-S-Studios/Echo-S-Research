# Verification notes — *Residual Return: Exact Learning Dynamics and Language over the Vector Substrate*

Paper: `papers/2026-06-residual-return-learning/residual_return_learning.tex`
Author: AceTheDactyl (@AceTheDactyl), Echo S Studios, June 2026.

The paper is the dynamic/linguistic companion to *The Vector Substrate*. It assembles an exact
streaming learner over a number field `K=Q(sqrt2+sqrt3)`, a cross-field growth mechanism (including a
non-disjoint compositum construction), an automatic detector-driven loop, and a residual-valued
"language" over the Clifford algebra `Cl(2,0) = M_2(R)`, unified by the golden-law keystone `x^2-x-1`.

All displayed numbers, matrices, polynomials, factorisations, and the SHA-256 witness digest are
**independently re-derived** here (sympy exact symbolic; mpmath dps>=40 for transcendental/Mahler
checks) from the paper's stated premises and compared to the paper's values. Trace-form Grams are
rebuilt from the regular representation (companion-matrix power traces), never transcribed.

**Final run:** `71 passed, 1 xfailed` (0 failed, 0 errors).

---

## Claims table

| Claim (paper location) | How tested | Status |
|---|---|---|
| Power-basis trace Gram G (eq Gpower) | Rebuild G_ij=Tr(θ^{i-1}θ^{j-1}) from trace(companion^k) | verified |
| G positive definite (totally real) | Sylvester minors > 0 and eigenvalues > 0 | verified |
| Exact projector idempotent P^2=P (Prop 2.1) | Build P=B(BᵀGB)⁻¹BᵀG over Q, assert P²=P, PB=B | verified |
| Worked episode w=2√6=θ²−5, coords (−5,0,1,0) (Ex 2.4) | Symbolic (√2+√3)²−5=2√6; G-orthogonality to col(B) | verified |
| Residual r=w, ‖r‖²_G=Tr(24)=96 (Ex 2.4) | Compute r=x−Px and rᵀGr; independent trace route 24·4 | verified |
| Seed minpoly x²−24 (Ex 2.4/3.3) | sympy minimal_polynomial(2√6) | verified |
| Witness digest `31f1f1e05ac9a35a` (Ex 2.6) | Recompute SHA-256("genesis"+canonical JSON)[0:16] from full body | verified |
| Bridge φ→x²−x−1 (Ex 3.3) | charpoly of ½I+½·rep(√5) | verified |
| Bridge 2√6→x²−24, 3→x−3, θ→x⁴−10x²+1 (Ex 3.3) | sympy minimal_polynomial | verified |
| Bridge ½ rejected (not algebraic integer) (Ex 3.3) | Monic minpoly x−½ ∉ Z[x] | verified |
| φ⊕φ vs C((x²−x−1)²): same charpoly/trace/Mahler, not similar (Rem 3.4) | Equal charpoly & trace 2; A²−A−I=0, B²−B−I≠0; shared Mahler φ² | verified |
| Jordan J₂⊕J₂ vs J₂⊕J₁⊕J₁: charpoly+minpoly insufficient (Rem 3.4) | Equal charpoly x⁴ & minpoly x²; rank 2≠1 ⇒ not similar | verified |
| Disjoint Grams G_K=diag(4,8,12,24), G_L=diag(2,14) (Ex 3.7) | Tr(c)=deg·c on square-root basis | verified |
| Kronecker Gram diag(8,56,16,112,24,168,48,336) (Ex 3.7) | G_K⊗G_L | verified |
| Out-of-field √7 score 56 (Ex 3.7) | Tr_{KL/Q}(7)=8·7 | verified |
| det relation 9216²·28⁴ = ∏diag (Rem 3.8) | det(G_K⊗G_L) vs (detG_K)²(detG_L)⁴ vs ∏ of 8 entries | verified |
| Non-disjoint: factor over Q(√2) (Ex 3.9) | Expand (x²−2√2x−1)(x²+2√2x−1)=x⁴−10x²+1 | verified |
| Operator poly x⁶−25x⁴+91x²−75=(x²−3)(x⁴−22x²+25) (Ex 3.9) | Squarefree part of Res_y(y²−2, m_β(x−y)) | verified |
| Selected factor m_θ=x⁴−22x²+25 = minpoly(2√2+√3) (eq mtheta) | sympy minimal_polynomial | verified |
| Spurious x²−3 from α=−√2 embedding (θ→√3) (Ex 3.9) | −√2+(√2+√3)=√3, √3²=3 | verified |
| β=−7/20·θ+1/20·θ³ reconstruction (Ex 3.9) | Symbolic expand equals √2+√3 | verified |
| Northcott gate (deg,height)=(2,24), Σcᵢ²=577 (Ex 4.2) | Coefficient arithmetic | verified |
| Landau certificate Mah²≤577 (x²−24), Mah≤8 (x²−7) (Ex 4.2) | mpmath Mahler = 24, 7; bounds hold | verified |
| Gate GROW/REJECT/STOP on integers (Ex 4.2) | Integer-comparison decision reproduces verdicts | verified |
| Fisher matrices for Q(√5), Q(√2,√3), Q(√2,√3,√7) (Ex 4.4) | Fisher=(1/n)(G−(1/n)ttᵀ); match displayed diags | verified |
| √5=2φ−1=(−1,2), ‖√5‖²_G=10=n·Fisher (Ex 4.4) | Quadratic forms of G and Fisher | verified |
| G=n·Fisher on trace-zero subspace (Ex 4.4) | Symbolic: diff = (tᵀr)²/n, vanishes when tᵀr=0 | verified |
| λ=2c identity, D_KL=(1/2c)‖r‖²_G (Thm 4.5) | Symbolic KL with Fisher=G/c; read-off λ=2c, λ|c=1=2, λ|c=n=2n | verified |
| Smyth μ_S=1.3247…, plastic root of x³−x−1 (Rem 4.6) | mpmath findroot, μ_S³−μ_S−1=0 | verified |
| Floors 0.5624, 2.2496, 4.4992 from 2c·log μ_S (Rem 4.6/4.10) | mpmath 2logμ_S, ×4, ×8 | verified |
| Certified GROW: 2log7≤3.89184≤56, 2log24≤6.35612≤96 (Ex 4.3b) | mpmath interval enclosures | verified |
| STOP: lattice noise 1/10 < 0.5624 (Ex 4.3b) | Rational comparison | verified |
| C-ladder eigenvalues (−1±√(1+4C))/2 (Prop 4.7) | charpoly R_C = x²+x−C, eigenvals | verified |
| Self-action spectrum {0,0,±√(1+4C)} (Prop 4.7) | 4×4 adjoint [R_C,·] eigenvalues | verified |
| 0-eigenspace = centraliser span{I,R_C} (Prop 4.7) | Nullspace of adjoint, dim 2, contains I & R_C | verified |
| Golden gap √5=φ−ψ at C=1 (Prop 4.7) | Symbolic + adjoint eigenvalue | verified |
| R_1 & keystone share self-action gap √5 (both disc 5) (Prop 4.7) | tr²−4det=5 for both; equal nonzero adjoint eigenvalues | verified |
| **"R_1 conjugate to keystone R" (Prop 4.7 aside)** | charpoly(R_1)=x²+x−1 ≠ x²−x−1=charpoly(R) | **xfail (flagged)** |
| Frame-shift gap √(1+4C)=√Mah (Def 4.8) | Mah(x²−(1+4C))=1+4C | verified |
| c=√(1+4C)/(2C), λ=√5, c=√5/2 at C=1 (Def 4.8) | Symbolic algebra | verified |
| Canon table gaps {√2,√3,√5}, c {2√2,√3,√5/2} (Table 4) | Evaluate on C∈{¼,½,1} | verified |
| Matrix iso is algebra hom (e₁²=e₂²=1,i²=−1,i=e₁e₂) (Sec 5.1) | mat(·) generator relations | verified |
| Keystone mat(R)=[[0,1],[1,1]], R²=R+I, det−1, τ=½ (Sec 5.2) | Matrix identities & scalar part | verified |
| Gate P₀=½(1+e₂) idempotent (Ex 5.3) | mat(P₀)=[[1,0],[0,0]], P₀²=P₀ | verified |
| Cayley–Hamilton X²−tr·X+det·I=0 (Sec 5.1) | On keystone + generic holding | verified |
| L(R)=5/2+e₁−½e₂≠0, R∉kerL (Prop 5.8) | L(R)=2R²−R=R+2I | verified |
| ker L=span{e₁+2e₂, i}, dim 2 (Thm 5.4) | Exact nullspace of 4×4 L | verified |
| Sylvester eigenvalues {√5,0,0,−√5} (Thm 5.4) | Eigenvalues of L matrix | verified |
| Disproof L(H(0,−1,1,0))=H(−3,0,0,0)≠0 (Rem 5.5) | Compute L; cl([[0,−1],[1,0]])=i∈kerL | verified |
| Commit projector eq(proj) exact matrix, ≠I (Thm 5.6) | K(KᵀK)⁻¹Kᵀ | verified |
| Generalisation: E₁, E₁+1 merge; i, 2i distinct (Thm 5.6/Ex 5.7) | commit(·) values | verified |
| Lexicon: 5 tokens → 4 entries (Ex 5.7) | Distinct committed residues | verified |
| Firewall counts 20+5=25, +2=27 (Table 6) | Arithmetic | verified |
| φ-keystone one object, 3 routes, Mah=φ (Prop 6.2) | companion=mat(R)=cl(½,1,−½,0); Mahler=φ | verified |

---

## VERIFIED (independently re-derived, exact where possible)

All 71 passing tests above. Highlights where the independent re-derivation is strongest:

- **Trace-form Grams** rebuilt from the regular representation (`trace(companion(m_θ)^k)`), not
  transcribed — the eq (Gpower) matrix `[[4,0,20,0],[0,20,0,196],[20,0,196,0],[0,196,0,1940]]`
  and the Fisher/Kronecker diagonals all fall out.
- **Non-disjoint witness** (the paper's headline "resolved open problem"): the degree-6 operator
  polynomial `x⁶−25x⁴+91x²−75` is reproduced as the *squarefree part of the resultant*
  `Res_y(y²−2, m_β(x−y))` (its distinct roots are the 6 distinct values α_i+β_j), it factors as
  `(x²−3)(x⁴−22x²+25)`, `x⁴−22x²+25` is verified to be `minpoly(2√2+√3)`, and the reconstruction
  `β = −7/20·θ + 1/20·θ³ = √2+√3` is confirmed symbolically. Fully consistent.
- **Witness digest** `31f1f1e05ac9a35a` reproduced bit-for-bit from the full canonical JSON body.
- **Language kernel & commit projector**: the 4×4 matrix of `L(X)=RX+XR−X` has exact rational
  nullspace `span{e₁+2e₂, i}` and the commit projector matches eq (proj) exactly, `≠ I`.
- **λ=2c** and the **self-action spectrum** `{0,0,±√(1+4C)}` verified symbolically for the whole
  ladder; the frame-shift `c=√(1+4C)/(2C)` and canon-table values `{2√2,√3,√5/2}` all check out.

## FAILED / FLAGGED for human review

**1 item — `test_R1_conjugate_to_keystone_literal` (xfail, strict).**

- **Location:** Proposition 4.7 (`prop:trifurcation`), parenthetical aside in the last sentence:
  *"…R_1 is conjugate to the keystone R=(0 1; 1 1) (charpoly=x²−x−1; both have discriminant 5,
  hence the same self-action gap)."*
- **Expected by paper:** `R_1` conjugate (similar) to the keystone `R`.
- **Reproduced:** **Not similar.** `R_1 = [[0,1],[1,-1]]` has characteristic polynomial `x²+x−1`;
  the keystone `R = [[0,1],[1,1]]` has `x²−x−1`. The characteristic polynomial is a similarity
  invariant, so two matrices with different char. polynomials are **not conjugate over any field**.
  (`R_1`'s eigenvalues are `{−φ,−ψ}`, so `R_1` is in fact similar to `−R`, not to `R`.)
- **Assessment — minor wording slip, not a substantive error.** The proposition's actual, load-bearing
  content — the self-action **spectrum** `{−√(1+4C),0,+√(1+4C)}`, the centraliser kernel, and at
  `C=1` the **gap √5 because both `R_1` and `R` have discriminant 5** — is fully correct and is
  independently verified (`test_R1_and_keystone_share_self_action_gap`,
  `test_ladder_selfaction_spectrum`, `test_golden_gap_is_sqrt5`). The word "conjugate" should read
  "shares the same self-action gap (both discriminant 5)" (or "conjugate to `−R`"). Flagged so a human
  can correct the phrasing; nothing downstream depends on the literal conjugacy.

No other value failed to reproduce from the paper's own premises.

## UNTESTABLE (documented, not mechanically checkable in this archive)

These are design/engineering guarantees, external-repo provenance, or cited theorems — not
mathematical derivations with a reproducible number. The *mathematical* core each rests on **is**
tested above.

1. **External suite counts** — L00M training `127/137`, full repo `544/574`, kira-language `121`,
   ZFP `74/74`, substrate probe `13`, companion probe `20`, calibration `4`. These live in the
   L00M / Plate-Matrices / kira-language repositories, not this paper archive.
2. **"Machine-verified by `test_X` (commit `Y`)" provenance** throughout — references to external
   test names and commit hashes; not reproducible here (we re-derive the underlying math instead).
3. **Witness-chain tamper-evidence** (design property: any edited field flips verification to False).
   We verify the digest *value*; the tamper mechanism is implementation behavior.
4. **Firewall runtime wiring** (only THEOREM/COMPUTED cross a wire). We verify only the arithmetic
   counts (25 wired, 27 bank); the dispatch behavior is implementation behavior.
5. **Lexicon content-hash ids** `SHA-256(canonical Fraction string)[0:16]` — the exact canonical
   string format is not specified in the paper, so the specific 16-hex ids cannot be reproduced
   (determinism/order-independence is an implementation property; the *values* are unspecified).
6. **Welford recurrence & variance-calibration gate** (eq calib `Σ M2_i ≤ τ n Σ mean_i²`) — an
   algorithmic design gate with no numeric instance given to reproduce.
7. **Čencov's theorem** (Fisher metric unique up to positive scale ⇒ `c` irreducible by invariance) —
   a cited external impossibility theorem, not a computation.
8. **Northcott finiteness** (admissible seed set finite) — a cited theorem; we verify the concrete
   Landau/height *certificates*, not the finiteness meta-claim.
9. **Design guarantees**: sole-mutator (confirm), propose-idempotence, one-way module isolation,
   float-rejection at intake, "no float in a decision". These are guarantees about the code; the
   mathematics they protect (exact capture criterion, residual return, exact projector) is tested.
10. **The KL expansion being "exact for the Gaussian location family"** — a modeling premise; we
    verify the resulting algebraic identity `D_KL=(1/2c)‖r‖²_G` given `Fisher=G/c`.

---

### How to run

```
py -m pytest tests\2026-06-residual-return-learning -v -p no:cacheprovider
```

## Correction applied 2026-07-04

Prop 4.7 aside corrected "R_1 conjugate to keystone R" -> "R_1 conjugate to -R"; test test_R1_conjugate_to_negative_keystone now asserts R_1 ~ -R (shared charpoly x^2+x-1) and passes. Shared self-action gap sqrt5 unchanged.
