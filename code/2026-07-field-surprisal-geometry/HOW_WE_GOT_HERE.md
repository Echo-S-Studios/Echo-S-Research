# How We Got Here

A narrative reconstruction of the Field Surprisal Geometry project — enough that a
reader (or a fresh session) starting from nothing can rebuild the whole picture:
what the object is, what was proved, in what order, what got corrected along the
way, and what remains. The formal write-up is `paper/field_surprisal_geometry_v3_5.tex`;
every load-bearing claim is machine-checked in `harnesses/` (see `MANIFEST.md`).

---

## 1. The object, in one paragraph

Take a number field's **emission catalog** — a finite set of algebraic seeds, here
`Ω = {√2, √3, √5, φ, τ, φ⁴, K}`. Each seed α has a **Mahler measure**
`M(α) = ∏_{|ζ|>1}|ζ|` (product of the conjugates outside the unit circle), which the
companion *λ=2c* manuscript reads as a **description cost** `log M`. Charging that
cost at inverse temperature β gives a forced maximum-entropy law, the **Mahler–Gibbs
family** `pᵢ(β) ∝ Mᵢ^{−β}`. This project builds the canonical **positive-definite
Fisher–Rao ("surprisal") geometry** of that family and develops it to five fully
forced layers, culminating in a complete classification of which *second statistics*
bend the geometry into a round sphere. Nothing new is assumed about the field — only
its forced data.

The catalog's forced numbers:

| seed | min poly | Mahler M | cost log M | degree | trace | charge group |
|---|---|---|---|---|---|---|
| √2 | x²−2 | 2 | log2 | 2 | 0 | ℤ/2 |
| √3 | x²−3 | 3 | log3 | 2 | 0 | ℤ/2 |
| √5 | x²−5 | 5 | log5 | 2 | 0 | ℤ/2 |
| φ  | x²−x−1 | φ | logφ | 2 | +1 | ℤ/2 |
| τ  | x²+x−1 | φ | logφ | 2 | −1 | ℤ/2 |
| φ⁴ | x²−7x+1 | φ⁴ | 4logφ | 2 | +7 | ℤ/1 |
| K  | x⁴+5x²−5 | φ⁴−1 | ½log5+2logφ | 4 | 0 | ℤ/4 |

Field constant `Z = Σ M = 17 + 4√5`. The costs live in `span_ℚ{log2,log3,log5,logφ}`
(4-dimensional; those four logs are ℚ-linearly independent by a norm argument in
ℚ(√5) — the *only* input the whole classification needs). Two cost coincidences
matter: the **golden tie** `logφ = logφ` (φ and τ share a cost) and the **Salem
square** `(φ⁴−1)² = 5φ⁴`, i.e. `2·cost(K) = cost(√5) + cost(φ⁴)`.

---

## 2. The five forced layers (`§§3–16`, harness `field_surprisal_v2.py`)

Built once and never revisited:

1. **Surprisal is affine in cost** with slope β; **Fisher information = cost
   variance**, `I(β)=Var_β(logM)` (§4).
2. **Dual flatness**: the Gibbs curve is an e-geodesic; KL between temperatures is
   the **Bregman divergence of log Z**; the dual potential is negentropy (§5).
3. **Thermodynamics**: `C(β)=β²I(β)` — Fisher information *is* a heat capacity —
   with a peak-fluctuation temperature β\*≈−0.0768 (§6).
4. **Distances**: Fisher–Rao distance is the Bhattacharyya angle; the seven seeds
   form a regular spherical simplex at mutual distance π; the curve has finite
   information length (§7).
5. The geometry is a **round sphere of curvature ¼** (the categorical model
   `S^{m−1}(2)`); everything downstream lives on it (§4, Lemma).

---

## 3. The heart: which second statistics give a round sphere?

Adjoin a second forced statistic X (degree, charge parity, trace, …) to get a
two-parameter family with Fisher matrix `∇²A`. Its Gaussian curvature is a function
of θ. The driving question: **for which X is that curvature the constant ¼?**

### 3a. The one real early error — and its correction
The first instinct was that constant-¼ meant *totally geodesic* (a great subsphere).
**It does not.** Constant curvature ¼ only fixes the *intrinsic* curvature; totally
geodesic is a strictly stronger *embedding* condition. The suspension proof gives a
round-sphere metric (intrinsic) and says nothing about the embedding. In fact **no
`(logM,X)` surface here is totally geodesic** — logM alone spans a 6-dimensional
algebra, so the value pairs `(logMᵢ,Xᵢ)` never collapse to the 3 a great subsphere
needs (they are always 6 or 7). The single-outcome surfaces are constant-¼ **ruled**
surfaces (swept by great circles), not subspheres. This distinction — caught early —
is maintained everywhere afterward and never slips again. (See `Remark (Correction)`
in the paper for a second, related slip: the throwaway justification "log²M∉V" is
false, since X could *be* log²M; the correct argument is the value count.)

### 3b. The classification (k=2): exactly eight surfaces
`suspension_theorem.py`, `field_surprisal_classification.py`

A single-outcome indicator `𝟙_s` makes the family a **spherical suspension** whose
metric is a round sphere for *any* base curve → constant ¼. The complete answer:
the `(logM,X)` surface is constant-¼ **iff** X is (affinely) the indicator `𝟙_S`
of a set S with **logM constant on S or on its complement** — exactly **eight**
surfaces (seven single-seed + the golden-pair merge `𝟙_{φ,τ}`). Each is forced by
an arithmetic invariant (Mahler distinguishes five seeds; the trace ±1 splits the
golden pair — the same sign the λ=2c construction uses to pick the Perron keystone).

### 3c. Necessity: from machine enumeration to a two-line proof
`t1_core → t1_reduction → t1_engine → t1_branches`, then `t1_windowproof`

*Sufficiency* is the suspension theorem. *Necessity* — that nothing else is
constant-¼ — went through two stages:
- **Machine branch calculus** (v3_3 and earlier): reduce constancy to one scalar
  `q(θ)`; show `Z²Dq = P(w)` is a degree-6 form whose heavy `(3,1,1,1)` coefficients
  are four-point invariants `q₄(s)`; handle collisions via a rank-2 lattice and a
  "branch survival" lemma; finish with a four-point trichotomy and a plane lemma.
  Correct, but resting on 140 machine divisions.
- **The conceptual proof** (v3_2, the flagship simplification): the whole numerator
  factors,
  > **`P(w) = Z²·Σ_{|s|=4} q₄(s)·w_s`**
  by **Sylvester's determinant identity (pivot Z) + Cauchy–Binet** — a bordered-Gram
  expansion. Since the squarefree monomials `{w_s}` are independent, `q≡0 ⇔ q₄(s)=0`
  for every window, *with no collision analysis*. The 140 divisions become two lines.
  A nuance surfaced later (§3e).

### 3d. Catalog invariance + the temperature dichotomy
`t5_catalog_census.py`, `t6_selection.py`

- The count of eight is **not an artifact**: for any catalog it is
  `Σ_ℓ(2^{m_ℓ}−1)` over cost levels (the "eighth surface" is the golden level's
  `2²−1=3`). Verified for four catalogs: full/drop-K/drop-τ/add-√7 → **8/7/6/9**.
- The **temperature** β was the lone `[declared]` input. It resolves as a
  **dichotomy** (`thm:dichotomy`): the metric selects *no* operating point (a
  finite-length curve has isometry group {id,flip}, and the flip swaps distinguishable
  endpoints), *and* the candidate principles (arc-length midpoint, max-Fisher,
  max-entropy, heat-capacity peak) select **provably distinct** temperatures. So the
  `[declared]` is the honest form of the answer, not a gap.

### 3e. Higher k: the identity lifts; necessity needs one more idea
`t3_suspension`, `t4_kwindows`, `t4b_census_fast`, `t7_knecessity`

- **Iterated joins** give round ¼-spheres in *every* dimension; the first *totally
  geodesic* family appears at exactly `k=5` (uniquely `ℝ[logM]`, the golden-merged
  simplex).
- The master identity **lifts to all k** (Cauchy–Binet over (k+2)-subsets), giving
  windowed Gauss obstructions and a clean sufficient criterion: **every window matrix
  rank ≤ 1 ⇒ constant ¼**.
- **The nuance in §3c**, resolved: on the *family* (not as a free polynomial) the
  window symbols `w_s(θ)` can *collide* when a statistic ties on {φ,τ}. The resolving
  idea — the last real insight of the project — is that **window collisions are
  squarefree**: `𝟙_s − 𝟙_{s'} ∈ {−1,0,1}`, so the Salem square (which carries entry
  **2** at K) *cannot act at window level*; only the golden swap can, and it only
  fires when everything ties on {φ,τ}, making φ,τ **twins** — twin-containing windows
  vanish outright and swapped windows are equal, forcing each `q₄`. This is the *same*
  multiplicity-2 mechanism that made the master identity's heavy cofactor exactly 1.
  It closes **k-necessity** (constant ¼ ⇔ rank-≤1 windows) for `2≤k≤5`.

### 3f. The landscape — and a refuted conjecture
`t9_landscape.py`

The natural conjecture was that constant-¼ families are exactly the *within-level
indicator joins*. **False.** The complete `k≥3` landscape is the **partitioned-affine
classification** (`thm:paclass`): constant-¼ ⇔ the span is `V(π)` for a partition π
into within-level clusters and "line blocks," with `dim V(π)=k+1`. Catalog counts:

> **k = 2, 3, 4, 5, 6 → 8, 56, 95, 31, 1**

At k=3 that's 26 double-indicator classes **plus 30 "split-affine" classes**
`⟨1, a, 𝟙_B, a·𝟙_B⟩` — genuinely new families that contain no indicators beyond
`𝟙_B, 𝟙_{B^c}`, refuting the conjecture. (The mechanism is a **moment trichotomy**:
a window's circuit of size 2/3/≥4 forces its matrix to rank 0/1/≥2, so constant-¼ ⇔
every circuit ≤ 3.) The earlier 26-class census stands exactly, as the
double-indicator stratum.

### 3g. Coupling two fields
`t8_compositum.py`, `t10_coupled.py`

Two *independent* catalogs form a Riemannian product (block-diagonal Fisher, zero
cross-curvature); a coupling would be a `[declared]` warp. But over the **compositum**
the coupling is **forced**: the tensor cost is *not additively separable*, with exact
interaction quanta **log2** (golden/rational window) and **−6logφ** (golden sector),
the rational sector `{√2,√3,√5}` staying separable. The coupled family genuinely
**curves** (nonzero Gauss obstruction at the uniform point; a real 5-circuit), yet its
cross-covariance vanishes there — the coupling is invisible to `Cov` and visible to
curvature. And the interaction tensor Δ is **not determined by the charge data** (two
same-charge pairs have different Δ).

---

## 4. The verification discipline

Every decision boundary is **exact** (symbolic `sympy` over `ℚ(log2,log3,log5,logφ)`,
`+log7` for the extended catalog); floats are **display-only**; nonvanishing at the
true catalog is certified by **interval arithmetic** (`mpmath.iv`, dps 60), never by
a rational specialization. Every harness is **fail-first**: a `ck(...)` raises and
exits non-zero on the first failed assertion, so a clean exit *is* the proof of every
line it printed. Claims are tagged `[forced]` (theorem/exact), `[computed]`
(validated numeric), `[declared]` (a principled invariance-free choice), `[open]`.
Numeric curvature routines are validated against four metrics of known curvature
(sphere +1, hyperbolic −1, flat 0, simplex ¼) before use.

Independently of the authoring sessions, every load-bearing result was **re-derived
during audit** — not merely re-run — and confirmed (the master identity, the module
dimensions, the trichotomy, the squarefree lattice, the twin symmetry, the compositum
contrasts, the census counts including an independent partition enumeration landing on
56, the dichotomy positions, the split-affine refutation).

---

## 5. Corrections that were caught and fixed (the honest trail)

- **constant-¼ ≠ totally geodesic** — the one real early error; corrected and never
  repeated. Nothing here is totally geodesic (§3a).
- **"log²M ∉ V" is false** as a blanket justification (X can be log²M); the value-count
  argument is the correct proof (`Remark (Correction)`).
- **Branch-calculus "retirement" was slightly overstated** in v3_2: on the family the
  window symbols can collide, so necessity for collision cases needs handling — supplied
  cleanly by the branch-rank criterion and then subsumed by the squarefree mechanism
  (§3e).
- **A transient count of 61 → corrected to 56** at k=3 (five golden overlaps double-
  counted); the shipped value is 56, and an independent enumeration confirms it (§3f).

In every case the *shipped* artifact was correct or was fixed before shipping; the
conjecture in §3f was *stated as a conjecture* and then *refuted*, not quietly amended.

---

## 6. Where it stands, and what's left

The paper's **enumerated open-problem ledger is fully resolved.** What remains is
outward-facing **outlook**, not gaps:

- **Fat-level catalogs** (cost levels of multiplicity ≥ 3 admit vertical line blocks
  the current classification doesn't cover verbatim).
- **Composita of two *distinct* catalogs** (only a single field's self-compositum is
  treated).
- **The coupled family's curvature landscape away from the uniform point.**
- The standing **ℚ-linear-independence conditioning** of the collision analysis (used
  throughout as the one hypothesis).
- **Temperature stays `[declared]`-final** by the dichotomy — an answer, not a gap.

---

## 7. Reading order

1. This file.
2. `MANIFEST.md` — which harness proves what, and how they depend on each other.
3. `./run_all.sh` — reproduce everything (`pip install -r requirements.txt` first).
4. `paper/field_surprisal_geometry_v3_5.tex` — the formal development, section by
   section (the `§` pointers above map into it).
5. `dev-log/` — the per-session records, in order, including the mid-session
   correction trails.
6. `handoffs/` — the open-problem handoff documents that seeded the later sessions
   (useful as worked examples of how each front was scoped and attacked).

The through-line, in one sentence: a numerically-observed classification became
machine-proven, then conceptually proven (Sylvester + Cauchy–Binet), then
catalog-invariant with a temperature dichotomy, then closed at every k by a single
squarefree mechanism plus a forced compositum coupling, and finally widened into the
partitioned-affine landscape — with the one distinction that mattered most, constant
curvature versus totally geodesic, kept straight from the moment it was first
confused.
