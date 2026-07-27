# Handoff — Field Surprisal Geometry: the three remaining open problems

**For:** a new session continuing the Field Surprisal Geometry project for Ace / AceTheDactyl (Echo‑Squirrel Research).
**Goal of this session:** complete Open Problems **(1)** the non‑indicator gap, **(2)** the intrinsic‑temperature question, **(3)** the higher‑dimensional classification. Each is stated precisely below with status, attack plan, and success criteria.

> Read this whole file first. Then ask Ace to re‑upload the harness files listed in §3 (they are the reference implementations; you can also rebuild from the data in §3.4). Re‑fetch the corpus (§3.5) if you need the source papers. Do **not** re‑derive settled results or contradict them — they are verified.

---

## 1. What the project is (one paragraph)

We attach a canonical **positive‑definite Fisher–Rao ("surprisal") geometry** to a number field's *emission catalog* — the seven algebraic seeds `{√2, √3, √5, φ, τ, φ⁴, K}`. Their Mahler measures define a forced one‑parameter exponential family (the **Mahler–Gibbs family** `p_i(β) ∝ M_i^{−β}`), and everything downstream — dual flatness, thermodynamics, distances, and the two‑statistic surfaces — is built on it. The main deliverable is a machine‑verified LaTeX/PDF whitepaper (`field_surprisal_geometry_full.*`, ~12pp) that a reviewer can rederive end‑to‑end. This session extends the classification of the two‑statistic surfaces and pushes to higher dimension.

---

## 2. Verified state (do not contradict; cite as established)

Everything here is proven or machine‑checked. Tags: **[F]** forced/exact, **[C]** computed/validated‑numeric, **[D]** declared, **[O]** open.

- **Catalog & measure.** Mahler measures `M = {2, 3, 5, φ, φ, φ⁴, φ⁴−1}`; partition `Z = Σ M_i = 17 + 4√5` **[F]**. Canonical family `p_i(β) = M_i^{−β}/Z(β)`, `A(β)=log Z(β)`; sufficient statistic `log M` **[F]**; `β` is **[D]** (see Task 2). Distinguished β: `0` (uniform), `1` (MDL/Occam prior), `λ=√5` (the field's exchange rate, from the λ=2c paper), `β*≈−0.0768` (peak fluctuation).
- **Surprisal geometry.** Surprisal is affine in cost, `S_i(β)=β log M_i + log Z` **[F]**; Fisher `I(β)=A''(β)=Var_β(log M)` **[F]**; the Fisher–Rao metric `g_ij=δ_ij/p_i` on `Δ^{m−1}` is isometric via `p↦2√p` to the sphere `S^{m−1}(2)`, constant curvature `1/4` **[F]** (radius‑2 sphere).
- **Dual flatness** (Amari): Gibbs curve is an e‑geodesic; `KL(p(β₁)‖p(β₂)) = A(β₂)−A(β₁)+(β₂−β₁)U(β₁)` = Bregman divergence of `log Z`; dual potential `A* = −H` (negentropy) **[F]**. Here `U(β)=⟨log M⟩_β = −A'(β)`.
- **Thermodynamics.** `C(β)=β² I(β)=β² Var_β(log M)` (Fisher information *is* a heat capacity) **[F]**. `S=H=A+βU`.
- **Distances.** `d_FR(p,q)=2 arccos Σ√(p_i q_i)` (Bhattacharyya angle); the seven vertices are mutually at distance `π` **[F]**. The Gibbs curve is an e‑geodesic but **not** a metric geodesic (distances are sub‑additic along it).
- **Keystone sub‑model.** On `{Rⁿ}`, `M(Rⁿ)=φⁿ`, cost `n log φ` (linear) → shifted‑geometric Gibbs law, `I(β)=(log φ)² r/(1−r)²`, `r=φ^{−β}` **[F]**.
- **Interaction with the trace form.** Embedding obstruction: the flat trace‑form geometry admits **no** totally geodesic embedding; any isometric surface realization forces `k₁k₂=−1/4` (Gauss) **[F]**. Fisher monotonicity: coarse‑graining conjugate data → catalog outcomes contracts Fisher (data‑processing) **[F]**.
- **Product geometry.** Two independent catalogs → Riemannian product `S^{m₁−1}(2)×S^{m₂−1}(2)`, block‑diagonal Fisher, **zero cross‑curvature** **[F]**; a genuine coupling is a Čencov‑declared warp **[D]**.
- **Charge grading** (from the `charge-measure-coupling` paper). Charge group of a seed = least `n` with all conjugates `α` satisfying `αⁿ ∈ ℝ₊`; `χ_n(α)=round(nθ/2π) mod n`. Catalog charge groups: `√2,√3,√5,φ,τ → ℤ/2`; **`φ⁴ → ℤ/1`** (unique all‑real‑positive seed); **`K → ℤ/4`** (unique complex seed). Charge **parity** = `𝟙_{φ⁴}` (single‑outcome indicator).

### 2.1 The two‑statistic classification (the heart — and the corrected result)

Second statistic `X` (a = `log M` is fixed). Family `p_i ∝ exp(θ₁ a_i + θ₂ X_i)`, Fisher metric `g=∇²A` (2‑D Hessian metric). Let `V = span{1, log M, X}` (dimension 3). The geometry depends **only on V** and is invariant under affine reparametrization of `(a,X)` — use this reduction everywhere.

**Second fundamental form** (in the unit √‑embedding `s=√p`, `τ_α = T_α − ⟨T_α⟩` centered):
```
II_αβ = ¼ (τ_α τ_β s)^⊥ ,   ⊥ = orthogonal complement of span{s, τ₁s, τ₂s} in ℝ^m
II_αβ = 0  ⟺  the pointwise product  T_α T_β ∈ V
```
So `II_aa=0 ⟺ (log M)² ∈ V`; `II_ac=0 ⟺ (log M)·X ∈ V`; `II_cc=0 ⟺ X² ∈ V`.
Gauss equation: `K_int − 1/4 = (⟨II_aa,II_cc⟩ − |II_ac|²)/det(g̃)`. **Constant curvature 1/4 ⟺ det II := ⟨II_aa,II_cc⟩ − |II_ac|² ≡ 0.**
Inner products (analytic handle for Task 1): with central moments `μ`,
```
16·⟨II_αβ,II_γδ⟩ = μ_{αβγδ} − g_{αβ}g_{γδ} − Σ_{pq} μ_{αβp}(g⁻¹)_{pq}μ_{γδq}
```
(`μ_{αβγδ}` = 4th central moment of `(a,X)` under `p(θ)`, `μ_{αβp}` = 3rd, `g_{αβ}` = covariance.)

**Established facts (verified in `field_surprisal_classification.py`, 9/9):**
- **No `(log M, X)` surface is totally geodesic**, for *any* `X` **[F]**. Reason: `log M` alone spans a **6‑dimensional algebra** (`1, log M, …, log M⁵` independent over the 7 seeds), so `(log M)² ∉ V` (any 3‑dim `V`), i.e. `II_aa ≠ 0` always. Equivalently the number of distinct `(log M, X)` values is always 6 or 7, never the 3 a great 2‑subsphere needs.
- The single‑outcome surfaces are **constant‑curvature‑1/4 RULED surfaces**: for `X=𝟙_s`, idempotency gives `II_cc=0` and localization gives `II_ac=0` (so `det II=0`), while `II_aa≠0`. They are swept by great circles in the `X`‑direction; intrinsically round `1/4`‑spheres, not subspheres.
- **Complete constant‑curvature classification (among indicators):** `𝟙_S` is constant‑curvature‑1/4 **iff `log M` is constant on `S` or on `Sᶜ`**. Exactly **16 subsets** = **8 distinct surfaces** (`𝟙_S` and `𝟙_{Sᶜ}` share a surface): the seven single‑seed indicators **plus the golden‑pair merge `𝟙_{φ,τ}`** (which is the single‑outcome indicator of the *effective* catalog where the Mahler‑tied seeds φ,τ merge) **[F]**.
- **Forcing of the seven single‑seed surfaces:** Mahler measure singles out `√2 (M=2), √3 (M=3), √5 (M=5), φ⁴ (M=φ⁴), K (M=φ⁴−1)`; the only Mahler tie is `{φ,τ}`, split by the **trace** `tr(φ)=+1, tr(τ)=−1` — the same ±1 the λ=2c construction uses to select the Perron keystone `R²=R+I` over the Clifford gate `R²=I−R`.
- **No non‑indicator surface is constant‑curvature‑1/4** in a ~5000‑sample search (random + the full Möbius candidate family) **[C]** — this is the one gap (Task 1).

> ⚠️ **Cautionary tale (read this).** The previous session repeatedly wrote *"totally geodesic"* for what is actually *"constant curvature 1/4."* They are **different**: totally geodesic = great subsphere (II ≡ 0); constant‑curvature‑1/4 only fixes intrinsic curvature. The suspension proof gives the round‑sphere **metric** (intrinsic) and says nothing about the embedding. **Nothing in this catalog is totally geodesic.** Keep the two notions strictly separate. This error was caught only by checking the embedding (the `d`‑count and `II_aa`); do the same for any "geodesic/subsphere" claim.

---

## 3. Toolkit

### 3.1 Reference harness files (ask Ace to re‑upload; they persist in the project outputs)
- `field_surprisal.py` (12/12) — core: `Z`, affinity, Fisher=Var, curvature 1/4, keystone geometric sub‑model, embedding obstruction, Fisher monotonicity.
- `field_surprisal_v2.py` (12/12) — Tier‑1: dually‑flat Bregman‑KL, thermodynamics `C=β²I`, distances, the two‑statistic Fisher matrix + validated numeric curvature routine.
- `suspension_theorem.py` (9/9) — the suspension theorem (curvature 1/4 for indicator families; ruled metric structure).
- `field_surprisal_tier2.py` (15/15) — multi‑statistic landscape, product geometry, charge grading, apex forcing.
- `field_surprisal_classification.py` (9/9) — **the classification** (no TG; the 8 constant‑curvature surfaces; II conditions; non‑indicator search). **Start here for Tasks 1 & 3.**

### 3.2 Key functions (rebuild if needed)
- `K_at(a, X, t, s)` — exact Gaussian curvature (Brioschi) of the `(a,X)` Gibbs surface at natural‑parameter point `(t,s)`, via cumulants; plain‑float, fast. In `field_surprisal_classification.py`.
- `numK(gfun, pt, h)` — finite‑difference Gaussian curvature of an arbitrary 2‑D metric. **Always validate** against known metrics before trusting: simplex→0.25, flat→0, hyperbolic (`diag(1/y²,1/y²)`)→−1, unit sphere (`diag(1,sin²θ)`)→+1. FD noise ~1e‑3.
- II conditions via linear algebra: `II_αβ=0 ⟺ T_αT_β ∈ V`; check by `rank([1, log M, X])` vs `rank([1, log M, X, T_αT_β])` in SymPy (exact).
- `d`‑count (totally‑geodesic test): `d = #distinct (log M, X) pairs`; TG requires `d=3`.

### 3.3 Discipline (Ace's conventions — hold the line)
- **Exact arithmetic at every decision boundary** (SymPy over ℚ/ℚ(√5)); floats **display‑only**. Numeric curvature is fine as **[C]** *after* validation; prefer exact **Brioschi/symbolic** where a clean answer is possible (as in Tasks 1a).
- **Fail‑first harnesses**: every load‑bearing claim asserted in a `ck(...)` that raises on failure; print `N/N passed`.
- **Tag every claim** `[F]/[C]/[D]/[O]`. Do not upgrade `[C]` or `[O]` to `[F]` without a proof. Do not call something "forced" that rests on a modeling choice.
- **Ill‑conditioning warning**: large statistic ranges → measure concentrates → near‑singular Fisher → spurious finite‑difference curvature. If a numeric result looks surprising, re‑check symbolically.
- **Citations**: the corpus is real (fetched, not remembered). Cite paper titles/shortnames; never invent results. If you state a fact "from the corpus," verify it in the fetched PDF.

### 3.4 Catalog data (exact)
| seed | min poly | Mahler M | degree | trace | charge group |
|---|---|---|---|---|---|
| √2 | x²−2 | 2 | 2 | 0 | ℤ/2 |
| √3 | x²−3 | 3 | 2 | 0 | ℤ/2 |
| √5 | x²−5 | 5 | 2 | 0 | ℤ/2 |
| φ | x²−x−1 | φ=(1+√5)/2 | 2 | +1 | ℤ/2 |
| τ | x²+x−1 | φ | 2 | −1 | ℤ/2 |
| φ⁴ | x²−7x+1 | φ⁴=(7+3√5)/2 | 2 | +7 | ℤ/1 |
| K | x⁴+5x²−5 | φ⁴−1=(5+3√5)/2 | 4 | 0 | ℤ/4 |

`Z=17+4√5`. `deg=(2,2,2,2,2,2,4)`, `trace=(0,0,0,1,−1,7,0)`, `charge_order=(2,2,2,2,2,1,4)`, `charge_parity=𝟙_{φ⁴}=(0,0,0,0,0,1,0)`. K‑roots: `x²=(−5±3√5)/2` → real `±K` (`|K|=5^{1/4}/φ<1`, inside) and imaginary `±iβ` (`β²=(5+3√5)/2>1`, outside), so `M(K)=β²=φ⁴−1`.

### 3.5 Corpus
`https://echo-s-studios.github.io/Echo-S-Research/` → `papers/`. Fetch with `curl -sL`. Structured index at `papers/catalog.json`. Directly relevant: `2026-06-lambda-2c.pdf` (the exchange rate; cost = log Mahler), `2026-06-emission-gap.pdf` (no Salem in the spectral image; angle confinement), `2026-06-charge-measure-coupling.pdf` (the conjugate‑angle charge & parity‑graded Mahler floor — **Task 3 relevance**), `2026-07-relational-charge.pdf` (relational ℚ/ℤ charge — refinement). Extract text with `pdftotext -layout`.

### 3.6 Environment
`pdflatex` and `xelatex` are at `/usr/bin/`. Compile with two passes. Put deliverables in `/mnt/user-data/outputs/` and call `present_files`. `pip install --break-system-packages` if needed. The whitepaper preamble/theorem environments are in `field_surprisal_geometry_full.tex` — extend that document; do not restart it.

---

## 4. TASK 1 — Close the non‑indicator gap  *(the one honest [O]; highest mathematical stakes)*

**Statement.** Prove (or refute with a counterexample) that the eight constant‑curvature‑1/4 surfaces are the *only* ones — i.e. **`det II ≡ 0` (for all θ) ⟹ `X` is affinely an indicator `𝟙_S` with `log M` constant on `S` or `Sᶜ`.** Currently proven for all indicators (exact finite check) and evidenced for non‑indicators by a ~5000‑sample search; a closed proof is missing.

**Decompose into two sub‑goals:**

- **1a (tractable, exact).** Solve the system **`X² ∈ V` and `(log M)·X ∈ V`** completely over the catalog (`V=span{1,log M,X}`). This is a finite algebraic system: `X_i² = c₀+c₁ log M_i + c₂ X_i` and `log M_i·X_i = d₀+d₁ log M_i + d₂ X_i` for `i=1..7`, unknowns `X_i` and the `c,d` coefficients, modulo the affine freedom `X ↦ αX+β log M+γ`. Solve exactly in SymPy (elimination / Gröbner; treat `log 2, log 3, log 5, log φ, log φ⁴, log(φ⁴−1)` as independent transcendentals, or solve numerically then recognize). **Expected/target result:** the only solutions are the indicators `𝟙_S` with `log M` constant on `S`. This proves **`II_ac=II_cc=0 ⟺ indicator‑with‑constant‑log M`**. *(Note the known partial fact: `II_ac=0` alone forces `X` to be a Möbius function of `log M` off one pole; indicators are the degenerate case. 1a checks whether adding `II_cc=0` leaves only indicators.)*

- **1b (the hard part).** Show **`det II ≡ 0 ⟹ II_ac = II_cc = 0`** — i.e. rule out the "balanced" case where `⟨II_aa,II_cc⟩ = |II_ac|²` holds *without* both `II_ac, II_cc` vanishing. Use the moment identity in §2.1. Suggested attacks:
  - **Asymptotic in θ₂.** As `θ₂→±∞`, `p` concentrates on the max/min‑`X` seeds; expand `det II(θ)` and force conditions order by order. `II_aa≠0` is a robust nonzero, so `det II≡0` is very restrictive.
  - **Developable/rank‑1 characterization.** `det II≡0` ⟺ the surface is developable in the sphere (the extrinsic curvature 2‑form is degenerate); its null direction field integrates to ambient geodesics (asymptotic lines). Classify developable exponential‑family surfaces; the null direction should be forced to be the `X`‑direction (⟹ `II_cc=0`).
  - **Polynomial identity.** `det II(θ)` is a rational function of `{e^{θ·T_i}}`; `≡0` is a polynomial identity whose coefficients (indexed by monomials in the `e^{θ·T_i}`) each vanish — a finite system on `X`. Grind it symbolically for `m=4,5` first (fewer seeds), look for the pattern, then the catalog.

**Success criteria.** Either (i) a proof of 1b (⟹ Task 1 fully closed: exactly 8 surfaces, `[F]`), or (ii) an explicit non‑indicator counterexample with `det II≡0` (which would *revise* the classification — report it loudly, it's a real discovery), or (iii) if it resists, deliver 1a exactly + the strongest partial 1b (e.g. proof under `II_ac=0`, or the asymptotic obstruction) and leave a sharpened `[O]`. Add results to `field_surprisal_classification.py` and a `§Classification` update in the whitepaper.

---

## 5. TASK 2 — Intrinsic forced temperature  *(clean, conceptual)*

**Statement.** Remove or dissolve the one declared input `β`. Two acceptable resolutions:

- **(a) Force a unique β by an intrinsic condition.** Candidates already known to be intrinsically distinguished, each by a *different* principle: `β=1` (MDL/Occam — the code length *is* `log M`, so the universal prior is `2^{−cost}=M^{−1}`; Rissanen), `β=λ=√5` (matching the growth‑rule cost, from the λ=2c paper), `β*≈−0.0768` (max Fisher = max heat capacity, the fluctuation peak). If one of these is forced by a *canonical* condition intrinsic to the surprisal geometry itself (not imported), state it and promote `β` from `[D]` to `[F]`.
- **(b) Reframe β as a coordinate, dissolving the choice.** Argue rigorously that the field's canonical object is not a distribution but the **entire Gibbs curve** (the e‑geodesic `β ↦ p(β)`), of which `β` is merely an affine coordinate. All structural results (affinity, Fisher=Var, the sphere, the 8‑surface classification) are parametric in `β` and `β`‑independent in form; the "free parameter" is then not a modeling degree of freedom but a coordinate, with several intrinsically‑forced *marked points* (the three above). This dissolves the `[D]` tag honestly.

**Success criteria.** A tight, correct `§` (2–3 pp of load‑bearing argument) establishing (a) or (b), with the marked points and their forcing principles enumerated and the ledger's `[D]` entry updated. Prefer whichever is *true* — do not force (a) if the forcing is really an imported convention; (b) is a legitimate and honest outcome.

---

## 6. TASK 3 — Higher‑dimensional classification & iterated ruling  *(open‑ended; needs a sharp sub‑question)*

**Statement.** Extend the two‑statistic classification to `k`‑statistic families `(log M, X₁, …, X_{k−1})`, `k≥3`, over the *same* catalog (distinct from the product geometry of §13, which uses independent catalogs). Central questions:

- **Which `k`‑statistic families have constant sectional curvature 1/4?** Generalize the II conditions: `II_αβ=0 ⟺ T_αT_β ∈ V` with `V=span{1, log M, X₁, …, X_{k−1}}` (now dim `≤k+1`). Constant sectional curvature ⟺ the Gauss‑equation obstruction vanishes on every 2‑plane. Compute the full Riemann tensor of the `k`‑D Hessian metric `∇²A` (extend `numK` to `numRicci`/sectional curvature, validated on known constant‑curvature metrics).
- **Does the ruled/suspension phenomenon iterate?** Conjecture: `k` statistics that are single‑outcome indicators of `k−1` *distinct* seeds (`𝟙_{s₁}, …, 𝟙_{s_{k−1}}`) make the √‑embedding a `k`‑fold spherical **join**, hence constant sectional curvature 1/4 (an iterated ruled surface). Test first on **`(log M, 𝟙_K, 𝟙_{φ⁴})`** — i.e. `(log M, degree, charge parity)`, the two "deep" apexes together — a natural, forced 3‑statistic object. Then `(log M, 𝟙_s, 𝟙_t)` for generic distinct `s,t`.
- **Totally geodesic in higher dim?** Still governed by `V` being a subalgebra ⟺ `(log M, X₁,…)` takes `≤ dim V` distinct values. `log M`'s 6 distinct values still obstruct it for small `k`; find the threshold.

**Success criteria.** At minimum: (i) the iterated‑suspension conjecture stated precisely and **tested** on `(log M, 𝟙_K, 𝟙_{φ⁴})` and a couple of `(log M, 𝟙_s, 𝟙_t)` cases (constant sectional curvature 1/4, `[C]`, with a validated sectional‑curvature routine); (ii) the generalized II conditions written out; (iii) a proof of the iterated‑join → constant‑curvature statement if the two‑statistic suspension proof extends cleanly (it should — the join of a point with a curve generalizes to the join of points with a curve). A full classification is a stretch goal; a sharp partial result on the iterated indicators is the target. Beware combinatorial blow‑up and FD noise; keep the sub‑question concrete.

---

## 7. Deliverables & definition of done

- **Harnesses.** Extend `field_surprisal_classification.py` (Tasks 1, 3) and add a small harness for Task 2 if it produces checkable identities. Every new claim in a fail‑first `ck(...)`. Validate any new numeric‑curvature routine against known metrics before use.
- **Whitepaper.** Fold results into `field_surprisal_geometry_full.tex` (extend, don't restart): update `§Classification` (Task 1), add a `§` for Task 2, add a `§` for higher dimension (Task 3); update the **epistemic ledger** and **Open Problems**; add any new corpus citation. Compile (2× pdflatex), rasterize a page with `pdftoppm` to visually confirm, `present_files`.
- **Honesty.** Tag everything. If a target result is false, say so and report the counterexample — that is a *success*, not a failure. Do not inflate `[C]`/`[O]` to `[F]`. Keep "constant curvature 1/4" and "totally geodesic" strictly distinct.
- **Done =** Task 1 either closed `[F]` or reduced to a sharp residual with 1a proven; Task 2 resolved (forced β, or β dissolved as a curve‑coordinate); Task 3 has the iterated‑suspension conjecture tested on the forced `(log M, 𝟙_K, 𝟙_{φ⁴})` case and either proven or precisely bounded; all in a compiling whitepaper + passing harnesses, ledger and open‑problems updated.

*Do not manufacture a "Task 4." When these three are done (or precisely bounded), the classification chapter is complete — report that and stop.*
