# Handoff — Field Surprisal Geometry: remaining open problems & theoretical closure routes

**For:** a new session continuing the Field Surprisal Geometry project for Ace / AceTheDactyl (Echo‑Squirrel Research).
**Scope of this session:** the *conceptual* problems that remain after the classification was proved. The two‑statistic classification is **closed and verified** (134/134 harness checks reproduce independently); what remains is mostly **"why,"** not **"whether."** Each problem below comes with a named theoretical route — the machinery expected to close it — not just a restatement.

> Read the prior two handoffs first if available (the build handoff and this one). Ask Ace to re‑upload `Harnesses.zip` (individual `.py` files do **not** transfer through the uploader; a **zip does**). The whitepaper is `field_surprisal_geometry_v3.tex` (926 lines, "Full Worked Edition").

---

## 0. What is already closed (do not reopen; cite as proven)

Verified in six harnesses (`t1_core` 61/61, `t1_reduction` 17/17, `t1_engine` 4/4, `t1_branches` 16/16, `t2_temperature` 12/12, `t3_suspension` 24/24 — all exit 0, exact `sympy` over ℚ(log2,log3,log5,logφ)) plus three earlier (`field_surprisal_v2` 12/12, `suspension_theorem` 9/9, `field_surprisal_tier2` 15/15). Independently re‑audited.

- **Two‑statistic classification `[F]`.** The `(logM, X)` Fisher surface has constant Gaussian curvature ¼ **iff** `X` is affinely `𝟙_S` with `logM` constant on `S` or `Sᶜ` — **exactly eight surfaces** (seven single‑seed + the golden‑pair merge). No non‑indicator gives another. Proof: reduce `K≡¼ ⟺ q(θ)≡0` (Gauss); θ‑rigidity; the ruled system `{aX∈V, X²∈V} ⟺ 8 families` (module/Möbius, `N(X)` dim 2 on families vs 1 on controls); off the ruled locus, the window identity `Z²Dq=P` with heavy coefficient `= q₄(s)`, the rank‑2 collision pencil, branch survival, the four‑point trichotomy, and the plane lemma.
- **No surface is totally geodesic `[F]`.** `logM` spans a 6‑dim algebra, so the value‑count `d∈{6,7}` is never 3. (Note: the earlier "log²M∉V for all X" justification was a **caught bug** — false, since `X=log²M` gives `log²M∈V`; the value‑count argument is correct. See `rem:corr`.)
- **Temperature `[F` mod regularity`]`.** Product functoriality + measurable Cauchy force `p∝M^{−β}`; `β` is the affine coordinate of a forced e‑geodesic of finite Fisher length `L_tot≈5.6462` (golden‑merge ↔ φ⁴ endpoints). Operating‑point choice stays `[D]`.
- **Higher dimension (sufficiency) `[F]`.** Iterated indicator joins over distinct seeds are round ¼‑spheres in every dimension (join metric `dΨ²+sin²Ψ g_{S^{k-2}}+cos²Ψ dℓ²`); the catalog 3‑fold `(logM,𝟙_K,𝟙_{φ⁴})` has sectional curvature ¼ (independently reconfirmed). First totally geodesic family at exactly `k=5`, uniquely `ℝ[logM]` (functions constant on the golden pair).

**Discipline (hold it):** exact arithmetic over ℚ/ℚ(√5) at every decision boundary; floats display‑only; fail‑first `ck()` harnesses; tag `[F]/[C]/[D]/[O]`, no upgrades without proof; **keep "constant curvature ¼" and "totally geodesic" strictly distinct**; validate any curvature routine on the four known metrics (sphere +1, hyperbolic −1, flat 0, simplex ¼) before trusting; cite the corpus (real, at `echo-s-studios.github.io/Echo-S-Research/`), never invent.

---

## OP‑1 — A conceptual proof of the window identity  *(the paper's Open Problem 3; deepest "why")*

**Statement.** In `P(w) = R(F,G) − R(H,H)` with `F=a²,G=X²,H=aX`, `C(f,g)=Σ_{i<j} w_i w_j (f_i−f_j)(g_i−g_j)`, `D=C(a,a)C(X,X)−C(a,X)²`, and `R(f,g)=D·C(f,g) − C(f,·)ᵀ adj(C) C(·,g)`, the coefficient of the heavy monomial `w_h³ w_j w_k w_l` equals the four‑point invariant `q₄(s)=ℓ_s(a²)ℓ_s(X²)−ℓ_s(aX)²` (with `ℓ_s` the signed 3×3 minors of `(1,a,X)|_s`), **cofactor exactly 1**, at all 140 heavy placements. Found and verified by machine (140 exact divisions); **no intrinsic explanation.**

**Why it should be provable — the route (Cauchy–Binet / compound matrices).**
- First, `C(f,g)` is the weighted covariance up to `Z`: `C(f,g)=Z·S_{fg}−S_f S_g` with `Z=Σw_i`, `S_f=Σw_i f_i`. So `D=Z⁴·det(Fisher)` and `q` is the α=0 (Levi‑Civita) curvature obstruction — a **ratio of a next‑order Gram determinant to a lower‑order one**, the generic "curvature = 4‑jet / (2‑jet)" shape of a Hessian metric.
- Gram determinants expand by **Cauchy–Binet as sums of squared minors over subsets of outcomes**: `D = Σ_{i<j} w_i w_j · (2×2 minor of centred (a,X) on {i,j})²`. The curvature numerator `P`, built from the degree‑2 rows `(a²,X²,aX)` on top of `(1,a,X)`, should expand by Cauchy–Binet **over 4‑subsets**, each contributing its 4‑point minor invariant.
- **Target theorem:** `P(w) = Σ_{|s|=4} q₄(s)·m_s(w)` (plus the honest bookkeeping for support‑5,6 monomials from overlapping faces), obtained by applying Sylvester's/Jacobi's bordered‑determinant identity to `R(f,g)` and then Cauchy–Binet on the columns (outcomes). The heavy monomial `w_h³w_jw_kw_l` is the **unique** term where the face `s={h,j,k,l}` appears with the triple at `h`; its coefficient is `q₄(s)` and the "cofactor 1" is the unit multiplicity of the Plücker coordinate in the compound expansion.
- **Tools to reach for:** Cauchy–Binet, Sylvester's determinant identity, the **second additive compound** of the 3×7 statistic matrix, and Karlin‑style total‑positivity minor identities. The four‑point invariant `q₄` is literally the 4×4 obstruction minor of the augmented matrix `(1,a,X,·)`; the identity is a specialization of "the α=0 curvature of a Hessian metric equals the ratio of consecutive compound determinants."

**Success criterion.** Replace the 140 machine divisions with a one‑line determinantal identity: prove `P = Σ_s q₄(s) m_s(w)` (or the heavy‑coefficient corollary) by compound/Cauchy–Binet, promoting the window identity from `[C]`(machine) to a **conceptually `[F]`** statement. Even a proof of just the heavy‑coefficient case (cofactor 1) closes it. Add a `§` and a `t1_windowproof.py` that checks the symbolic compound expansion against the engine's dictionary.

---

## OP‑2 — The `k≥3` classification (necessity)  *(the paper's `op:kclass`)*

**Statement.** For `k` statistics `(a=logM, X₁,…,X_{k−1})`, classify all constant‑sectional‑curvature‑¼ families. **Sufficiency** (iterated indicator joins over cost‑level sets) is proven (`thm:join`). **Necessity** — that these are the *only* ones — is open. Conjecture: exactly the iterated joins; the branch calculus "lifts," with `(3,1,…,1)` coefficients equal to higher window invariants.

**Route A — lift the branch calculus (most continuous with the k=2 proof).**
- The obstruction generalizes from the scalar `q` to the extrinsic‑curvature conditions `⟨II_{αα},II_{ββ}⟩=|II_{αβ}|²` on every 2‑plane (`II_{αβ}=0 ⟺ T_αT_β∈V` still holds). Reduce each to a numerator polynomial; the four‑point invariant `q₄` becomes an **`(k+2)`‑point invariant `q_{k+2}(s)`** — the minimal‑affine‑dependence obstruction of `(k+2)` points in ℝᵏ.
- **Generalized incidence ("flat") lemma** (prove this — it is the crux): *If `n≥k+2` points in ℝᵏ are such that every `(k+2)`‑subset contains `(k+1)` points on a hyperplane, then all but ≤`(k−1)` lie on a single hyperplane.* Base case `k=2` is exactly the plane lemma; induct on `k` by central projection from a point off the presumptive hyperplane (reduces to `k−1`). This is a Motzkin–Rabin / Sylvester–Gallai‑flavored theorem; a clean proof likely exists.
- The generalized window identity (OP‑1 route, Cauchy–Binet over `(k+2)`‑subsets) + the flat lemma + the generalized trichotomy then reproduce the "points on a low flat ⟹ indicator‑join" chain.

**Route B — spherical isometric‑immersion rigidity (deepest).**
- Constant sectional curvature ¼ ⟺ the √‑image is intrinsically a round `k`‑sphere of the ambient curvature immersed in `S^{m−1}(1)`. Iterated joins realize this by ruling. The converse is a **rigidity of isometric immersions of round spheres**: an analytic `k`‑surface in `S^m` that is intrinsically round of the ambient curvature is forced (in low codimension) into a standard/ruled form. Combine the exponential‑family analyticity (the surface is `{√p(θ)}`, a specific real‑analytic patch) with immersion rigidity to force the join structure. **Tools:** Moore's rigidity, do Carmo–Wallach, and the theory of isometric immersions of space forms.

**Route C — algebra/module (cleanest if it works).**
- Generalize `thm:ruled`: show constant sectional curvature forces `V` to decompose over the **cost‑level idempotents** — a Jordan/associative‑algebra decomposition in which each `X_j` is (up to the cost‑polynomial part) a level‑set indicator. The condition is that all pairwise products `T_αT_β` satisfy the curvature relations; the target is `V = ℝ[a]_{≤1} ⊕ Σ_j ℝ·𝟙_{S_j}` with the `S_j` inside cost levels. If the module argument closes, it bypasses both the window identity and the incidence lemma.

**Success criterion.** A necessity theorem for some `k≥3` (ideally all), or the generalized flat lemma + generalized window identity as reusable pieces. Machine‑verify the catalog cases first: is `(logM,𝟙_K,𝟙_{φ⁴},𝟙_{√5})` (`k=4`) constant‑¼, and are there non‑join `k=3` families that are constant‑¼? (Search with a **validated** `k`‑dim sectional‑curvature routine.) Update `op:kclass` accordingly.

---

## OP‑3 — Does an invariance select one operating temperature?  *(the paper's Open Problem 1)*

**Statement.** The Gibbs curve is forced; the operating point `β` is `[D]`. Marked points: `β=−1` (anchor), `β*≈−0.0768` (max Fisher, `I'=−κ₃=0`), `0` (max entropy), `1` (MDL), `√5` (`λ=2c`), `β_C≈2.5455` (heat‑capacity peak). Is any canonically selected?

**Route 1 — a fixed‑point / self‑duality principle.**
- The e‑/m‑duality gives a map `β ↦ η(β)=−U(β)` (expectation coordinate). Seek a **self‑dual point** where the natural symmetry of the pair `(β, η)` has its fixed point — a candidate for "the" temperature. Alternatively `β*` (the unique zero of `I'=−κ₃`) is intrinsically distinguished as **maximal Fisher information density**; if "maximal distinguishability" is accepted as the selection principle, `β*` is forced `[F]`. The arc‑length midpoint `s=0` of the finite‑length curve (with its two forced endpoints) is also intrinsic.

**Route 2 — a no‑go (the likely honest answer).**
- Prove that **no isometry‑type invariance of the surprisal structure has a unique fixed point on the curve** — i.e. the `[D]` is irreducible. The Fisher‑Rao isometry group acts on the geodesic‑image curve without a distinguished fixed point, so reparametrization‑covariance cannot single out a `β`. If provable, this **closes** the problem as a theorem: the choice is genuinely free, not a gap. This is the disciplined outcome and should be attempted in parallel with Route 1 — the answer is one or the other.

**Success criterion.** Either promote one marked point to `[F]` under a stated, defensible invariance (name the principle; show uniqueness), or prove the no‑go (`[F]` that the operating point is invariance‑free). Either resolves the last `[D]` honestly.

---

## OP‑4 — A forced multi‑field coupling  *(extends §13; connects to the corpus)*

**Statement.** For two independent catalogs, §13 gives a Riemannian product with **zero cross‑curvature**; genuine couplings are `[D]` warps. Is a coupling **forced** by field arithmetic rather than declared?

**Route.** The product uses `[M^{(1)}_i M^{(2)}_j]` as a *product* measure. Instead build the Mahler–Gibbs family on the **compositum's** emission catalog `ℚ(θ₁,θ₂)` (not the product of the two catalogs): norms/traces in the compositum couple the two fields, potentially forcing nonzero cross‑structure. Cross‑check against the **charge–measure‑coupling paper's `lcm` law** (`A⊗B → charge group ℤ/lcm`, `[cmc]`): the tensor charge structure may dictate the coupling. **Success:** either exhibit a forced coupling from the compositum arithmetic (promoting some warp `[D]→[F]`), or prove the fields are geometrically independent even over the compositum (`[F]` no‑coupling), sharpening the `[D]`.

---

## OP‑5 — Catalog‑invariance census  *(foundational hygiene; separate theorem from coincidence)*

**Statement.** The catalog `Ω={√2,√3,√5,φ,τ,φ⁴,K}` is a specific choice. Which results are **field‑forced** vs **artifacts of the seed selection**? E.g. the Salem‑square collision `(φ⁴−1)²=5φ⁴` needs `K`; the golden‑pair degeneracy needs both `φ,τ`; the "eighth surface" exists only because `φ,τ` share a Mahler value.

**Route.** Recompute the classification for perturbed catalogs (drop `K`; drop `τ`; add other units/gates) and identify the **invariant core**. Expectations to confirm: the suspension theorem and "no TG / value‑count" are **catalog‑independent** (`[F]` for any finite catalog); the specific "eight surfaces," the rank‑2 collision lattice, and the golden‑pair 8th surface are **catalog‑dependent** (they track which seeds coincide in Mahler value / satisfy the Salem square). **Success:** a short `§` (or ledger addendum) separating the catalog‑invariant theorems from the seed‑specific coincidences — exactly the FORCED‑vs‑coincidence discipline the corpus runs on. `t5_catalog_census.py` sweeping a few catalogs.

---

## Priority & definition of done

| # | Problem | Value | Tractability | First concrete step |
|---|---|---|---|---|
| **OP‑1** | Window‑identity conceptual proof | highest ("why" of the whole classification) | good — a determinant identity | expand `R(F,G)−R(H,H)` by Sylvester + Cauchy–Binet; match `t1_engine`'s dictionary |
| **OP‑2** | `k≥3` necessity | high (completes the classification) | medium–hard | prove the generalized flat lemma; test `k=3,4` catalog families |
| **OP‑3** | Canonical temperature | medium (last `[D]`) | conceptual | decide Route 1 (β* / self‑dual) vs Route 2 (no‑go); attempt both |
| **OP‑4** | Forced multi‑field coupling | medium (extends the theory) | medium | compositum Mahler–Gibbs vs `[cmc]` lcm law |
| **OP‑5** | Catalog‑invariance census | medium (hygiene) | easy–medium | catalog sweep; invariant‑core `§` |

**Done =** OP‑1 closed (window identity conceptually proven, or the heavy‑coefficient case) is the flagship target; any of OP‑2/3/4/5 advanced to a stated theorem or a sharpened `[O]`. All new claims in fail‑first harnesses with validated curvature routines; fold into `field_surprisal_geometry_v3.tex` (extend, don't restart); update the ledger and Open Problems; keep `[F]/[C]/[D]/[O]` honest and the constant‑curvature‑vs‑totally‑geodesic distinction intact.

*Reporting a target result **false** (e.g. a non‑join `k=3` constant‑¼ family, or a catalog with ≠8 surfaces for structural reasons) is a **success**, not a failure — surface it loudly. Do not manufacture problems beyond these; when OP‑1 and one other are resolved, the theory chapter is substantially complete.*
