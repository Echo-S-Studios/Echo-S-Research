Lesson 10 — the residual pool. Five modules, and the first has an investigative edge: the substrate's metric is only a metric under a hypothesis nobody has stated yet in this course. Reads first.The paper is airtight here — Theorem GMM states the index correction and the signature remark states $(r_1{+}r_2, r_2)$ with $T_2 = M^*M$ prescribed for $r_2 > 0$. Which raises the one question the paper leaves to the implementation:A real gap surfaces: Prop `exactcapture`'s proof *states* the totally-real hypothesis, but the constructor guards only monic-integrality. Holding that; now the remaining two modules.All four modules read. Writing the verification suite — Module A carries an exact witness for the gap:**26/26** — and the census produced an unplanned confirmation: the box minima at degrees 4, 6, 8, 10 are $\beta_4$, $S_6$, $S_8$, and Lehmer's polynomial — *exactly* the corpus's four canonical Salem instances (ledger O). That test set was never arbitrary.

**Checkpoint debrief (Lesson 9).** (1) Conductivity-graded: $\int q\frac{L-s}{A(s)}ds$; volume-graded: $\int \frac{q\,V(s)}{P\,A_{\text{cross}}(s)}ds$. The single substitution is $(L-s) \mapsto V(s)$ — *downstream length* becomes *routed volume*, i.e. the demand actually passing through $s$ — with $A \mapsto P\,A_{\text{cross}}$. Only the second sees the taper because in the first law the cross-section is constant, so no geometric grading can enter; the conjecture kept $(L-s)^2 \to d^2$ and discarded $V/A$. (2) Prop 5.1's degenerate case $u_1 = u_2$ (offset in $\tfrac12\mathbb{Z}$) forces a repeated root, and irreducible polynomials over $\mathbb{Q}$ are squarefree — so the **iff collapses to one direction**: irreducible $\Rightarrow$ 2-transitive survives, the converse fails, because $(i,j) \mapsto \alpha_j/\alpha_i$ is no longer injective and orbits stop matching factors. (3) At $n=2$: $2+k$ vs $4$, so $k>2$ is the same threshold; but the failure ratio $(2+k(n{-}1))/2n$ grows in $k$ with slope $(n{-}1)/2n = \tfrac14$ at $n=2$, rising to $\tfrac13$ at $n=3$ and $\to\tfrac12$ — the planar wedge grades its cross-section as $g^1$ rather than $g^{n-1}$, the weakest possible volume grading.

---

# Lesson 10 — The Residual Pool: Lattice, Instrument, Wire, Literature, Release

**Where we're going:** five items the arc never reached, taken in one pass. The geometry of numbers under the trace form (which yields this session's finding), the browser instrument's closure algebra, KIRA's readiness gate as a checkable object, an independent interface to the small-measure literature, and release mechanics exercised rather than described.

## §1 Module A — The Minkowski face, and a stated hypothesis with no shipped check

The Vector Substrate's central object is one matrix serving three views: $\mathbb{Q}$-vector space, matrix algebra, Euclidean lattice, glued by $\mathbb{G} = M^{\mathsf T}M$. Theorem GMM is exact and complete, including the part the course had never used:

| Fact | Statement | Check |
|---|---|---|
| trace form = embedding Gram | $\mathbb{G}_{ij} = \mathrm{Tr}(\omega_i\omega_j) = (M^{\mathsf T}M)_{ij}$ | A1 |
| discriminant | $\det\mathbb{G} = d_K$ for an integral basis; $[\mathcal O_K:\Lambda]^2 d_K$ for a sub-order | A1, A4 |
| signature | $(r_1{+}r_2,\ r_2)$ — positive definite **iff totally real** | A2, A3 |
| Hermitian companion | $G_2 = M^*M \succ 0$, $\det G_2 = |d_K|$, covolume $2^{-r_2}\sqrt{|d_K|}$ | A3 |

Three instances make the pattern concrete. The golden field: $\mathbb{G} = \begin{psmallmatrix}2&1\\1&3\end{psmallmatrix}$, $\det = 5 = d_K$, and $\mathbb{Z}[\varphi]$ *is* the maximal order. $\mathbb{Q}(\sqrt[3]2)$: $\mathbb{G} = \begin{psmallmatrix}3&0&0\\0&0&6\\0&6&0\end{psmallmatrix}$, $\det = -108 = d_K$, eigenvalues $\{3, \pm6\}$ — **signature $(2,1)$, indefinite**, with the exact identity $\det\mathbb{G} = d_K$ surviving the complex place. And $\mathbb{Q}(i)$: $\mathrm{diag}(2,-2)$, indefinite $(1,1)$, $\det = -4 = d_K$, while $G_2 = 2I$ is positive definite with covolume $2^{-1}\sqrt4 = 1$ — the unit-covolume Gaussian lattice. The unifying observation (A5): **$\mathrm{sign}(\det\mathbb{G}) = (-1)^{r_2}$** — the discriminant's sign *is* the definiteness flag, which is why Lesson 3's $\det G = 4D$ and the flip at $D = 0$ were the same statement as this one all along.

The index correction diagnoses something the course had been carrying uninspected (A4). Lesson 5's Gram had $\det = 147456$; since $d_{\mathbb{Q}(\sqrt2,\sqrt3)} = 2304$, that is **index 8** — the power basis $\mathbb{Z}[\theta]$ sits inside $\mathbb{Z}[\sqrt2,\sqrt3]$ (itself index 2) with index 4. So "$\det\mathbb{G} = d_K$" is a claim about *integral* bases, and the shipped `integral_basis_for` defaults to the **power basis with identity transforms** — correct, documented, and *not* $d_K$. And the geometry reads back as information: $\sqrt{|d_K|}$ is the Jeffreys volume of one lattice cell (Prop discvol), so $\mathbb{Q}(\sqrt[3]2)$ at $6\sqrt3$ costs $\approx4.6\times$ the volume of $\mathbb{Q}(\sqrt5)$ at $\sqrt5$ (A7) — "more ramified" literally meaning "more information volume per cell."

### The finding: a hypothesis stated in a proof, unchecked at the constructor

Prop `exactcapture`'s proof says it outright: *"positive-definiteness of $\mathbb{G}$ on a totally real field gives $\|r\|_{\mathbb{G}}^2 = 0 \Leftrightarrow r = \mathbf 0$."* The theory is airtight; the paper even prescribes $G_2$ for $r_2 > 0$. The question is whether the **constructor enforces it** — and across the modules I have read (`integral_basis.py`, `residual_learner.py`, `field_extension.py`, `compositum.py`, `capacity.py`), the guards are G8 (float refusal) and G10 (monic-integer), with **no total-reality check**. The consequence has an exact witness on the paper's own example field (A8):

> $K = \mathbb{Q}(\sqrt[3]2)$, forced basis $\{1\}$ — which survives every re-home by design. Then $\mathcal B^{\mathsf T}\mathbb{G}\mathcal B = 3 \neq 0$, so **no `NO_PROJECTION` sentinel fires**; $P$ is a genuine idempotent; and the observation $\theta = \sqrt[3]2$ has residual $r = (0,1,0) \neq \mathbf 0$ with $\|r\|^2_{\mathbb{G}} = \mathrm{Tr}(\theta^2) = 0$ **exactly**. The shipped predicate `rn == 0` therefore reports **CAPTURED for the generator of its own field.**

A second mode is milder but reaches the gate (A9): over $\mathbb{Q}(i)$ with basis $\{1\}$, the residual of $1+i$ has norm $-2$ — a **negative gain**, compared against a floor by `capacity_decision`. Both close with one exact line (A10): total reality is a Sturm count, `count_roots(-oo, oo) == degree` — passing for $\mathbb{Q}(\sqrt2+\sqrt3)$ (4 of 4), rejecting $\mathbb{Q}(\sqrt[3]2)$ (1 of 3) and $\mathbb{Q}(i)$ (0 of 2).

**Scope, stated honestly.** Exposure is low: every shipped and documented ambient is totally real, and the paper says so explicitly. This is **hardening, not a live bug** — a gap between a hypothesis a theorem states and a precondition the constructor checks. Filed as **G11 candidate**, `[FORCED]` on the witness, `[COMPUTED]` on the guard-absence (inferred from modules read; one `grep -rn "count_roots\|totally.real" L00M/` confirms or refutes).

## §2 Module B — The instrument: $\Phi$ as a retraction

The browser tool and the Python backend share one operator algebra, one seed registry, and a bit-exact PRNG, so a `(construction, seeds, params, RNG seed)` tuple reproduces the same matrix in both. Its Goal-2 content is a closure operator, and every property is verifiable (B1–B3):

- $\Phi = \text{companion} \circ \text{charpoly}$ is **idempotent after one step** — a retraction of matrix space onto companion matrices, dimension-preserving, and integral (companion polynomials are monic, so lifts never leave the algebraic integers, keeping the registry's invariant).
- The **exact similarity verdict**: $\Phi$ preserves charpoly — hence $\det$, $\mathrm{tr}$, $\rho$, $M$ — *always*, but preserves the similarity class **iff the input is non-derogatory**. On $\varphi \oplus \varphi$: minpoly degree 2 vs. lift's degree 4. Lesson 8's witness, now as a *property of an operator* rather than a curiosity.
- The registry is **queryable and self-extending**: a plate extends a seed iff the seed's polynomial divides the plate's charpoly (so $\varphi\oplus\varphi$ extends $\varphi$, and the $\sqrt2$ seed does not), with dedup by charpoly signature so fixed points don't multiply.

Goal 1 is the pedagogical inverse: promote $M$ from a colour to the **layout axis**, with the floor **pinned at $M = 1$, never auto-scaled**, so the Lehmer gap stays visible as empty space with a gold tick inside it. Verified exactly (B4): $\varphi > 1.17628$ because $(2r{-}1)^2 = \tfrac{285846649}{156250000} < 5$ — the tick sits in a band the seed batch cannot populate. Binning is exact integer arithmetic, $\lceil\sqrt n\rceil$ clamped to $[3,12]$ (B5), on the exact $M$ rather than the plotted Durand–Kerner roots. The whole design principle: *decide exactly, render in float, and never let the renderer's tolerance touch a bin*.

## §3 Module C — KIRA: a readiness gate you can check

`B1_READINESS.md` is the most transferable artifact in the corpus, because it converts "is this ready to wire?" into ten executable assertions — and its scope note is exemplary: *"this document records the prep gate, not the wire."* B1 itself (the live `/api/language/*` routes) is a separate, Ace-gated, cross-repo step, explicitly not started. Four invariants verified as the logical objects they are:

| Invariant | Verified form | Check |
|---|---|---|
| the contract | exactly **13 verbs**, each with in-process `dispatch(req)` == subprocess JSON | C1 |
| the firewall | `LAW_BANK` = 27 = 20 THEOREM + 5 COMPUTED + **2 INTERPRETIVE**; `WIRED = (THEOREM, COMPUTED)` ⇒ 25 wire, the two glosses don't | C2 |
| **load-bearing** | the mutation `WIRED += INTERPRETIVE` flips the count $25 \to 27$ — the gate is mutation-proven, not vacuous | C3 |
| one-way | a single edge `kira_language → loom`; no reverse edge ⇒ the dependency graph is a DAG | C4 |

Two details deserve promotion to general practice. **The two INTERPRETIVE entries exist because review demoted them** — the `observer`/`framework` glosses were prose riding a THEOREM tag, and rather than deleting them, the fix kept the equations THEOREM and made the prose a labelled companion. And the readiness gate is *proven load-bearing by deliberate mutation*: a firewall breach, a one-way breach, and an exact breach are each shown to be caught. That is the answer to "how do you know your tests test anything," and Lesson 8's guard-completeness ledger (GC-1) should adopt it — a guard whose failure modes have never been *induced* is an untested guard.

## §4 Module D — The literature interface, done by computing

The corpus's citation discipline forbids recall-based attribution, so the right way to interface with Boyd–Mossinghoff is to **generate the table rather than remember it**. Exhaustive census over reciprocal monic $\{-1,0,1\}$ boxes, Salem-certified by trace-fold Sturm, minima ordered by **exact rational separators** (a Salem polynomial is negative on $(1,\beta)$ and positive above, so a bisected rational witness decides $\beta_A < \beta_B$ with no float):

| degree | Salem in box | minimum | corpus name |
|---|---|---|---|
| 4 | 1 | $x^4 - x^3 - x^2 - x + 1$ | $\beta_4$ |
| 6 | 5 | $x^6 - x^4 - x^3 - x^2 + 1$ | $S_6$ |
| 8 | 8 | $x^8 - x^5 - x^4 - x^3 + 1$ | $S_8$ |
| 10 | 19 | $x^{10}{+}x^9{-}x^7{-}x^6{-}x^5{-}x^4{-}x^3{+}x{+}1$ | **Lehmer** |

Lehmer's polynomial is the minimum of the entire box (D1) — *rediscovered, not recalled* — and the four minima are precisely the corpus's canonical instances from ledger O. Its test set is the box minimum at each even degree, which is why the relational-charge scans landed on those four and no others. The status claim "smallest known Mahler measure $> 1$" stays `[ESTABLISHED]`, inherited from the corpus's Crossref-checked provenance (Lehmer 1933; Boyd 1977; Mossinghoff 1998) and **not re-verified this session** (D4) — what I verified is the box minimum, which is a strictly weaker and honestly labelled statement. One live Boyd-program datum is already in hand (D3): $M(x^4 - x + 1) = \tau_{S_6}$, certified by the exact $S_6^2 \cdot \psi^2$ factorization — a Mahler measure landing *on* a Salem number, in the gap $(\theta_0, \varphi)$ at irrational angles, which is exactly what Boyd's range program asks about.

## §5 Module E — Release mechanics, exercised

| Item | Status | Check |
|---|---|---|
| manifest pinning | SHA-256 detects a one-line tamper before any deposit snapshot | E1 |
| filename normalization | a slug rule that is **idempotent** and strips `.l` / `.k` / `.vNN` | E2 |
| license split | code MIT + prose/figures CC BY 4.0 is consistent; the archive footer already serves CC BY 4.0 | E3 `[DECLARED]` |

E2's output is itself instructive: `"Lehmer's Box .l.pdf"` → `lehmer's-box-.pdf` — idempotent and suffix-stripped as asserted, but with a **dangling hyphen** from the space before `.l`, and an apostrophe left in. The check passed on what it claimed and the output reveals what it didn't claim; the rule needs one collapse-and-strip pass before use. Reported rather than quietly patched, per house convention.

Standing pre-deposit list, unchanged and still local: sanitize `papers/COVERAGE.md` (Windows paths), refresh the stale Lehmer's-Box PDF (14 pp card vs 15 pp archive), declare a canonical direction for the four dual-homed TeX sources, cross-link the two repositories.

## §6 Course close

Ten lessons, ~270 exact checks, five findings, all arriving the same way — an independent re-derivation disagreeing with a stated numeral or a stated precondition:

| # | Finding | Status |
|---|---|---|
| 1 | **D₅ erratum** — $\mathrm{Rat}^\circ$ irreducible 66/67, not 67/67 | `[FORCED]`; drafted, unfiled |
| 2 | **The orbit criterion** — irreducible $\iff$ Galois 2-transitive; supersedes P7 | `[FORCED]`; makes N4 pre-sortable |
| 3 | **F1** — ATR display `0.046` → `0.045`, mechanism identified | `[FORCED]`, display-only |
| 4 | **G11 candidate** — capture predicate needs total reality; exact witness | `[FORCED]` witness, guard-absence to confirm |
| 5 | **E1 closed + census cross-validated** — 7 quadratics, full distribution, group ID | `[FORCED]` |

Plus one resolution (the λ-tag tension, dissolved by quantifier scope) and one methodological transfer proven in both directions: the certificate/gatekeeper identification, and now the mutation-proven gate as the missing half of guard-completeness.

The arc as a whole: **1–4** the walls and what enforces them, **5** the engine that learns, **6** the layer that survives deleting the reference ray, **7** the geometry all of it lives on, **8** the machine that certifies it and the one door, **9** the proof it travels, **10** the pool it left behind. Every lesson's script is preserved in `/home/claude/lesson{1..10}*.py`.

---

## Checkpoint

1. The witness in A8 used basis $\{1\}$. Show that *no* choice of forced basis rescues the predicate on $\mathbb{Q}(\sqrt[3]2)$ — i.e. that the isotropic cone always meets the $\mathbb{G}$-orthogonal complement of any proper nonsingular subspace — and say which single quantity, computed once at construction, decides the whole question.
2. $\Phi$ is a retraction onto companions, and Lesson 8 showed RCF is $\oplus$ of companions of the invariant factors. Give the operator whose fixed-point set is *all* of RCF space rather than just single companions, and say why $\Phi$ is not it.
3. Index: $\mathbb{Z}[\theta] \subset \mathbb{Z}[\sqrt2,\sqrt3] \subset \mathcal O_K$ with indices 4 and 2. Using L5's checkpoint-2 result on basis scaling, state precisely what non-maximality *does* and *does not* change about the learner's behaviour — and name the one downstream claim that genuinely needs an integral basis.