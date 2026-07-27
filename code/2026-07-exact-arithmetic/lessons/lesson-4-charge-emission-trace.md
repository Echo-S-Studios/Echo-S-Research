Lesson 4(a) it is, built to hand Lesson 5(b) its objects pre-verified. Three sources to read first: `capacity.py` itself (pipelines, project knowledge), and the two archive papers that carry the parity floor and the dichotomy — both fetchable as sealed PDFs from the archive site.Falling back to the landing page (user-provided URL, always valid) to re-register the PDF links:25/25. (One housekeeping note in the spirit of the discipline: the script contains a harmless dead placeholder line before the embeddings matrix is actually built — the executed rows are the four genuine embeddings, as the diag output confirms. Flagged, no effect.) Now the lesson.

**Checkpoint debrief (Lesson 3).** (1) $c(-\tfrac14) = \sqrt{0}/(2\cdot(-\tfrac14)) = 0$, so $\lambda = 0$: the *gate* reading is that the pinned value is the collapsed spectral gap — information exchanges at rate zero; the *flip* reading is that $0$ is the crossing itself. One number carrying both events, exactly as the co-location table said. (2) $M(x^2-1) = 1$ — roots $\pm1$ on the circle, Kronecker's trivial locus. The $D=1$ rung is the floor itself: $\mathbb{Q}(\sqrt1) = \mathbb{Q}$, no extension, gap $\pm1$ made of units. The ladder's first *contentful* rung is $D = 2$ because content begins where $M > 1$ begins. (3) $\mathrm{ad}_R(I) = R - R = 0$ and $\mathrm{ad}_R(R) = R^2 - R^2 = 0$ for any square $R$, so the kernel contains $\mathrm{span}\{I, R\}$, forcing $x^2 \mid \mathrm{charpoly}(\mathrm{ad}_R)$ always. That argument is theorem-grade — valid for all $R$ symbolically — strictly stronger than E4b's two instances. The asymmetry rule cuts one way: instances never promote, but this isn't an instance.

---

# Lesson 4 — The Capacity Gate, the Parity-Graded Floor, and How the Dichotomy Assembles

**Where we're going:** the gate that *decides* growth (`capacity.py`, exactly as shipped), then the generalization of the charge from $\mathbb{Z}/4$ to all finite cyclic groups, then the principal coupling theorem — the parity-graded floor — then how the no-Salem dichotomy assembles from it, ending at the one door through which everything excluded re-enters. Every claim below carries the *paper's own tag*, because this paper's erratum trail (v1→v4) is itself the best tagging lesson in the corpus.

## §1 The gate in code

`capacity_decision` in `residual-return-verification/L00M/training/capacity.py` is **two layers**, and the layering *is* the exactness discipline:

| Layer | Decides | Over | Status |
|---|---|---|---|
| Northcott admissibility (always on) | REJECT vs. proceed: $\deg \leq D_{\max}$ **and** coeff-height $\leq H_{\max}$ | $\mathbb{Z}$ only | rests on theorems: Northcott 1949 finiteness, Landau $M(p) \leq \|p\|_2$ `[ESTABLISHED]` |
| Information threshold (opt-in, default **off** byte-for-byte) | STOP vs. GROW: exact-`Fraction` gain $\|r\|_G^2$ vs. the Lesson-3 floor $2c\log\mu$ | gain over $\mathbb{Q}$; cost by certified interval | the $\lambda = 2c$ machinery |

The G8 rule in the source: *the float Mahler is never consulted* — admissibility is degree + integer coefficient height, with Landau's inequality making coeff-height $\leq H_{\max}$ a **sufficient certified** bound for $M \leq H_{\max}$ (verified on instances: $\varphi \leq \sqrt3$ decided via $\varphi^2 \leq 3$; $M(L) < 2 \leq 3 = \|L\|_2$, A2). The design doc is explicit that this is a *conservative over-approximation* — it can reject a high-coefficient, low-$M$ seed — a named, accepted bias. My replica of the decision table reproduces the shipped demo's ACT 4 (A3): the $2\sqrt6$ seed $x^2 - 24$ with coeff-height $24 > H_{\max} = 10$ → **REJECT**, pure integer comparison; a real defect under budget → **GROW**; zero defect → **STOP**. Three verdicts, all decided over $\mathbb{Z}/\mathbb{Q}$.

Now the branch that made me stop when reading the source: the **reciprocity branch** `[RULING 1]`. Non-reciprocal candidates get the $\mu_S$ floor — Smyth's theorem, positive, unconditional. Reciprocal candidates get **Dobrowolski — which is vacuous at degree 2** (A4: $\log\log 2 < 0$, the bound drops below 1, floor $\to 0$), and the Lehmer value is opt-in heuristic only, never default. Read that again: *Lesson 1's open problem is compiled into the control flow.* The code cannot cite a reciprocal floor because none is proven; the branch structure is the epistemic status. One staleness observation for your ledger: `A3_DESIGN.md`'s honest ledger still lists F4 ($G \propto$ Fisher) as OPEN, but Lesson 3's trace-form Lemma `[FORCED]` resolved it (on trace-zero, $1 \in \mathrm{col}(B)$) and `CHANGES.md` records O1/O4 closed — the design doc predates the closure. Doc-sync candidate.

Also planted for Lesson 5: `effective_degree` lets compositum growth reuse the *same* gate whether the extension is disjoint or not — the budget check is disjointness-independent. Hold that.

## §2 The charge, generalized

Lesson 2's $\mathbb{Z}/4$ charge is the $n = 4$ instance of the paper's Character II: an object is **charge-admissible** if every conjugate satisfies $\alpha^n \in \mathbb{R}_{>0}$ for some $n$ — all arguments on the lattice $(2\pi/n)\mathbb{Z}$ — and the charge group is $\mathbb{Z}/n$ for the least such $n$. Structure, all `[forced]`:

- **Cyclicity** (Thm 3.1): charge groups are exactly the finite cyclic groups — finite subgroups of $\mathbb{Q}/\mathbb{Z}$ are cyclic; $\mathbb{Z}/p \times \mathbb{Z}/p$ is unrealizable. The phase character is a rank-one invariant.
- **Realizability** (Thm 3.2): $x^n - 2$ realizes $\mathbb{Z}/n$ with all charges and $M = 2$ exactly (B1).
- **CRT via Adams** (Thm 3.3): $\psi^{n/p^e}$ projects onto the $p$-primary component — your Lesson-2 checkpoint ($\psi^3$ as the unit $3 \bmod 4$) was the first taste of $\psi$ acting as multiplication on charges.
- **The lcm law** (Thm 4.3/4.4): $\otimes$ composes charge groups to $\mathbb{Z}/\mathrm{lcm}$; composition is lossless **iff coprime**. Engine instance replicated exactly (B8): $C(x^3{-}2) \otimes C(x^4{-}2)$ has charpoly $x^{12} - 128$, charge $\mathbb{Z}/12$, $M = 128 = 2^7$ — an off-circle case where factored happens to equal tropical, which is precisely the class of case that let the v1–v3 tensor-law error survive testing.

## §3 The parity-graded floor — and the erratum that teaches the tags

Define $\mu(n) = \inf\{M(O) : O$ charge-admissible, group $\mathbb{Z}/n$, $M > 1\}$. **Theorem 5.1**: $\mu(n) = \varphi$ for even $n$, $2$ for odd $n$. But the paper's v1→v2 erratum is the real content — this is your remembered "correct parity-graded floor" in its native habitat — and the corrected strength split must be preserved verbatim:

| Statement | Tag (paper's own) |
|---|---|
| $\mu(\text{even}) \leq \varphi$ — the witness attains it | `[forced]` (constructive) |
| $\mu(\text{even}) \geq \varphi$ | `[plausible]` (inherits the universal bound's status) |
| odd forced lower bound: $M \geq \mu_S$ (non-reciprocal, Smyth) | `[forced]` |
| $\mu(\text{odd}) = 2$ | `[forced]` for $\mathbb{Z}/3$; `[computed]` in general |
| universal $(1,\varphi)$ gap | `[forced]` for quadratics only; `[plausible]` in general |

v1 read the universal floor as forced and thereby implied a forced $\varphi$ for odd charge — "that tag was too strong." The forced odd floor is $\mu_S$. This is the asymmetry rule enforced against the corpus's own headline.

**The even witness** is a gift (B2): $q_k(x) = x^{2k} + x^k - 1$ — substitute $y = x^k$ and you get $y^2 + y - 1$: **Lesson 3's golden gate polynomial**, verbatim. The $y = -\varphi$ branch contributes $k$ roots of modulus $\varphi^{1/k}$ at *odd* multiples of $\pi/k$, product exactly $\varphi$; the $y = 1/\varphi$ branch sits inside (B2c). So $M(q_k) = \varphi$ with charge group $\mathbb{Z}/2k$, all charges. At $k = 2$: $x^4 + x^2 - 1$, measure-bearing pair $\pm i\sqrt\varphi$ on the imaginary axis, full $\mathbb{Z}/4$ charge, $M = \varphi$ exactly (B2b) — the even-sector floor of Lesson 2, now *attained*.

**Why parity?** One arithmetic fact (Lemma 5.2, `[forced]`, B3): the golden conjugate $\varphi' = -1/\varphi < 0$ sits at argument $\pi$, and $\pi \in (2\pi/n)\mathbb{Z}$ **iff $n$ is even**. The floor-attaining structure must host $\varphi'$ on the lattice; odd lattices have no $\pi$-ray. The entire dichotomy is the single parity bit *"does the charge lattice contain $\pi$."* And the odd side's reciprocal obstruction is an old friend (Lemma 2.6, B4): a real reciprocal unit pair with integer trace has $r + r^{-1} \geq 3$, minimum at $r = \varphi^2$ — attained by $x^2 - 3x + 1$. That is E1's missing polynomial making its *fourth* appearance (L1 checkpoint → L2 $\psi^2R$ → here as the extremal reciprocal unit).

## §4 How the dichotomy assembles

**$\mathbb{Z}/3$, closed elementarily** `[forced]` (B5): a lattice cubic factors over $\mathbb{R}$ as $(x-\alpha)(x^2 + \rho x + \rho^2)$, and the coefficient identity $c_1 = \rho\, c_2$ holds *identically* (verified symbolically). So $c_2 \neq 0 \Rightarrow \rho = c_1/c_2 \in \mathbb{Q} \Rightarrow$ reducible; irreducible forces $c_2 = c_1 = 0$, i.e. $x^3 - m$, and $\mu(3) = 2$ at $x^3 - 2$. The closure works because $2\cos120° = -1 \in \mathbb{Z}$. Bonus exact exclusion with *no angle ever computed* (B5b): $x^3 - x - 1$ has $c_2 = 0$ but $c_1 = -1 \neq 0$ — off the lattice, by pure coefficient algebra.

**The keystone that guards the gap** (B6): $\beta_4 = x^4 - x^3 - x^2 - x + 1$, the minimal degree-4 Salem number. Reciprocal, irreducible, trace-fold $t^2 - t - 3$ with Sturm counts one-in/one-out — Salem structure derived by Lesson 1's own method — and $M(\beta_4) \in (1.7, 1.8) \subset (\varphi, 2)$ certified. It sits *inside* the odd gap, and it is charge-inadmissible (irrational conjugate angles, Lemma 6.1). **Theorem 6.2**: the odd gap $(\varphi, 2)$ *is* the no-Salem closure — the measures that would fill it are realized only by the objects the charge excludes. `[forced]` for $\mathbb{Z}/3$; `[computed]/[plausible]` for all degrees (17 screen candidates, every one an irrational-angle Salem/Pisot on exact reclassification).

**$\mathbb{Z}/5$ — where elementary dies** (§7, sourced from the companion z5 note). Pentagon cosines (Prop 7.1): $2\cos72° = \varphi - 1$ and $2\cos144° = -\varphi$ are the two roots of — again — $x^2 + x - 1$, Galois-conjugate under $\sqrt5 \mapsto -\sqrt5$. A $72°$ pair *drags in* its $144°$ partner so the $\sqrt5$-parts cancel; that forced coupling is what $\mathbb{Z}/3$ never needs and what kills the coefficient trick. What survives, each verified:

- **Degree-four sector closed at $\varphi^4$** `[forced]` (B7): the Galois-coupled form bottoms out at $(k,m) = (3,1)$, and I verified the minimizer *by exact cancellation*: $(x^2 - \varphi x + \varphi^4)(x^2 + \varphi^{-1}x + \varphi^{-4}) = x^4 - x^3 + 6x^2 + 4x + 1$ on the nose, $M = \varphi^4$.
- **Pure powers** `[forced]`: $x^5 - m$ realizes $M = m$, non-reciprocal, so $\mu(5) = 2$ at $x^5 - 2$ — held at `[computed]` equality.
- **The residual, named** (Prop 7.6/Rem 7.7): non-reciprocal charge-$\mathbb{Z}/5$ objects in $[\mu_S, 2)$ — empty over the window, unpromotable. Realification is circular ($\psi^5$ gives only $2^{1/5}$), Smyth stops at $\mu_S$, classification runs out at degree $\geq 5$. The paper names it: **Lehmer's problem restricted to the pentagon lattice — Schur–Siegel–Smyth trace-problem territory** — matching your standing framing of the dichotomy as an SSS-class trace problem. And it cites the asymmetry rule *in its own prose* as the reason the four forced sub-results plus a clean window are "the honest maximum."

## §5 The commutator door

Every result above governs the **abelian** semiring $(\oplus, \otimes, \psi^n)$. Proposition 9.1 `[forced]` marks the exact boundary: by Albert–Muckenhoupt + Laffey–Reams `[ESTABLISHED]`, every trace-zero integer matrix is a commutator — and the door is wide enough for everything. Verified witness (C1): $L(x)(x-1) = x^{11} - x^9 - x^8 + x^3 + x^2 - 1$, whose $x^{10}$ coefficient is $0$, so its companion is trace-zero, hence a commutator — with $M = M(L) \in (1, \varphi)$ certified (C1b) and charge group $\perp$ (C1c, via Lesson 1's $\gcd(L, x^4-1) = 1$). **Below the floor and off every lattice simultaneously.** Consequences: the floor, gap, and charge structure are operation-relative; no single theorem governs every free commutator; the no-Salem closure excludes Salems *by the restriction to the abelian charge-admissible semiring*, not by the ambient matrix algebra. The paper's closing line is the program's thesis: *confinement is exactly the choice to stay abelian.* And your pipeline's `INVALID_CLOSURE` guard firing on a Lehmer carrier is this proposition as executable code — the tripwire at the door, working as designed — with "guard-completeness as the achievable residue of the free-commutator front" as the standing open item.

## §6 The bridge — Lesson 5's objects, pre-verified

The gate is static; Lesson 5 is the loop that feeds it. Every object it needs is now on the table with a green check:

1. **The residual** (D1): captured basis $\{1\}$ in $\mathbb{Q}(\varphi)$, $G = \binom{2\ 1}{1\ 3}$; the $G$-orthogonal projection of $\varphi$ gives $r = (-\tfrac12, 1)$ with $\|r\|_G^2 = \tfrac52$ — an exact `Fraction`, exactly the *gain* the gate consumes. Gate at $c = 1$: $\tfrac52 > 2\log\varphi$ certified via $\varphi^2 < e$ → **GROW**.
2. **Capture $\iff r = 0$** (D2): adjoin $\varphi$; then $\varphi^2 = 1 + \varphi$ lies in the span and its residual is *exactly* zero. Capture is an algebraic identity, not a tolerance.
3. **Field growth through the Kronecker Gram** (D3): $G_{\mathbb{Q}(\sqrt2,\sqrt3)} = G_{\mathbb{Q}(\sqrt2)} \otimes G_{\mathbb{Q}(\sqrt3)}$ — verified as an exact $4\times4$ identity via the four embeddings, $\mathrm{diag}(4,12,8,24)$, with $\det = (\det G_A)^2(\det G_B)^2 = 9216$ (D3b): discriminants compose, Lesson 2's $\otimes$ acting on the *metric* itself.
4. **Disjointness-independence**: `effective_degree` means the gate's budget check is identical whether the compositum is the full Kronecker degree or smaller — the decision decouples from the construction.

So the loop Lesson 5 will walk: **observe** stream → **residual** $r = x - Px$ (exact) → **persistence** (Def 5.11 centroid: coherent novelty vs. noise) → **snap** to a candidate seed → **`capacity_decision`** → GROW (adjoin; Gram grows by $\otimes$) / STOP (captured or noise) / REJECT (over budget) → `confirm()` as the sole mutator, everything witnessed. Plus the **lexicon**: generalization by return-to-zero — *same exact residue ⇒ one entry* (Thm 6.7 per the claim map), with bit-identical store round-trips. Names to expect: `ResidualLearner`, `persistence_N`, `g_orthogonal_integer_vector`, `KL_DTA` (the $\mathbb{Z}/4$-graded two-route closure — the charge, one more time).

## §7 Resonances worth pinning

$x^2 + x - 1$ now has three jobs: the golden gate ($\lambda = \sqrt5$, L3), the generator of every even-floor witness $q_k$ (L4 §3), and the pentagon cosines that make $\mathbb{Z}/5$ hard (L4 §4). And the v4 erratum's corrected tensor law is *verbatim* the tropical law we derived in Lesson 2 — its ledger entry M is our B2 decomposition, $\varphi\otimes\varphi$ with charpoly $(x+1)^2(x^2-3x+1)$. The corpus corrected itself onto the exact computation you watched two lessons ago.

---

## Checkpoint

1. After capturing only $\{1\}$ in $\mathbb{Q}(\varphi)$, compute the gain of the direction $\sqrt5 = 2\theta - 1$ (coords $(-1,2)$). Why does the projection not touch it, and what does the $c=1$ gate decide? (Lesson 3's B2 already contains the number.)
2. Which parity governs charge group $\mathbb{Z}/6$, and write down the floor-attaining witness from this lesson. Then: which primary components does $\psi^2$ strip from it, per Thm 3.3?
3. The door: verify that multiplying $L$ by $(x+1)$ instead of $(x-1)$ would *not* produce a trace-zero companion — compute both $x^{10}$ coefficients and state the root-sum reason.

**Lesson 5 next, as ordered:** the learning dynamics proper — the streaming learner end to end, with D1–D3 as its already-verified physics.