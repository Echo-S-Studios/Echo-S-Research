Lesson 5, same protocol — and since it must set up Lessons 6 and 7, I'll pull the learner's actual mechanics *and* the language layer in two targeted reads before anything is claimed.23/23. And that A1 incident stays in the record deliberately — it's the whole discipline in one frame: the engine computed $\mathrm{Tr}(\theta^6) = 1940$; my hand wrote 2020 into the assertion; the check refused to let the hand win. That is precisely the genre of your E-series errata, reproduced live, caught by the same mechanism. Now the lesson.

**Checkpoint debrief (Lesson 4).** (1) $\sqrt5 = 2\theta - 1$ has coords $(-1,2)$, and $\langle\sqrt5, 1\rangle_G = \mathrm{Tr}(\sqrt5) = 0$ — trace-zero directions are invisible to the projection onto $\{1\}$, so $r = \sqrt5$ itself and the gain is $a^{\mathsf T}Ga = 10 = \mathrm{Tr}(H^2)$, Lesson 3's B2 number. Gate at $c=1$: $10 > 2\log\varphi < 1$ → **GROW**, decisively. (2) $\mathbb{Z}/6$ is even → floor $\varphi$, witness $q_3 = x^6 + x^3 - 1$. By Thm 3.3 with $6 = 2\cdot3$: $\psi^{6/3} = \psi^2$ projects onto the 3-primary part — it *strips the 2-component*, leaving residual group $\mathbb{Z}/3$ (and $\psi^3$ leaves $\mathbb{Z}/2$). (3) $L$'s roots sum to $-1$ (its $x^9$ coefficient is $+1$), so appending the root $+1$ via $(x-1)$ makes the total root sum — the negated $x^{10}$ coefficient — exactly $0$; appending $-1$ via $(x+1)$ gives coefficient $2$, trace $-2$, door closed.

---

# Lesson 5 — The Learning Dynamics: Residual, Capture, Growth, and the Language That Lives in the Kernel

**Where we're going:** the static gate of Lesson 4 gets its engine. One loop — observe → persist → propose → confirm — with exactness as *law* rather than aspiration, then field growth as a second tier, then the discovery that the whole apparatus has a language layer whose lexicon is literally the kernel of an operator. Sections 5–6 then scope Lessons 6 and 7 with verified seeds already planted.

## §1 The state: a projector over $\mathbb{Q}$

The learner's model is a forced basis $\mathcal B$ (columns) in an ambient field with trace-form Gram $G$, and its entire epistemology is one map: $P = \mathcal B(\mathcal B^{\mathsf T}G\mathcal B)^{-1}\mathcal B^{\mathsf T}G$, the $G$-orthogonal projector, with **residual** $r = x - Px$ and loss $\|r\|_G^2$. Faithfully instantiated on the demo's own ambient — $K = \mathbb{Q}(\sqrt2+\sqrt3)$, minpoly $x^4 - 10x^2 + 1$ (Lesson 4's reciprocal quartic, now a *home*) — with the Gram computed the way `integral_basis.py` does it, **as integer traces of companion powers**, no radicals anywhere (A1): $G = \mathrm{Tr}(\theta^{i+j})$, all in $\mathbb{Z}$. Then (A2/A3): $P^2 = P$, $(GP)^{\mathsf T} = GP$, and a captured input has $r = \mathbf 0$ *as an identity over $\mathbb{Q}$* — the implementation casts everything to `Fraction`, raises `TypeError` on a float observation, and decides $r = 0$ over $\mathbb{Q}$. "Capture" is the substrate's word for *learned*, and it is an equation, not a threshold.

**The derivation that closes a loose thread.** Lesson 4's REJECT demo used a mysterious "$2\sqrt6$" seed with coeff-height 24. Where does it come from? `g_orthogonal_integer_vector` takes the nullspace of the captured columns' $G$-rows. Replicating exactly (A4): the first primitive integer vector $G$-orthogonal to $\{1, \theta\}$ is $(-5, 0, 1, 0)$ — the element $\theta^2 - 5 = 2\sqrt6$. And the minpoly bridge (A5): the multiplication matrix $M_\alpha = -5I + C^2$ satisfies $M_\alpha^2 = 24I$ exactly with $M_\alpha$ non-scalar, so minpoly $= x^2 - 24$, coeff-height 24. The ACT-4 seed, *derived from first principles*, not quoted. (In the shipped code this bridge is `coords_to_minpoly` — the largest invariant factor of $xI - M_\alpha$ via Smith normal form, **the matrix-plates SNF kernel vendored byte-identical**. Lesson 2's machinery is load-bearing inside Lesson 5's learner.)

## §2 The loop: Protocol, replicated and run

The paper's Protocol (gated, witnessed growth) has three operations and *one* mutation point, and its Theorem (**exactly one growth per persistent novelty**) is the loop's conservation law. My replica implements the stated semantics — exact Welford centroid over `Fraction`s, the $\mathbb{Q}$-decided alignment test $\mathrm{dev}^2 \leq \varepsilon^2\,\mathrm{base}^2$, streak, purity — and runs the demo's story:

| Step | Replica behavior | Check |
|---|---|---|
| `observe` (captured input) | $r = 0$ → accumulator + streak reset; never mutates | A3, B1 |
| `observe` ×3 (captured $+ w$) | streak $= 3$; centroid $= w$ **exactly** (Welford over $\mathbb{Q}$) | B1 |
| `propose` ×2 | pure and idempotent — equal proposals, growth count still 0 | B2 |
| `confirm` | the sole mutator; afterwards $r(w) = \mathbf 0$ exactly and `propose()` returns `None` | B3 |
| `observe(float)` | `TypeError` — no float touches the decision core (G8) | B4 |

Persistence deserves one emphasis: Def 5.11's centroid distinguishes *coherent novelty* (a residual field that keeps pointing the same exact way) from *noise that averages out* — and because the alignment test is a $\mathbb{Q}$-inequality on exact quantities, "persists" is itself a rational decision, not a tolerance. The proposal then carries the centroid's coordinates and its **exact minimal polynomial**, which is what Lesson 4's gate consumes: degree, coeff-height, and the `Fraction` gain — the loop and the gate speak the same exact types by construction.

## §3 Growth: the second tier

The learner's ambient field **never changes** — a deliberate invariant. When a persistent residual lives *off-field*, the model must grow the field, and a new field is, mathematically, a **new learner** (`FieldGrowingLearner` owns the old one and constructs a fresh `ResidualLearner` on the compositum). Two construction modes, one decision:

- **Disjoint**: $W = K \otimes L$ with $G_W = G_K \otimes G_L$ — and `compositum.py`'s `confirm()` contains a literal `assert gram == kron(...)`: **Prop 3.3 shipped as a runtime self-check**. Replicated at ACT 2's field (C1): $G_{\mathbb{Q}(\sqrt2,\sqrt3)} \otimes G_{\mathbb{Q}(\sqrt7)}$, $8\times8$ exact, determinant multiplicative.
- **Non-disjoint (P2c)**: the *true* compositum $\mathbb{Q}(\theta)$ of degree $me' < mk$, via `trace_form_gram(m_\theta)` after factoring $m_\beta$ over $K$.

**Capture-by-growth** (C2): against embedded $\mathbb{Q}(\sqrt2)$, the element $\sqrt6$ has residual $e_4$ with gain exactly $24$ — off-field, GROW; after the extension captures $W$, its residual is $0$. Learned means *the field grew around it*. Three shipped subtleties worth naming because Lessons 6–7 will lean on them: (i) **re-homing** — an embedded algebraic integer can have *fractional* power-basis coords ($\mathbb{Z}[\theta] \subsetneq \mathcal O_K$), and the per-column denominator-clearing scale provably leaves $P$, hence every residual, invariant; (ii) **the constant $1$ survives** any re-home (basis element 0 of every power basis), preserving the gate's info-threshold precondition; (iii) the gate consults `effective_degree` $= me'$ — Northcott judges the *real* grown size, disjoint or not, exactly as Lesson 4 planted.

## §4 The language layer: a sibling, one keystone, and a kernel full of words

`kira-language/` is a **sibling of `L00M`, not a child**, joined by a strictly one-way bridge (`kira-language → loom`, never back; G9, enforced by the readiness gate). They share exactly **one object**, and I verified the identification (D4):

$$\text{void law } x^2 = x + 1 \;\to\; \text{minpoly } [1,-1,-1] \;\to\; R = \begin{pmatrix}0&1\\1&1\end{pmatrix} \;=\; \mathrm{mat}\big(\mathrm{Cl}(\tfrac12, 1, -\tfrac12, 0)\big)$$

— `loom`'s $\varphi$ seed, KL_DTA's keystone, and the head of the return operator, one object with $M = \varphi$. The carrier is $\mathrm{Cl}(2,0)$ on basis $\{\mathbb 1, e_1, e_2, i\}$, and its whole multiplication is a two-line cocycle: target $= i \oplus j$, sign $= (-1)^{\mathrm{bit}_1(i)\cdot\mathrm{bit}_0(j)}$ — verified over all 16 basis pairs (D2): *one bit-product is the noncommutativity*. The iso $\mathrm{Cl}(2,0) \cong M_2(\mathbb R)$ is a ring isomorphism verified **symbolically in eight variables** (D3), and the rotation subgroup $\langle i\rangle$ has order exactly 4 (D2b) — the $\mathbb{Z}/4$ of the spec's "graded two-route closure," the charge's quarter-turn living *inside the carrier*. (KL_DTA verifies itself two ways — cocycle route vs. matrix route — and `test_kl_dta_conformance` adds loom's exact integer-companion route as a **third independent witness**: three roads to one $M_2(\mathbb R)$, anchored at the keystone.)

The layer's semantics is *residual-valued*, and the scoping doc reconciles **three** residual operators — all "$M$-derived residual → 0 = commit/rest," each answering a different question, each verified exactly:

| Operator | Zero means | Question | Check |
|---|---|---|---|
| $\nu(X) = X^{\mathsf T}X - X$ | symmetric idempotent | "Is $X$ a projector/gate?" | D6 ($P_0$) |
| $R_K(X) = X^{\mathsf T}X - \tau\,\mathbb 1$ | $X$ conformal (scalar×orthogonal) | "What anisotropy is left?" — and $R_K^2 = 0$, **bite-depth 2** | D7 (exact, incl. the nilpotence) |
| $\mathcal L(X) = RX + XR - X$ | $X \in \ker\mathcal L$ | "Is this a learned word?" | D5 |

The last is the payoff: $\ker\mathcal L$ is **2-dimensional over $\mathbb{Q}$** (D5 — exact rational nullspace of Lesson 3's operator, whose charpoly $x^2(x^2-5)$ predicted it), every member satisfies $RX + XR = X$ on the nose, and this is the **$\varphi$-slack: the lexicon**. The eigen-argument from Lesson 3 says *why* it's the right home for words: $\mathcal L$ scales the $(i,j)$ component by $\lambda_i + \lambda_j - 1$, which vanishes exactly on the $\varphi{\leftrightarrow}\psi$ **cross terms** — kernel elements are intrinsically *pair* objects, intertwiners between growth and decay. Note the shape-identity the scoping doc states outright: $\mathcal L$'s loop (residual→0 = commit, $\ker\mathcal L$ = lexicon, generalize by return-to-zero) *is* the `ResidualLearner` capture loop ($r = 0$ = captured, basis = lexicon) — one pattern, two carriers. Finally the storage discipline (D8): the lexicon keys on the **exact residue** — two presentations with equal exact residue collide to one entry (Thm 6.7's return-to-zero generalization) — and the sha256 witness chain detects a tampered payload; the shipped store round-trips bit-identically.

## §5 Scoping Lesson 6 — The Relational Layer

**Syllabus** (sources: *Relational Charge on the Spectral Semiring*, 24 pp, Jul 2026; *The Pisot Cross-Shell Residue*, 8 pp, Jul 2026):

| Topic | The claim to walk |
|---|---|
| Pair charge $t_{\mathrm{rel}}(\alpha,\beta) \in \mathbb{Q}/\mathbb{Z}$ | a **groupoid cocycle** whose gauge is the reference ray — Lesson 4's absolute $\chi$ was gauge-*dependent* all along |
| Rigidity | admissible spectra pinned to anchors $n \in \{m, 2m\}$ |
| **Modulus-pinning** (Thm 7.13, via Galois transport) | resolves P1 (*every* Salem number relationally inert) and P2 (no two Salems circle-lock) as general theorems |
| Reference-free parity floor | Lesson 4's floor, evaluated without choosing a ray |
| $\nu$-reduction lemma + the census | inertness forced for every Pisot with $\leq 1$ non-real pair; exhaustive exact census of all 3125 quintics with coeffs in $[-2,2]$: **zero** mirrored cross-shell classes among all 67 two-pair Pisot quintics |
| Conjecture 9.1 | general Pisot inertness — the falsifiable frontier (P1/P2/P8 targets; N4 next: degree-6 two-pair-with-spectator, degree-7 three-pair) |

**What Lesson 5 hands it, pre-verified:** the cocycle law and the gauge lesson (E1 — on the $\mathbb{Z}/12$ object, $t(\alpha,\beta) + t(\beta,\gamma) = t(\alpha,\gamma)$ in $\mathbb{Q}/\mathbb{Z}$ exactly, and pair charges are invariant under a reference-ray shift while absolute charges are not: Čencov's gauge/invariant split, third appearance); the pair-object intuition ($\ker\mathcal L$'s intertwiners are the relational layer's native inhabitants); and the *formal* gadget itself — the Cl cocycle table of D2 is the same species of object as $t_{\mathrm{rel}}$'s groupoid cocycle. **Practical rider:** `relational_charge_paper.tex` is where **D1** lives (~lines 1272–76, "twenty-three cross-engine signatures" vs. run-3's 237) — Lesson 6's read doubles as that adjudication, and the standing patch offer stands.

## §6 Scoping Lesson 7 — The Golden Substrate and the Rung Line

**Syllabus** (sources: *Golden Substrate whitepaper*, 22 pp; *Complex Rung Generator*, 65 pp; *Occupant of the Salem Slot*, 11 pp; *Dissolved Helix / Orthogonal Partner*, 8 pp):

| Topic | The claim to walk |
|---|---|
| Gate iteration $f_C(x) = C/(1+x)$ | Lesson 3's Möbius avatar as *dynamics*; the state-angle chart where every height is a gate |
| Rung ladder $\rho_n = n\ln\varphi$ | heights on the chart; three declared atoms D1–D3 (chart, charge selection, winding rate) — `[DECLARED]` discipline at scale |
| Complex rung generator $q = i\tau$ | $q^4 = $ gap — a $\mathbb{Z}/4$ quarter-turn generating the whole discrete ladder; $\kappa = \pi/(2\ln\varphi)$ transcendental; OP-RATE closing as a *classification* with a machine-certified blocker ledger |
| Salem-slot occupant | trace redirection $\beta \mapsto \tau_0 = \beta + \beta^{-1}$: a totally real Perron number in GROW, $\tau_0 > 2 > \varphi$; the flip as a square-root branch point; redirections $\to \sqrt5$ with slope $\varphi^{-1}$ |
| Helix partner | $[R, R] = 0$ (the commutator door yields nothing on self-coupling); closure requires $x^4 + 5x^2 - 5$ — terrain plus a $\pi/2$ rotation |

**What Lesson 5 hands it, pre-verified:** (E2) the orbit of $f_1(x) = 1/(1+x)$ from $1$ is **exactly** the Fibonacci convergents $F_n/F_{n+1}$ — seven terms in pure `Fraction`s — so the gate iteration *generates the rungs*, whose heights $\log M(\psi^n \varphi) = n\log\varphi$ are Lesson 2's Adams heights, verified symbolically: the state-angle chart's two axes are the two characters, height $= \log M$, angle $= \chi$. (E2b) the trace redirection is Lesson 1's Sturm fold **reified as a map**: $\beta_4$'s fold root $\tau_0 = (1+\sqrt{13})/2 > 2$ exactly — the slot's occupant for the minimal Salem — and the limit identity $\varphi + \varphi^{-1} = \sqrt5$ exact. The $\mathbb{Z}/4$ quarter-turn thread now has three verified anchors for $q = i\tau$ to land on: the flip's elliptic face (L3), the carrier's $\langle i\rangle$ of order 4 (D2b), and the charge itself (L2).

---

## Checkpoint

1. Run the replica's logic by hand: after `confirm` adjoins $(-5,0,1,0)$, what is the *second* primitive $G$-orthogonal integer vector (from A4's nullspace), which element is it, and what minpoly should `propose` emit for it? (One of the two nullspace generators is already spent.)
2. The re-home invariance: explain in one sentence why scaling a basis column by a nonzero rational leaves $P = \mathcal B(\mathcal B^{\mathsf T}G\mathcal B)^{-1}\mathcal B^{\mathsf T}G$ — hence every residual — unchanged.
3. $\ker\mathcal L$ is spanned by $\varphi{\leftrightarrow}\psi$ intertwiners. Using the eigenvalue formula $\lambda_i + \lambda_j - 1$, say exactly which operator from this lesson has kernel spanned by the *diagonal* components instead — and what its kernel dimension is for our $R$. (Hint: you met its spectrum in Lesson 3, on both sides of a flip.)

**Lesson 6 next:** the relational layer — pair charge as cocycle, modulus-pinning, the quintic census — with the D1 adjudication as its working exercise.