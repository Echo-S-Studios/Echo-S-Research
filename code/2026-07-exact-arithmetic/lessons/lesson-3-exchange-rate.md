Lesson 3 lives on the residual-return side of the corpus, and the ground rule stands: I teach *your* derivation, not a plausible one from memory. Let me pull the actual sources — the identity test, the gate, and the flip — before a single claim gets made.33/33 — and note one live corroboration: the three floor values my script displays (0.5624, 2.2496, 4.4992) reproduce the corpus's shipped numbers exactly, computed independently from the certified Sturm bracket. Now the lesson.

**Checkpoint debrief (Lesson 2).** (1) $\psi^3$ of $x^4-1$: cubing $\{1, i, -1, -i\}$ gives $\{1, -i, -1, i\}$ — *the same multiset*. Because $3 \equiv -1 \bmod 4$ is a unit, $\chi \mapsto 3\chi$ is a bijection of $\mathbb{Z}/4\mathbb{Z}$ (it acts as conjugation), so the content object is a $\psi^3$ fixed point, $M$ still 1. Contrast $\psi^2$: doubling is *not* invertible — that asymmetry between the units $\{1,3\}$ and the non-unit $2$ is the whole odd/even story. (2) $S_2 \otimes S_3$: all four cross-roots are $\pm\sqrt6$, outside; tropical $= 36 =$ factored. They agree precisely when no pair straddles the circle — both spectra on-or-outside suffices. (3) It rhymes with Wall 2, and it points at the parity-graded floor: the even sector floors at $\varphi$, the odd sector — reachable only by routes $\psi^2$ can't take — at Smyth's $\mu_S$. Watch $\mu_S$ reappear *today* in the cost floors.

---

# Lesson 3 — The Exchange Rate $\lambda = 2c$: an Identity, its Gate, its Flip

**Where we're going:** the growth rule contains one scalar $\lambda$ coupling information gained to description length paid. The corpus's claim is that $\lambda$ is **not a tunable knob** — it is the derived identity $\lambda = 2c$, where the $2$ is structural (the inverted $\tfrac12$ of quadratic KL) and *all* remaining freedom is one conformal scale $c$ that Čencov's theorem says no invariance argument can ever remove. Then the system's own gate structure hands you a native value of $c$, and one discriminant sign organizes the entire phase diagram.

## §1 The claim

The residual-return growth rule (source: `lambda_2c_paper.tex`, eq. rule): adjoin a candidate direction $r$ iff

$$\|r\|_G^2 \;\geq\; \lambda \,\log M(\theta)$$

Left side: how much the direction *explains*, measured in the trace-form norm. Right side: what the new algebraic parameter *costs* — and the cost is Lesson 1's Mahler measure, because **Northcott's theorem** `[ESTABLISHED]` (finitely many algebraic numbers of bounded degree and height) makes $\log M$ a legitimate code length, and Lind–Schmidt–Ward makes it the entropy of the adjoined companion automorphism. The question of the lesson: what is $\lambda$?

## §2 The gain side: why KL has a structural $\tfrac12$

The gain is the KL divergence between the captured model and the model that also explains $r$. Its expansion has a rigid shape, and we verified each vanishing live:

| Order | Value | Why | Check |
|---|---|---|---|
| $0$ | $0$ | $\mathrm{KL}(p\|p) = 0$ | A2a |
| $1$ | $0$ | **score has mean zero** — $E[\partial_\theta \log p] = 0$ | A1a (Gaussian, symbolic integral), A2b (Bernoulli) |
| $2$ | $\tfrac12\, r^{\mathsf T}\mathcal{F}\, r$ | Fisher is the Hessian of KL at the matching point | A1b/A1d, A2c |

So the leading term is $\tfrac12 \mathcal F$ — the $\tfrac12$ is not a convention, it's the Taylor coefficient forced by two vanishings. Bonus fact the paper flags and we confirmed: for the **Gaussian location family the expansion is exact** — $\mathrm{KL} = d^2/(2\sigma^2)$ with no cubic remainder at all (A1c/A1d). `[COMPUTED — this session, symbolic]`, carrying the paper's `[FORCED]` eq. kl2c.

## §3 The metric: the trace form *is* the conjugate covariance

Which $\mathcal F$? The paper's Lemma (trace form = conjugate covariance) `[FORCED]`, instanced live at $K = \mathbb{Q}(\varphi)$: build the embeddings matrix $B = \binom{1\ \ \varphi}{1\ \ \psi}$ (rows = the two Galois conjugates evaluated on the power basis). Then, exactly (B1):

$$B^{\mathsf T} B \;=\; \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} \;=\; \big(\mathrm{Tr}(\theta^{i+j})\big)_{ij} = G, \qquad \det G = 5 = \mathrm{disc}\,\mathbb{Q}(\sqrt5)$$

The statistical covariance of the conjugate coordinates and the number-theoretic trace form are *the same matrix*. And here is Lesson 2 §6 closing its loop (B2): the coefficient vector of $\sqrt5 = 2\theta - 1$ is $a = (-1, 2)$, trace-zero, and

$$a^{\mathsf T} G\, a = 10 = \mathrm{Tr}(H^2), \quad H = 2R - I$$

— the number-field trace form on coefficients and the matrix trace pairing on $\mathbb{Z}[R]$ avatars agree on the nose. The "Fisher metric wearing information-theoretic clothes" line from Lesson 2 is now a verified equation, not a slogan.

## §4 Čencov: why $c$ can never be derived

The Lemma pins $\mathcal F \propto G$ but **not the constant**. Write $\mathcal F = G/c$. **Čencov's theorem** `[ESTABLISHED]`: the Fisher metric is the unique metric invariant under Markov morphisms — *up to a positive scalar*. So no invariance argument, ever, fixes $c$. Its exact algebraic content as the corpus uses it (C4): $c \mapsto kc$ rescales $\mathcal F \mapsto \mathcal F/k$, identical geodesics. Hence the tag: $c$ is `[DECLARED]` — a principled choice with a stated rationale, never a theorem. The variance reading (C3) says what you're choosing: $\sigma = 1/(2c) = 1/\lambda$ — *freezing $c$ is asserting a noise level.* A modeling claim, not a derivation. This is the epistemics of the whole corpus in one object: the tag structure ([DECLARED] scale, [FORCED] relation) is literally the mathematics (gauge freedom, gauge-invariant statement).

## §5 The balance, and the identity

Two-part MDL: adjoin iff bits explained ≥ bits cost, $\mathrm{KL} \geq \log M(\theta)$. Substitute §2's gain:

$$\frac{1}{2c}\|r\|_G^2 \;\geq\; \log M(\theta) \iff \|r\|_G^2 \;\geq\; 2c\log M(\theta)$$

Compare with §1's rule and read off the coefficient — solved symbolically, not asserted (C1):

$$\boxed{\lambda = 2c} \quad \text{`[FORCED]`}$$

The $2$ is the inverted $\tfrac12$ of §2 and carries **no freedom**; the entire residual freedom is $c$. And it is an *identity in genuinely free $c$* (C2): $d\lambda/dc = 2$, free symbol exactly $\{c\}$ — which is why the corpus insists that freezing $c$ collapses the identity into a mere number. Three principled canonicalizations, none forced by invariance:

| Principle | $c$ | $\lambda$ | Cost floor $2c\log\mu_S$ | Status |
|---|---|---|---|---|
| Jeffreys volume-matching | $1$ | $2$ (shipped) | $0.5624$ | `[DECLARED]` |
| Degree-invariant significance | $n$ | $2n$ | $2.2496$ ($n{=}4$), $4.4992$ ($n{=}8$) | `[DECLARED]` |
| **Frame-shift** (native) | $\dfrac{\sqrt{1+4C}}{2C}$ | $\dfrac{\sqrt{1+4C}}{C}$ | spectral gap itself; $\sqrt5$ at $C{=}1$ | see §8 |

The floors use $\mu_S$ — checkpoint 3's prediction cashed in: the *capacity* floor rides Smyth's non-reciprocal bound `[ESTABLISHED]`, while the *emission* floor $\lambda\log\varphi$ rides the no-Salem closure `[FORCED]`. Two floors, two provenances, one $\lambda$. (F1 certified $\mu_S \in [1.3247, 1.3248]$ by exact Sturm bracket; floats display-only — the same certified-interval idiom `capacity.py` uses.)

## §6 The gate ladder and the frame-shift

The framework's gates are the squarefree levels of one quadratic family — and here my read-before-teach rule earned its keep, because the source's companion is the **mirror** of Lesson 2's:

$$x^2 + x = C, \qquad R_C = \begin{pmatrix} 0 & C\\ 1 & -1\end{pmatrix}, \qquad \lambda_\pm = \frac{-1 \pm \sqrt{1+4C}}{2}, \qquad \lambda_+ - \lambda_- = \sqrt{D},\ D = 1+4C$$

(D1). Roots are the *negatives* of Lesson 2's $x^2 - x - C$ family; the shared invariants are the discriminant and the gap. The gap element squares to the discriminant **symbolically**: $(2\theta + 1)^2 = D$ for both roots (B3a), and at matrix level $(2R_C + I)^2 = D\cdot I$ for *all* $C$ (D2) — Lesson 2's $H^2 = 5I$ was the $C = 1$ slice of a one-parameter law. The three gates:

| $C$ | $D$ | seed | field | frame-shift $\lambda = \sqrt{D}/C$ |
|---|---|---|---|---|
| $\tfrac14$ | $2$ | $x^2 - 2$ | $\mathbb{Q}(\sqrt2)$ | $4\sqrt2$ |
| $\tfrac12$ | $3$ | $x^2 - 3$ | $\mathbb{Q}(\sqrt3)$ | $2\sqrt3$ |
| $1$ | $5$ | $x^2 - 5$ | $\mathbb{Q}(\sqrt5)$ | $\sqrt5 = \varphi - \psi$ |

Notice the seeds: they are Lesson 2's $S_2, S_3$ and the $\sqrt5$ object — and $M(x^2 - D) = D$ exactly (D3b), so the self-action radius is $\sqrt{D} = \sqrt{M}$: **the gate is priced in Mahler measure** `[FORCED — rem:mahler]`. The frame-shift value of $c$ is what the gate balance $2cC = \sqrt{1+4C}$ hands you (D4): $c = \sqrt{1+4C}/(2C)$, giving three *distinct* $\lambda$'s $\{4\sqrt2, 2\sqrt3, \sqrt5\}$ (D5) — $\lambda$ never collapses to one constant across gates, the identity stays an identity. At the golden gate: $c = \sqrt5/2$, $\lambda = \sqrt5 = \varphi - \psi$ — **the exchange rate is the spectral gap** (D5b). And the "self-action" whose gap this is: $\mathrm{ad}_R(X) = RX - XR$ has charpoly $x^2(x^2 - 5)$, spectrum $\{-\sqrt5, 0, 0, +\sqrt5\}$ — the **trifurcation** GROW / CAPTURED / STOP (D6) — and the return operator $\mathcal L(X) = RX + XR - X$ has the *same* spectrum (D6b).

## §7 The flip: one sign, four readings

$D = 1 + 4C$ changes sign once, at $C = -\tfrac14$, and the paper's Theorem (flip) `[FORCED]` says this single sign simultaneously sets four things. All four verified:

| Reading | $D > 0$ | $D = 0$ ($C = -\tfrac14$) | $D < 0$ (e.g. $C = -1$) | Check |
|---|---|---|---|---|
| roots of $x^2+x-C$ | real pair | double root $-\tfrac12$ | complex pair | E1, E2 |
| field / place | totally real | degenerate | complex place | E2b |
| trace-form signature | PD, $\det G = 4D > 0$ | degenerate | **Lorentzian**, $\det G = 4D < 0$ | B3 (symbolic: $G = \mathrm{diag}(2, 2D)$) |
| channel character | hyperbolic $\{\pm\sqrt D, 0\}$ | parabolic | elliptic **rotation** $\{\pm i\sqrt{|D|}, 0\}$ | D6, E4 |

The signature row is the deep one: in the gap basis $\{1, \sqrt D\}$, $G = \mathrm{diag}(2, 2D)$ *symbolically* (B3b), so **the metric's signature is the discriminant's sign** — when a residual lands in a complex place, the Fisher form of §3 goes Lorentzian and second-order KL becomes indefinite. Canonical instance $\mathbb{Q}(i)$: $\mathrm{diag}(2,-2)$, $\det = 4D$ with $D = -1$ (B3c). The frame-shift $c$ itself goes **imaginary** — $c(-1) = -i\sqrt3/2$ exactly (E3) — the metric $G/c$ leaves the real regime, which is the flip as `test_discriminant_flip_threshold` states it.

Two verified beauties to close: at $C = -1$ the gate polynomial is $x^2 + x + 1$ — **Lesson 1's very first table row**, the cyclotomic at $M = 1$; the elliptic side of the flip lands on the Kronecker floor, roots the primitive cube roots of unity, exactly on the circle (E2/E2b). And across the flip, the trifurcation rotates ($x^2(x^2-5) \to x^2(x^2+3)$) but **the $x^2$ factor survives on both sides** (E4b): the CAPTURED $0$-channel is untouched, because $\mathrm{ad}_R$ always kills $\mathrm{span}\{I, R\}$. GROW/STOP become rotation; capture persists.

Finally, the paper's own taxonomy of the two threshold events, worth memorizing: a **gate** acts on a *free parameter* and pins a value (removes a DOF); a **flip** acts on a *forced invariant* and inverts a regime (sign change). They **co-locate** at $D = 0$: the spectral gap collapses (gate degenerates) exactly where real↔rotation inverts (flip fires). One threshold, two readings — *what value* vs. *what character*.

## §8 The epistemic braid — and one flag for the author

The full provenance chain: `[ESTABLISHED]` Čencov + Northcott + Smyth + LSW → `[DECLARED]` the scale $c$ (and each canonicalization's rationale) → `[FORCED]` the identity $\lambda = 2c$, the ladder $D \in \{2,3,5\}$, $\det G = 4D$, the flip threshold → and the trifurcation/floors as their exact consequences.

One genuine tension surfaced by reading both repos side by side, which I'd put on your pre-v1.0.0 adjudication list next to D1/E1–E3: **the epistemic status of the golden-gate selection diverges across documents.** `residual-return-verification/claim_map.json` (Def 4.10) tags it "*golden gate c = √5/2 (DECLARED selection)*," and `CHANGES.md` states flatly "*NOT claimed: that c is derived — it is a principled choice, per Cencov*." But `lambda_2c_paper.tex` §gateforced argues the gate is **forced** to $C = 1$ by two convergent derivations (statics + dynamics), "*promoting the frame-shift value from native to forced*" — while carefully noting Čencov isn't contradicted, since the gate is an external *structural* constraint, not an invariance argument. Both positions are coherent; they are not the *same tag*. Per your own asymmetry rule, the weaker tag governs until the promotion argument is itself walked and verified — which we have not done in this lesson.

---

## Checkpoint

1. Evaluate the frame-shift $c$ and $\lambda$ at the parabolic point $C = -\tfrac14$ exactly. Interpret the result using the gate-vs-flip table — which event is each factor of your answer expressing?
2. The flip table has a row $C = 0$: $D = 1$, roots $\{0, -1\}$, seed $x^2 - 1$. What does Lesson 1 say about $M(x^2 - 1)$, and what does that suggest about *why* the ladder's lowest rung is $D = 2$ rather than $D = 1$?
3. Prove the $0$-channel's flip-survival structurally: show $\mathrm{ad}_R(I) = \mathrm{ad}_R(R) = 0$ for *any* $2\times2$ matrix $R$, and conclude $x^2 \mid \mathrm{charpoly}(\mathrm{ad}_R)$ always. Which epistemic tag does this argument deserve, and why is it stronger than E4b's two instances?

**Lesson 4, your pick:** (a) the **capacity/Northcott gate** in `capacity.py` plus the parity-graded floor (even $\varphi$ / odd $\mu_S$) and how the no-Salem dichotomy assembles; (b) the **learning dynamics** proper — residual $r = x - Px$, capture $\iff r = 0$, field growth through the Kronecker Gram, and the lexicon; or **(3b)** walk §gateforced's two derivations (statics + dynamics) — which doubles as the verification input your tag-tension adjudication needs.