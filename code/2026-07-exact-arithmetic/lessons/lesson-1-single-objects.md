Good — let's do this properly. Before I teach anything, let me pull the corpus's own pedagogical documents so I teach *your* derivation chain, not a generic textbook one.16/16 in exact arithmetic. (One cosmetic slip on my side: the C6 note string contains a stray draft comment — the *decision* was exact and correct, φ² − φ = 1 > 0; the label was sloppy. Noted, owned, math stood.) Now the lesson, with every claim tagged and pinned to the check that just ran.

---

# Lesson 1 — The Confinement Story: from one definition to the two walls

**Where we're going:** one definition (Mahler measure) → one classical floor (Kronecker) → one open question (Lehmer) → one structural split (Smyth / Salem) → the corpus's move (build an emitter, prove two walls). Everything below is either **derived live** (check IDs from the run above) or **imported with a tag** naming the file that carries it. That asymmetry — derive what you can, pin what you import — *is* the method being taught.

## §1 The object: one definition, three values

For a monic integer polynomial $P(x) = \prod_{i=1}^d (x - \alpha_i)$:

$$M(P) = \prod_{i=1}^{d} \max(1, |\alpha_i|)$$

Read it as: **only roots outside the unit circle contribute; everything on or inside costs a factor of 1.** Three computations you should be able to reproduce by hand:

| Polynomial | Root picture | $M$ | Verified |
|---|---|---|---|
| $x^2 + x + 1$ | both roots *on* the circle ($e^{\pm 2\pi i/3}$) | $1$ | C1c `[COMPUTED]` |
| $x^2 - x - 1$ | $\varphi \approx 1.618$ outside; $-1/\varphi$ inside | $\varphi$ (exactly $\tfrac{1+\sqrt5}{2}$) | C1a `[COMPUTED]` |
| $x^2 - 2$ | $\pm\sqrt2$, both outside | $\sqrt2\cdot\sqrt2 = 2$ | C1b `[COMPUTED]` |

Note the middle row: only **one** root escapes the circle, so $M$ is that single root, $\varphi$. Hold onto that polynomial — the entire corpus is built on it.

## §2 The floor at 1 — completely understood

Since every factor is $\geq 1$, $M(P) \geq 1$ always. **Kronecker (1857)** `[ESTABLISHED]` settles equality: $M(P) = 1$ iff every root is zero or a root of unity, i.e., $P$ is cyclotomics times a power of $x$. So $M = 1$ is the *trivial locus* — fully classified. Every open question lives strictly above 1.

## §3 Lehmer's question — and a derivation, not an assertion

**Lehmer (1933)** `[OPEN]`: is there $c > 1$ with $M(P) \geq c$ whenever $M(P) > 1$? Or can measures crowd down onto 1? His record-holder, unbeaten for 93 years:

$$L(x) = x^{10}+x^9-x^7-x^6-x^5-x^4-x^3+x+1, \qquad M(L) = 1.17628\ldots$$

Textbooks *state* that $L$ has one root outside the circle. We just **derived** it, using the house Sturm method — this is worth internalizing because it's the exact-decision idiom the pipelines use everywhere:

1. **Reciprocity is a polynomial identity**, not a numerical observation: $x^{10}L(1/x) = L(x)$ exactly (C2a). So roots pair as $\alpha \leftrightarrow 1/\alpha$.
2. **Fold the symmetry out.** Because $L$ is reciprocal of even degree, $L(z)/z^5$ is a polynomial in $t = z + 1/z$. Using the recurrence $p_k = t\,p_{k-1} - p_{k-2}$ for $z^k + z^{-k}$, we computed it exactly (C2c):
$$Q(t) = t^5 + t^4 - 5t^3 - 5t^2 + 4t + 3$$
3. **The circle becomes an interval.** $|z| = 1 \iff t = z + 1/z \in [-2, 2]$ (real). So counting on-circle roots of $L$ = counting real roots of $Q$ in $(-2,2)$ — a **Sturm count**, decided over $\mathbb{Q}$, no floats. Result (C2d/C2e): $Q$ has 5 real roots — **4 inside $(-2,2)$, 1 in $(2,\infty)$**.
4. **Unfold.** Each interior root ⇒ a conjugate pair on the circle: $8$ on-circle roots. The one exterior trace root ⇒ a real pair $(\beta, 1/\beta)$ with $\beta > 1$. Total accounted: $8 + 2 = 10 = \deg L$. Hence **exactly one root outside the circle**, and $M(L) = \beta \approx 1.17628$ (float for display only; C2f pins $\beta \in (1,2)$ exactly).

You have just derived that $L$ is a **Salem polynomial**: irreducible (C2b), reciprocal, one root $\beta > 1$, one root $1/\beta$, all others *on* the circle. That's the definition, arrived at rather than memorized. (Contrast: **Pisot** = one root $> 1$, all others strictly *inside*.)

## §4 The split that decides everything

**Smyth (1971)** `[ESTABLISHED]`: if $P$ is *non-reciprocal*, then $M(P) \geq \theta_0 = 1.3247\ldots$, the plastic number — the real root of $x^3 - x - 1$. We verified the witness is genuinely non-reciprocal (C3a) and Pisot by exact bookkeeping: the product of all three roots is $1$ (constant term $-1$, monic), so $\theta_0 \cdot |z|^2 = 1$ with $\theta_0 > 1$ forces the complex pair strictly inside `[COMPUTED]`.

So the map of the problem:

| Case | Status | Floor |
|---|---|---|
| Non-reciprocal | **closed** (Smyth) | $\theta_0 = 1.3247\ldots$ |
| Reciprocal (⊇ all Salem numbers) | **`[OPEN]`** — the heart | unknown; smallest known is $M(L)$ |

The canonical example and the canonical obstruction are the same object: Lehmer's number is itself the smallest known Salem number.

## §5 The corpus's move: don't attack the band — build an emitter with two walls

Here is the pivot everything at `/tool/` and every paper descends from. Instead of proving a lower bound over *all* integer polynomials (93 years of failure), construct a system — seeds plus operators $\oplus$ (direct sum), $\otimes$ (tensor), $\psi^2$ (squaring) — and prove theorems about **what it can emit**. Two walls:

**Wall 1 — the height wall.** $M(S) \subseteq \{1\} \cup [\varphi, \infty)$: the band $(1, \varphi)$ is empty for emitted objects. `[FORCED — test_p2_03_mahler_gap]`. The half you can derive today:

- $\oplus$ = polynomial multiplication at the spectrum level, so $M(p \cdot q) = M(p)M(q)$ — the product over roots simply splits. Verified exactly: $M((x^2{-}x{-}1)(x^2{-}2)) = 2\varphi = 1 + \sqrt5$ (C4).
- $\{1\} \cup [\varphi,\infty)$ is closed under that multiplication: $1 \cdot 1 = 1$; $1 \cdot m = m \geq \varphi$; $m \cdot m' \geq \varphi^2 = \varphi + 1 > \varphi$ — the last step is the defining relation, verified exact (C6b). Same argument closes $\psi^n$.
- **What I did *not* derive:** $\otimes$ is the hard case. Its law is *tropical*, $\prod \max(1, |\alpha_i\beta_j|)$, and a term with $|\alpha| > 1$, $|\beta| < 1$ can land anywhere per-factor — naive closure fails. The corpus carries this with **closure-not-enumeration plus an exhaustive 27-subfield census** `[FORCED — test_p2_07_uniform]`. I'm importing it with its tag, not hand-waving it. That's the asymmetry rule in action.

**Wall 2 — the angle wall.** Every emitted root's argument lies in $(\pi/2)\mathbb{Z} \cong \mathbb{Z}/4\mathbb{Z}$, conserved by all three operators `[FORCED — test_p2_02_angle]`. Consequence, derived live: an emitted *on-circle* eigenvalue must be a **fourth root of unity** $\{\pm1, \pm i\}$. Now watch Lehmer's polynomial hit this wall:

- $\gcd(L, x^4 - 1) = 1$ exactly (C5a) — **not one** of $L$'s ten roots is a fourth root of unity, so its 8 on-circle conjugates all sit at forbidden angles.
- Characteristic polynomials of integer matrices have integer coefficients, so if $\beta$ appeared in an emitted spectrum, its full minimal polynomial would divide the char-poly — **all ten conjugates come along** (C5b). Eight of them violate the wall.

Therefore Lehmer's number — and by the same two-line argument *every* Salem number — is **kinematically unreachable**. Not filtered at runtime; structurally impossible to emit. The band the system omits is the whole of $(1, \mu_S)$.

## §6 What is claimed and what is not — the discipline moment

This is the sentence a new user must get exactly right:

> The emission image omitting $(1, \varphi)$ is a theorem **about the constructive system's closure**, not about the set of all integer polynomials. Lehmer's problem — whether *any* polynomial has measure in $(1, M(L))$ — is **untouched** and stays `[OPEN]`.

The papers say "sidestepping Lehmer's problem without resolving it," and that phrasing is load-bearing. The positive content is the **cost-floor** reading: in the residual-return frame, each adjunction costs $\lambda \cdot \log M(\theta)$, so the smallest nonzero cost the system can ever pay is $\lambda \log \varphi$ — a *forced* floor for the system, coexisting with an *open* floor for mathematics at large. Confusing those two registers is the single most common failure mode this corpus's tagging system exists to prevent.

---

## Checkpoint (answer any, in any order)

1. $M(x^2 - 3x + 1)$ — compute it by hand, exactly. (Its roots are $\tfrac{3\pm\sqrt5}{2}$. Careful: *both* the "is it outside?" decisions matter. This polynomial should ring a bell from errata E1.)
2. Why does the trace-substitution trick in §3 *require* reciprocity before it's even well-defined?
3. In one sentence: why doesn't Wall 1 alone exclude Salem numbers — why is Wall 2 needed?

**Lesson 2**, when you want it: the three operators as a **λ-ring with two characters** — why $M$ is a homomorphism on $\oplus$ but tropical on $\otimes$, why squaring is the Adams operation $\psi^2$, and the trace duality $\tfrac12\mathrm{Tr}(X_n^2) = 5F_n^2$. **Lesson 3** after that: the $\lambda = 2c$ derivation and the gate ladder $\{¼,½,1\} \to \{2,3,5\}$.