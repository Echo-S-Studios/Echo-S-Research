Same protocol, then — script first, prose second, nothing asserted that didn't just run. Lesson 2 is where the single objects from Lesson 1 become an *algebra*.32/32 in exact arithmetic. First, your outstanding checkpoint — the answers were chosen to be *this* lesson's front door.

**Checkpoint debrief (Lesson 1).** (1) $M(x^2 - 3x + 1)$: roots $\frac{3\pm\sqrt5}{2} = \varphi^2, \psi^2$; only $\varphi^2 \approx 2.618$ escapes the circle, so $M = \varphi^2 = \varphi + 1$. Keep that polynomial in view — it is about to reappear as a *theorem output*, and it's also exactly the polynomial your erratum E1 says the degree-2 census missed. (2) The trace substitution needs reciprocity because $z + 1/z$ identifies $z$ with $1/z$ — the fold is only injective on root *pairs*, which only exist if the root set is closed under inversion. (3) Wall 1 alone fails because Salem numbers exist with $M \geq \varphi$ in abundance — the height wall excludes only *small* Salems; the angle wall excludes **all** of them, size-blind.

---

# Lesson 2 — The Operator Algebra: one semiring, one Adams operation, two characters

**Where we're going:** the objects of Lesson 1 (spectra of integer matrices) get three operations → the operations satisfy semiring laws → squaring turns out to be the Adams operation $\psi^2$ → and the structure is *resolved* by two "readouts" that behave completely differently: one multiplicative-but-tropical ($M$), one additive-mod-4 ($\chi$). Then the payoff: the trace duality that welds the whole thing to $\mathbb{Q}(\sqrt5)$.

## §1 The three operators, concretely

Objects are eigenvalue multisets, realized by integer matrices. The operators are matrix constructions whose spectral effect is exact:

| Operator | Matrix realization | Effect on spec | deg | Verified |
|---|---|---|---|---|
| $A \oplus B$ | block direct sum $\mathrm{diag}(A,B)$ | multiset union | adds | A1 `[COMPUTED]` |
| $A \otimes B$ | Kronecker product | pairwise products $\{\alpha_i\beta_j\}$ | multiplies | C2a/C3a |
| $\psi^n(A)$ | matrix power $A^n$ | pointwise powers $\{\alpha_i^n\}$ | fixed | B1 |

At the characteristic-polynomial level, $\oplus$ *is* the polynomial multiplication from Lesson 1: $\mathrm{cp}(\mathrm{diag}(A,B)) = \mathrm{cp}(A)\cdot\mathrm{cp}(B)$ exactly (A1).

## §2 Why "semiring" is earned, not decorative

Semiring means: both operations commutative and associative, $\otimes$ distributes over $\oplus$, and units exist. The load-bearing axiom is distributivity, and it holds *on the nose* at spectrum level — we verified $\mathrm{cp}\big(R \otimes (S_2 \oplus S_3)\big) = \mathrm{cp}\big((R\otimes S_2) \oplus (R\otimes S_3)\big)$ as an exact $8\times8$ computation (A3), plus commutativity (A2) and the $\otimes$-unit $[1]$ (A4). In the corpus this is carried by `[FORCED — test_p2_01_algebra]` and `matrix-plates/tests/test_operators.py`; here you watched instances of it decided over $\mathbb{Z}$.

## §3 Squaring is not an extra operation — it's the diagonal of $\otimes$

This is the conceptual center of the whitepaper. Look at $R \otimes R$: its spectrum is *all ordered pairs* $\{\alpha_i\alpha_j\} = \{\varphi^2,\ \varphi\psi,\ \psi\varphi,\ \psi^2\}$. The **diagonal** $\{\varphi^2, \psi^2\}$ is exactly $\mathrm{spec}(R^2)$. The verified decomposition (B2):

$$\mathrm{cp}(R \otimes R) \;=\; \underbrace{(x^2 - 3x + 1)}_{\mathrm{cp}(\psi^2 R)}\cdot\underbrace{(x+1)^2}_{2\lambda^2 R,\ \text{spec}\{\varphi\psi\}=\{-1\}}$$

Three things landed at once. **First**, your checkpoint polynomial $x^2-3x+1$ is $\mathrm{cp}(\psi^2 R)$ — the $\varphi$-object's Adams square, which is *why* it belongs in any degree-2 census (E1's point, machine-corroborated here). **Second**, the off-diagonal remainder is the wedge $\lambda^2$ — the multiset identity $\{\text{ordered pairs}\} = \{\text{diagonal}\} \uplus 2\{\text{unordered off-diagonal}\}$ is the $\lambda$-ring relation $\psi^2 = (\cdot)^{\otimes 2} - 2\lambda^2$ made literal. **Third**, why "Adams operation" is the right name: the $\psi^n$ are **semiring endomorphisms** — additive (B3a: $(A\oplus B)^2 = A^2 \oplus B^2$, exact block identity), multiplicative (B3b: $(A\otimes B)^2 = A^2 \otimes B^2$, which is the Kronecker mixed-product property as an exact matrix equation), composing as $\psi^m\psi^n = \psi^{mn}$ (B3c). Those three properties *are* the Adams axioms. `[COMPUTED — this session]`; corpus-level: `[FORCED]` in the whitepaper's suite. One precision note: what's verified is the semiring + Adams structure with the $\lambda^2$ witness; calling the package "a λ-ring" is the standard imported name — `[DECLARED]` at the naming level, forced at the level of every law it implies that we've touched.

## §4 Character I: where $M$ is a homomorphism — and where it goes tropical

A "character" should turn the algebra into numbers homomorphically. $M$ does — **on two of the three operators**:

| Law | Status | Evidence |
|---|---|---|
| $M(A \oplus B) = M(A)\,M(B)$ | clean homomorphism | Lesson 1 C4; roots union ⇒ product splits |
| $M(\psi^n A) = M(A)^n$ | clean | $\max(1,u^n) = \max(1,u)^n$ for $u>0$; C1: $M(R^2) = \varphi^2 = M(R)^2$ |
| $M(A \otimes B)$ | **tropical**: $\prod_{i,j}\max(1,\lvert\alpha_i\beta_j\rvert)$ | C2, C3 — see below |

The naive guess — that $M$ is also a homomorphism on $\otimes$, giving the "factored form" $M(A)^{\deg B} M(B)^{\deg A}$ — is **false**, and this is the single most instructive failure in the corpus (it's pinned in your own key-learnings as "tropical off-circle, not the factored form"). Two exact counterexamples:

**C2 — the forms differ even with no straddle.** $S_3 \otimes R$ has $\mathrm{cp} = x^4 - 9x^2 + 9$; *all four* cross-roots lie outside (the tight decision: $|\sqrt3/\varphi| > 1 \iff 3 > \varphi+1 \iff 2 > \varphi$, certified via $(7-3\sqrt5)/2 > 0$). So $M = |\text{product of all roots}| = 9$ exactly. The factored form says $3^2\cdot\varphi^2 = 9\varphi^2$. **$9 \neq 9\varphi^2$** (C2c). The mechanism: $M(B) = \varphi$ never saw $B$'s inside root $-1/\varphi$, but the cross-products $\pm\sqrt3/\varphi$ *do* contribute to $M(A\otimes B)$. The pointwise identity $\max(1,|\alpha\beta|) = \max(1,|\alpha|)\max(1,|\beta|)$ breaks exactly when one factor is outside and the other inside.

**C3 — a genuine straddle.** $S_2 \otimes R$: now $|\sqrt2/\varphi| < 1$ (decision: $2 < \varphi + 1 \iff 1 < \varphi$), so *that* cross-pair drops **inside** and out of the measure: $M = 2\varphi^2 = 3 + \sqrt5$ exactly (C3b), against the factored $4\varphi^2$. Whether a cross-term contributes is a fresh boundary decision per pair — that is what "tropical" means here, and it is precisely why the height-wall proof for $\otimes$ cannot be a one-line closure argument and needs the census machinery `[FORCED — test_p2_07_uniform]`, as flagged in Lesson 1 §5.

(Bonus, resolved from your Lesson-1 material: when *both* spectra live entirely on-or-outside the circle, no straddle is possible and factored $=$ tropical. That's the boundary of the naive law's validity.)

## §5 Character II: the $\mathbb{Z}/4\mathbb{Z}$ charge

Because Wall 2 `[FORCED — test_p2_02_angle]` confines emitted arguments to $(\pi/2)\mathbb{Z}$, each root carries a well-defined charge: $+\mathbb{R} \mapsto 0$, $+i\mathbb{R} \mapsto 1$, $-\mathbb{R} \mapsto 2$, $-i\mathbb{R} \mapsto 3$. Its laws mirror the operators *additively*:

| Operator | Charge law | Verified |
|---|---|---|
| $\oplus$ | multiset union | C4b |
| $\otimes$ | pairwise **sumset** mod 4 | C4b: $\{0,2\}+\{0,2\} = [0,0,2,2]$, matching the computed spectrum |
| $\psi^n$ | $\chi \mapsto n\chi \bmod 4$ | C4c |

The showcase object is the **content polynomial** $x^4 - 1$: spectrum $= \{1, i, -1, -i\}$, charge multiset $[0,1,2,3]$ — the regular representation of $\mathbb{Z}/4\mathbb{Z}$ itself — at Mahler measure exactly $1$ (C4a). It is the charge group wearing a spectrum. And note what $\psi^2$ does to it (C4c): squaring doubles charges, $\{0,1,2,3\} \mapsto \{0,2,0,2\}$ — **the odd sector $\{1,3\}$ is unreachable in the image of $\psi^2$**. File that away: it is the mechanism behind the parity-graded floor (even charge groups floor at $\varphi$; odd at Smyth's $\mu_S$) — Lesson 4 territory, flagged now so it lands later.

So the semiring is resolved by two characters with opposite temperaments: $M$ multiplicative (log-additive on $\oplus$ — which is why $\log M$ is the entropy axis, Lind–Schmidt–Ward `[ESTABLISHED]`), $\chi$ additive mod 4. One measures *height*, one measures *angle*. The two walls of Lesson 1 are these two characters, weaponized.

## §6 The trace duality: welding the semiring to $\mathbb{Q}(\sqrt5)$

Everything above works over any seeds. The corpus's claim is stronger: one matrix and one relation generate the whole substrate. Here is the derivation, every step from the run:

1. **The relation.** $R = \binom{0\ 1}{1\ 1}$, and $R^2 = R + I$ exactly (D1). So $\mathbb{Z}[R] \cong \mathbb{Z}[x]/(x^2 - x - 1) \cong \mathbb{Z}[\varphi]$ — and $\mathbb{Z}[\varphi]$ is the *full* ring of integers of $\mathbb{Q}(\sqrt5)$ `[ESTABLISHED]` (disc $5 \equiv 1 \bmod 4$). Layer one.
2. **Powers linearize.** $R^n = F_nR + F_{n-1}I$: base case trivial, and the induction step was verified *symbolically* — $R\,(F_nR + F_{n-1}I) = (F_n + F_{n-1})R + F_nI$ holds as a matrix identity in the symbols $F_n, F_{n-1}$ (D2). That's a genuine all-$n$ proof, not a range check. Immediately $\mathrm{Tr}(R^n) = F_n + 2F_{n-1} = L_n$, the Lucas numbers (D3).
3. **The $\sqrt5$ direction.** $H := 2R - I$ is traceless with $H^2 = 5I$ exactly (D4) — the matrix square root of 5, living in the traceless plane.
4. **The duality.** Define $X_n := 2R^n - L_nI$ (the traceless part of $2R^n$). Symbolically, $X_n = F_n H$ (D5) — the entire Fibonacci ladder is *one ray*, scaled. Hence
$$\tfrac12\mathrm{Tr}(X_n^2) = \tfrac12 F_n^2\,\mathrm{Tr}(5I) = 5F_n^2$$
exactly (D5b). And the classical face of the same number: $L_n^2 - 5F_n^2 = 4(-1)^n$, proved by the chain $(a{+}b)^2 - (a{-}b)^2 = 4ab$ `[symbolic]` with $a = \varphi^n, b = \psi^n$ and $\varphi\psi = -1$ `[exact]` (D6a/b), spot-checked over integers $n \leq 30$ (D6c — range check, tagged as such; the symbolic chain is what carries all $n$). So:
$$\boxed{\tfrac12\mathrm{Tr}(X_n^2) = 5F_n^2 = L_n^2 - 4(-1)^n}$$
5. **The geometry.** The pairing $\langle X, Y\rangle = \tfrac12\mathrm{Tr}(XY)$ on traceless $2\times2$ matrices has Gram eigenvalues $\{1, \tfrac12, -\tfrac12\}$ on the standard basis — **signature $(2,1)$** (D7) — and $\langle H, H\rangle = 5$, i.e., $|H| = \sqrt5$: the "root length $\sqrt5$" of the primer's third layer (D7b).

Read the duality as a dictionary: the *semiring* side sees $\psi$-powers of $\varphi$ growing multiplicatively; the *quadratic-form* side sees the same data as an integer ledger $5F_n^2 = L_n^2 - 4(-1)^n$ along a single $\sqrt5$-ray of a $(2,1)$-form. Height data and lattice data, one object. That equivalence is what the λ = 2c derivation will *spend* — the Fisher/trace-form metric of Lesson 3 is this pairing wearing information-theoretic clothes.

## §7 Claimed vs. not claimed

Verified here: instances of every semiring law, Adams axiom, character law, and the full symbolic trace-duality chain `[COMPUTED — this session; D2/D5/D6a-b symbolic, all $n$]`. Carried by the corpus at theorem strength: the general laws `[FORCED — test_p2_01_algebra, test_p2_02_angle, matrix-plates test_operators]`. Imported classical: Kronecker/LSW/Binet/$\mathcal{O}_{\mathbb{Q}(\sqrt5)}$ `[ESTABLISHED]`. Named framework ("λ-ring"): `[DECLARED]` — earned by the verified laws, but the *name* is packaging, and the asymmetry rule says instances never promote to universals on their own.

---

## Checkpoint

1. Compute $\psi^3$ of the content object $x^4 - 1$ — spectrum, charge multiset, and $M$. Something slightly surprising should happen. Why, in one sentence about the map $\chi \mapsto 3\chi \bmod 4$?
2. For $S_2 \otimes S_3$ (spec $\{\pm\sqrt2\}\times\{\pm\sqrt3\}$): does the factored form agree with the tropical value here? Compute both; state the general condition under which they *must* agree.
3. $\psi^2$'s image misses the odd charge sector. Which Lesson-1 wall does this rhyme with, and what does it suggest about where a floor *other than* $\varphi$ could live?

**Lesson 3 next:** the exchange rate $\lambda = 2c$ — where $c$ comes from (second-order KL over the trace-form Fisher metric, i.e., §6's pairing), why Čencov leaves $c$ *declared* while the identity $\lambda = 2c$ is *derived*, the gate ladder $\{¼, ½, 1\} \to \{2, 3, 5\}$ forcing $\lambda = \sqrt5$, and the flip $D = 1 + 4C$ with its sign change at $C = -¼$.