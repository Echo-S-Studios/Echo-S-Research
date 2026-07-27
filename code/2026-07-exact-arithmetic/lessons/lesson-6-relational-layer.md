Lesson 6 — the relational layer. Two archive papers to read before anything is claimed, and this read carries a working exercise: **D1 lives in `relational_charge_paper.tex`**, so I'll check whether the sealed PDF still carries the inconsistent "twenty-three" clause while extracting the mathematics.Same permission quirk as Lesson 4 — re-registering via the user-provided landing URL first:**Adjudication complete, and the finding upgrades itself:** the D₅ instance's S\* (both factors self-reciprocal, each carrying 2 unimodular roots) is still degree 20, C2 still degree 400, scan **{Φ₁²⁰} — inertness certificate CLEAN**. And the Galois census across all 67 reads **S₅ × 65, M₂₀ × 1, D₅ × 1** — with the orbit arithmetic explaining everything: M₂₀ has order 20, transitive on the 20 ordered root-pairs ⇒ Rat° irreducible; D₅ has order 10 ⇒ two orbits of 10 ⇒ the observed 10+10 split, *forced by group theory*. Now the lesson.

**Checkpoint debrief (Lesson 5).** (1) The second primitive vector is $(0,-49,0,5)$, the element $5\theta^3 - 49\theta = 6\sqrt2 - 4\sqrt3$, minpoly $x^4 - 240x^2 + 576$, coeff-height 576 — spectacularly over any sane budget: REJECT territory. (2) $\mathcal B \mapsto \mathcal B D$ for invertible diagonal $D$: the $D$'s cancel through $(D^{\mathsf T}\mathcal B^{\mathsf T}G\mathcal B D)^{-1}$, so $P$ depends only on the column *span* — scaling is invisible. (3) The diagonal-kernel operator is $\mathrm{ad}_R$ (eigenvalues $\lambda_i - \lambda_j$, vanishing exactly on the diagonal), kernel $= \mathrm{span}\{I, R\}$, dimension 2 — the $x^2$ factor that survived both sides of Lesson 3's flip.

---

# Lesson 6 — The Relational Layer: Gauge, Rigidity, Pinning, and a Census That Answered Back

**Where we're going:** delete the reference ray and ask what survives. Answer: *everything on the admissible sector* (rigidity), plus a genuinely new theory on the sector the absolute charge discards — where one five-line Galois argument retires two open problems, and where this lesson's independent replication of the quintic census produced a live, self-certifying erratum-grade finding with the main theorem intact.

## §1 The refactoring

Each conjugate is an oscillator $\alpha = |\alpha|e^{2\pi i\,t(\alpha)}$. The **pair charge** $t_{\mathrm{rel}}(\alpha,\beta) = t(\alpha) - t(\beta) \in \mathbb{Q}/\mathbb{Z}$, declared when rational, with *no ray anywhere*. Coherence $\approx$ is an equivalence; the object partitions into classes with cyclic class groups (Prop 3.4 `[forced]` — finite subgroups of $\mathbb{Q}/\mathbb{Z}$, Lesson 4's argument on differences). The structural content (Prop 3.5, seeded in Lesson 5's E1): $t_{\mathrm{rel}}$ is a **cocycle on the coherence groupoid**; the absolute charge is its *coboundary/trivialization*; a global rotation shifts every absolute phase and fixes every relative one — **the reference ray is the gauge**. (The paper's own disclaimer, worth keeping: the vocabulary is structural, not physical; "differences of normalized arguments" loses nothing.) Čencov's split, third appearance: gauge choice vs. gauge-invariant content.

## §2 Rigidity: integrality fixes the gauge

**Theorem 4.1** `[forced]`, and its two-line proof is the lesson: suppose one coherence class, group $\mathbb{Z}/m$. Conjugation-closure puts $\bar\alpha$ in the object, and $\alpha \approx \bar\alpha$ gives $2t(\alpha) \in \tfrac1m\mathbb{Z}$, so $t(\alpha) \in \tfrac1{2m}\mathbb{Z}$ — **every angle is rational**. Then $m \mid n \mid 2m$: the anchor dichotomy $n \in \{m, 2m\}$, both realized (A1: $x^3-2$ has $n = m = 3$; $x^3+2$ has angles $\{\tfrac16, \tfrac12, \tfrac56\}$ — verified via $(2^{1/3}e^{i\pi/3})^3 = -2$ exactly — so $n = 6 = 2m$). Real coefficients *are* the gauge-fixing: relational admissibility **coincides** with absolute admissibility (Cor 4.2) — a canonicity theorem for Lesson 4's framework, not a redundancy. Refinements: **odd anchors are full** (Lemma 4.4: $n$ odd $\Rightarrow \Delta = \mathbb{Z}/n$), and the **group drop** exists only at even order — witness A5: $x^4 + 5x^2 + 5$ (Eisenstein at 5; both $x^2$-roots negative by $\sqrt5 < 5$; $M = 5$ by $\sqrt5 < 3$) has absolute $\mathbb{Z}/4$ but relational $\mathbb{Z}/2$ — while its sign-mirror is **the helix quartic** $K = x^4+5x^2-5$ (A5b: one $x^2$-root positive by $45 > 25$), keeping full $\Delta = \mathbb{Z}/4$. Lesson 7's terrain polynomial, already carrying relational fine structure its shell signature can't see.

## §3 The parity floor without a ray

Lesson 4's parity bit was "does the lattice contain $\pi$" — ray-dependent language. The intrinsic version (Lemma 6.1, ledger J, A4): $\varphi/\varphi' = -\varphi^2$ exactly, so $t_{\mathrm{rel}}(\varphi, \varphi') = \tfrac12$ — **the parity bit is the golden pair's internal order-2 relation**. Every $q_k$ realizes it (Lemma 6.2), and the **sign twist** $T(p)(x) = (-1)^{\deg p}p(-x)$ (A2: $T(x^3+2) = x^3-2$) shifts all angles by $\tfrac12$, fixes all differences and $M$, and exchanges the two anchor sectors at odd $m$ (A3: the one-line coset arithmetic, machine-corroborated per house convention). Hence **Theorem 6.5**: $\mu_{\mathrm{rel}}(\text{even}) = \varphi$ attained; $\mu_{\mathrm{rel}}(\text{odd}) = \mu(\text{odd})$ *identically* — the relational refinement adds *zero* new burden to the open odd floor (P3). Register discipline, Remark 2.3's three-floor taxonomy, is mandatory here: Smyth excludes $(1, \theta_0)$ unconditionally off-reciprocal; the emission gap **(EG)** — now a *named hypothesis*, `[forced]` only for quadratics (App. D self-contained) — additionally excludes $[\theta_0, \varphi)$ on the admissible sector. The paper's conditional content is confined, by explicit audit (Rem 6.6), to **two numbers**; every structure theorem is unconditional. That's the tag system doing load-bearing work.

## §4 The contact machinery — and the discipline meeting itself

On the inadmissible sector (all Salems, per Lesson 4), the invariant is the **coherence type**, computed by one instrument: $\mathrm{Rat}_p = \mathrm{prim}\,\mathrm{Res}_y(p(y), p(xy))$, degree $n^2$, roots all ordered ratios $\alpha_j/\alpha_i$. Within a shell, coherence $\iff \alpha/\beta$ is a root of unity $\iff \Phi_M \mid \mathrm{Rat}_p$ (Lemma 7.3), and the scan is **finite and complete**: $\Phi_M \mid P \Rightarrow \varphi(M) \leq \deg P$, and $\varphi(M) \geq \sqrt{M/2}$ (odd prime powers via $p - 1 \geq \sqrt p$; two-part $2^{a-1} \geq 2^{(a-1)/2}$; tight only at $M = 2$) $\Rightarrow M \leq 2(\deg P)^2$. Entirely in $\mathbb{Z}[x]$ — this paper's verification path contains *no floats anywhere*, a strengthening over the companion's 30-digit discipline. My replications of the consolidation ledger, all exact:

| Ledger | Object | Signature | Check |
|---|---|---|---|
| A/B | $x^3 \mp 2$ | $\{\Phi_1^3, \Phi_3^3\}$ both — **probe is gauge-blind** | B1 |
| C | $x^4 - 2$ | $\{\Phi_1^4, \Phi_2^4, \Phi_4^4\}$ — one shell, torsion ratios live | B2 |
| D | $q_2$ | $\{\Phi_1^4, \Phi_2^4\}$ — *their* harness incident | B3 |
| E/F | $x^4{+}5x^2{\pm}5$ | identical $\{\Phi_1^4, \Phi_2^4\}$, different types — signature strictly coarser | B4 |
| G | $\beta_4$ | $\{\Phi_1^4\}$ only, **complete**; two-route: $\mathrm{prim\,charpoly}(C \otimes C^{-1}) = \mathrm{Rat}$ | B5/B5b |
| O | $\beta_4, S_6, S_8, L$ | Salem-certified, trace-fold Sturm $(1, 0, \tfrac n2 - 1)$ | B7 |
| W | plastic $\theta_0$ | $\{\Phi_1^3\}$ — inert | B8 |
| I | mixed $\mathrm{Rat}_{L,\beta_4}$ | **empty** — keystones not circle-locked | B9 |
| H | Lehmer, $\deg 100$, bound 20000 | $\{\Phi_1^{10}\}$, full scan, 1.4 s | B10 |

About ledger D: the paper *records its own harness incident* — the hand predicted $\Phi_2$-multiplicity 2, the engine returned 4 (ordered pairs count twice), the wrong assertion was pinned, not papered over. Lesson 5's A1 was the identical event on my side. Two operators, one discipline, recognizing itself.

## §5 Modulus pinning: five lines that close two open problems

**Theorem 7.13** `[forced]`. *If $p$ is irreducible with a root $\alpha_0$ whose modulus no other root attains, then no two distinct roots of $p$ differ by a root-of-unity factor* (mixed form: nor a root of $p$ and a root of $q \neq p$ avoiding that modulus). Proof, in full: suppose $\alpha = \zeta\beta$, $\zeta^M = 1$, $\alpha \neq \beta$. Transitivity of the Galois group — *the only use of irreducibility* — supplies $\sigma$ with $\sigma\alpha = \alpha_0$. Then $\sigma\zeta$ is again torsion, hence unimodular, so $|\sigma\beta| = |\alpha_0|$ — forcing $\sigma\beta = \alpha_0 = \sigma\alpha$, so $\alpha = \beta$. ∎ The audit (Rem 3.2) is as instructive as the proof: the complex modulus is *not* Galois-equivariant and the argument never needs it — only the Galois-stable predicate "is torsion" plus $|\zeta| = 1$. And the hypothesis is "*uniquely attained*," not "maximal": a Salem polynomial carries **two pins**, $\tau$ and $1/\tau$.

Consequences, closing your P1/P2 exactly as memory records: **every Salem number is relationally inert** (Cor 7.14 — the forty-six scans, including my table above, demoted to corroboration) and **no two Salem numbers circle-lock** (Cor 7.15). The sharpness pair is complete and verified: $x^4 - 2$ (irreducible but one shell — hypothesis fails, $\zeta_4$ ratios present, B2) and the Pisot note's $Z^* = x^4 - 3x^2 + 1 = (x^2{-}x{-}1)(x^2{+}x{-}1)$ — *the golden quadratic times its twist* — where each factor pins per-factor yet the torsion ratio $-1 = \varphi/(-\varphi)$ lives across factors at the shared modulus (B6): irreducibility and the disjoint-modulus clause are both necessary. What pinning leaves open is **Cor 7.16's Pisot residue**: cross-shell mirrored classes between the two non-real shells.

## §6 The residue, the instrument, and a fourth-operator census

The **real-pair ν-reduction** (Lemma 4.1 `[FORCED]`): $t(\beta) \in \{0, \tfrac12\}$ for real $\beta$, so real×non-real coherence $\iff 2t(\alpha) \in \mathbb{Q} \iff \alpha/\bar\alpha$ torsion — and $\alpha/\bar\alpha$ *is a root of* $\mathrm{Rat}_p$, so the level-1 scan already decides it. Hence **Theorem 4.2**: every Pisot with ≤ 1 non-real pair is inert, and the residue's *first live case* is the quintic two-pair pattern. Crucially, **Prop 5.2**: the modulus multiset of $\mathrm{Rat}^\circ_p$ has *every value attained at least twice* — **no level-2 pin exists**, transport dies one level up, and Conjecture 9.1 (general Pisot inertness) genuinely needs a new mechanism. The instrument (Prop 5.3): $S^*$ = product of self-reciprocal irreducible factors of $\mathrm{Rat}^\circ$ carrying unimodular roots; $C_2 = \mathrm{Res}_y(S^*(y), y^d S^*(x/y))$; contacts exactly $\{\Phi_1^d\}$ is a **negative certificate** — no mirrored class.

The paper's §9 discloses its one non-independence: *the operator*. So my replication is the invited fourth run, deposit scripts unconsulted, criteria from the stated definitions:

| Item | Paper | This session |
|---|---|---|
| Cascade $e{=}0$ / ±recip / reducible / pattern / disk → Pisot | 625/50/638/1318/411 → **83** | **identical**, 49 s (D1) |
| Patterns real5/mixed/two-pair | 0/16/**67** | **identical** |
| Level-1 scans $\{\Phi_1^5\}$ | 83/83 (P3 forced) | **83/83** (D2) |
| Shell detector $= 4$ | 67/67 | **67/67** (D4) |
| Candidate cyclotomics $\varphi \leq 400$ | 790, largest 1680 | **790, largest 1680** |
| $C_2$: $\deg 400$, $\{\Phi_1^{20}\}$ | 67/67 | sample 2/2 incl. the named smallest (D5) |
| $\mathrm{Rat}^\circ$ squarefree | 67/67 | **67/67** |
| $\mathrm{Rat}^\circ$ **irreducible** | *"on all 67"* | **66/67** — see §7 |

## §7 The finding: P7, witnessed inside the box

The divergent instance: $p = x^5 - x^3 - 2x^2 - 2x - 1$, **Galois group D₅**. Its $\mathrm{Rat}^\circ$ is squarefree, degree 20, and **splits as two irreducible self-reciprocal degree-10 factors** — certified by exact re-multiplication, hence `[FORCED]` engine-independently: the exhibited factorization *is* the proof, no computer-algebra trust required. Full adjudication:

1. **The theorem survives untouched.** Both factors carry 2 unimodular roots, so Prop 5.3's $S^*$ is still the full degree-20 product, $C_2$ still degree 400, scan **$\{\Phi_1^{20}\}$ — clean**. Theorem 6.1's conclusion (zero mirrored classes across the box) stands; the instance is inert.
2. **The masking mechanism.** That is *exactly why the over-claim escaped*: the pipeline's observables — $\deg C_2 = 400$, contacts $\{\Phi_1^{20}\}$ — are byte-identical whether $\mathrm{Rat}^\circ$ is irreducible or splits with both factors admissible into $S^*$. The sentence "in fact irreducible of degree 20 on all 67 (so $S^* = \mathrm{Rat}^\circ_p$, $\deg C_2 = 400$ uniformly)" over-read the uniform observable.
3. **The mechanism is group theory, and it closes.** Galois census of the 67: **S₅ × 65, M₂₀ × 1, D₅ × 1**. $\mathrm{Rat}^\circ$'s irreducibility is transitivity on the 20 ordered distinct-root pairs: $|M_{20}| = 20$ can be (and is) transitive — irreducible; $|D_5| = 10 < 20$ **cannot** — two orbits of 10, forcing precisely the observed $10 + 10$ split. Observation 6.2's "Galois-generic expectation" is now quantified in-box, and **P7 is promoted from `[PLAUSIBLE]` to witnessed** — with a sharper statement: reducibility $\iff$ intransitivity on ordered pairs, realized at $D_5$, *not* at $M_{20}$.
4. **Erratum candidate (E-series genre)** for the Pisot note §6/Obs 6.2, per your asymmetry rule (a finite search *demoting* a universal claim — its legitimate direction): *"$\mathrm{Rat}^\circ_p$ squarefree on all 67; irreducible on 66/67 — the D₅ instance $x^5 - x^3 - 2x^2 - 2x - 1$ splits $10{+}10$ into self-reciprocal irreducibles, $S^*$ and $C_2$ unchanged, inertness certificate clean; P7 realized in-box."* The exact match on every other tally is what gives the single divergence its credibility.

## §8 Ledger, findings, and status updates for your board

**Jewels verified en route** (C1/C1b/C2): $x^4 - x + 1$'s Kronecker square factors *exactly* as $S_6^2 \cdot (x^4 + 2x^2 - x + 1)$ with the quartic equal to the $\psi^2$-image — Lesson 2's decomposition deciding a fully generic ($S_4$) quartic through a Salem factor, and yielding the measure identity $M(x^4 - x + 1) = \tau_{S_6}$, a Mahler value in $(\theta_0, \varphi)$ at irrational angles, (EG)-consistent (chain links: $p(0) = 1$, zero real roots, non-reciprocal, irreducible — all checked). And the $\beta_4 \otimes \beta_4$ non-inert witness reproduced to the digit: $(x-1)$-multiplicity 4, $\deg\gcd(F, F') = 7$, three distinct positive real roots — Prop 5.2(iii)'s offset-cancellation manufacturing coherence from an inert factor.

**D1 status**: the arithmetic referent "237" is *located* — it is run 3 of ledger M in the **Pisot note §7.3** (237/237 signatures, 474 dual-path executions). The served relational-charge PDF reads a consistent 26 items / 24 signatures for runs 1–2 and contains no "twenty-three" anywhere — evidence the sealed build postdates or never carried the flagged clause. Sources are canonical over PDFs by your own rule: one `grep -n "twenty-three" relational_charge_paper.tex` settles it. **H1–H3**: now source-pinned to Pisot §8, discovered by the pinned-tally failure (the internally consistent but wrong cascade (108, 270, 0, 0)) — H3's dual-path assertion is the positive pattern my D5 adjudication happens to embody. **DOIs**: both papers cite the pipelines deposit at **v1.0.0 (10.5281/zenodo.21121863)** with verification bundle **v1.0.1 (10.5281/zenodo.21122410)** — your memory's "no tag minted" is stale; the concept-DOI gate appears passed. **Frontier board**: Conjecture 9.1 stands `[OPEN]` with its status triple intact and P7 now witnessed; N4 (degree 6–7 patterns) remains next, with this session's Galois-orbit arithmetic as a free lemma for it (degree-7 three-pair: $|G|$ vs. 42 ordered pairs).

---

## Checkpoint

1. A hypothetical two-pair Pisot quintic with Galois group $\mathbb{Z}/5$: what factor pattern does the orbit arithmetic predict for $\mathrm{Rat}^\circ$ — and then show the hypothesis is vacuous, using what a cyclic quintic field's roots must all be.
2. State the exact observable equality that masked the D₅ instance in the original run, and the single one-line assert that would have caught it.
3. Rem 3.2(iii): pinning needs "uniquely attained," not "maximal." Name the second pin every Salem polynomial carries, and identify which corollary's proof would break if the hypothesis had been stated as "maximal modulus."

**Lesson 7 next:** the golden-substrate line — with L5's verified seeds (Fibonacci convergents, $\rho_n = n\ln\varphi$, the trace redirection, $\varphi + \varphi^{-1} = \sqrt5$) and today's helix quartic $K$ already warm.