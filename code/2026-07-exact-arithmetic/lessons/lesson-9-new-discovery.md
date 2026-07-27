Lesson 9 — the domain jump, then the filing lab. Two unread papers in a domain the spectral machinery can't help with, so I'll read both fully before any claim. Re-registering first:Both papers read in full. Now the filing-lab evidence — D1 first, since it's a one-search adjudication against the canonical source:Now E1–E3 and the D1 clause against canonical sources — two targeted searches:One more evidence step — the canonical TeX for D1, via the repo tree the landing page links:Both papers read, filing evidence gathered (with one blocked path to report). Now the verification suite — Part A/B exact witnesses, Part C the orbit criterion:The symbolic-exponent integral is tripping sympy's meijerg path — reformulating that block with the exponent as a single symbol:**32/32.** And Part C produced something better than a filing: the D₅ anomaly turned into a theorem. First the debrief.

**Checkpoint debrief (Lesson 8).** (1) For Salem $R$, factor $R(\varphi) = \prod(\varphi - \alpha_i)$: each circle pair contributes $|\varphi - \lambda|^2 > 0$, and $\varphi - 1/\beta > 0$ since $1/\beta < 1 < \varphi$ — so *every factor but one is forced positive* and $\mathrm{sign}\,R(\varphi) = \mathrm{sign}(\varphi - \beta)$. A second real root above 1 would break the single-sign read and could mis-assign `FORCED_ABOVE_FLOOR` vs `INVALID_CLOSURE`; the Sturm count `above == 1` (exactly one trace-down root in $(2,\infty)$) already excludes it, backed by `n_real == m`. (2) $r_1 = 8$: $g^{-1}cg = c$ requires $g$ in the centralizer of the reflection, order $4 \times 4 = 16$, giving $16/2 = 8$ cosets — so $K$ is $(16, 8, 4)$, i.e. $(2k,k)$ with $k = 4$, confirmed as C1's last row today. (3) The eight mixed products $(\pm\rho)(\pm i\beta) = \pm i\rho\beta$ with $\rho^2\beta^2 = \frac{3\sqrt5-5}{2}\cdot\frac{5+3\sqrt5}{2} = \frac{45-25}{4} = 5$ — so $\rho\beta = \sqrt5$ exactly and they *are* the roots of $(x^2+5)^4$; the grow factors split $\{\beta^2, -\rho^2\}$ and $\{\rho^2, -\beta^2\}$. It's $x^2+5$ rather than $x^2+1$ because charge is the *argument* $\pm\pi/2$ while the modulus stays free — $x^4-1$ records where the charges *would* sit on the circle; the algebra realizes 1 and 3 at modulus $\sqrt5$.

---

# Lesson 9 — The Domain Jump, and the Filing Lab

**Where we're going:** the course's meta-claim has been *the method is the subject*. Today it gets tested where no Mahler measure can help — a second-order elliptic PDE and a biology-facing protocol — and then the course converts its own findings into filed work. Verdict up front: the method transfers completely, one structural identification (gatekeeper ≡ certificate) is exact rather than analogical, and Part C closes E1, sharpens the D₅ finding into a theorem, and produces one new display-level defect.

## §1 Part A: exactness in an analytic domain

The setting is the mixed-boundary torsion problem $-\nabla\!\cdot(\mathcal A\nabla u) = q$, Dirichlet on a supply portion $\Sigma$, natural (no-flux) on the remainder $\Gamma$, with the discipline stated in one word: **ratio-free** — every constant depends on $P_{\max}$ alone, never on the ellipticity ratio, volume, or regularity moduli. That is the same move as `[FORCED]`-vs-`[COMPUTED]`: a declared budget of what a constant may depend on, enforced throughout.

| Result | Content | Check |
|---|---|---|
| Thm 2.1 | exact 1-D law $u = q\int_0^x \frac{L-s}{A}\,ds$; far-end $\geq qL^2/2P_{\max}$, **equality iff $A = P_{\max}$ a.e.** | A1 (equality case exact), A2 (a strict instance: $A = P/(1{+}s/L)$ gives $\tfrac23 qL^2/P > \tfrac12$) |
| Thm 3.3 | ball-mean bound with $n(n{+}2)$, via same-operator comparison + rigidity anti-monotonicity | A4 (radial mean $= 2r^2/(n{+}2)$ symbolic; $15$ at $n{=}3$), A5 |
| Prop 4.1 / Thm 4.2 | $\mathcal A\nabla f = 2(z{-}x_0)$, $\nabla\!\cdot(\mathcal A\nabla f) = 2n$, $f \geq |z{-}x_0|^2/P_{\max}$; star bound $u(x_0) \geq \frac{q}{2n}\inf_\Sigma f$ | A7 (anisotropic $\mathcal A$, symbolic in 3 vars), A9 (eigenbasis sum termwise) |
| Prop 4.3(iii) | equality on ellipsoids: $u = \frac{q}{2n}(r^2 - f)$ | A8 |

Two things are worth naming. The **method of proof is the guard's method**: Lemma 3.1's comparison and every bound below it run by testing the weak form against a Stampacchia truncation $(w-m)^+$ and concluding "$\nabla\varphi \equiv 0$, and a constant with zero trace on $\Sigma$ is zero." That is a *certificate pattern* — construct an explicit admissible witness, test, conclude — structurally identical to Lesson 8's `validate_closure`: build the trace-down, test the straddle, read one sign. And **Corollary 4.4 is the honest half**: distance to the *full* boundary always controls $u$ from below. What fails is distance to the *supply* portion, and only past a sharply located hypothesis.

## §2 The closure: an exactly solvable family, a rational witness

The falsification is built, not searched. On the truncated horn $g(x) = \varepsilon h(x/h)^{k/2}$, one explicit field $v = a(h^2 - x^2) - b|y|^2$ with $a = q/(P(2 + k(n{-}1)))$, $b = ak/2$ does three things at once, all verified:

- $-P\Delta v = q$ **identically** (A10) — the coefficient $a$ is precisely what makes the constant come out;
- the lateral-wall flux vanishes **exactly** (A11) — because $xg' = \tfrac k2 g$ for the power profile *and* $b = ak/2$; two matched choices, one cancellation;
- the end-face flux has the right sign and $\sup v = ah^2(1-\delta^2)$ (A12).

Comparison then gives $0 \le u \le v + b(\varepsilon h)^2$, hence $\frac{P_{\max}\|u\|_\infty}{qR^2} \le \beta = \frac{1-\delta^2 + k\varepsilon^2/2}{(2+k(n-1))(1-\delta)^2}$, and the witness lands (A13): at $n=3$, $k=4$, $\delta = \tfrac3{20}$, $\varepsilon = \tfrac15$,

$$\beta = \frac{423/400}{2890/400} = \frac{423}{2890} < \frac16, \qquad \frac16 - \frac{423}{2890} = \frac{88}{4335}.$$

A bounded Lipschitz domain, constant scalar coefficient, flat Dirichlet mouth — every hypothesis of the conjecture satisfied, the constant beaten by 12.2%, **in exact rational arithmetic with a named deficit**. And the infimum is 0 (A14: $\varepsilon^2 = \delta = 1/k$ drives $\beta \to 0$), so all three metric variants fail together — the nesting $d_E \le d_\Omega \le \sqrt{P_{\max}}d_{\mathcal A}$ collapsing on this family (A19). This is *precisely* Lesson 7's move — an exactly solvable family closing a conjecture by driving an infimum, the redirections' $\sqrt5$ limit in another domain — and Lesson 6's asymmetry rule in its legitimate direction: a finite explicit witness **demoting** a universal claim.

## §3 The cone, and the mechanism

The dichotomy is sharp, and the transition is an identity (A15–A17): $2 + k(n{-}1) = 2n$ **identically at $k = 2$ in every dimension**, so the straight cone is pinched at exactly the conjectured constant $1/(2n)$; $k \le 2$ satisfies the star condition at every axis point and the bound holds; $k > 2$ makes the wall integrand negative for all $x > x_*k/(k{-}2)$ — solved exactly — and failure becomes unbounded. The paper's scope note is exemplary: **necessity across all geometries is not asserted**. A sharp *sufficient* hypothesis, with the failure mode located, and no overclaim.

The mechanism (Rem 7.1, A18) is the lesson's intellectual payoff: there are **two different one-dimensional balances**. The conductivity-graded one (Thm 2.1) has sharp constant $1/2 = 1/(2n)|_{n=1}$; the volume-graded one — sink $q$ per unit *volume* in a tube of cross-section $A_{\text{cross}}$ — obeys $V(s)/A_{\text{cross}}(s) = s/(m{+}1)$ with $m = k(n{-}1)/2$, giving $\frac qP\int_0^h \frac{V}{A_{\text{cross}}} = \frac{qh^2}{P(2+k(n-1))} = v(0,0)$ **exactly**. The conjecture extrapolated the first balance: it kept the distance and **discarded the volume grading $V/A$**. The horns are the geometries where the two functionals separate without bound. Internal consistency is also exact (A20): on a horn every inscribed ball has radius $\le \varepsilon h$, so the mean bound asserts only $1/375$ of $qh^2/P$ against an actual scale of $1/10$ — the two results *never touch*.

## §4 Part B: the unit split, and one object read two ways

The ATR paper's first act is a **units audit**, and it is the exact-arithmetic instinct applied to dimensional bookkeeping. Two tube integrals exist: geometric resistance $R_{\text{geom}} = \int \frac{V}{PA_{\text{cross}}}$ (geometry only) and predicted concentration drop $\Delta C_{\text{tube}} = \int \frac{Q_\downarrow}{PA_{\text{cross}}}$ (demand already inside). They coincide only under uniform demand — verified both ways: $\Delta C_{\text{tube}} = qR_{\text{geom}}$ exactly when $Q_\downarrow = qV$ (B1), and **genuinely failing otherwise** (B2 — density $2q_0s/h$ with the same mean gives $q_0h^2/3$ against $q_0h^2/2$). Multiplying the demand-weighted object by $q$ again would double-count the load; naming the two objects separately removes the ambiguity at the source. Same disease, same cure as Lesson 6's absolute-vs-relational charge.

A structural observation the papers state in two places without joining (B3): the "tube quantity" falling and the "distance-only failure factor" climbing are **reciprocal readings of one number**, $1/(2+k(n{-}1))$ vs $2+k(n{-}1)$ — $\tfrac16 \leftrightarrow 6$ at $k=2$, $\tfrac1{22} \leftrightarrow 22$ at $k=10$. And the monotonicity is *forced*, not sampled (B4): $\partial_k[1/(2+k(n{-}1))] = -(n{-}1)/(2+k(n{-}1))^2 < 0$ for all $n \geq 2$. The maintenance number's correction is also exactly one-directional (B5): dividing by the margin $C_\Sigma - C_{\text{crit}}$ rather than $C_\Sigma$ can only *tighten* the viability call.

## §5 The gatekeeper is the certificate

This is the identification the scoping predicted, and it holds exactly rather than by analogy:

| Lesson 8, spectral | Lesson 9, empirical |
|---|---|
| the **theorem** (angle invariant ⇒ no Salem constructed) | the **theorem** (horn family ⇒ distance gives no lower bound) |
| the **guard**: an exact runtime verdict on each produced object | the **gatekeeper**: solve (1) on the actual geometry, ask whether $R_{\text{ATR}}$ tracks the exact $u$ |
| `FORCED` / `FORCED_ABOVE_FLOOR` / `INVALID_CLOSURE` — a scoped verdict set, *not* the strong theorem | "testing a computational surrogate inspired by the exact horn functional" — the licensed register, *never* "we transported the theorem into tissue" |
| the scope disclaimer: floor certificate ≠ no-Salem theorem | the scope disclaimer: exact on its stated geometry; a surrogate until the gate is passed |
| GC ledger: guard-completeness `[OPEN]` | §13's three-way split: forced / declared constructions / open |

Both are the same epistemic object: **a verdict plus an exact witness plus a stated scope**, standing between a theorem and a claim. Three further features are worth importing back into the spectral corpus. The **ablation ladder** — six rungs adding one information source at a time, with the decomposition itself reported as a result independent of whether the surrogate is ever adopted — is a *tag-resolution instrument*: it answers "which ingredient carries the signal" (here: the restriction × demand interaction, not distance, restriction, or demand alone). The **"why a high correlation is not enough"** section is the same lesson as Lesson 8's derogatory witness: a surrogate at Pearson 0.995 still inverts 20 matched pairs, exactly as equal characteristic polynomials still fail to imply similarity — *the aggregate invariant is not the complete invariant*. And the FEP demotion, with the honest note that muscle may break both models, is `[POSITED]` discipline under empirical risk. Note carefully: the numerical claims ($\rho_u = 1.83$; the $1.66\times$/$1.63\times$/$1.6\%$ triple; the correlations $0.00/{-}0.85/{+}0.23/{+}0.92$; collinearity $0.369$) are PDE-solver outputs — **`[COMPUTED, harness]`, not verified here**, since replicating them requires the solver, not exact arithmetic. I did not check them and do not assert them.

## §6 The rounding audit: one pass, one fail

Two display-level values got interval audits, and the pair is instructive.

**Passes** (B6): the eikonal tail asymmetry "$3.95\times$" is attainable from unrounded percentiles rounding to $5.4°$ and $21.2°$ — window $[3.881, 3.972]$ contains 3.95. Internally consistent.

**Fails** (B7) — **new finding, display-level**: the tube quantity at $k = 10$, $n = 3$ is exactly $1/22 = 0.045454\ldots$, which rounds to **0.045**; the ATR paper displays **0.046** in three places (§2.1 text, Table 1, §14 C4). Mechanism identified: $0.167 \times \tfrac6{22} = 0.045545 \to 0.046$ — the value was scaled from the *already-rounded* 0.167 rather than computed from $1/22$. Genre and weight: identical to Lesson 7's near-miss guard — no decision rests on it, the paper's own claims (monotone decrease, negative $k$-derivative) are exact and unaffected, and the correct display is a one-character fix. Filed below as **F1**.

## §7 Part C: the filing lab

### E1 — adjudicated, and closed on the pipelines side

The referent was in the compendium, not the TeX: the subfield census table. My Lesson 8 group-theoretic census, extended today to the **full signature distribution**, matches the corpus table row for row (C1):

| $(\deg, r_1, r_2)$ | corpus | this session (subgroup lattice) |
|---|---|---|
| $(1,1,0)$ | 1 | 1 |
| $(2,2,0)$ | **7** | **7** |
| $(4,2,1)$ | 4 | 4 |
| $(4,4,0)$ | 7 | 7 |
| $(8,4,2)$ | 6 | 6 |
| $(8,8,0)$ | 1 | 1 |
| $(16,8,4)$ | 1 | 1 |
| total | 27 | 27 |

— and the group itself is confirmed independently (C2): element orders $\{1{:}1, 2{:}23, 4{:}8\}$, matching the corpus's PARI identification `SmallGroup(32,46)` $= C_2 \times C_2 \times D_4$. **Verdict**: the degree-2 row is 7, correctly; the seven quadratic subfields are $\mathbb{Q}(\sqrt d)$, $d \in \{2,3,5,6,10,15,30\}$, and $x^2 - 3x + 1$ (roots $\varphi^{\pm2}$, discriminant 5) generates $\mathbb{Q}(\sqrt5)$ — **exactly the member whose omission yields 6** (C3). Two notes: the served pipelines compendium already shows the corrected 7, so E1 is not an outstanding defect there; and the omitted row was the *golden* one — the census had dropped the field generated by its own keystone.

### D1, E2, E3 — status, with an access boundary reported

I attempted the canonical grep. `project_knowledge_search` does not surface `relational_charge_paper.tex`, `163`/`160`, or the degree-7 Pisot display — consistent with those living in the **archive** repo, while project knowledge indexes the **pipelines** repo; and the direct route (`github.com/.../tree/main/papers`) returned `ROBOTS_DISALLOWED`. So:

| Item | Evidence held | Remaining action |
|---|---|---|
| **D1** | served PDF §9 reads "all **twenty-six** items… **twenty-four** contact signatures"; no "twenty-three" anywhere; the 237 is Pisot §7.3 run-3 | one local `grep -n "twenty-three" papers/relational_charge_paper.tex` — expected empty |
| **E2, E3** | target numerals absent from every indexed source, both sessions | same local grep + the v1.0.0 deposit diff |
| **E1** | **closed** on the pipelines side, count independently confirmed | optional: confirm the deposit TeX also reads 7 |

### The D₅ finding, promoted: an iff-criterion

Lesson 6 explained the $10{+}10$ split by orbit arithmetic. Today it becomes a clean statement, verified:

> **Lemma (orbit criterion).** If $\mathrm{Rat}^\circ_p$ is squarefree, then $\mathrm{Rat}^\circ_p$ is **irreducible over $\mathbb{Q}$ iff $\mathrm{Gal}(p)$ acts 2-transitively on the roots of $p$** — irreducibility is transitivity on the $n(n{-}1)$ ordered distinct root pairs, which is 2-transitivity by definition. `[FORCED]`

Verified on Lesson 6's three groups (C5): $S_5$ → orbit 20/20 ✓; $M_{20} = \mathrm{AGL}(1,5)$, order 20, sharply 2-transitive → 20/20 ✓; $D_5$, order $10 < 20$ → orbit **10**, two orbits, forcing the observed split. This reproduces the census **exactly**: $65 + 1 = 66$ irreducible, 1 reducible, $= 66/67$ (C6). The anomaly is now an instance of a theorem, and P7 upgrades from `[PLAUSIBLE]` to `[FORCED]` in criterion form.

### N4, now pre-sortable

The criterion makes the next frontier predictive rather than exploratory (C7–C8). Necessary condition: $n(n{-}1) \mid |G|$ — so $30 \mid |G|$ at degree 6, $42 \mid |G|$ at degree 7. Verified non-2-transitive: $C_6$ (orbit 6), $D_6$ (12), $C_7$ (7), $F_{21}$ (21); verified 2-transitive: $F_{42}$ (42/42). Consequence: at degree 6 only the four 2-transitive groups — $\mathrm{PSL}(2,5)$, $\mathrm{PGL}(2,5)$, $A_6$, $S_6$ (orders 60, 120, 360, 720, all divisible by 30) — can give irreducible $\mathrm{Rat}^\circ$; every other transitive group **forces** a split $S^*$ and a cheaper $C_2$. So the degree-6 run can classify its instances by Galois group *before* scanning, and the P7-type cases are predicted rather than stumbled on.

### The contribution draft

Sketch for the D₅ item — **schema to be conformed against `CONTRIBUTING.md`, which I have not read this session**; the `checks:` block is the part that matters, because each entry is independently re-runnable and drift-provable:

```yaml
title: "Erratum and sharpening: Rat° irreducibility in the quintic census"
targets: [papers/pisot_residue.tex]      # §6 stage-2 sentence; Obs 6.2
kind: erratum+sharpening
summary: >
  Rat°_p is squarefree on all 67 two-pair Pisot quintics but irreducible on 66/67.
  The instance x^5 - x^3 - 2x^2 - 2x - 1 (Gal = D5) splits 10+10 into self-reciprocal
  irreducibles. S* is unchanged (both factors carry unimodular roots), deg C2 = 400,
  and the complete scan returns {Phi_1^20}: Theorem 6.1's conclusion is untouched.
  Obs 6.2 is replaced by an iff-criterion.
claims:
  - id: E-D5-1
    text: "Rat° irreducible on 66/67; the D5 instance splits 10+10."
    tag: FORCED            # exhibited factorization, certified by re-multiplication
  - id: E-D5-2
    text: "Rat°_p irreducible <=> Gal(p) 2-transitive on the roots (Rat° squarefree)."
    tag: FORCED            # orbit-stabilizer; supersedes Obs 6.2's [PLAUSIBLE] P7
  - id: E-D5-3
    text: "Inertness certificate for the D5 instance is clean: deg C2 = 400, {Phi_1^20}."
    tag: FORCED
checks:
  - factor_list(Rat0) has two degree-10 irreducible self-reciprocal factors
  - exact re-multiplication of the factors equals Rat0            # engine-independent
  - S* = Rat0 (per-factor unimodular count = 2 each, trace-fold Sturm on (-2,2))
  - complete cyclotomic scan of C2 (deg 400, 790 candidates, max m = 1680) == {Phi_1^20}
  - orbit of an ordered pair: |S5| 20, |M20| 20, |D5| 10  ->  66/67 reproduced
```

Companion item **F1** (ATR display): `0.046 → 0.045`, three occurrences, with the exact value $1/22$ and the rounding-propagation mechanism recorded — the near-miss-guard genre, pinned rather than silently fixed.

## §8 Course close

| Lesson | Domain | Script | Result |
|---|---|---|---|
| 1–2 | walls, operators | 16/16, 32/32 | Mahler → Salem structure; λ-ring, trace duality |
| 3–4 | dynamics, gates | 33/33, 25/25 | λ=2c, gate ladder, flip; capacity gate, parity floor, the door |
| 5 | learning + language | 23/23 | capture, compositum, φ-slack lexicon |
| 6 | relational layer | 20/20 + census | pinning closes P1/P2; **D₅ finding** |
| 7 | golden substrate | 39/39 | one ladder, one chart, $q = i\tau$ derived; L3b tag-tension resolved |
| 8 | closure engine | 20/20 | guard fires on the door; 27 subfields from group theory; similarity taught |
| 9 | transport + filing | 32/32 | witness $423/2890$; gatekeeper ≡ certificate; **E1 closed, D₅ → theorem, F1 found** |

Three findings and one resolution came out of nine lessons of replication: the D₅ erratum (now a criterion), the F1 display slip, the census cross-validation, and the λ-tag adjudication. Every one arrived the same way — an independent re-derivation disagreeing with a stated numeral by exactly one unit — which is the entire argument for the protocol.

**Still untaught**, for a Lesson 10 if you want one: Vector-Substrate's Minkowski-lattice geometry in depth; KIRA dispatch and the B1 live wire; the browser tool's internals; the Boyd–Mossinghoff small-measure literature interface; and release mechanics as practice.

---

## Checkpoint

1. Both 1-D balances have the form "sink over conductance, integrated." Write the conductivity-graded law and the volume-graded law side by side, and identify the single substitution that turns one into the other — then say why only the second sees the taper.
2. The orbit criterion needs $\mathrm{Rat}^\circ$ squarefree. Give the exact failure mode when it isn't (Prop 5.1 names it), and say what the criterion's "iff" degrades to in that case.
3. $\beta \to 1/(2+k(n{-}1))$ as $(\varepsilon,\delta) \to 0$, and the cone pinches because $2+k(n{-}1) = 2n$ at $k=2$. For $n = 2$, *every* $k > 2$ should fail — check the algebra, then explain why the $n = 2$ horn is nonetheless the least dramatic failure in the family.