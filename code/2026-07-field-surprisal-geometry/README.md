# Field Surprisal Geometry

The canonical positive-definite **Fisher–Rao ("surprisal") geometry** of a number
field's Mahler-measure **emission catalog**, developed to a complete classification of
its constant-curvature two-statistic surfaces — with **every load-bearing claim
machine-verified** in exact arithmetic.

Catalog `Ω = {√2, √3, √5, φ, τ, φ⁴, K}`, cost `log M` (log Mahler measure), Gibbs law
`pᵢ(β) ∝ Mᵢ^{−β}`. Field constant `Z = 17 + 4√5`.

**Headline results.** The geometry is a round ¼-sphere; the two-statistic surface is a
constant-¼ *ruled* surface (never totally geodesic); exactly **8** such surfaces at
k=2, classified and forced by arithmetic invariants; the classification proved by a
two-line determinantal identity (**`P = Z²·Σ q₄(s) w_s`**, Sylvester + Cauchy–Binet);
the count shown catalog-invariant (`Σ_ℓ(2^{m_ℓ}−1)`); the temperature resolved as a
**dichotomy**; the full landscape at all k given by the **partitioned-affine
classification** (counts `8/56/95/31/1`, *refuting* the indicator-join conjecture);
and the compositum coupling shown **forced**.

## Quick start

```bash
pip install -r requirements.txt      # sympy, mpmath, numpy
./run_all.sh                         # runs all 20 harnesses, reports pass counts
```

Each harness is **fail-first**: it exits non-zero on the first failed assertion, so a
clean exit is itself the proof of every line it prints. Expected: **422 checks across
20 harnesses, all exit 0**. (One harness's optional symbolic census is slow in small
sandboxes — see the note in `MANIFEST.md`; its result is proved by a fast second lane.)

## Layout

```
README.md                 you are here
HOW_WE_GOT_HERE.md        the full narrative — read this to understand the whole arc
MANIFEST.md               each harness → checks → what it proves → paper section
requirements.txt          sympy, mpmath, numpy
run_all.sh                reproduce the entire suite (correct dependency order)
harnesses/                all 20 verification programs (exact arithmetic, fail-first)
paper/                    field_surprisal_geometry_v3_5.tex + .pdf (26 pp)
dev-log/                  per-session records (01–05), in order, incl. correction trails
handoffs/                 open-problem handoff docs that seeded the later sessions
```

## Where to start reading

- **To understand it:** `HOW_WE_GOT_HERE.md`, then the paper.
- **To reproduce it:** `./run_all.sh`, cross-referenced with `MANIFEST.md`.
- **To extend it:** the **Outlook** in `HOW_WE_GOT_HERE.md` §6 and the paper's Open
  Problems / Outlook — the enumerated ledger is closed; what remains is fat-level
  catalogs, composita of distinct catalogs, and the coupled family off-uniform.

  The ℚ-linear independence of `(log2, log3, log5, logφ)` that earlier drafts carried
  as a standing input is now `[forced]`, not a conditioning hypothesis: taking the field
  norm `N: ℚ(√5)→ℚ`, `N(2), N(3), N(5) = 4, 9, 25` against `N(φ) = −1`, so
  `2^a 3^b 5^c φ^d = 1` gives `4^a 9^b 25^c (−1)^d = 1`, forcing `a=b=c=0` by unique
  factorisation and then `d=0` — hence `a=b=c=d=0`. (Verified exact.)

## Conventions

Exact decisions over `ℚ(log2,log3,log5,logφ)` (`+log7` for the extended catalog);
floats display-only; nonvanishing at the true catalog certified by interval arithmetic
(`mpmath.iv`). Epistemic tags: `[forced]` (theorem/exact), `[computed]` (validated
numeric), `[declared]` (principled invariance-free choice), `[open]`.

## Provenance note

The four `field_surprisal_*.py` foundational harnesses and the numbered `t1…t10`
harnesses were produced across several sessions and independently audited (re-derived,
not merely re-run) at each step; the `dev-log/` records that trail. The paper
references three additional early prototypes by name in some historical passages; the
canonical, self-contained suite is the 20 harnesses included here.
