# Code -- *The Generative Content of a Conserved Emptiness*

Producer scripts for `papers/2026-06-generative-emptiness/generative_emptiness.tex`
(*Kinematic Voids as Superselection Generators: the Salem Slot and the Five
Objects Its Charge Produces*, AceTheDactyl / Echo S Studios).

Each script re-derives one of the paper's results from the emission-algebra
premises (Def. 2.1 seeds `phi = x^2-x-1`, `K = x^4+5x^2-5`, and the three
operators `(x)` tensor, `( )^2` squaring, `(+)` direct sum) and **writes a
machine-readable artifact** into `data/2026-06-generative-emptiness/`.

These are **producers**: functions plus a `main()` that *emits* data. They are
independent of, and do not import from, the verifier suite under
`tests/2026-06-generative-emptiness/` (which instead *asserts* against the
paper's stated values). The shared engine `ge_core.py` re-implements the
operators, the Mahler measure, and the Z/4Z angle charge from scratch.

Requirements: `py` launcher with `sympy`, `mpmath` (data); `matplotlib` is not
needed, but `pymupdf` (`import fitz`) + MiKTeX XeLaTeX are needed for figures.

## Scripts

| Script | Run command | Paper result produced | Output |
|---|---|---|---|
| `ge_core.py` | `py code/2026-06-generative-emptiness/ge_core.py` | Shared engine (operators, Mahler measure, Z/4Z charge; Def. 2.1). Self-check only, emits no data. | (console) |
| `object1_grading.py` | `py code/2026-06-generative-emptiness/object1_grading.py` | **Object I** -- the Z/4Z grading (Thm 2.2): charge multisets, operator action add/double/union, Cayley closure, orbit stays on-lattice, K irreducible. | `data/.../object1_grading.json` |
| `object2_content.py` | `py code/2026-06-generative-emptiness/object2_content.py` | **Object II** -- content polynomial `x^4-1 = Phi_1 Phi_2 Phi_4` (Prop 3.1, Rem 3.2): root lattice, realised `+-1`, K's place `+-i*beta` (beta=2.4195), full Z/4Z in charge / Z/2Z on circle. | `data/.../object2_content.json` |
| `object3_gap.py` | `py code/2026-06-generative-emptiness/object3_gap.py` | **Object III** -- the gap / cost floor `phi` (Prop 4.1, Cor 4.2, Prop 7.1(3)): quadratic enumeration, empty band `(1,phi)`, floor at `x^2-x-1`, first realised values, smallest Perron, floor-set closure. | `data/.../object3_gap.json`, `data/.../object3_gap_quadratics.csv` |
| `object4_normalform.py` | `py code/2026-06-generative-emptiness/object4_normalform.py` | **Object IV** -- graded normal form (Prop 5.1): the three diagnostic factorizations with sector classification and grow measures; `M(G)=M(P)`; the `a*b=sqrt5` imaginary sector. | `data/.../object4_normalform.json` |
| `object5_current.py` | `py code/2026-06-generative-emptiness/object5_current.py` | **Object V** -- conserved current (Prop 6.1, 6.2): two-generation orbit (block-diagonal, clean radial growth) and the seven listed measures `phi,phi^2,phi^4,46.98,76.63,122.99,8049.92`. | `data/.../object5_measures.json`, `data/.../object5_orbit.csv` |
| `minimality.py` | `py code/2026-06-generative-emptiness/minimality.py` | **Minimality chain** (Prop 7.1) + scope/ledger identities: ternary lock `d^2-d+1=3 iff d=2`, spectrum `{-sqrt5,0,+sqrt5}`, `(pi/2)Z ~= Z/4Z`, `Q(5^{1/4})` contains `i`, `sqrt5 = phi + phi^{-1}`. | `data/.../minimality.json` |
| `make_figures.py` | `py code/2026-06-generative-emptiness/make_figures.py` | Renders the paper's one `tikzpicture` (Figure 1, `fig:charge`) to PNG+PDF using the paper's own preamble. | `figures/.../figure1.{png,pdf}` |

## Reproduce everything

```
py code/2026-06-generative-emptiness/object1_grading.py
py code/2026-06-generative-emptiness/object2_content.py
py code/2026-06-generative-emptiness/object3_gap.py
py code/2026-06-generative-emptiness/object4_normalform.py
py code/2026-06-generative-emptiness/object5_current.py
py code/2026-06-generative-emptiness/minimality.py
py code/2026-06-generative-emptiness/make_figures.py
```
