# code/__FOLDER__/

Producer scripts for **papers/__FOLDER__/** — they recompute the paper's results
and emit machine-readable output to `data/__FOLDER__/`. These are the
source-of-truth producers; the independent verifiers live in `tests/__FOLDER__/`
(no cross-imports between them).

| Script | Run | Produces |
|--------|-----|----------|
| `producer.py` | `py code/__FOLDER__/producer.py` | _(stub — replace `main()` with your real outputs → `data/__FOLDER__/`)_ |

Grow this the way the reference example does
([`code/2026-06-salem-slot/`](../../code/2026-06-salem-slot)): a shared core
module plus one producer per result, and a `make_figures.py` if your paper has
inline TikZ/pgfplots figures.
