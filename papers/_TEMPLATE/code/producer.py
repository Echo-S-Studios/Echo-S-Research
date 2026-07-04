"""Producer stub for code/__FOLDER__/ (source-of-truth for papers/__FOLDER__/).

Producers COMPUTE and EMIT machine-readable output to data/__FOLDER__/ — they do
NOT assert (that is the job of the independent verifiers in tests/__FOLDER__/).
Keep code/ and tests/ independent (no cross-imports). Stamp every output with its
provenance. Mirror the reference example code/2026-06-salem-slot/ (a shared core
module + per-result producer scripts + optional make_figures.py).

Run:  py code/__FOLDER__/producer.py
"""

import json
import os

PAPER_TEX = "papers/__FOLDER__/__SHORTNAME__.tex"
GENERATED_BY = "code/__FOLDER__/producer.py"
DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data",
    "__FOLDER__",
)


def write_json(name, payload):
    """Emit a JSON artifact with provenance fields at the top."""
    os.makedirs(DATA_DIR, exist_ok=True)
    out = {"_source_paper": PAPER_TEX, "_generated_by": GENERATED_BY}
    out.update(payload)
    with open(os.path.join(DATA_DIR, name), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


def main():
    # TODO: compute your paper's results and emit them, e.g.:
    #   write_json("my_result.json", {"claim": "…", "value": 42})
    # Until you replace this, the producer emits nothing, so the drift CI stays
    # green (no committed data yet to compare against).
    print("code/__FOLDER__/producer.py: stub — no outputs yet. Replace main().")


if __name__ == "__main__":
    main()
